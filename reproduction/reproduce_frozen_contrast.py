"""Reproduce the paper's frozen-model provenance contrast end to end.

The script starts from public, hash-verified CATH data and a pinned ProteinMPNN
revision. It recomputes logits, constructs target-derived and backbone-derived
channels at approximately matched marginal decodability, fits the injection
coefficient on the training population, and evaluates the held-out contrast.

    fetch + hash-verify CATH  ->  clone ProteinMPNN  ->  download released
    weights  ->  recompute logits for the 928-chain test population  ->
    rebuild target/backbone codebooks  ->  chain bootstrap  ->  compare against
    the recorded values.

By default, generated data and the external checkout live under `.cache/` next
to this script. Set `SURFACE_PROVENANCE_CACHE` to use another cache directory.
The full run takes approximately 15--25 minutes on CPU.

    python3 reproduction/reproduce_frozen_contrast.py
    python3 reproduction/reproduce_frozen_contrast.py --smoke-test
"""
import os, sys, json, hashlib, subprocess, urllib.request, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.environ.get("SURFACE_PROVENANCE_CACHE", os.path.join(HERE, ".cache"))
WORK = os.path.join(CACHE, "external")
DATA = os.path.join(CACHE, "data")
PROTEINMPNN_COMMIT = "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57"
PROTEINMPNN_WEIGHT_SHA256 = "c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd"
ALPHA = "ACDEFGHIKLMNPQRSTVWY"; AAI = {a: i for i, a in enumerate(ALPHA)}
KD = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"G":-0.4,"T":-0.7,"S":-0.8,
      "W":-0.9,"Y":-1.3,"P":-1.6,"H":-3.2,"E":-3.5,"Q":-3.5,"D":-3.5,"N":-3.5,"K":-3.9,"R":-4.5}
MAXLEN, MINLEN = 260, 40
EXPECT = {"chain_set.jsonl": "1944d21b975a11543b1b92fe906a0bb4a9a276a52ca3d493e5aadafd09f464c1",
          "chain_set_splits.json": "8e9a587a50c7f6c026e4ed00f6c1c30b106100f36f7a01de47542bdfc060adc2"}
# Values reported in the paper.
TARGETS = {"base_pmpnn": 46.07, "contrast_cath_pmpnn": 12.49, "backbone_gain": 0.00}
ok = fail = 0


def chk(label, got, want, tol):
    global ok, fail
    good = abs(float(got) - float(want)) <= tol
    print(f"  {'OK  ' if good else 'FAIL'} {label:44s} got {float(got):8.3f}  recorded {float(want):8.3f}  (tol {tol})")
    ok += good; fail += (not good)


def ensure_data():
    os.makedirs(DATA, exist_ok=True)
    base = "https://people.csail.mit.edu/ingraham/graph-protein-design/data/cath/"
    for fn, want in EXPECT.items():
        p = os.path.join(DATA, fn)
        if not os.path.exists(p):
            print(f"  fetching {fn} ...", flush=True)
            partial = p + ".partial"
            try:
                urllib.request.urlretrieve(base + fn, partial)
                os.replace(partial, p)
            except Exception as e:
                print(f"  fetch failed: {e}")
                return False
            finally:
                if os.path.exists(partial):
                    os.remove(partial)
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        if h != want:
            print(f"  HASH MISMATCH {fn}: {h[:16]} != {want[:16]}"); return False
        print(f"  {fn} sha256 verified")
    return True


def ensure_pmpnn():
    repo = os.path.join(WORK, "ProteinMPNN")
    if not os.path.exists(repo):
        os.makedirs(WORK, exist_ok=True)
        print("  cloning ProteinMPNN ...", flush=True)
        steps = [
            ["git", "init", "-q", repo],
            ["git", "-C", repo, "remote", "add", "origin",
             "https://github.com/dauparas/ProteinMPNN.git"],
            ["git", "-C", repo, "fetch", "-q", "--depth", "1", "origin",
             PROTEINMPNN_COMMIT],
            ["git", "-C", repo, "checkout", "-q", "--detach", "FETCH_HEAD"],
        ]
        for cmd in steps:
            r = subprocess.run(cmd, capture_output=True)
            if r.returncode != 0:
                print("  checkout failed:", r.stderr.decode()[:300]); return None
    head = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    if head != PROTEINMPNN_COMMIT:
        print(f"  commit mismatch: {head} != {PROTEINMPNN_COMMIT}")
        return None
    w = os.path.join(repo, "vanilla_model_weights", "v_48_020.pt")
    if not os.path.exists(w): print("  weights missing in clone"); return None
    weight_hash = hashlib.sha256(open(w, "rb").read()).hexdigest()
    if weight_hash != PROTEINMPNN_WEIGHT_SHA256:
        print(f"  weight hash mismatch: {weight_hash}")
        return None
    print(f"  ProteinMPNN commit and weights verified ({head[:12]})")
    return repo


def backbone_desc(N, CA, C, O):
    n = len(CA)
    def dih(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1n = b1 / (np.linalg.norm(b1, axis=-1, keepdims=True) + 1e-8)
        v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
        w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
        return np.degrees(np.arctan2((np.cross(b1n, v) * w).sum(-1), (v * w).sum(-1)))
    phi = np.zeros(n); psi = np.zeros(n)
    if n > 1:
        phi[1:] = dih(C[:-1], N[1:], CA[1:], C[1:]); psi[:-1] = dih(N[:-1], CA[:-1], C[:-1], N[1:])
    d = np.linalg.norm(CA[:, None, :] - CA[None, :, :], axis=-1); eye = np.eye(n, dtype=bool)
    dens = [((d < r) & ~eye).sum(1).astype(float) for r in (8., 10., 12., 14.)]
    b = CA - N; c = C - CA; a = np.cross(b, c)
    cb = -0.58273431 * a + 0.56802827 * b - 0.54067466 * c
    cb /= np.linalg.norm(cb, axis=-1, keepdims=True) + 1e-8
    proj = ((CA[None, :, :] - CA[:, None, :]) * cb[:, None, :]).sum(-1)
    near = (d < 13.0) & ~eye
    return np.stack([np.sin(np.radians(phi)), np.cos(np.radians(phi)),
                     np.sin(np.radians(psi)), np.cos(np.radians(psi)), *dens,
                     (near & (proj > 0)).sum(1).astype(float),
                     (near & (proj <= 0)).sum(1).astype(float)], 1)


def logP(codes, y, ncode):
    M = np.full((ncode, 20), 0.5); np.add.at(M, (codes, y), 1.0)
    return np.log(M / M.sum(1, keepdims=True))


def stack_eval(Gtr, ytr, ctr, Gte, yte, cte, nc):
    lp = logP(ctr, ytr, nc); best = (0, -1)
    for lam in [0, .25, .5, 1, 2, 4, 8, 16, 32, 64]:
        a = ((Gtr + lam * lp[ctr]).argmax(1) == ytr).mean()
        if a > best[1]: best = (lam, a)
    return ((Gte + best[0] * lp[cte]).argmax(1) == yte)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chains", type=int, default=0,
                    help="limit held-out chains; 0 uses the full population")
    ap.add_argument("--smoke-test", action="store_true",
                    help="run four training and two test chains; checks execution, not paper values")
    a = ap.parse_args()
    print("FROZEN-MODEL PROVENANCE CONTRAST\n")
    if not ensure_data(): sys.exit(1)
    repo = ensure_pmpnn()
    if repo is None: sys.exit(1)
    sys.path.insert(0, repo)
    import torch
    from protein_mpnn_utils import ProteinMPNN, tied_featurize
    dev = torch.device("cpu")
    ck = torch.load(os.path.join(repo, "vanilla_model_weights", "v_48_020.pt"),
                    map_location=dev, weights_only=False)
    mdl = ProteinMPNN(ca_only=False, num_letters=21, node_features=128, edge_features=128,
                      hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                      augment_eps=0.0, k_neighbors=ck["num_edges"])
    mdl.load_state_dict(ck["model_state_dict"]); mdl.eval()

    sp = json.load(open(f"{DATA}/chain_set_splits.json"))
    te, tr = set(sp["test"]), set(sp["train"])
    def collect(keep, cap):
        G, Y, D, C = [], [], [], []; k = 0
        with torch.no_grad():
            for line in open(f"{DATA}/chain_set.jsonl"):
                e = json.loads(line)
                if e["name"] not in keep: continue
                s = e["seq"]
                if len(s) > MAXLEN or len(s) < MINLEN or (set(s) - set(ALPHA)): continue
                ca = np.asarray(e["coords"]["CA"], float); m = ~np.isnan(ca[:, 0])
                co = {x: np.asarray(e["coords"][x], float)[m] for x in ("N", "CA", "C", "O")}
                s2 = "".join(np.array(list(s))[m])
                if len(s2) < MINLEN or np.isnan(np.concatenate(list(co.values()))).any(): continue
                try:
                    b = [{"seq_chain_A": s2, "coords_chain_A":
                          {f"{x}_chain_A": co[x].tolist() for x in ("N", "CA", "C", "O")},
                          "name": e["name"], "num_of_chains": 1, "seq": s2}]
                    o = tied_featurize(b, dev, {e["name"]: (["A"], [])}, None, None, None, None, None,
                                       ca_only=False)
                    lg = mdl(o[0], o[1], o[2], o[4] * o[10], o[12], o[5], torch.randn(o[4].shape))
                    G.append(lg[0].numpy().astype(np.float32)[:, :20])
                    Y.append(np.array([AAI[c] for c in s2]))
                    D.append(backbone_desc(co["N"], co["CA"], co["C"], co["O"]))
                    C.append(np.full(len(s2), k)); k += 1
                except Exception:
                    continue
                if cap and k >= cap: break
                if k % 100 == 0 and k: print(f"    {k} chains", flush=True)
        return [np.concatenate(x) for x in (G, Y, D, C)], k
    train_cap = 4 if a.smoke_test else 593
    test_cap = 2 if a.smoke_test else a.chains
    print("\n  recomputing TRAIN logits (fit population) ...", flush=True)
    (Gtr, ytr, Dtr, _), ktr = collect(tr, train_cap)
    print(f"  train: {ktr} chains")
    print("  recomputing TEST logits ...", flush=True)
    (Gte, yte, Dte, Cte), kte = collect(te, test_cap)
    print(f"  test: {kte} chains, {len(yte)} residues\n")

    base_corr = (Gte.argmax(1) == yte); base = base_corr.mean() * 100
    if not a.smoke_test:
        chk("frozen ProteinMPNN base (%)", base, TARGETS["base_pmpnn"],
            1.0 if a.chains else 0.3)

    from sklearn.cluster import MiniBatchKMeans
    mu, sg = Dtr.mean(0), Dtr.std(0) + 1e-8
    km = MiniBatchKMeans(16, random_state=0, n_init=3, batch_size=4096).fit((Dtr - mu) / sg)
    rs = np.random.RandomState(201); mp = {x: int(v) for x, v in zip(ALPHA, rs.randint(0, 2, 20))}
    tk = sorted(set(mp.values())); tm = {c: i for i, c in enumerate(tk)}
    tgt = stack_eval(Gtr, ytr, np.array([tm[mp[ALPHA[v]]] for v in ytr]), Gte, yte,
                     np.array([tm[mp[ALPHA[v]]] for v in yte]), len(tk))
    bbn = stack_eval(Gtr, ytr, km.predict((Dtr - mu) / sg), Gte, yte,
                     km.predict((Dte - mu) / sg), 16)
    gb = (bbn.mean() - base_corr.mean()) * 100
    con = (tgt.mean() - bbn.mean()) * 100
    if a.smoke_test:
        finite = np.isfinite([base, gb, con]).all()
        print(f"\n  smoke outputs: base={base:.2f}%, backbone gain={gb:+.2f} pp, "
              f"target-backbone={con:+.2f} pp")
        print("\nSMOKE TEST PASSED" if finite else "\nSMOKE TEST FAILED")
        sys.exit(0 if finite else 1)
    chk("backbone-derived gain (pp)", gb, TARGETS["backbone_gain"], 0.5)
    chk("target - backbone contrast (pp)", con, TARGETS["contrast_cath_pmpnn"],
        1.5 if a.chains else 0.6)

    rng = np.random.RandomState(0); uq = np.unique(Cte)
    idx = {c: np.where(Cte == c)[0] for c in uq}
    d = []
    for _ in range(1000):
        sampled = rng.choice(uq, len(uq), True)
        ii = np.concatenate([idx[c] for c in sampled])
        d.append((tgt[ii].mean() - bbn[ii].mean()) * 100)
    lo, hi = np.percentile(d, [2.5, 97.5])
    print(f"\n  chain bootstrap 95% CI on the contrast: [{lo:+.2f}, {hi:+.2f}]")
    print(f"\n{ok} OK, {fail} FAIL")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
