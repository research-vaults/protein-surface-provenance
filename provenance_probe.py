#!/usr/bin/env python3
"""provenance_probe — audit dependence on substituted target labels.

VERSION 0.3.0. A with/without-channel score does not identify how the channel
was constructed. `intervene` compares supplied channel values after labels-only
substitution; `screen` measures association and cannot establish provenance.

The caller must hold every other preprocessing input, correspondence and random
state fixed. This is a software-input intervention, not a physical mutation.
Any exact channel change is evidence of target dependence conditional on that
contract. Unchanged values do not certify absence of dependence: the chosen
substitution may preserve the channel's equivalence classes. Full nonempty geom
fields are required for a controlled verdict, and only the supplied fields can
be checked. This tool cannot establish completeness of the caller's fixture.

Usage: provenance_probe.py intervene --original A.jsonl --substituted B.jsonl
       provenance_probe.py screen --train T.jsonl --test E.jsonl [--pred P.txt]
Intervention row: {"chain":"id", "channel":[v,...], "geom":[g,...]}
Screen row: {"seq":"ACDE...", "channel":[v,...], "resid":[i,...]}
Licence: BSD-3-Clause. Standard library only.
"""
import argparse, collections, difflib, json, sys

__version__ = "0.3.0"


def keyof(v):
    return json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else v


def load(path):
    with open(path) as f:
        return [json.loads(l) for l in f if l.strip()]


# ================================================================ intervene ===
def intervene(orig_path, sub_path):
    """Exact comparison of supplied values; scientific verdict is conditional."""
    A, B = load(orig_path), load(sub_path)
    if len(A) != len(B):
        sys.exit(f"chain count differs: {len(A)} vs {len(B)}")
    if not A:
        sys.exit("empty intervention fixture")
    chain_ids = [r.get("chain") for r in A]
    if any(x is None for x in chain_ids) or len(set(chain_ids)) != len(chain_ids):
        sys.exit("chain IDs must be present and unique")
    geometry_chains = 0
    moved = total = 0
    geom_moved = geom_total = 0
    per_chain = []
    for a, b in zip(A, B):
        if a.get("chain") != b.get("chain"):
            sys.exit(f"chain id mismatch: {a.get('chain')!r} vs {b.get('chain')!r}")
        ca, cb = a["channel"], b["channel"]
        if len(ca) != len(cb):
            sys.exit(f"chain {a.get('chain')}: channel length differs "
                     f"({len(ca)} vs {len(cb)}) -- the structures are not matched")
        if not ca:
            sys.exit(f"chain {a.get('chain')}: empty channel")
        m = sum(keyof(x) != keyof(y) for x, y in zip(ca, cb))
        moved += m; total += len(ca)
        per_chain.append(m / len(ca) if ca else 0.0)
        ga, gb = a.get("geom"), b.get("geom")
        if ga is not None and gb is not None:
            if len(ga) != len(gb):
                sys.exit(f"chain {a.get('chain')}: geometry length differs")
            geometry_chains += bool(ga)
            geom_moved += sum(keyof(x) != keyof(y) for x, y in zip(ga, gb))
            geom_total += len(ga)
    per_chain.sort()
    n = len(per_chain)
    out = dict(mode="intervene", exact=True,
               n_chains=n, n_values=total,
               channel_moved_frac=moved / total if total else 0.0,
               channel_moved_mean_per_chain=sum(per_chain) / n if n else 0.0,
               n_changed_values=moved,
               geometry_checked=geometry_chains == n,
               n_geometry_chains=geometry_chains,
               geometry_coverage=geometry_chains / n,
               geometry_moved_frac=(geom_moved / geom_total) if geom_total else None,
               n_geometry_values=geom_total)
    if geom_moved:
        out["verdict"] = "INVALID_INTERVENTION"
    elif geometry_chains != n:
        out["verdict"] = "INCONCLUSIVE_UNVERIFIED_CONTROL"
    elif moved:
        out["verdict"] = "TARGET_DEPENDENCE_DETECTED"
    else:
        out["verdict"] = "NO_DEPENDENCE_DETECTED"
    out["interpretation"] = (
        "Conditional on labels being the only changed input, fixed random state and "
        "stable correspondence. Geometry coverage refers only to supplied fields. "
        "No finite unchanged substitution certifies C=f(B), and target dependence "
        "alone does not imply positive I(Y;C|B) or use by a trained model.")
    return out


# =================================================================== screen ===
def screen_association(rows):
    """Observational only. Reports how strongly channel value is associated with
    residue identity. NOT a substitution rate and NOT a provenance verdict."""
    pair = collections.defaultdict(collections.Counter)
    for r in rows:
        for v, i in zip(r["channel"], r["resid"]):
            if 0 <= i < len(r["seq"]):
                pair[r["seq"][i]][keyof(v)] += 1
    if not pair:
        return dict(association=0.0, n_classes=0, note="no (value, residue) pairs")
    n_a = {a: sum(c.values()) for a, c in pair.items()}
    N = sum(n_a.values())
    classes = {k for c in pair.values() for k in c}
    P = {a: {k: pair[a][k] / n_a[a] for k in classes} for a in pair}
    pa = {a: n_a[a] / N for a in pair}
    same = sum(pa[a] * sum(P[a][k] ** 2 for k in classes) for a in pair)
    num = den = 0.0
    for a in pair:
        for b in pair:
            if a == b:
                continue
            w = pa[a] * pa[b]
            num += w * sum(P[a][k] * P[b][k] for k in classes)
            den += w
    diff = num / den if den else 0.0
    return dict(association=max(same - diff, 0.0),
                p_same_class_same_residue=round(same, 4),
                p_same_class_different_residues=round(diff, 4),
                n_classes=len(classes), n_residue_types=len(pair))


def standalone_ceiling(train, test):
    cls = collections.defaultdict(collections.Counter)
    unc = collections.Counter()
    for r in train:
        named = set()
        for v, i in zip(r["channel"], r["resid"]):
            if 0 <= i < len(r["seq"]):
                cls[keyof(v)][r["seq"][i]] += 1; named.add(i)
        for i, ch in enumerate(r["seq"]):
            if i not in named:
                unc[ch] += 1
    rule = {k: c.most_common(1)[0][0] for k, c in cls.items()}
    fb = unc.most_common(1)[0][0] if unc else next(iter(rule.values()))
    corr = tot = cov = 0
    per = []
    for r in test:
        named = {i: keyof(v) for v, i in zip(r["channel"], r["resid"])
                 if 0 <= i < len(r["seq"])}
        c = sum((rule.get(named[i], fb) if i in named else fb) == ch
                for i, ch in enumerate(r["seq"]))
        corr += c; tot += len(r["seq"]); cov += len(named); per.append(c / len(r["seq"]))
    per.sort(); n = len(per)
    return dict(micro=corr / tot, mean_per_chain=sum(per) / n,
                n_classes=len(rule), coverage=cov / tot, n_chains=n, n_residues=tot)


def margin(test, pred_path, ceiling):
    preds = [l.strip() for l in open(pred_path)]
    if len(preds) != len(test):
        sys.exit(f"--pred has {len(preds)} lines, test has {len(test)} chains")
    corr = tot = 0
    for r, p in zip(test, preds):
        sm = difflib.SequenceMatcher(None, p, r["seq"], autojunk=False)
        corr += sum(bl.size for bl in sm.get_matching_blocks()); tot += len(r["seq"])
    return dict(model_micro=corr / tot, ceiling_micro=ceiling["micro"],
                margin_pp=100 * (corr / tot - ceiling["micro"]))


# ===================================================================== cli ====
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"provenance_probe {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("intervene", help="Exact comparison after controlled label substitution")
    pi.add_argument("--original", required=True)
    pi.add_argument("--substituted", required=True)
    pi.add_argument("--json")

    ps = sub.add_parser("screen", help="OBSERVATIONAL screen; cannot identify provenance")
    ps.add_argument("--train", required=True)
    ps.add_argument("--test", required=True)
    ps.add_argument("--pred")
    ps.add_argument("--json")

    a = ap.parse_args()
    print("=" * 70)
    print(f"provenance probe {__version__}")
    print("=" * 70)

    if a.cmd == "intervene":
        r = intervene(a.original, a.substituted)
        print(f"\nMODE: intervene  —  exact supplied-value comparison, {r['n_values']:,} values "
              f"over {r['n_chains']} chains\n")
        print(f"  channel values that moved    {100*r['channel_moved_frac']:7.2f}%")
        if r["geometry_checked"]:
            print(f"  geometry values that moved   {100*r['geometry_moved_frac']:7.2f}%"
                  f"   ({r['n_geometry_values']:,} checked)")
        else:
            print("  geometry values that moved       n/a   (no `geom` field supplied;")
            print("                                          complete invariance was not verified)")
        print()
        print(f"  VERDICT: {r['verdict']}")
        print(f"  Geometry coverage: {r['n_geometry_chains']}/{r['n_chains']} chains")
        print(f"  {r['interpretation']}")
        rep = r

    else:
        train, test = load(a.train), load(a.test)
        assoc = screen_association(train)
        ceil = standalone_ceiling(train, test)
        rep = dict(mode="screen", exact=False,
                   observational_association=assoc, standalone_ceiling=ceil,
                   metric_note="standalone_ceiling is a legacy key for train-fitted lookup accuracy, not a population Bayes bound")
        print("\nMODE: screen  —  OBSERVATIONAL. This cannot identify provenance.\n")
        print(f"  value/residue association    {100*assoc['association']:7.2f}%")
        print(f"    P(same class | same residue)        {assoc['p_same_class_same_residue']:.4f}")
        print(f"    P(same class | different residues)  {assoc['p_same_class_different_residues']:.4f}")
        print(f"\n  train-fitted lookup          {100*ceil['micro']:7.2f}%  micro"
              f"   ({100*ceil['mean_per_chain']:.2f}% mean per chain)")
        print(f"    distinct classes                    {ceil['n_classes']}")
        print(f"    coverage                            {100*ceil['coverage']:.2f}% of residues named")
        if a.pred:
            m = margin(test, a.pred, ceil); rep["margin"] = m
            print(f"\n  model                        {100*m['model_micro']:7.2f}%")
            print(f"  margin over channel          {m['margin_pp']:+7.2f} pp")
            if m["margin_pp"] <= 0:
                print("    The trained model does not beat a parameter-free lookup on its")
                print("    own input channel.")
        print("\n  NO PROVENANCE VERDICT IS ISSUED.")
        if assoc["association"] > 0.5:
            print("  High association. This is consistent with a target-derived channel AND")
            print("  with a design-valid channel that happens to correlate with residue")
            print("  identity. Build the intervention fixture and run `intervene`.")
        else:
            print("  Low association. A target-derived channel is unlikely but not excluded.")
            print("  Inspect the construction and run controlled substitutions.")

    if a.json:
        json.dump(rep, open(a.json, "w"), indent=2)
        print(f"\nreport written to {a.json}")
    print()


if __name__ == "__main__":
    main()
