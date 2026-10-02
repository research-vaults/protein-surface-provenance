# Provenance-aware evaluation of a computed model input

A with/without-channel score measures performance under a specified learner. It
does not identify how the channel was constructed. The checks below separate
construction provenance, fixed-model reliance, conditional information, and
finite-learner utility.

## 1. State the admissible input contract

Define the information available at deployment. For inverse folding, let `B`
be the supplied backbone, `Y` the unknown target sequence, and `C` a computed
channel. State whether a valid implementation requires `C=f(B)` or may use
native residue identities for the scientific task at hand.

The same channel can be legitimate for interaction analysis, where native
chemistry is observed, and inadmissible for designing an unknown sequence. The
contract—not the feature name—determines the provenance question.

## 2. Substitute labels while clamping every other input

Change `Y`, hold the coordinates, correspondence, preprocessing configuration,
and random state fixed, then recompute `C`.

```bash
python3 provenance_probe.py intervene \
  --original original.jsonl \
  --substituted substituted.jsonl
```

Report both channel and geometry movement. Any channel change supports target
dependence under this software-input intervention if geometry and the other
controls remain fixed. Geometry movement invalidates that interpretation.

In the audited implementations, fixed-coordinate substitutions changed 89.12%
of SurfDesign's inspected chemical-channel values and 90.12% of SurfPro's,
while checked geometry changed 0.00%. These are implementation-specific
measurements, not a prevalence estimate.

An unchanged result is only `NO_DEPENDENCE_DETECTED`. Non-injective codebooks
can preserve a chosen substitution, and finite tests cannot certify `C=f(B)`.

## 3. Measure simple channel-only predictability

Fit a train-only lookup from channel class to target label and evaluate it on
held-out data. Keep the scorer and denominator explicit.

The audited SurfDesign hydropathy partition has a channel-only analytic bound
of 85.89% on its stated population. For SurfPro, the train-fitted lookup reaches
69.72% under the coverage scorer and 54.67% under the official scorer. Those
numbers are not interchangeable and neither identifies provenance by itself.

## 4. Test whether fixed weights follow a specified channel change

Removing or corrupting a channel can cause generic out-of-distribution failure.
Prefer a targeted counterfactual: change the channel according to a prespecified
code transformation and test whether predictions move toward the corresponding
alternative target rather than merely losing native accuracy.

For the released SurfPro checkpoint and two tested code-table permutations, the
mean intended-minus-other target preference is +35.84 percentage points
[34.89, 36.74]. This supports code-specific reliance under those interventions;
it does not rule out every distribution-shift mechanism or establish a result
for all checkpoints.

## 5. Separate conditional information from finite-learner utility

A design-valid channel `C=f(B)` satisfies `I(Y;C|B)=0` but may still help a
bounded learner by presenting the backbone in a more usable representation. A
target-derived channel can be redundant given `B`. Therefore:

- do not infer conditional information from standalone decoding accuracy;
- do not infer provenance from matched marginal accuracy;
- report model class, training data, budget, seeds, scorer, and tuning rule for
  every utility contrast; and
- use clean joint training with competitive controls when claiming a useful
  design-valid representation.

The exact example in `results/identification_example/` demonstrates these
distinctions without sampling or fitting.

## 6. Include positive, negative, and invalid controls

At minimum, exercise the diagnostic on:

- a known target-dependent construction;
- a known backbone-derived construction that can remain correlated with the
  target;
- rare movement;
- a class-preserving substitution;
- missing or partial geometry controls; and
- changed geometry, mismatched IDs, and empty inputs.

The offline regression suite implements these cases using synthetic data.

## Reporting minimum

Report the input contract, construction code/version, intervention, clamped
fields, geometry coverage, movement rates, population, scorer, uncertainty,
fixed-model response, finite-learner configuration, and prohibited inferences.
`REPORT_TEMPLATE.md` supplies a compact structure.

## What this protocol does not establish

It does not show that surfaces are useless, certify that an unchanged channel
is design-valid, identify physical mutation effects, estimate prevalence across
methods, or replace competitive clean training. It separates the questions so
that each claim is supported by the appropriate evidence.
