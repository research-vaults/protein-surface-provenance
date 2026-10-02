# Auditing a computed input: seven operations, four questions

The task's inference-time contract comes first. State which information is available when the target is unknown. The following checks support different conclusions; the list is neither a minimal identification theorem nor a universally necessary protocol.

| Operation | Question answered | Limit |
|---|---|---|
| 1. Substitute target labels; hold all other preprocessing inputs and randomness fixed | Does construction depend on those labels? | Movement supports dependence; invariance under one substitution does not certify independence. Check complete geometry and correspondence. |
| 2. Fit a simple channel-only decoder on train; test it on held-out chains | How well can this decoder read the target? | Lookup accuracy is not generally the Bayes ceiling or conditional information. Analytic bounds require an explicit known code and target distribution. |
| 3. Compare channels with similar standalone accuracy | Do channels of similar marginal decodability have similar utility for this learner? | Does not match cardinality, complementary errors, geometry or learnability. No isolated provenance effect follows. |
| 4. Compare native, constant, shuffled, predicted and counterfactual inputs | What does this fixed model rely on? | Neutralisation alone can cause distribution shift; wrong-target controls make a targeted-redirection claim more specific. |
| 5. Use frozen published models on a common population and scorer | What incremental gain does the chosen fusion rule obtain? | The trained fusion rule is still a finite learner, and zero gain does not establish zero conditional information. |
| 6. Jointly train admissible backbone-only and surface arms, with capacity and alignment controls | Does this design-valid representation help at the tested budget? | A weak-model null cannot rule out stronger meshes or training. Training is not a minutes-long, GPU-free check. |
| 7. Remove, restore and repair an internal pathway | Where is reliance expressed, and can it be changed? | Mechanistic evidence is conditional on the tested models and interventions. |

## Three concrete comparisons from the paper

- Source substitution: 89.12% chemistry movement in SurfDesign and 90.12% in SurfPro, with 0% movement in the respective checked geometry. The standalone utility and fixed-model reliance must still be measured separately.
- SurfPro, same 887 chains: train-fitted lookup versus released model is 54.67% versus 57.79% on the official scorer (model margin +3.12 pp [2.46, 3.78]); on a shared coverage scorer it is 69.72% versus 70.33% (margin +0.61 pp [−0.05, 1.28]). Do not mix scorers or populations.
- Three-seed matched retraining: changing meanings while keeping the original 18 classes changes recovery by at most 0.342 pp; one random 18-class partition improves by 3.347 pp [2.393, 4.301]. Same cardinality is not the same partition or information. These results concern sequence recovery, not biochemical function.

A useful report records the contract, source intervention, simple rival, same-scorer margin, dependence of the trained predictor, and scope of clean joint training. Conditional mutual information requires its own argument or estimator; none of these lookup or payoff numbers is that quantity. Surfaces may still help a model when they are computed entirely from admissible inputs.
