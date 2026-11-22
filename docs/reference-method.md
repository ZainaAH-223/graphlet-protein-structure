# Reference method: GRAFENE

Faisal FE, Newaz K, Chaney JL, Li J, Emrich SJ, Clark PL, Milenković T (2017).
GRAFENE: Graphlet-based alignment-free network approach integrates 3D structural
and sequence (residue order) data to improve protein structural comparison.
Sci Rep 7:14890. doi:10.1038/s41598-017-14411-y (+ author correction 2020,
doi:10.1038/s41598-020-69025-8)
Software: https://academicweb.nd.edu/\~cone/PSN/

This project follows GRAFENE at every core step. Gaps (marked TODO) are to be
filled from the paper; anything not from GRAFENE is marked with its own source.

| Pipeline step | GRAFENE's choice | This project | Status |
|---|---|---|---|
| Network nodes | Residues | Same | ✅ |
| Contact definition | 4 options: any heavy atom ≤4, 5 or 6 Å; Cα ≤7.5 Å (GR-Align default) | Cα ≤7.5 Å (main); one heavy-atom definition for cross-check vs their software | TODO: pick heavy-atom cutoff |
| Network filters | TODO (software suggests ≥101 nodes, diameter ≥6, 1 component) | Same as GRAFENE | TODO: confirm in paper |
| Dataset | 9,509 CATH + 11,451 SCOP domains, 35 sets | Small CATH subset | TODO: CATH version, categories, group size (software: ≥30) |
| Graphlets | Whole-graphlet counts: 8 (3–4 nodes) or 29 (3–5 nodes); not orbits | Same | TODO: counting tool |
| Order-aware | Ordered graphlets (one version per distinct node ordering) | Same | TODO: how they count them |
| Long-range | K: count a graphlet only if all edges join residues ≥K apart (K=1 = standard) | Same | TODO: which K; what the 2020 correction changed |
| Normalisation | Count / total graphlet count of the network | Same | ✅ |
| Dimensionality reduction | PCA; fewest PCs (≥2) explaining ≥90% variance | Same | ✅ |
| Distance | 1 − cosine similarity | Same | ✅ |
| Evaluation | Rank all pairs by distance; AUPR + AUROC for same- vs different-category pairs; unsupervised | Same | ✅ |
| Simple-stats baseline | Other network measures (software: average degree, assortativity, …) | Same | TODO: exact list |
| Aligner comparison | GR-Align, TM-align, DaliLite: AUPR/AUROC + runtime | TM-align on a subset | ✅ |
| Biological check | α domains enriched in dense graphlets (G2, G8, G29); β in path-like (G1, G3, G9) | Same, plus DSSP | ✅ |

## Not from GRAFENE (own additions, each with a source)
- Sensitivity to cutoff and protein length (course brief)
- Disorder extension on AlphaFold models: pLDDT bands (Jumper et al. 2021;
  Tunyasuvunakool et al. 2021), DisProt, conditionally folded IDRs (Alderson et al. 2023)
