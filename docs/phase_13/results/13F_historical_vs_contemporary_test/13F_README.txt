PHASE 13F — HISTORICAL VS CONTEMPORARY TEST

PRIMARY QUESTION
Do C3 ballooning assemblages show stronger relative contemporary-environment
tracking and weaker relative historical-boundary structure than N0
non-ballooning assemblages?

PRIMARY RESPONSE
C3 minus N0 Simpson replacement dissimilarity.

SECONDARY RESPONSE
C3 minus N0 Jaccard total dissimilarity.

PRIMARY CELL SET
Only cells with at least one C3 genus and at least one N0 genus are retained.
This creates complete matched dissimilarity matrices and avoids interpreting
empty assemblages as Simpson replacement.

MODEL
delta ~ z(log1p geographic distance) + z(contemporary environmental distance)
        + B01 La Paz strict crossing + B03 Vizcaino strict crossing

INFERENCE
Freedman-Lane-style residual permutation with simultaneous row/column
(cell-label) permutation. Pairwise rows are NOT treated as iid observations.
Permutation p-values are the inferential p-values.

PERMUTATIONS
4999
Seed: 20260718

DIRECTIONAL EXPECTATIONS
- Environment coefficient > 0: C3 is relatively more associated with modern
  environmental dissimilarity than N0.
- Boundary coefficient < 0: C3 is relatively less discontinuous across that
  historical boundary than N0.

All reported coefficient permutation tests are two-sided. Direction is
interpreted from the sign of the observed coefficient.

IMPORTANT
The primary historical hypothesis uses only B01 (La Paz) and B03
(Vizcaino/mid-peninsular), frozen before Phase 13 outcome testing.
B02 Loreto and B04 northern transition are not included in the primary model.
