PHASE 13E — PAIRED MASTER TABLE

This step joins four pre-locked information layers on exactly the same unordered
25-km cell pairs:

1. Historical boundary separation:
   Primary = B01 Isthmus of La Paz + B03 Vizcaino/mid-peninsular.
   Secondary = B02 Loreto + B04 northern climatic transition.
2. Contemporary environmental distance from 13C.
3. Great-circle geographic distance between cell centroids.
4. Paired C3 and N0 community dissimilarities from 13D.

Primary historical coding:
- strict_cross requires one cell south and one north of the frozen boundary.
- For broad primary zones, a cell inside the zone does NOT create a strict
  crossing. A separate touches_or_crosses variable retains those cases.

No statistical inference is performed in 13E.
No boundary location, environmental predictor, or pair is selected based on
C3/N0 outcome magnitude.

The paired-valid table is the primary input for 13F and contains only pairs for
which both C3 and N0 Jaccard/Simpson values are defined.
