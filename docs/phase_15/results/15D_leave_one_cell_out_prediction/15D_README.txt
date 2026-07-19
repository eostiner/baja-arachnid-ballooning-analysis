PHASE 15D — LEAVE-ONE-CELL-OUT PREDICTIVE VALIDATION
========================================================
Classification: PREDICTIVE_IMPROVEMENT_POSITIVE_BUT_UNCERTAIN
Held-out cells: 12
Cells where stress model had lower CRPS: 9
Mean CRPS improvement (null minus stress): 0.035160
Exact cell sign-flip p: 0.292969 (4096 patterns)

This is a prediction check, not another test fitted to the same observations. Each
cell is withheld in turn, and the stress model is compared with a null model containing
the same sampling and time controls. Positive CRPS improvement means the stress model
predicted the held-out cell more accurately.
