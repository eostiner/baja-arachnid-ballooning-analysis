# Step 10I centerpiece map quality control

The centerpiece map is produced by:

```text
scripts/mapping/step_10I/10I_build_final_map_and_ratios_v7.py
```

It combines two independently standardized summaries:

- ecoregion donuts: equal-cell C3 and N0 richness using the common Step 10C ecoregion cell count;
- latitude-band inset: median C3:N0 ratios and 2.5th–97.5th percentile equal-cell resampling intervals from Step 10G.

The script reads the actual `equal_cells` and `iterations` values from the output tables. It will therefore label a fresh rerun correctly if the common cell limit changes.

## Required visual checks after a rerun

1. The title states C3 = D1–D3 versus fixed N0 and D4 excluded.
2. Hatched ecoregions have fewer than the formal eight-cell threshold and receive no donut.
3. All four a priori test latitudes (24°, 26°, 28°, and 30°N) are shown as dashed lines.
4. The latitude-band panel contains all five bands in north-to-south display order.
5. The legend says “a priori test latitude,” not “published break.”
6. The lower note states that no adjacent-band ratio shift survived BH correction when the rerun table supports that statement.
7. PDF and SVG exports retain editable vector text and geometry.

Compare the result with `docs/reference_figures/Figure_3_reference_2026-07-16.*`. Values may change after a deliberate input revision; the reference is principally for layout and annotation regression. Font rendering and antialiasing may differ across operating systems.
