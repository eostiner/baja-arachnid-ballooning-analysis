# Phase 15G v0.5 — editable SVG text

The statistical content is unchanged from the frozen v0.4 figure.

## SVG export convention

- `svg.fonttype = none` preserves labels as SVG text rather than glyph outlines.
- Every Matplotlib `Text` artist is grouped under its original `text_*` group.
- Multiline strings emitted by Matplotlib as separate sibling text elements are collapsed
  into one `<text>` element containing one positioned `<tspan>` per line.
- Each `text_*` group receives an `inkscape:label` containing the human-readable phrase.
- The exporter fails if any `text_*` group contains outlined glyph paths or more than one
  SVG `<text>` object after normalization.

This makes phrases and multiline labels substantially easier to select and edit in
Inkscape, Adobe Illustrator, and Affinity Designer. The SVG references DejaVu Sans and
therefore requires that font, or a deliberate font substitution, on the editing system.
