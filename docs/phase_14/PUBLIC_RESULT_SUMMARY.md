# Phase 14 public result summary

"
        "Phase 14 tested whether ballooning-capable assemblages showed greater temporal replacement as recent environmental stress worsened. "
        "The analysis retained 15 adjacent-period comparisons across 12 repeatedly sampled 25-km cells. Eleven of 12 frozen primary or sensitivity estimates were positive, but none of their frequentist 95% intervals excluded zero. "
        "The result is therefore classified as **positive but uncertain across primary models** and remains exploratory, geographically limited, observational, and non-causal.

"
        "The selected public outputs preserve the eligibility audit, paired temporal turnover, original and extended stress models, sampling-continuity and trait-threshold sensitivities, and the final evidence synthesis. Raw Earth Engine tables are regenerated locally and are not committed.
",
        encoding="utf-8",
    )
    p15 = repo / "docs" / "phase_15" / "PUBLIC_RESULT_SUMMARY.md"
    p15.write_text(
        The primary Bayesian paired hierarchical measurement-error model estimated a standardized recent-stress effect of **0.200426**, with a 95% credible interval of **0.022920 to 0.358910** and `Pr(beta > 0 | data) = 0.989100`. The positive direction remained stable across the frozen prior and observation-uncertainty sensitivities. However, held-out prediction improved only modestly and was not decisive (`exact sign-flip p = 0.292969`). The frozen synthesis is therefore **positive but uncertain**: H3 is plausible but unconfirmed.

Step 15G combines the original frozen spatial C3-N0 turnover contrast with the exploratory temporal result. The spatial and temporal panels answer distinct questions and must not be interpreted as one demonstrated causal pathway.
