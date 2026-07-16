from __future__ import annotations

import importlib.util
import ast
import csv
import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, relative_path: str):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class PublicationCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.step11g = load_module(
            "step11g",
            "scripts/core/step_11G/11G_equal_cell_trait_partitioned_richness_turnover.py",
        )
        cls.normalize_traits = load_module(
            "normalize_traits",
            "tools/normalize_trait_table.py",
        )

    def test_explicit_trait_tiers(self):
        expected = {
            "D1": "D1",
            "D2": "D2",
            "D3": "D3",
            "D4": "D4",
            "N0": "N0",
            "C3": "C3",
            "non-ballooning": "N0",
            "D4 excluded": "D4",
        }
        for value, result in expected.items():
            with self.subTest(value=value):
                self.assertEqual(self.step11g.parse_evidence_class(value), result)

    def test_simpson_replacement_formula(self):
        # A={a,b,c}; B={b,c,d,e}: shared=2, unique=1 and 2, beta_sim=1/(2+1).
        a = 0b00111
        b = 0b11011
        metrics = self.step11g.beta_metrics(a, b)
        self.assertAlmostEqual(metrics["simpson_turnover"], 1 / 3)

    def test_rarefaction_axis_uses_mean_expected(self):
        text = (ROOT / "scripts/core/step_11F/11F_equal_cell_rarefaction.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("Mean expected genus richness", text)
        self.assertNotIn('ax.set_ylabel("Median genus richness"', text)

    def test_authoritative_two_column_trait_schema(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            trait_dir = root / "ANALYSIS_READY_INPUTS/03_trait_tables"
            trait_dir.mkdir(parents=True)
            path = trait_dir / "07_reviewed_genus_trait_lookup_final.csv"
            rows = [
                {"genus": "Directus", "exclusive_tier": "D1", "primary_C3_group": "Ballooning"},
                {"genus": "Juvenilis", "exclusive_tier": "D3", "primary_C3_group": "C3"},
                {"genus": "Broadus", "exclusive_tier": "D4", "primary_C3_group": "Excluded (D4)"},
                {"genus": "Groundus", "exclusive_tier": "", "primary_C3_group": "Non-ballooning"},
            ]
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            result = self.normalize_traits.normalize_project_traits(root, write=True)
            self.assertEqual(result.counts["C3"], 2)
            self.assertEqual(result.counts["N0"], 1)
            self.assertEqual(result.counts["D4_excluded"], 1)
            with result.output.open(encoding="utf-8", newline="") as stream:
                normalized_rows = list(csv.DictReader(stream))
            groundus = next(row for row in normalized_rows if row["genus"] == "Groundus")
            self.assertEqual(groundus["evidence_class"], "N0")
            self.assertEqual(groundus["analysis_class"], "N0")


    def test_step11b_band_metrics_calls_are_complete(self):
        path = ROOT / "scripts/core/step_11B/11B_equal_cell_resampling.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        required = {
            "genus_mask",
            "selected_indices",
            "cell_masks",
            "balloon_mask",
            "n0_mask",
            "d4_mask",
            "classified_mask",
        }
        calls = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "band_metrics"
        ]
        self.assertGreaterEqual(len(calls), 2)
        for call in calls:
            supplied = {kw.arg for kw in call.keywords if kw.arg is not None}
            self.assertTrue(
                required.issubset(supplied),
                f"Incomplete band_metrics() call on line {call.lineno}: "
                f"missing {sorted(required - supplied)}",
            )

    def test_step11f_allows_missing_confidence_field(self):
        step11f = load_module(
            "step11f",
            "scripts/core/step_11F/11F_equal_cell_rarefaction.py",
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "traits.csv"
            rows = [
                {"genus": "Alpha", "evidence_class": "D1"},
                {"genus": "Beta", "evidence_class": "N0"},
            ]
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            result = step11f.load_low_confidence_mask(path, ["Alpha", "Beta"])
            self.assertEqual(result["low_mask"], 0)
            self.assertEqual(result["keep_mask"].bit_count(), 2)
            self.assertEqual(result["confidence_counts"], {"UNSPECIFIED": 2})
            self.assertIsNone(result["confidence_field"])


    def test_authoritative_workbook_trait_schema(self):
        from openpyxl import Workbook
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            path = root / "USE_GOOD_BalloonID_Baja_Arachnid_GenusSpecies_Long_D1_D4_AUTHORITATIVE.xlsx"
            wb = Workbook()
            ws = wb.active
            ws.title = "Genus_Trait_Master_267"
            ws.append(["genus", "order_final", "family_final", "final_tier_for_current_build", "primary_C3_group"])
            ws.append(["Directus", "Araneae", "Testidae", "D1", "Ballooning (C3)"])
            ws.append(["Groundus", "Scorpiones", "Groundidae", None, "Non-ballooning (N0)"])
            ws.append(["Broadus", "Araneae", "Broadidae", "D4", "Excluded (D4)"])
            wb.save(path)
            result = self.normalize_traits.normalize_project_traits(root, write=True)
            self.assertEqual(result.source, path)
            self.assertEqual(result.counts["C3"], 1)
            self.assertEqual(result.counts["N0"], 1)
            self.assertEqual(result.counts["D4_excluded"], 1)

    def test_rerun_compatibility_patches_present(self):
        checks = {
            "scripts/biogeography/step_10B/10B_build_cell_ecoregion_crosswalk.R": [
                '"ecoregion_label", drop = FALSE',
                'covered <- if (nrow(z))',
            ],
            "scripts/environment/step_12B/12B_install_and_audit_GEE_exports.R": [
                "auditing existing installed rasters instead",
                "same_file <- identical",
            ],
            "scripts/environment/step_12C/12C_build_cell_environment_model_table.R": [
                "primary_candidate_environment_complete",
                "primary_candidate_denominators_positive",
            ],
            "scripts/environment/step_12D/12D_screen_predictors_compare_models.R": [
                "candidate_sets_consistent_with_step12C",
            ],
            "scripts/environment/step_12E/12E_test_coefficient_stability_and_reduced_models.R": [
                "candidate_sets_consistent_with_step12C",
            ],
        }
        for relative_path, required_text in checks.items():
            text = (ROOT / relative_path).read_text(encoding="utf-8")
            for fragment in required_text:
                self.assertIn(fragment, text, f"Missing {fragment!r} in {relative_path}")
        combined = "\n".join((ROOT / path).read_text(encoding="utf-8") for path in [
            "scripts/environment/step_12D/12D_screen_predictors_compare_models.R",
            "scripts/environment/step_12E/12E_test_coefficient_stability_and_reduced_models.R",
        ])
        self.assertNotIn("expected_counts <- c(primary = 195L", combined)

    def test_step10i_support_module_is_present(self):
        self.assertTrue(
            (ROOT / "scripts/mapping/step_10I/_base_v2.py").is_file()
        )


if __name__ == "__main__":
    unittest.main()
