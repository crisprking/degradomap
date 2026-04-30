"""Minimal smoke tests."""

def test_imports():
    from degradomap import e3_list, structures, depmap, analysis, cli
    assert hasattr(e3_list, "build_e3_list")
    assert hasattr(structures, "compute_features")
    assert hasattr(depmap, "essentiality_summary")
    assert hasattr(analysis, "run_experiment")

def test_protac_genes_list():
    from degradomap.e3_list import PROTAC_E3_GENES
    assert "CRBN" in PROTAC_E3_GENES
    assert "VHL" in PROTAC_E3_GENES
    assert len(PROTAC_E3_GENES) >= 12

def test_feature_lists():
    from degradomap.analysis import STRUCTURAL_FEATURES, BIOLOGICAL_FEATURES, ALL_FEATURES
    assert len(STRUCTURAL_FEATURES) == 10
    assert len(BIOLOGICAL_FEATURES) == 4
    assert len(ALL_FEATURES) == 14
