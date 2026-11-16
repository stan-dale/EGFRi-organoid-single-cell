"""Paths, dataset registry, colour palettes and marker resources."""

import numpy as np
import pytest

from src import config
from src.markers import (EEC_markers, cell_type_markers, filter_marker_dict,
                         load_cellassign_matrix)
from src.palette import celltype_palette, eec_palette, get_color, normalize_celltype_name


def test_output_dirs_created_on_import():
    for d in config.OUTPUT_DIRS:
        assert d.is_dir(), d


def test_registry_paths_live_under_data_dir():
    for key, entry in config.DATASETS.items():
        for stage in ("raw", "clustered", "labelled"):
            assert config.DATA_DIR in entry[stage].parents, (key, stage)
        assert entry["condition"] in {"Dz", "Lapa", "AS"}


@pytest.mark.parametrize("name", [p.stem for p in config.CELLASSIGN_MARKERS_DIR.glob("*.csv")])
def test_cellassign_matrices_are_binary(name):
    m = load_cellassign_matrix(name)
    assert m.shape[0] > 0 and m.shape[1] > 0
    assert set(np.unique(m.to_numpy())) <= {0, 1}
    assert not m.index.duplicated().any(), "duplicate genes in marker matrix"


def test_aliases_map_to_canonical_colours():
    assert normalize_celltype_name("Inflammed secretory PCs") == "Inflamed Secretory PCs"
    assert normalize_celltype_name("Secretory progenitor cells") == "Secretory PCs"
    assert get_color("Proliferating progenitor cells") == get_color("Proliferating PCs")
    assert get_color("not a cell type") == "#808080"


def test_every_marker_cell_type_has_a_colour():
    for ct in cell_type_markers:
        assert ct in celltype_palette, ct
    for ct in EEC_markers:
        if ct in {"Pre EECs", "T4 cells"}:  # marker-only groups, never used as labels
            continue
        assert ct in eec_palette, ct


def test_filter_marker_dict_drops_absent_genes(adata):
    filtered, missing = filter_marker_dict(adata, {"ISCs": ["LGR5", "NOT_A_GENE"], "X": "ALSO_ABSENT"})
    assert filtered == {"ISCs": ["LGR5"]}
    assert missing == {"ISCs": ["NOT_A_GENE"], "X": ["ALSO_ABSENT"]}
