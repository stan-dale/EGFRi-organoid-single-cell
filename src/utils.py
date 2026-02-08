"""
Backwards-compatibility shim.

All contents have been moved to dedicated modules:
  - markers.py   → cell_type_markers, EEC_markers, cell_cycle_markers
  - palette.py   → celltype_palette

New code should import from those modules directly.
"""

from .markers import cell_type_markers, EEC_markers, cell_cycle_markers
from .palette import celltype_palette
