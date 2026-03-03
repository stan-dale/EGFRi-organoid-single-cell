#!/usr/bin/env python3
"""
Fix hardcoded paths in Jupyter notebooks to use src.config.

Handles:
1. Removing os.chdir() calls and prefixing subsequent reads with config DIRs
2. Replacing hardcoded /Users/stanleydale/... paths with config variables
3. Adding/updating config imports as needed
4. Removing stale sys.path.append('../src')
"""

import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"
ABS_PREFIX = "/Users/stanleydale/user_generated/breault-lab/single-cell"

# Skip these directory prefixes (relative to notebooks/)
SKIP_PREFIXES = ("dge/archive", "cellrank")


def get_rel_prefix(nb_path):
    """Return '../..' style prefix to reach PROJECT_ROOT from notebook dir."""
    rel = nb_path.parent.relative_to(NOTEBOOKS_DIR)
    depth = len(rel.parts) + 1          # +1 for notebooks/ itself
    return "/".join([".."] * depth)


def fix_notebook(nb_path):
    with open(nb_path) as f:
        nb = json.load(f)

    original = json.dumps(nb)
    rel_prefix = get_rel_prefix(nb_path)

    # ── state ──
    dir_state = None            # None | 'data' | 'analysis'
    config_imports_needed = set()

    # ── locate existing config import ──
    cfg_cell_idx = cfg_line_idx = None
    has_sys_path_setup = False
    for ci, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        for li, line in enumerate(cell["source"]):
            if "from src.config import" in line:
                cfg_cell_idx, cfg_line_idx = ci, li
            if "sys.path.insert" in line and "resolve" in line:
                has_sys_path_setup = True

    # ── process each code cell ──
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue

        new_lines = []
        for line in cell["source"]:

            # ---- detect os.chdir and update dir_state BEFORE removing ----
            chdir_m = re.search(r"os\.chdir\(['\"](.+?)['\"]\)", line)
            if chdir_m:
                p = chdir_m.group(1)
                if "data" in p and "analysis" not in p:
                    dir_state = "data"
                    config_imports_needed.add("DATA_DIR")
                elif "analysis" in p or "manual_labelled" in p:
                    dir_state = "analysis"
                    config_imports_needed.add("ANALYSIS_DIR")

                # Remove the os.chdir line entirely
                # But keep any other code on the same line (rare)
                stripped = line.strip()
                if stripped == chdir_m.group(0) or stripped == chdir_m.group(0) + "\n":
                    continue          # drop the whole line
                else:
                    # multi-statement line – just remove the chdir part
                    line = line.replace(chdir_m.group(0), "").lstrip("\n")
                    if not line.strip():
                        continue

            # ---- remove stale sys.path.append('../src') ----
            if re.match(r"\s*sys\.path\.append\(['\"]\.\.\/src['\"]\)\s*$", line):
                continue

            # ---- hardcoded absolute paths ----
            if ABS_PREFIX in line:
                line = _fix_abs(line, config_imports_needed)

            # ---- relative reads after os.chdir ----
            if dir_state == "data":
                line = _fix_rel_data(line, config_imports_needed)
            elif dir_state == "analysis":
                line = _fix_rel_analysis(line, config_imports_needed)

            new_lines.append(line)

        cell["source"] = new_lines

    # ── add / update config imports ──
    if config_imports_needed:
        _ensure_imports(nb, config_imports_needed, cfg_cell_idx,
                        cfg_line_idx, has_sys_path_setup, rel_prefix)

    new_json = json.dumps(nb)
    if new_json == original:
        return False

    with open(nb_path, "w") as f:
        json.dump(nb, f, indent=1)
        f.write("\n")
    return True


# ═══════════════════════════  helpers  ═══════════════════════════


def _fix_abs(line, needs):
    """Replace hard-coded /Users/stanleydale/… paths."""

    # sc.settings.figdir  (with or without Path(...))
    m = re.search(
        r'sc\.settings\.figdir\s*=\s*(?:Path\()?["\']'
        + re.escape(ABS_PREFIX)
        + r"/figures/?([^\"']*)[\"'](?:\))?",
        line,
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("FIGURES_DIR")
        repl = f'str(FIGURES_DIR / "{sub}")' if sub else "str(FIGURES_DIR)"
        return re.sub(
            r'sc\.settings\.figdir\s*=\s*(?:Path\()?["\']'
            + re.escape(ABS_PREFIX)
            + r"/figures/?[^\"']*[\"'](?:\))?",
            f"sc.settings.figdir = {repl}",
            line,
        )

    # intermediate_directory = '…/analysis/…'
    m = re.search(
        r"intermediate_directory\s*=\s*['\"]"
        + re.escape(ABS_PREFIX)
        + r"/analysis/?([^'\"]*)['\"]",
        line,
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("ANALYSIS_DIR")
        val = f'str(ANALYSIS_DIR / "{sub}")' if sub else "str(ANALYSIS_DIR)"
        return re.sub(
            r"intermediate_directory\s*=\s*['\"]"
            + re.escape(ABS_PREFIX)
            + r"/analysis/?[^'\"]*['\"]",
            f"intermediate_directory = {val}",
            line,
        )

    # pd.read_csv('/…/utilities/…')
    m = re.search(
        r"pd\.read_csv\(['\"]"
        + re.escape(ABS_PREFIX)
        + r"/utilities/([^'\"]+)['\"]",
        line,
    )
    if m:
        fname = m.group(1)
        needs.add("UTILITIES_DIR")
        return re.sub(
            r"pd\.read_csv\(['\"]"
            + re.escape(ABS_PREFIX)
            + r"/utilities/[^'\"]+['\"]",
            f'pd.read_csv(UTILITIES_DIR / "{fname}"',
            line,
        )

    # Catch-all: any remaining '/Users/stanleydale/.../figures/…'
    m = re.search(
        r"['\"]" + re.escape(ABS_PREFIX) + r"/figures/([^'\"]*)['\"]", line
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("FIGURES_DIR")
        return re.sub(
            r"['\"]" + re.escape(ABS_PREFIX) + r"/figures/[^'\"]*['\"]",
            f'str(FIGURES_DIR / "{sub}")',
            line,
        )

    # Catch-all: '/Users/stanleydale/.../analysis/…'
    m = re.search(
        r"['\"]" + re.escape(ABS_PREFIX) + r"/analysis/?([^'\"]*)['\"]", line
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("ANALYSIS_DIR")
        val = f'str(ANALYSIS_DIR / "{sub}")' if sub else "str(ANALYSIS_DIR)"
        return re.sub(
            r"['\"]" + re.escape(ABS_PREFIX) + r"/analysis/?[^'\"]*['\"]",
            val,
            line,
        )

    # Catch-all: '/Users/stanleydale/.../data/…'
    m = re.search(
        r"['\"]" + re.escape(ABS_PREFIX) + r"/data/?([^'\"]*)['\"]", line
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("DATA_DIR")
        val = f'str(DATA_DIR / "{sub}")' if sub else "str(DATA_DIR)"
        return re.sub(
            r"['\"]" + re.escape(ABS_PREFIX) + r"/data/?[^'\"]*['\"]",
            val,
            line,
        )

    # Catch-all: '/Users/stanleydale/.../utilities/…'
    m = re.search(
        r"['\"]" + re.escape(ABS_PREFIX) + r"/utilities/([^'\"]*)['\"]", line
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("UTILITIES_DIR")
        return re.sub(
            r"['\"]" + re.escape(ABS_PREFIX) + r"/utilities/[^'\"]*['\"]",
            f'str(UTILITIES_DIR / "{sub}")',
            line,
        )

    # Catch-all: any other '/Users/stanleydale/…/single-cell/…'
    m = re.search(
        r"['\"]" + re.escape(ABS_PREFIX) + r"/?([^'\"]*)['\"]", line
    )
    if m:
        sub = m.group(1).strip("/")
        needs.add("PROJECT_ROOT")
        val = f'str(PROJECT_ROOT / "{sub}")' if sub else "str(PROJECT_ROOT)"
        return re.sub(
            r"['\"]" + re.escape(ABS_PREFIX) + r"/?[^'\"]*['\"]",
            val,
            line,
        )

    return line


def _fix_rel_data(line, needs):
    """Prefix relative reads that follow an os.chdir to data/."""
    m = re.search(r"sc\.read_h5ad\(['\"]([^'\"/][^'\"]*\.h5ad)['\"]", line)
    if m:
        fname = m.group(1)
        needs.add("DATA_DIR")
        return re.sub(
            r"sc\.read_h5ad\(['\"][^'\"/][^'\"]*\.h5ad['\"]",
            f'sc.read_h5ad(DATA_DIR / "{fname}"',
            line,
        )
    return line


def _fix_rel_analysis(line, needs):
    """Prefix relative reads/writes that follow an os.chdir to analysis/."""
    # sc.read_h5ad
    m = re.search(r"sc\.read_h5ad\(['\"]([^'\"/][^'\"]*\.h5ad)['\"]", line)
    if m:
        fname = m.group(1)
        needs.add("ANALYSIS_DIR")
        return re.sub(
            r"sc\.read_h5ad\(['\"][^'\"/][^'\"]*\.h5ad['\"]",
            f'sc.read_h5ad(ANALYSIS_DIR / "{fname}"',
            line,
        )

    # sc.write('rel/path', adata)
    m = re.search(r"sc\.write\(['\"]([^'\"/][^'\"]*)['\"]", line)
    if m:
        fpath = m.group(1)
        needs.add("ANALYSIS_DIR")
        return re.sub(
            r"sc\.write\(['\"][^'\"/][^'\"]*['\"]",
            f'sc.write(str(ANALYSIS_DIR / "{fpath}")',
            line,
        )

    # .to_csv('rel/path')
    m = re.search(r"\.to_csv\(['\"]([^'\"/][^'\"]*)['\"]", line)
    if m:
        fpath = m.group(1)
        needs.add("ANALYSIS_DIR")
        return re.sub(
            r"\.to_csv\(['\"][^'\"/][^'\"]*['\"]",
            f'.to_csv(str(ANALYSIS_DIR / "{fpath}")',
            line,
        )

    # pd.read_csv('rel/path')
    m = re.search(r"pd\.read_csv\(['\"]([^'\"/][^'\"]*)['\"]", line)
    if m:
        fpath = m.group(1)
        needs.add("ANALYSIS_DIR")
        return re.sub(
            r"pd\.read_csv\(['\"][^'\"/][^'\"]*['\"]",
            f'pd.read_csv(str(ANALYSIS_DIR / "{fpath}")',
            line,
        )

    return line


def _ensure_imports(nb, needs, cfg_ci, cfg_li, has_sys_path, rel_prefix):
    """Add or update the src.config import block."""
    imports_str = ", ".join(sorted(needs))

    if cfg_ci is not None:
        # Update existing import line
        cell = nb["cells"][cfg_ci]
        old = cell["source"][cfg_li]
        m = re.search(r"from src\.config import (.+)", old)
        if m:
            existing = {s.strip().rstrip("\n") for s in m.group(1).split(",")}
            merged = existing | needs
            if merged != existing:
                cell["source"][cfg_li] = (
                    "from src.config import "
                    + ", ".join(sorted(merged))
                    + "\n"
                )
        return

    # No existing import – find the first code cell with imports
    for ci, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code" or not cell["source"]:
            continue
        src = "".join(cell["source"])
        if "import " not in src:
            continue

        additions = []
        if "import sys" not in src and "import os, sys" not in src:
            additions.append("import sys\n")
        if "from pathlib import Path" not in src:
            additions.append("from pathlib import Path\n")
        if not has_sys_path:
            additions.append(
                f'sys.path.insert(0, str(Path("{rel_prefix}").resolve()))\n'
            )
        additions.append(f"from src.config import {imports_str}\n")

        cell["source"].extend(["\n"] + additions)
        return

    # Fallback: insert a new cell at the top
    new_cell = {
        "cell_type": "code",
        "metadata": {},
        "source": [
            "import sys\n",
            "from pathlib import Path\n",
            f'sys.path.insert(0, str(Path("{rel_prefix}").resolve()))\n',
            f"from src.config import {imports_str}\n",
        ],
        "outputs": [],
        "execution_count": None,
    }
    # Insert after the first markdown cell (if any) or at position 0
    insert_at = 0
    for ci, cell in enumerate(nb["cells"]):
        if cell["cell_type"] == "markdown":
            insert_at = ci + 1
            break
    nb["cells"].insert(insert_at, new_cell)


# ═══════════════════════════  main  ═══════════════════════════


def main():
    notebooks = sorted(NOTEBOOKS_DIR.rglob("*.ipynb"))
    fixed = []
    skipped = []

    for nb_path in notebooks:
        if ".ipynb_checkpoints" in str(nb_path):
            continue
        if nb_path.is_dir():
            continue
        rel = str(nb_path.relative_to(NOTEBOOKS_DIR))
        if any(rel.startswith(s) for s in SKIP_PREFIXES):
            skipped.append(rel)
            continue

        if fix_notebook(nb_path):
            fixed.append(rel)
            print(f"  FIXED  {rel}")
        else:
            print(f"  ok     {rel}")

    print(f"\n{'='*50}")
    print(f"Fixed:   {len(fixed)} notebooks")
    print(f"Skipped: {len(skipped)} (archive)")
    for f in fixed:
        print(f"  - {f}")


if __name__ == "__main__":
    main()
