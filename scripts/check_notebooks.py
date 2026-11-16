#!/usr/bin/env python
"""Static checks that the notebooks can run end to end on a fresh clone.

Needs no data and no scientific stack (stdlib + PyYAML only). It reads each
notebook's code and checks, without executing anything:

  json        the notebook file parses
  manifest    every notebook is listed in pipeline.yaml (run or not_in_pipeline)
  import      every ``from src.x import y`` names a module and object that exist,
              and nothing imports a src module by its bare name (``from vis import``)
  abspath     no hard-coded machine-specific paths (/Users/..., /Volumes/...)
  relpath     no reads/writes relative to the working directory, no os.chdir()
  undefined   no names used without being defined (needs pyflakes; skipped if absent)
  resource    every file referenced under RESOURCES_DIR is tracked in resources/
  flow        every h5ad a pipeline step reads is raw input or written by an
              earlier step, following the order in pipeline.yaml
  registry    every path in src.config.DATASETS is written by its dataset's steps

Usage:
    python scripts/check_notebooks.py          # report, exit 1 on any issue
"""

from __future__ import annotations

import ast
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
NB_DIR = ROOT / "notebooks"
SRC_DIR = ROOT / "src"
sys.path.insert(0, str(ROOT))

from src import config  # noqa: E402  (pure paths, no heavy imports)

ABSPATH_RE = re.compile(r"""['"](/(?:Users|Volumes|home)/[^'"]*|~/[^'"]*)['"]""")
READ_FUNCS = {"read_h5ad", "read", "read_csv"}
WRITE_FUNCS = {"write_h5ad", "write", "to_csv"}
DATA_SUFFIXES = (".h5ad", ".csv")


@dataclass(frozen=True)
class Issue:
    notebook: str
    kind: str
    message: str

    def __str__(self) -> str:
        return f"[{self.kind}] {self.notebook}: {self.message}"


# ── Manifest ──────────────────────────────────────────────────────────────────

def load_manifest(path: Path = ROOT / "pipeline.yaml") -> dict:
    return yaml.safe_load(path.read_text())


def pipeline_order(manifest: dict) -> list[tuple[str, str | None]]:
    """(notebook path relative to notebooks/, dataset key or None) in run order."""
    order = []
    for key, ds in manifest["datasets"].items():
        for step in ds["steps"]:
            order.append((f"{ds['folder']}/{step['nb']}", key))
    for step in manifest.get("cross_dataset", []):
        order.append((step["nb"], None))
    return order


def is_excluded(rel: str, manifest: dict) -> bool:
    return any(rel == e or (e.endswith("/") and rel.startswith(e))
               for e in manifest.get("not_in_pipeline", []))


# ── Notebook parsing ──────────────────────────────────────────────────────────

def code_cells(nb_path: Path) -> list[str]:
    nb = json.loads(nb_path.read_text())
    cells = []
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        src = "".join(c["source"])
        # Drop IPython magics / shell escapes so the cell parses as Python.
        src = "\n".join("" if l.lstrip().startswith(("%", "!")) else l
                        for l in src.splitlines())
        cells.append(src)
    return cells


def parse(src: str) -> ast.Module | None:
    try:
        return ast.parse(src)
    except SyntaxError:
        return None


# ── src module introspection ──────────────────────────────────────────────────

def module_names(path: Path) -> set[str]:
    """Top-level names a module defines (functions, classes, assignments, imports)."""
    tree = ast.parse(path.read_text())
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            for t in node.targets:
                names |= {n.id for n in ast.walk(t) if isinstance(n, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            names |= {(a.asname or a.name).split(".")[0] for a in node.names}
    return names


SRC_MODULES = {p.stem: p for p in SRC_DIR.glob("*.py") if p.stem != "__init__"}
# Names that once lived in src/ under a bare module name, plus current stems.
LOCAL_MODULE_NAMES = set(SRC_MODULES) | {"labels", "utils", "vis", "config"}


def check_imports(rel: str, trees: list[ast.Module]) -> list[Issue]:
    issues = []
    for tree in trees:
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                top = node.module.split(".")[0]
                if top == "src":
                    sub = node.module.split(".", 1)[1] if "." in node.module else None
                    if sub is None:  # from src import figure_style
                        for a in node.names:
                            if a.name not in SRC_MODULES:
                                issues.append(Issue(rel, "import", f"src/{a.name}.py does not exist"))
                        continue
                    if sub not in SRC_MODULES:
                        issues.append(Issue(rel, "import", f"src/{sub}.py does not exist "
                                                           f"(needed for: {', '.join(a.name for a in node.names)})"))
                        continue
                    defined = module_names(SRC_MODULES[sub])
                    for a in node.names:
                        if a.name != "*" and a.name not in defined:
                            issues.append(Issue(rel, "import", f"src.{sub} has no '{a.name}'"))
                elif top in LOCAL_MODULE_NAMES and node.level == 0:
                    issues.append(Issue(rel, "import", f"bare 'from {node.module} import ...' "
                                                       f"only works with a CWD hack; use 'from src.{node.module} import'"))
            elif isinstance(node, ast.Import):
                for a in node.names:
                    if a.name in LOCAL_MODULE_NAMES:
                        issues.append(Issue(rel, "import", f"bare 'import {a.name}'; use 'from src import {a.name}'"))
    return issues


# ── Path evaluation ───────────────────────────────────────────────────────────

def base_namespace() -> dict:
    ns = {k: getattr(config, k) for k in dir(config) if k.isupper()}
    ns.update(Path=Path, str=str)
    return ns


def try_eval(expr: ast.expr, ns: dict):
    try:
        return eval(compile(ast.Expression(expr), "<nb>", "eval"), {"__builtins__": {}}, ns)
    except Exception:
        return None


def call_name(node: ast.Call) -> str | None:
    f = node.func
    if isinstance(f, ast.Attribute):
        return f.attr
    if isinstance(f, ast.Name):
        return f.id
    return None


def file_ops(trees: list[ast.Module]) -> list[tuple[str, object]]:
    """Ordered (op, evaluated path) for every read/write whose path can be resolved.

    op is "read", "write", or "unresolved-<op>" for paths that are string
    literals relative to the working directory.
    """
    ns = base_namespace()
    ops = []
    for tree in trees:
        # Visit statements in order so assignments are known before use.
        for stmt in tree.body:
            for node in ast.walk(stmt):
                if not isinstance(node, ast.Call) or not node.args:
                    continue
                name = call_name(node)
                if name in READ_FUNCS:
                    op = "read"
                elif name in WRITE_FUNCS:
                    op = "write"
                else:
                    continue
                arg = node.args[0]
                val = try_eval(arg, ns)
                if isinstance(val, (str, Path)) and str(val).endswith(DATA_SUFFIXES):
                    ops.append((op, Path(val)))
            if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1 \
                    and isinstance(stmt.targets[0], ast.Name):
                val = try_eval(stmt.value, ns)
                if isinstance(val, (str, Path)):
                    ns[stmt.targets[0].id] = val
    return ops


def rel_to_root(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


# ── Checks ────────────────────────────────────────────────────────────────────

def undefined_names(rel: str, src: str) -> list[Issue]:
    """Names used but never defined anywhere in the notebook (hidden kernel state).

    Uses pyflakes when it is installed; otherwise this check is skipped.
    """
    try:
        from pyflakes.api import check
        from pyflakes.reporter import Reporter
    except ImportError:
        return []
    import io
    out = io.StringIO()
    check(src, rel, Reporter(out, io.StringIO()))
    names = sorted({line.rsplit("undefined name ", 1)[1].strip("'")
                    for line in out.getvalue().splitlines() if "undefined name" in line})
    return [Issue(rel, "undefined", f"'{n}' is used but never defined in this notebook "
                                    f"(it only worked with leftover kernel state)") for n in names]


def check_notebook(rel: str, path: Path) -> tuple[list[Issue], list[tuple[str, Path]]]:
    try:
        cells = code_cells(path)
    except (json.JSONDecodeError, KeyError) as e:
        return [Issue(rel, "json", f"not a valid notebook ({e.__class__.__name__})")], []

    issues: list[Issue] = []
    trees = [t for t in map(parse, cells) if t is not None]
    issues += check_imports(rel, trees)

    src = "\n".join(cells)
    for m in sorted(set(ABSPATH_RE.findall(src))):
        issues.append(Issue(rel, "abspath", f"hard-coded path {m!r}"))

    for tree in trees:
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and call_name(node) == "chdir":
                issues.append(Issue(rel, "relpath", "os.chdir() makes paths depend on where "
                                                    "the notebook is run from; use src.config paths"))

    issues += undefined_names(rel, src)

    ops = file_ops(trees)
    for op, p in ops:
        if not p.is_absolute():
            issues.append(Issue(rel, "relpath", f"{op}s {str(p)!r} relative to the working directory"))
        elif op == "read" and config.RESOURCES_DIR in p.parents and not p.exists():
            issues.append(Issue(rel, "resource", f"{rel_to_root(p)} is not in the repository"))

    # Resource paths that are loaded via np.load etc. rather than a read function.
    for tree in trees:
        for node in ast.walk(tree):
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
                val = try_eval(node, base_namespace())
                if isinstance(val, Path) and config.RESOURCES_DIR in val.parents \
                        and val.suffix and not val.exists():
                    issues.append(Issue(rel, "resource", f"{rel_to_root(val)} is not in the repository"))
    return issues, ops


def check_all() -> list[Issue]:
    manifest = load_manifest()
    order = pipeline_order(manifest)
    in_pipeline = {nb for nb, _ in order}
    issues: list[Issue] = []

    for nb, _ in order:
        if not (NB_DIR / nb).exists():
            issues.append(Issue(nb, "manifest", "listed in pipeline.yaml but missing"))

    ops_by_nb: dict[str, list] = {}
    for path in sorted(NB_DIR.rglob("*.ipynb")):
        if ".ipynb_checkpoints" in path.parts:
            continue
        rel = str(path.relative_to(NB_DIR))
        listed = rel in in_pipeline or is_excluded(rel, manifest)
        if not listed:
            issues.append(Issue(rel, "manifest", "not listed in pipeline.yaml"))
        if is_excluded(rel, manifest) and rel.startswith(("dge/archive/", "exploratory/")):
            continue  # frozen history; not held to these checks
        nb_issues, ops = check_notebook(rel, path)
        issues += nb_issues
        ops_by_nb[rel] = ops

    # Flow: walk the pipeline in order, tracking which files exist.
    data_root = config.DATA_DIR
    produced: set[Path] = set()
    written_by_dataset: dict[str, set[Path]] = {}
    for nb, ds in order:
        for op, p in ops_by_nb.get(nb, []):
            if not p.is_absolute() or data_root not in p.parents or p.suffix != ".h5ad":
                continue
            if op == "read" and config.RAW_DIR not in p.parents and p not in produced:
                issues.append(Issue(nb, "flow", f"reads {rel_to_root(p)}, which no earlier step writes"))
            elif op == "write":
                produced.add(p)
                if ds:
                    written_by_dataset.setdefault(ds, set()).add(p)

    for key, entry in config.DATASETS.items():
        written = written_by_dataset.get(key, set())
        for stage in ("clustered", "labelled"):
            if entry.get(stage) and entry[stage] not in written:
                issues.append(Issue(f"src/config.py DATASETS[{key!r}]", "registry",
                                    f"{stage} = {rel_to_root(entry[stage])} is not written by any {key} step"))
    return issues


def main() -> int:
    issues = check_all()
    if not issues:
        print("All notebook checks passed.")
        return 0
    by_kind: dict[str, list[Issue]] = {}
    for i in issues:
        by_kind.setdefault(i.kind, []).append(i)
    for kind, items in by_kind.items():
        print(f"\n{kind} ({len(items)})")
        for i in items:
            print(f"  {i.notebook}: {i.message}")
    print(f"\n{len(issues)} issue(s)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
