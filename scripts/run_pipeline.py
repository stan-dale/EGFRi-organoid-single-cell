#!/usr/bin/env python
"""Run the analysis notebooks in order, as listed in pipeline.yaml.

Examples
--------
    python scripts/run_pipeline.py --list                  # show the plan
    python scripts/run_pipeline.py D2_DZ                   # one dataset, all steps
    python scripts/run_pipeline.py D2_DZ D2_Lapa --dry-run # print what would run
    python scripts/run_pipeline.py D10_Lapa --from 3-unsupervised_assignment.ipynb
    python scripts/run_pipeline.py --all                   # every dataset, then DGE + GSEA
    python scripts/run_pipeline.py --step dge              # every DGE comparison
    python scripts/run_pipeline.py --step dge --comparison d2_lapa_iscs

Each notebook runs from its own folder with the Jupyter kernel named after its
conda env in pipeline.yaml (single-cell-env, decoupler-env, pydeseq-env). If you
use one environment for everything, pass --kernel <name> to override.

Executed copies are written to data/executed-notebooks/ (gitignored), so the
tracked notebooks are left untouched unless you pass --inplace.

Requires papermill (included in single-cell-env).
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.config import DATA_DIR, DATASETS, DGE_OUTPUT_DIR, NOTEBOOKS_DIR  # noqa: E402

EXECUTED_DIR = DATA_DIR / "executed-notebooks"


@dataclass
class Job:
    notebook: Path               # absolute path to the source notebook
    kernel: str
    label: str                   # what to print
    output_name: str             # file name for the executed copy
    parameters: dict | None = None
    needs: list[Path] | None = None   # files that must exist before running


def load_manifest() -> dict:
    return yaml.safe_load((ROOT / "pipeline.yaml").read_text())


def dataset_jobs(manifest: dict, key: str, kernel: str | None) -> list[Job]:
    ds = manifest["datasets"][key]
    raw = DATASETS[key]["raw"]
    jobs = []
    for i, step in enumerate(ds["steps"]):
        nb = NOTEBOOKS_DIR / ds["folder"] / step["nb"]
        jobs.append(Job(
            notebook=nb,
            kernel=kernel or manifest["envs"][step["env"]],
            label=f"{key}: {step['nb']}",
            output_name=step["nb"],
            needs=[raw] if i == 0 else None,
        ))
    return jobs


def comparison_keys(only: list[str] | None) -> list[str]:
    cfg = yaml.safe_load((NOTEBOOKS_DIR / "dge" / "dge_comparisons.yaml").read_text())
    keys = list(cfg)
    if only:
        unknown = sorted(set(only) - set(keys))
        if unknown:
            sys.exit(f"Unknown comparison(s): {', '.join(unknown)}.\nAvailable: {', '.join(keys)}")
        keys = [k for k in keys if k in only]
    return keys


def cross_jobs(manifest: dict, step_name: str, comparisons: list[str] | None,
               kernel: str | None) -> list[Job]:
    """DGE or GSEA: one job per comparison in dge_comparisons.yaml."""
    step = next(s for s in manifest["cross_dataset"] if Path(s["nb"]).parts[0] == step_name)
    nb = NOTEBOOKS_DIR / step["nb"]
    cfg = yaml.safe_load((NOTEBOOKS_DIR / "dge" / "dge_comparisons.yaml").read_text())
    jobs = []
    for key in comparison_keys(comparisons):
        if step_name == "dge":
            needs = [DATASETS[d]["labelled"] for d in cfg[key]["datasets"]]
        else:
            needs = [DGE_OUTPUT_DIR / f"{key}_deseq2_results.csv"]
        jobs.append(Job(
            notebook=nb,
            kernel=kernel or manifest["envs"][step["env"]],
            label=f"{step_name}: {key}",
            output_name=f"{nb.stem}_{key}.ipynb",
            parameters={"COMPARISON": key},
            needs=needs,
        ))
    return jobs


def print_plan(manifest: dict) -> None:
    print("Datasets (run in this order within each dataset):")
    for key, ds in manifest["datasets"].items():
        raw = DATASETS[key]["raw"]
        status = "raw input present" if raw.exists() else f"needs {raw.relative_to(ROOT)}"
        print(f"\n  {key}  [notebooks/{ds['folder']}]  ({status})")
        for step in ds["steps"]:
            print(f"    {step['nb']:<52} {manifest['envs'][step['env']]}")
    print("\nThen, across datasets (one run per comparison in dge/dge_comparisons.yaml):")
    for step in manifest["cross_dataset"]:
        print(f"    {step['nb']:<52} {manifest['envs'][step['env']]}")
    print(f"\nComparisons: {', '.join(comparison_keys(None))}")


def run(jobs: list[Job], *, dry_run: bool, inplace: bool, keep_going: bool) -> int:
    if dry_run:
        for j in jobs:
            params = f"  {j.parameters}" if j.parameters else ""
            print(f"would run  {j.label:<62} kernel={j.kernel}{params}")
        return 0

    import papermill as pm  # imported here so --list/--dry-run work without it

    failures = []
    for n, job in enumerate(jobs, 1):
        missing = [p for p in (job.needs or []) if not p.exists()]
        if missing:
            msg = "missing input(s):\n" + "\n".join(f"      {p.relative_to(ROOT)}" for p in missing)
            print(f"[{n}/{len(jobs)}] SKIP {job.label}: {msg}")
            failures.append((job, msg))
            if not keep_going:
                break
            continue

        out = job.notebook if inplace else (
            EXECUTED_DIR / job.notebook.parent.relative_to(NOTEBOOKS_DIR) / job.output_name)
        out.parent.mkdir(parents=True, exist_ok=True)
        print(f"[{n}/{len(jobs)}] {job.label} (kernel {job.kernel}) ...", flush=True)
        t0 = time.time()
        try:
            pm.execute_notebook(
                str(job.notebook), str(out),
                parameters=job.parameters or {},
                kernel_name=job.kernel,
                cwd=str(job.notebook.parent),
                progress_bar=False,
            )
        except Exception as e:  # papermill raises PapermillExecutionError with the cell
            print(f"    FAILED after {time.time() - t0:.0f}s: {str(e).strip().splitlines()[-1]}")
            print(f"    executed notebook with the traceback: {out}")
            failures.append((job, str(e)))
            if not keep_going:
                break
            continue
        print(f"    done in {time.time() - t0:.0f}s -> {out.relative_to(ROOT)}")

    if failures:
        print(f"\n{len(failures)} step(s) did not complete.")
        return 1
    print("\nAll steps completed.")
    return 0


def main(argv=None) -> int:
    manifest = load_manifest()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("datasets", nargs="*", metavar="DATASET",
                   help=f"dataset keys: {', '.join(manifest['datasets'])}")
    p.add_argument("--all", action="store_true", help="every dataset, then DGE and GSEA")
    p.add_argument("--step", choices=["dge", "gsea"], help="run only a cross-dataset step")
    p.add_argument("--comparison", action="append", help="limit --step to these comparison keys")
    p.add_argument("--from", dest="start", metavar="NOTEBOOK",
                   help="with one dataset: start at this notebook (resume after a failure)")
    p.add_argument("--kernel", help="use this Jupyter kernel for every notebook")
    p.add_argument("--list", action="store_true", help="print the pipeline and exit")
    p.add_argument("--dry-run", action="store_true", help="print the jobs without running them")
    p.add_argument("--inplace", action="store_true",
                   help="overwrite the tracked notebooks with their executed versions")
    p.add_argument("--keep-going", action="store_true", help="continue after a failed step")
    args = p.parse_args(argv)

    if args.list:
        print_plan(manifest)
        return 0

    unknown = [d for d in args.datasets if d not in manifest["datasets"]]
    if unknown:
        p.error(f"unknown dataset(s) {unknown}; choose from {list(manifest['datasets'])}")
    if not (args.datasets or args.all or args.step):
        p.error("name dataset(s), or pass --all, --step or --list")

    jobs: list[Job] = []
    for key in (manifest["datasets"] if args.all else args.datasets):
        jobs += dataset_jobs(manifest, key, args.kernel)
    if args.start:
        if len(args.datasets) != 1:
            p.error("--from needs exactly one dataset")
        names = [j.output_name for j in jobs]
        if args.start not in names:
            p.error(f"{args.start} is not a step of {args.datasets[0]}: {names}")
        jobs = jobs[names.index(args.start):]
        jobs[0].needs = None
    steps = ["dge", "gsea"] if args.all else ([args.step] if args.step else [])
    for s in steps:
        jobs += cross_jobs(manifest, s, args.comparison, args.kernel)

    return run(jobs, dry_run=args.dry_run, inplace=args.inplace, keep_going=args.keep_going)


if __name__ == "__main__":
    sys.exit(main())
