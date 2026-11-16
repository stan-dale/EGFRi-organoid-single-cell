"""Static checks that every pipeline notebook can run on a fresh clone.

See scripts/check_notebooks.py for what is checked. Runs without data.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from check_notebooks import check_all  # noqa: E402

# Issues that need files from outside the repository. Each entry is
# (notebook, kind, substring of the message). Remove an entry once fixed;
# test_known_issues_still_present fails if an entry has gone stale.
KNOWN_ISSUES = [
    ("d10-lapa/4-eec_subanalysis.ipynb", "import", "src/labels.py does not exist"),
    ("d10-lapa/4-eec_subanalysis.ipynb", "resource", "d10_eec_canonical_embedding.npz"),
]


def _is_known(issue):
    return any(issue.notebook == nb and issue.kind == kind and text in issue.message
               for nb, kind, text in KNOWN_ISSUES)


@pytest.fixture(scope="module")
def issues():
    return check_all()


def test_no_new_notebook_issues(issues):
    new = [str(i) for i in issues if not _is_known(i)]
    assert not new, "\n" + "\n".join(new)


@pytest.mark.parametrize("known", KNOWN_ISSUES, ids=lambda k: k[2])
def test_known_issues_still_present(issues, known):
    nb, kind, text = known
    assert any(i.notebook == nb and i.kind == kind and text in i.message for i in issues), (
        f"{known} is fixed: remove it from KNOWN_ISSUES")
