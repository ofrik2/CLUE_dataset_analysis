import os
import sys
from pathlib import Path
import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--pathway-files",
        action="append",
        default=[],
        help="One or more pathway mapping file paths (comma-separated or repeated). Example: --pathway-files=/tmp/pathways.csv"
    )
    parser.addoption(
        "--signature-files",
        action="append",
        default=[],
        help="One or more signature file paths (comma-separated or repeated). Example: --signature-files=/tmp/signature.csv"
    )

@pytest.fixture(scope="session", autouse=True)
def add_src_to_path():
    """
    Ensure the project's 'src' directory is on sys.path so tests can import modules in src/.
    """
    repo_root = Path(__file__).resolve().parents[1]  # tests/ -> project root
    src_dir = repo_root / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))
    yield

def _collect_paths(option_values, env_var_name):
    raw = []
    for entry in option_values:
        if not entry:
            continue
        raw.extend([p.strip() for p in entry.split(",") if p.strip()])
    if not raw:
        env = os.environ.get(env_var_name, "")
        if env:
            raw.extend([p.strip() for p in env.split(",") if p.strip()])
    return [Path(p) for p in raw if p]

@pytest.fixture(scope="session")
def pathway_files(request):
    paths = _collect_paths(request.config.getoption("--pathway-files"), "PATHWAY_FILES")
    if not paths:
        pytest.skip("No pathway files provided. Use --pathway-files or set PATHWAY_FILES to run pathway IO tests.")
    return paths

@pytest.fixture(scope="session")
def signature_files(request):
    paths = _collect_paths(request.config.getoption("--signature-files"), "SIGNATURE_FILES")
    if not paths:
        pytest.skip("No signature files provided. Use --signature-files or set SIGNATURE_FILES to run signature IO tests.")
    return paths
