"""Load and fetch the numerical baseline specified in baseline.json."""
import argparse
import json
from pathlib import Path
import subprocess

CONFIG = Path(__file__).with_name('baseline.json')
FETCHED_REF = 'refs/tags/hicola-regression-baseline'


def load_baseline(path=CONFIG):
    return json.loads(Path(path).read_text())


def prepare_baseline(repository, path=CONFIG):
    config = load_baseline(path)
    subprocess.run(['git', '-C', str(repository), 'fetch', '--no-tags',
                    config['repository'], config['ref']], check=True)
    actual = subprocess.check_output(
        ['git', '-C', str(repository), 'rev-parse', 'FETCH_HEAD^{commit}'],
        text=True).strip()
    if actual != config['commit']:
        raise ValueError(f"Baseline ref resolves to {actual}, expected {config['commit']}")
    # A named ref keeps the verified commit available to fresh local clones.
    subprocess.run(['git', '-C', str(repository), 'update-ref', FETCHED_REF, actual], check=True)
    return actual


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path.cwd())
    args = parser.parse_args()
    print(prepare_baseline(args.repository))
