import os
import re
import subprocess
import sys
import tomllib
import unittest
from importlib.metadata import version
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPO_ROOT / 'HiCOLA'
# A module-level parse_args() call marks a file that is run directly as a script
SCRIPT_PATTERN = re.compile(r'^\S.*\.parse_args\(', re.MULTILINE)


def find_scripts():
    """Return the package files that parse command line arguments when run directly."""
    return sorted(
        path for path in PACKAGE_ROOT.rglob('*.py')
        if SCRIPT_PATTERN.search(path.read_text(encoding='utf-8'))
    )


def find_entry_points():
    """Return the console scripts declared in pyproject.toml as {name: 'module:function'}."""
    with open(REPO_ROOT / 'pyproject.toml', 'rb') as handle:
        return tomllib.load(handle).get('project', {}).get('scripts', {})


def find_tools():
    """Return (label, command) pairs for every command line tool in the package."""
    tools = [(str(path.relative_to(REPO_ROOT)), [sys.executable, str(path)])
             for path in find_scripts()]
    for name, target in sorted(find_entry_points().items()):
        module, function = target.split(':')
        code = (f'import sys; sys.argv[0] = {name!r}; '
                f'from {module} import {function}; sys.exit({function}())')
        tools.append((name, [sys.executable, '-c', code]))
    return tools


def run_tool(command, *args):
    env = dict(os.environ, MPLBACKEND='Agg')
    return subprocess.run(
        [*command, *args],
        cwd=REPO_ROOT, env=env, capture_output=True, text=True, timeout=120,
    )


class CliHelpTests(unittest.TestCase):
    def test_tools_found(self):
        self.assertTrue(find_scripts())
        self.assertTrue(find_entry_points())

    def test_help(self):
        for label, command in find_tools():
            with self.subTest(tool=label):
                result = run_tool(command, '--help')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)

    @unittest.skip('--version is not implemented for the command line tools yet')
    def test_version(self):
        for label, command in find_tools():
            with self.subTest(tool=label):
                result = run_tool(command, '--version')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(version('Hi_COLA'), result.stdout)


if __name__ == '__main__':
    unittest.main()
