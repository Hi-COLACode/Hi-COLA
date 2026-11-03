"""Check pinned baseline resolution and availability to fresh clones."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from baseline import prepare_baseline


class BaselineTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root/'source'
        self.target = self.root/'target'
        self.run_git('init', '-b', 'main', str(self.source))
        self.run_git('-C', str(self.source), '-c', 'user.name=Baseline test',
                     '-c', 'user.email=baseline@example.invalid', 'commit',
                     '--allow-empty', '-m', 'Baseline')
        self.commit = self.run_git('-C', str(self.source), 'rev-parse', 'HEAD').strip()
        self.run_git('-C', str(self.source), 'tag', 'baseline-v1')
        self.run_git('init', '--bare', str(self.target))
        self.config = self.root/'baseline.json'
        self.write_config(self.commit)

    def run_git(self, *args):
        return subprocess.check_output(['git', *args], text=True, stderr=subprocess.STDOUT)

    def write_config(self, commit):
        self.config.write_text(json.dumps(dict(repository=str(self.source),
            ref='refs/tags/baseline-v1', commit=commit)))

    def test_baseline_is_available_to_a_fresh_clone(self):
        self.assertEqual(prepare_baseline(self.target, self.config), self.commit)
        checkout = self.root/'checkout'
        self.run_git('clone', '--no-local', '--no-checkout', str(self.target), str(checkout))
        actual = self.run_git('-C', str(checkout), 'rev-parse', self.commit+'^{commit}').strip()
        self.assertEqual(actual, self.commit)

    def test_ref_pointing_at_an_unexpected_commit_is_rejected(self):
        self.write_config('0'*40)
        with self.assertRaisesRegex(ValueError, 'expected'):
            prepare_baseline(self.target, self.config)
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_git('-C', str(self.target), 'rev-parse', '--verify',
                         'refs/tags/hicola-regression-baseline')
