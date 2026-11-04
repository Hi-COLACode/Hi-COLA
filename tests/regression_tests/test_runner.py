"""Exercise checkout resolution with a real local Git transport."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from run import clone


class CheckoutTests(unittest.TestCase):
    def test_named_candidate_branch_resolves_after_fresh_clone(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'repository'
            subprocess.run(['git','init','-b','main',str(source)],check=True,capture_output=True)
            subprocess.run(['git','-C',str(source),'-c','user.name=Regression test',
                            '-c','user.email=regression@example.invalid','commit',
                            '--allow-empty','-m','Fixture'],check=True,capture_output=True)
            subprocess.run(['git','-C',str(source),'branch','candidate'],check=True)
            expected = subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
            result = clone(str(source),'candidate',root/'checkout')
            self.assertEqual(result['commit'],expected)
            self.assertEqual(result['submodules'],'')


if __name__=='__main__':
    unittest.main()
