"""Check that every command line tool starts and answers --help."""
import unittest

from tools import find_entry_points, find_scripts, find_tools, run_tool


class HelpTests(unittest.TestCase):
    def test_tools_found(self):
        self.assertTrue(find_scripts())
        self.assertTrue(find_entry_points())

    def test_help(self):
        for label, command in find_tools():
            with self.subTest(tool=label):
                result = run_tool(command, '--help')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)


if __name__ == '__main__':
    unittest.main()
