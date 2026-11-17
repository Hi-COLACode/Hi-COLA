"""Check that every command line tool reports the package version with --version.

The tools do not implement --version yet, so there is nothing to check.
"""
import unittest


class VersionTests(unittest.TestCase):
    pass


if __name__ == '__main__':
    unittest.main()
