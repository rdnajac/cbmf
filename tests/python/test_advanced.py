from .context import module

import unittest


class AdvancedTestSuite(unittest.TestCase):
    """Advanced test cases."""

    def test_thoughts(self):
        self.assertIsNone(module.hmm())


if __name__ == '__main__':
    unittest.main()
