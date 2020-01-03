import unittest

import test_Sequence
import test_Primer
import test_NestedLoop

import sys

""" 
HOW TO ADD MODULES
First, import the test module following the convention above.
Then, use suite.addTests in suite() using the same pattern.
"""


def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromModule(test_Sequence))
    suite.addTests(loader.loadTestsFromModule(test_Primer))
    suite.addTests(loader.loadTestsFromModule(test_NestedLoop))
    return suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite())