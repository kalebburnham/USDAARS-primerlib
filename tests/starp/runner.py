import unittest

import test_AmasFactory
import test_AmasPrimer
import test_Parsers
import test_Sequence
import test_Snp

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromModule(test_AmasFactory))
    suite.addTests(loader.loadTestsFromModule(test_AmasPrimer))
    suite.addTests(loader.loadTestsFromModule(test_Parsers))
    suite.addTests(loader.loadTestsFromModule(test_Sequence))
    suite.addTests(loader.loadTestsFromModule(test_Snp))
    return suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite())