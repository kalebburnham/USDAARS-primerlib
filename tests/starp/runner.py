import unittest

import test_AmasFactory
import test_AmasPrimer
import test_BindingSites
import test_PairwiseParser
import test_Parsers
import test_Sequence
import test_SingleBlastParser
import test_Snp
import test_XmlParser

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromModule(test_AmasFactory))
    suite.addTests(loader.loadTestsFromModule(test_AmasPrimer))
    suite.addTests(loader.loadTestsFromModule(test_BindingSites))
    suite.addTests(loader.loadTestsFromModule(test_PairwiseParser))
    suite.addTests(loader.loadTestsFromModule(test_Parsers))
    suite.addTests(loader.loadTestsFromModule(test_Sequence))
    suite.addTests(loader.loadTestsFromModule(test_SingleBlastParser))
    suite.addTests(loader.loadTestsFromModule(test_Snp))
    suite.addTests(loader.loadTestsFromModule(test_XmlParser))
    return suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite())