import json
import unittest

from nestedloop.nestedloop import NestedLoop

class TestNestedLoop(unittest.TestCase):

    def test_one_pair(self):
        """ one_pair data files are constructed so that there will be
        exactly one possibility for f and r primers. Here we check
        that NL found them as expected and attributes them correctly. """
        datafile = 'data/one_pair.xml'
        with open(datafile) as f:
            xml = f.read()
        inputdata = 'data/one_pair.json'
        with open(inputdata) as f:
            data = json.load(f)

        nl = NestedLoop(data['sequence'],
                        (data['tm_min'], data['tm_opt'], data['tm_max']),
                        data['f_from'],
                        data['f_to'],
                        data['r_from'],
                        data['r_to'],
                        data['pcr_min'],
                        data['pcr_max'],
                        data['num_to_return'],
                        data['custom_forward_primer'],
                        data['custom_reverse_primer'],
                        xml)

        nl.run()
        self.assertEqual(nl.pairs[0].reverse_primer.sequence, data['r']['sequence'])
        self.assertEqual(nl.pairs[0].forward_primer.sequence, data['f']['sequence'])
        

        self.assertEqual(nl.pairs[0].forward_primer.start, data['f']['start'])
        self.assertEqual(nl.pairs[0].forward_primer.end, data['f']['end'])
        self.assertEqual(nl.pairs[0].forward_primer.__len__(), data['f']['len'])
        self.assertEqual(nl.pairs[0].forward_primer.gc, data['f']['gc'])
        self.assertEqual(nl.pairs[0].forward_primer.complementary_score, 
                         data['f']['self_complementary'])

        self.assertEqual(nl.pairs[0].reverse_primer.start, data['r']['start'])
        self.assertEqual(nl.pairs[0].reverse_primer.end, data['r']['end'])
        self.assertEqual(nl.pairs[0].reverse_primer.__len__(), data['r']['len'])
        # Asssert almost equal because this data is a float. 
        self.assertAlmostEqual(nl.pairs[0].reverse_primer.gc, data['r']['gc'])
        self.assertEqual(nl.pairs[0].reverse_primer.complementary_score, 
                         data['r']['self_complementary'])


if __name__ == "__main__":
    unittest.main()

