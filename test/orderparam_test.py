import unittest

import numpy as np

import orderparam

class TestOrderParam(unittest.TestCase):

    def setUp(self):
        # two particles
        self.positions = np.array([[0.5, 0.5, 0.5],
                                   [1.0, 1.0, 1.0]])

        # (most of) the actual values of the parameters don't matter
        # for the test, it only matters that they are defined.  note
        # however that these parameters mean that both particles are
        # identified as being fcc by LD criterion!
        self.params = {'npartot': 2,
                       'stillsep': 1,
                       'nparsurf': 0,
                       'zperiodic': True,
                       'q6link': 0.5,
                       'q6numlinks': 1,
                       'usenearest': True,
                       'lboxx': 3.0,
                       'lboxy': 3.0,
                       'lboxz': 3.0
                       }

    def test_allfracld_cpp(self):
        allfracs = orderparam.allfracld_cpp(self.positions, self.params)
        print allfracs
        pass

    def test_allfracldtf_cpp(self):
        pass

    def test_fracld_cpp(self):
        fracld = orderparam.fracld_cpp(self.positions, self.params)
        self.assertEqual(fracld, 1.0)

    def test_fractf_cpp(self):
        fractf = orderparam.fractf_cpp(self.positions, self.params)
        self.assertEqual(fractf, 0.0)

    def test_nclusld_cpp(self):
        nld = orderparam.nclustf_cpp(self.positions, self.params)
        self.assertEqual(nld, 0.0)

    def test_nclustf_cpp(self):
        ntf = orderparam.nclustf_cpp(self.positions, self.params)
        self.assertEqual(ntf, 0.0)

    def test_q6global_cpp(self):
        q6 = orderparam.q6global_cpp(self.positions, self.params)
        print q6
#        self.assertEqual(q6, 0.0)
        

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestOrderParam)
    unittest.TextTestRunner(verbosity=2).run(suite)
#    unittest.main()
