import unittest

import numpy as np

import orderparam

class TestOrderParam(unittest.TestCase):
    """Test order parameters on a simple configuration with only two particles."""

    def setUp(self):
        # two particles
        self.positions = np.array([[0.5, 0.5, 0.5],
                                   [1.0, 1.0, 1.0]])

        # (most of) the actual values of the parameters don't matter
        # for the test, it only matters that they are defined.
        self.params = {'npartot': 2,
                       'stillsep': 0.5,
                       'nparsurf': 0,
                       'zperiodic': True,
                       'q6link': 0.5,
                       'q6numlinks': 1,
                       'usenearest': False,
                       'lboxx': 3.0,
                       'lboxy': 3.0,
                       'lboxz': 3.0
                       }

        # most of the functions in orderparam.py should return (0.0,)
        self.zerotuple = (0.0,)

    def test_allfracld_cpp(self):
        allfracs = orderparam.allfracld_cpp(self.positions, self.params)
        # all the system should be liquid
        self.assertEqual(allfracs, (0.0, 0.0, 0.0, 1.0, 0.0, 0.0))

    def test_allfracldtf_cpp(self):
        allfracs = orderparam.allfracldtf_cpp(self.positions, self.params)
        # all the system should be liquid
        self.assertEqual(allfracs, (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0))

    def test_fracld_cpp(self):
        fracld = orderparam.fracld_cpp(self.positions, self.params)
        self.assertEqual(fracld, self.zerotuple)

    def test_fractf_cpp(self):
        fractf = orderparam.fractf_cpp(self.positions, self.params)
        self.assertEqual(fractf, self.zerotuple)

    def test_nclusbcld_cpp(self):
        npoly = orderparam.nclusbcld_cpp(self.positions, self.params)
        self.assertEqual(npoly, (0, 0))

    def test_nclusld_cpp(self):
        nld = orderparam.nclustf_cpp(self.positions, self.params)
        self.assertEqual(nld, self.zerotuple)

    def test_nclustf_cpp(self):
        ntf = orderparam.nclustf_cpp(self.positions, self.params)
        self.assertEqual(ntf, self.zerotuple)

    def test_q6global_cpp(self):
        q6 = orderparam.q6global_cpp(self.positions, self.params)
        self.assertEqual(q6, self.zerotuple)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestOrderParam)
    unittest.TextTestRunner(verbosity=2).run(suite)
