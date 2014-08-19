import unittest

import writeoutput

class TestWriteOutput(unittest.TestCase):
    """Test writing of output files."""

    def setUp(self):
        pass


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestOrderParam)
    unittest.TextTestRunner(verbosity=2).run(suite)
