import unittest
import rosalindFunctions as rf

#Most basic approach
#assert rf.independent_alleles(kth_gen=2, n_min=1) == 0.68359375

class TestRosalindFunctions(unittest.TestCase):
    def test_independent_alleles(self):
        result = rf.independent_alleles(kth_gen=2, n_min=1)
        self.assertAlmostEqual(result, 0.68359375)

if __name__ == '__main__':
    unittest.main()