import unittest
import os
from Code.maxcfl import modify_cfl

# FILE: test_maxcfl.py


class TestModifyCFL(unittest.TestCase):

    def setUp(self):
        # Create a temporary input file for testing
        self.test_filename = 'test_input.txt'
        with open(self.test_filename, 'w') as f:
            f.write("CaseName\n")
            f.write("1.4 287.0\n")
            f.write("0.4 0.5 0.0001\n")
            f.write("5000\n")
            f.write("53 37\n")

    def tearDown(self):
        # Remove the temporary input file after tests
        if os.path.exists(self.test_filename):
            os.remove(self.test_filename)

    def test_modify_cfl(self):
        new_cfl = 1.0
        modify_cfl(self.test_filename, new_cfl)
        
        with open(self.test_filename, 'r') as f:
            lines = f.readlines()
        
        # Check if the CFL value is correctly modified
        self.assertEqual(lines[2].split()[0], str(new_cfl))

if __name__ == '__main__':
    unittest.main()