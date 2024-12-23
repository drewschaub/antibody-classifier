import unittest
from src.antibody_classifier.data_io import load_structures_from_directory
import os

class TestDataIO(unittest.TestCase):
    def test_load_structures_from_directory(self):
        # Suppose you have a small test PDB file in test_data/
        test_dir = "tests/test_data"
        structures = load_structures_from_directory(test_dir)
        self.assertGreater(len(structures), 0, "Should load at least one structure")

if __name__ == '__main__':
    unittest.main()
