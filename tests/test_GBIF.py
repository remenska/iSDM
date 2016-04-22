import unittest
from iSDM.species import GBIFSpecies
import pandas as pd

class TestGBIF(unittest.TestCase):

    def test_GBIF(self):
      test_species = GBIFSpecies(name_species="Etheostoma_blennioides")
      test_species.find_species_occurrences()
      self.assertEqual(test_species.ID, 2382397)
      self.assertEqual(test_species.name_species, "Etheostoma_blennioides")
      self.assertIsInstance(test_species.data_full, pd.DataFrame)

if __name__ == '__main__':
    unittest.main()

