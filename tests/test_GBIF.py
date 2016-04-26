import unittest
from iSDM.species import GBIFSpecies
import pandas as pd


class TestGBIF(unittest.TestCase):

    def test_GBIF_DATA_api(self):
        test_species = GBIFSpecies(name_species="Etheostoma_blennioides")
        test_species.find_species_occurrences()
        self.assertEqual(test_species.ID, 2382397)
        self.assertEqual(test_species.name_species, "Etheostoma_blennioides")
        self.assertIsInstance(test_species.data_full, pd.DataFrame)

        test_species1 = GBIFSpecies(name_species="Some_nonsense")
        with self.assertRaises(ValueError):
            test_species1.find_species_occurrences()

    def test_GBIF_data_CSV(self):
        test_species = GBIFSpecies(name_species="Etheostoma_blennioides")
        test_species.load_csv("./data/GBIF.csv")
        self.assertEqual(test_species.ID, 2382397)
        self.assertIsInstance(test_species.data_full, pd.DataFrame)
        self.assertIsNotNone(test_species.data_full)

if __name__ == '__main__':
    unittest.main()
