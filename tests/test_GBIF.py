import unittest
from iSDM.species import GBIFSpecies
import pandas as pd
import geopandas as gp
from shapely.geometry import Point


class TestGBIF(unittest.TestCase):

    def setUp(self):
        self.test_species = GBIFSpecies(name_species="Pseudecheneis crassicauda")
        self.test_species1 = GBIFSpecies(name_species="Some_nonsense")
        self.test_species2 = GBIFSpecies(name_species="Etheostoma_blennioides")

    def test_GBIF_find_species_occurrences(self):
        self.test_species.find_species_occurrences()
        self.assertEqual(self.test_species.ID, 2341467)
        self.assertEqual(self.test_species.name_species, "Pseudecheneis crassicauda")
        self.assertIsInstance(self.test_species.data_full, pd.DataFrame)

        # with self.assertRaises(ValueError):
        #    test_species1.find_species_occurrences()
        data_empty = self.test_species1.find_species_occurrences()
        self.assertIsInstance(data_empty, pd.DataFrame)
        self.assertEqual(data_empty.empty, True)

    def test_GBIF_load_csv(self):
        self.test_species2.load_csv("./data/GBIF.csv")
        self.assertEqual(self.test_species2.ID, 2382397)
        self.assertIsInstance(self.test_species2.data_full, pd.DataFrame)
        self.assertIsNotNone(self.test_species2.data_full)

    def tearDown(self):
        self.test_species = None
        self.test_species1 = None
        self.test_species2 = None

    def test_GBIF_geometrize(self):
        self.test_species.find_species_occurrences()
        self.test_species.geometrize()
        self.assertIsInstance(self.test_species.data_full, gp.GeoDataFrame)
        self.assertIsNotNone(self.test_species.data_full.geometry)
        self.assertIsInstance(self.test_species.data_full.geometry, gp.geoseries.GeoSeries)
        self.assertIsInstance(self.test_species.data_full.geometry.iat[0], Point)
        self.assertIsNotNone(self.test_species.data_full.crs)

        self.test_species.find_species_occurrences()
        number_species = self.test_species.data_full.shape[0]
        self.test_species.geometrize(dropna=False)
        self.assertEqual(number_species, self.test_species.data_full.shape[0])

if __name__ == '__main__':
    unittest.main()
