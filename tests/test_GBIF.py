import unittest
from iSDM.species import GBIFSpecies
import pandas as pd
import geopandas as gp
from shapely.geometry import Point
import numpy as np

class TestGBIF(unittest.TestCase):

    def setUp(self):
        self.test_species = GBIFSpecies(name_species="Pseudecheneis crassicauda")
        self.test_species1 = GBIFSpecies(name_species="Some_nonsense")
        self.test_species2 = GBIFSpecies(name_species="Etheostoma blennioides")

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

    def test_GBIF_geometrize(self):
        self.test_species.find_species_occurrences()
        self.test_species.geometrize()
        self.assertIsInstance(self.test_species.data_full, gp.GeoDataFrame)
        self.assertIsNotNone(self.test_species.data_full.geometry)
        self.assertIsInstance(self.test_species.data_full.geometry, gp.geoseries.GeoSeries)
        self.assertIsInstance(self.test_species.data_full.geometry.iat[0], Point)
        self.assertIsNotNone(self.test_species.data_full.crs)

        self.test_species.find_species_occurrences()
        number_point_records = self.test_species.data_full.shape[0]
        self.test_species.geometrize(dropna=False)
        self.assertEqual(number_point_records, self.test_species.data_full.shape[0])

    def test_GBIF_rasterize(self):
        # with self.assertRaises(AttributeError):
            # self.test_species.rasterize()
        # self.test_species2.load_csv("./data/GBIF.csv")
        self.test_species2.find_species_occurrences()
        pixel_size = 0.2
        number_point_records = self.test_species2.data_full.shape[0]
        # self.test_species.load_shapefile("./data/fish/selection/acrocheilus_alutaceus/acrocheilus_alutaceus.shp")
        result = self.test_species2.rasterize(pixel_size=pixel_size, raster_file="./data/fish/tmp.tif", all_touched=True)
        transform = self.test_species2.raster_affine
        self.assertEqual(result.shape, (np.abs(transform.yoff) * (2 / pixel_size), np.abs(transform.xoff) * (2 / pixel_size)))
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(set(np.unique(result)), {0, 1})
        self.assertGreater(number_point_records, np.sum(result))
        result1 = self.test_species2.rasterize(pixel_size=0.5, no_data_value=55, default_value=11)
        self.assertEqual(set(np.unique(result1)), {55, 11})

    def tearDown(self):
        del self.test_species
        del self.test_species1
        del self.test_species2

if __name__ == '__main__':
    unittest.main()
