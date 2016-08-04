import unittest
from iSDM.species import IUCNSpecies
# import pandas as pd
import geopandas as gp
from shapely.geometry import MultiPolygon
from rasterio.transform import Affine
import numpy as np


class TestIUCN(unittest.TestCase):

    def setUp(self):
        self.test_species = IUCNSpecies(name_species="Acrocheilus alutaceus")

    def test_IUCN_load_shapefile(self):
        with self.assertRaises(TypeError):
            self.test_species.load_shapefile()
            with self.assertRaises(OSError):
                self.test_species.load_shapefile("./nonexistent_folder/file.shp")
        self.test_species.load_shapefile('./data/fish/selection/acrocheilus_alutaceus/acrocheilus_alutaceus.shp')
        self.assertIsInstance(self.test_species.data_full, gp.GeoDataFrame)
        self.assertIsNotNone(self.test_species.data_full.geometry)
        self.assertIsInstance(self.test_species.data_full.geometry, gp.geoseries.GeoSeries)
        self.assertIsInstance(self.test_species.data_full.geometry.iat[0], MultiPolygon)

    def test_IUCN_find_species_occurrences(self):
        with self.assertRaises(AttributeError):
            self.test_species.find_species_occurrences()
        self.test_species.load_shapefile("./data/fish/selection/acrocheilus_alutaceus/acrocheilus_alutaceus.shp")
        self.test_species.find_species_occurrences()
        self.assertEqual(self.test_species.ID, 201940.0)
        self.assertIsNotNone(self.test_species.data_full.binomial)
        with self.assertRaises(ValueError):
            self.test_species.find_species_occurrences(name_species="Some name")

    def test_IUCN_save_shapefile(self):
        self.test_species.load_shapefile("./data/fish/selection/acrocheilus_alutaceus/acrocheilus_alutaceus.shp")
        with self.assertRaises(AttributeError):
            self.save_shapefile()
            self.save_shapefile(overwrite=False)

    def test_IUCN_rasterize(self):
        with self.assertRaises(AttributeError):
            self.test_species.rasterize()
        self.test_species.load_shapefile("./data/fish/selection/acrocheilus_alutaceus/acrocheilus_alutaceus.shp")
        result = self.test_species.rasterize(pixel_size=0.5, raster_file="./data/fish/tmp.tif")
        self.assertEqual(result.shape, (360, 720))
        self.assertIsNotNone(self.test_species.raster_file)
        self.assertIsInstance(self.test_species.raster_affine, Affine)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(set(np.unique(result)), {0, 1})
        result1 = self.test_species.rasterize(pixel_size=0.5, raster_file="./data/fish/tmp.tif", no_data_value=55, default_value=11)
        self.assertEqual(set(np.unique(result1)), {55, 11})


if __name__ == '__main__':
    unittest.main()
