import unittest
from iSDM.environment import ClimateLayer
from iSDM.environment import ContinentsLayer
from iSDM.environment import Source
# from iSDM.species import IUCNSpecies
# import pandas as pd
import geopandas as gp
from shapely.geometry import Polygon
from rasterio.transform import Affine
import numpy as np
# import logging
# import sys

# logger = logging.getLogger()
# logger.level = logging.DEBUG
# stream_handler = logging.StreamHandler(sys.stdout)
# logger.addHandler(stream_handler)


class TestsEnvironment(unittest.TestCase):

    def setUp(self):
        self.climate_layer = ClimateLayer(file_path="./data/watertemp/max_wt_2000.tif")
        self.climate_layer_bad = ClimateLayer()
        self.biomes_layer = ClimateLayer(file_path="./data/rebioms/w001001.adf")
        continents = ContinentsLayer(file_path="./data/continents/continent.shp", source=Source.ARCGIS)
        continents.load_data()
        continents_rasters = continents.rasterize(raster_file="./data/continents/continents_raster.tif", pixel_size=0.5)
        continents_rasters[0] = continents_rasters[0] + continents_rasters[2]   # combine Europe and Asia
        continents_rasters[0][continents_rasters[0] > 1] = 1
        self.continents_rasters = np.delete(continents_rasters, 2, 0)

    def test_RasterEnvironmentalLayer_load_data(self):
        with self.assertRaises(AttributeError):
            self.climate_layer_bad.load_data()
        self.climate_layer.load_data()
        self.assertIsNotNone(self.climate_layer.raster_reader)
        self.assertEqual(self.climate_layer.resolution, (0.5, 0.5))
        self.assertIsInstance(self.climate_layer.raster_affine, Affine)
        self.assertIsInstance(self.climate_layer.metadata, dict)

    def test_RasterEnvironmentalLayer_pixel_to_world_coordinates(self):
        with self.assertRaises(AttributeError):
            self.climate_layer_bad.pixel_to_world_coordinates()
            self.climate_layer.read(1)
        self.climate_layer.load_data()
        band = self.climate_layer.read(1)
        self.assertEqual(band.shape, (360, 720))
        coordinates = self.climate_layer.pixel_to_world_coordinates()
        self.assertIsInstance(coordinates, tuple)
        self.assertIsInstance(coordinates[0], np.ndarray)
        self.assertEqual(len(coordinates[0]), 259200)
        self.assertEqual(coordinates[0][0], 89.75)

    def test_RasterEnvironmentalLayer_polygonize(self):
        with self.assertRaises(AttributeError):
            self.climate_layer.polygonize()
        self.climate_layer.load_data()
        df_polygons = self.climate_layer.polygonize()
        self.assertIsInstance(df_polygons, gp.GeoDataFrame)
        self.assertIsNotNone(df_polygons.geometry)
        self.assertIsInstance(df_polygons.geometry.iat[0], Polygon)

    def test_RasterEnvironmentalLayer_close_dataset(self):
        with self.assertRaises(AttributeError):
            self.climate_layer.close_dataset()
        self.climate_layer.load_data()
        self.assertFalse(self.climate_layer.raster_reader.closed)
        self.climate_layer.close_dataset()
        self.assertTrue(self.climate_layer.raster_reader.closed)
        self.assertIsNotNone(self.climate_layer.raster_reader)

    def test_RasterEnvironmentalLayer_read(self):
        with self.assertRaises(AttributeError):
            self.climate_layer.read(1)
        self.climate_layer.load_data()
        band = self.climate_layer.read(1)
        self.assertEqual(band.shape, (360, 720))
        self.assertIsInstance(band, np.ndarray)
        with self.assertRaises(IndexError):
            self.climate_layer.read(2)

    def test_RasterEnvironmentalLayer_reproject(self):
        self.climate_layer.load_data()
        original_resolution = self.climate_layer.resolution
        self.climate_layer.reproject(destination_file="./data/tmp.tif", resolution=(original_resolution[0] * 2, original_resolution[1] * 2))
        self.climate_layer.load_data("./data/tmp.tif")
        self.assertEqual(original_resolution, (self.climate_layer.resolution[0] / 2, self.climate_layer.resolution[1] / 2))

    def test_RasterEnvironmentalLayer_sample_pseudo_absences(self):
        self.biomes_layer.load_data()
        some_species = np.ones_like(self.biomes_layer.read(1))
        some_species[0][0] = 0  # only one pixel set to zero, species covers entire range
        pixels_to_sample_from, sampled_pixels = self.biomes_layer.sample_pseudo_absences(species_raster_data=some_species)
        self.assertIsNone(sampled_pixels)

        # set half of the pixels to 0, now the species covers about half of the map
        for index in range(int(some_species.shape[0] / 2)):
            some_species[index] = 0
        pixels_to_sample_from, sampled_pixels = self.biomes_layer.sample_pseudo_absences(species_raster_data=some_species)
        self.assertIsNotNone(sampled_pixels)
        self.assertIsInstance(pixels_to_sample_from, np.ndarray)
        self.assertIsInstance(sampled_pixels, np.ndarray)
        self.assertEqual(pixels_to_sample_from.shape, self.biomes_layer.read(1).shape)
        self.assertEqual(sampled_pixels.shape, self.biomes_layer.read(1).shape)

        pixels_to_sample_from_1, sampled_pixels_1 = self.biomes_layer.sample_pseudo_absences(species_raster_data=some_species,
                                                                                             continents_raster_data=self.continents_rasters)
        self.assertIsNotNone(sampled_pixels_1)
        self.assertIsInstance(pixels_to_sample_from_1, np.ndarray)
        self.assertIsInstance(sampled_pixels_1, np.ndarray)       
        self.assertGreater(np.sum(pixels_to_sample_from), np.sum(pixels_to_sample_from_1))

    def tearDown(self):
        del self.climate_layer

if __name__ == '__main__':
    # logging.getLogger( "iSDM.environment" ).setLevel( logging.DEBUG )
    unittest.main()
