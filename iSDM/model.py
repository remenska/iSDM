
"""
Another interesting module

:synopsis: Another useful module indeed, model.

.. moduleauthor:: Daniela Remenska <remenska@gmail.com>

"""
from enum import Enum
from iSDM.environment import RasterEnvironmentalLayer, VectorEnvironmentalLayer
import pandas as pd
import numpy as np

import logging
logger = logging.getLogger('iSDM.model')
logger.setLevel(logging.DEBUG)


class Algorithm(Enum):
    GAM = 1
    GLM = 2
    MAXENT = 3


class Evaluation(Enum):
    ROC = 1
    KAPPA = 2
    TSS = 3
    AUC = 4


class Model(object):
    def __init__(self, pixel_size=0.5, **kwargs):
        logger.info("Preparing a base dataframe...")
        x_min, y_min, x_max, y_max = -180, -90, 180, 90  # global
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        self.base_layer = RasterEnvironmentalLayer()
        logger.info("Base layer: Computing world coordinates...")
        all_coordinates = self.base_layer.pixel_to_world_coordinates(raster_data=np.zeros((y_res, x_res)), filter_no_data_value=False)
        self.base_dataframe = pd.DataFrame([all_coordinates[0], all_coordinates[1]]).T
        self.base_dataframe.columns = ['decimallatitude', 'decimallongitude']
        self.base_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)

    def cross_validate(self, percentage, random_seed):
        pass

    def evaluate_performance(self, method=Evaluation.ROC, **kwargs):
        self.method = method

    def fit(self):
        pass

    def add_environmental_layer(self, layer, discard_threshold=None):
        logger.info("Loading environmental layer from %s " % layer.file_path)
        if isinstance(layer, RasterEnvironmentalLayer):
            layer_reader = layer.load_data()
            layer_data = layer_reader.read(1)
            if discard_threshold is not None:
                logger.info("Discarding values below %s " % discard_threshold)
                layer_data[layer_data < discard_threshold] = 0
            logger.info("Computing world coordinates...")
            layer_coordinates = layer.pixel_to_world_coordinates(raster_data=layer_data,
                                                                 filter_no_data_value=True)

            logger.info("Constructing dataframe for %s  ..." % layer.name_layer)
            layer_dataframe = pd.DataFrame([layer_coordinates[0], layer_coordinates[1]]).T
            layer_dataframe.columns = ['decimallatitude', 'decimallongitude']
            flattened_data = layer_data.reshape(np.product(layer_data.shape))
            if discard_threshold is not None:
                logger.info("Ignoring values below %s " % discard_threshold)
                layer_dataframe[layer.name_layer] = flattened_data[flattened_data != discard_threshold]
            else:
                layer_dataframe[layer.name_layer] = flattened_data
            layer_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
            logger.info("Shape of layer_dataframe: % s " % (layer_dataframe.shape, ))
            logger.info("Finished constructing dataframe for %s ..." % layer.name_layer)
            logger.info("Merging with base dataframe...")
            self.base_dataframe = pd.merge(self.base_dataframe, layer_dataframe,
                                           how='left',
                                           left_index=True,
                                           right_index=True)
            logger.info("Shape of base_dataframe: %s " % (self.base_dataframe.shape, ))
            del layer_dataframe
        elif isinstance(layer, VectorEnvironmentalLayer):
            layer.load_data()
            if not hasattr(layer, 'pixel_size'):
                logger.error("Please provide a pixel size to use for rasterizing the vector layer.")
                return
            if not hasattr(layer, 'raster_file'):
                logger.error("Please provide a target raster_file location.")
                return
            if not hasattr(layer, 'classifier_column'):
                logger.info("No classifier column specified; will rasterize everything in one band...")
                layer.set_classifier(classifier_column=None)
            raster_data = layer.rasterize(raster_file=layer.get_raster_file(), pixel_size=layer.get_pixel_size(), classifier_column=layer.get_classifier())
            logger.info("Raster data has shape: %s " % (raster_data.shape,))
            self.base_dataframe[layer.name_layer] = np.nan
            for idx, band in enumerate(raster_data):
                logger.info("Computing world coordinates for band number %s " % idx)
                band_coordinates = self.base_layer.pixel_to_world_coordinates(raster_data=band)
                logger.info("Constructing dataframe for band number %s " % idx)
                band_dataframe = pd.DataFrame([band_coordinates[0], band_coordinates[1]]).T
                band_dataframe.columns = ['decimallatitude', 'decimallongitude']
                band_dataframe[layer.name_layer] = idx + 1
                band_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
                logger.info("Finished constructing dataframe for band number %s " % idx)
                logger.info("Merging with base dataframe...")
                self.base_dataframe.update(band_dataframe)  # much fastter than combine_first
                logger.info("Shape of base_dataframe: %s " % (self.base_dataframe.shape, ))
                logger.info("Finished processing band number %s " % idx)
                del band_dataframe

    def get_base_dataframe(self):
        return self.base_dataframe
