
"""
A module for preparing all environmental layers used for fitting a model, into a single dataframe, of a convenient format.

      .. moduleauthor:: Daniela Remenska <remenska@gmail.com>

"""
from enum import Enum
from iSDM.environment import RasterEnvironmentalLayer, VectorEnvironmentalLayer
import pandas as pd
import numpy as np
from rasterio.transform import Affine
import gc
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
    """
    Model

    A class for encapsulating the data preparation functionality for fitting a model. Makes sure all layers are at a consistent
    resolution, global scale. Constructs the "base" dataframe and allows adding environmental layers data to this dataframe incrementally.

    :ivar pixel_size: The size of the pixel in degrees, i.e., the resolution to use for rasterizing.

    :ivar x_res: The number of pixels in the global-scale raster, along the x-axis.

    :ivar y_res: The number of pixels in the global-scale raster, along the x-axis.

    :ivar base_layer: A "base" raster environmental layer to use for deducing the world coordinates used as index in the base dataframe.
    All burned pixels (value==1) will be converted to global-scale latitude/longitude coordinates (the index).

    :vartype base_layer: RasterEnvironmentalLayer

    :ivar base_dataframe: A base dataframe containing latitude/longitude columns as an index. Further environmental layers added
    to the model will be merged with this base dataframe. Therefore, it is expected to contain all the necessary latitude/longitude combinations
    that may show up in any further environmental layers. If raster_data is not provided as a base(see below), then all world coordinates
    (at a particular resolution) are taken into account.

    :vartype base_dataframe: pandas.DataFrame
    """
    def __init__(self, pixel_size, raster_data=None, **kwargs):
        logger.info("Preparing a base dataframe...")
        x_min, y_min, x_max, y_max = -180, -90, 180, 90  # global
        self.pixel_size = pixel_size
        self.x_res = int((x_max - x_min) / self.pixel_size)
        self.y_res = int((y_max - y_min) / self.pixel_size)
        self.base_layer = RasterEnvironmentalLayer()
        self.base_layer.raster_affine = Affine(pixel_size, 0.0, x_min, 0.0, -pixel_size, y_max)
        logger.info("Base layer: Computing world coordinates...")
        if raster_data is not None:
            all_coordinates = self.base_layer.pixel_to_world_coordinates(raster_data=raster_data)
        else:
            all_coordinates = self.base_layer.pixel_to_world_coordinates(raster_data=np.ones((self.y_res, self.x_res), dtype=np.uint8))
        # self.base_dataframe = pd.DataFrame([all_coordinates[0], all_coordinates[1]]).T
        # self.base_dataframe.columns = ['decimallatitude', 'decimallongitude']
        self.base_dataframe = pd.DataFrame(columns=['decimallatitude', 'decimallongitude'])
        self.base_dataframe['decimallatitude'] = np.array(all_coordinates[0])
        self.base_dataframe['decimallongitude'] = np.array(all_coordinates[1])
        self.base_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
        del all_coordinates
        gc.collect()

    def cross_validate(self, percentage, random_seed):
        # Needs to be implemented
        pass

    def evaluate_performance(self, method=Evaluation.ROC, **kwargs):
        self.method = method

    def fit(self):
        pass

    def add_environmental_layer(self, layer, discard_threshold=None, discard_nodata_value=True):
        """
        Adds an environmental layer (Raster or Vector) by converting its values to latitude/longitude decimals, and merging with the
        base dataframe. If the added layer is already raster, the first band is expected to contain all the (pixel) data values.
        The resolution of the provided raster needs to match the Model resolution. The pixels are then converted to world coordinates,
        and merged with the base dataframe. If the environmental layer is vector, then it is first rasterized, and the logic proceeds
        the same way as with a raster layer.

        :param EnvironmentalLayer layer: The environmental layer to be added to the model dataframe.

        :param int discard_threshold: Optional pixel value to use for discarding layer pixels below a certain value, before adding the layer.
        Default is `pickle`. Another possibility is "msgpack", as it has shown as 10% more efficient in terms of time and memory, \
        for the type of data we are dealing with.

        :param bool discard_nodata_value: Optionally filter out "nodata" pixel values from the raster, when converting the layer
        pixels to world coordinates.

        """
        logger.info("Loading environmental layer from %s " % layer.file_path)
        if isinstance(layer, RasterEnvironmentalLayer):
            layer_reader = layer.load_data()
            layer_data = layer_reader.read(1)
            if layer_data.shape != (self.y_res, self.x_res):
                logger.error("The layer is not at the proper resolution! Layer shape:%s " % (layer_data.shape, ))
                return
            logger.info("Computing world coordinates...")
            if discard_nodata_value:
                logger.info("Filtering out no_data pixels.")
                layer_data = np.where(layer_data != layer_reader.nodata, layer_data, np.nan)
            if discard_threshold is not None:
                logger.info("Discarding values below %s " % discard_threshold)
                layer_data = np.where(layer_data > discard_threshold, layer_data, np.nan)

            layer_coordinates = layer.pixel_to_world_coordinates(raster_data=layer_data,
                                                                 filter_no_data_value=discard_nodata_value,
                                                                 no_data_value=layer_reader.nodata)

            logger.info("Constructing dataframe for %s  ..." % layer.name_layer)
            # layer_dataframe = pd.DataFrame([layer_coordinates[0], layer_coordinates[1]]).T
            # layer_dataframe.columns = ['decimallatitude', 'decimallongitude']
            layer_dataframe = pd.DataFrame(columns=['decimallatitude', 'decimallongitude'])
            layer_dataframe['decimallatitude'] = np.array(layer_coordinates[0])
            layer_dataframe['decimallongitude'] = np.array(layer_coordinates[1])
            del layer_coordinates
            gc.collect()
            flattened_data = layer_data.reshape(np.product(layer_data.shape))
            # if discard_threshold:
            #     logger.info("Ignoring values below %s " % discard_threshold)
            #     flattened_data = flattened_data[flattened_data > discard_threshold]
            layer_dataframe[layer.name_layer] = flattened_data[~np.isnan(flattened_data)]
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
            gc.collect()
        elif isinstance(layer, VectorEnvironmentalLayer):
            layer.load_data()
            # if not hasattr(layer, 'pixel_size'):
            #     logger.error("Please provide a pixel size to use for rasterizing the vector layer.")
            #     return
            if not hasattr(layer, 'raster_file'):
                logger.error("Please provide a target raster_file location.")
                return
            if not hasattr(layer, 'classifier_column'):
                logger.info("No classifier column specified; will rasterize everything in one band...")
                layer.set_classifier(classifier_column=None)
            raster_data = layer.rasterize(raster_file=layer.get_raster_file(), pixel_size=self.pixel_size, classifier_column=layer.get_classifier())
            logger.info("Raster data has shape: %s " % (raster_data.shape,))
            self.base_dataframe[layer.name_layer] = np.nan
            for idx, band in enumerate(raster_data):
                logger.info("Computing world coordinates for band number %s " % idx)
                band_coordinates = self.base_layer.pixel_to_world_coordinates(raster_data=band)
                logger.info("Constructing dataframe for band number %s " % idx)
                # band_dataframe = pd.DataFrame([band_coordinates[0], band_coordinates[1]]).T
                # band_dataframe.columns = ['decimallatitude', 'decimallongitude']
                band_dataframe = pd.DataFrame(columns=['decimallatitude', 'decimallongitude'])
                band_dataframe['decimallatitude'] = np.array(band_coordinates[0])
                band_dataframe['decimallongitude'] = np.array(band_coordinates[1])
                band_dataframe[layer.name_layer] = idx + 1
                band_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
                logger.info("Finished constructing dataframe for band number %s " % idx)
                logger.info("Merging with base dataframe...")
                self.base_dataframe.update(band_dataframe)  # much fastter than combine_first
                logger.info("Shape of base_dataframe: %s " % (self.base_dataframe.shape, ))
                logger.info("Finished processing band number %s " % idx)
                del band_dataframe
                gc.collect()

    def get_base_dataframe(self):
        return self.base_dataframe
