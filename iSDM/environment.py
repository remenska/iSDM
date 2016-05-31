
import logging
from enum import Enum
import rasterio
import pprint
from rasterio.warp import calculate_default_transform, RESAMPLING
from geopandas import GeoSeries
from iSDM.species import IUCNSpecies
from affine import Affine
import rasterio.features
import numpy as np


logger = logging.getLogger('iSDM.environment')
logger.setLevel(logging.DEBUG)


class Source(Enum):
    WORLDCLIM = 1
    GLOBE = 2


class EnvironmentalLayer(object):
    def __init__(self, source=None, file_path=None, **kwargs):
        # you want to be able to agregate at a different resolution
        # and back/forth, right?
        # self.resolution = kwargs['resolution']
        if source:
            if not isinstance(source, Source):
                raise AttributeError("The source can only be one of the following: %s " % list(Source.__members__))
            self.source = source
        if file_path:
            self.file_path = file_path

    def load_dataset(self, file_path=None):
        pass

    def set_source(self, source):
        if not isinstance(source, Source):
            raise AttributeError("The source can only be one of the following: %s " % list(Source.__members__))
        self.source = source

    def get_source(self):
        return self.source.name


class ClimateLayer(EnvironmentalLayer):
    def __init__(self, source=None, file_path=None, **kwargs):
        EnvironmentalLayer.__init__(self, source, file_path, **kwargs)

    def load_dataset(self, file_path=None):
        if file_path:
            self.file_path = file_path

        if not self.file_path:
            raise AttributeError("Please provide a file_path argument to load the data from.")

        logger.info("Loading data from %s " % self.file_path)
        climate_data = rasterio.open(self.file_path, 'r')
        self.metadata = climate_data.meta
        self.resolution = climate_data.res
        self.bounds = climate_data.bounds
        pp = pprint.PrettyPrinter(depth=5)

        logger.info("Metadata: %s " % pp.pformat(self.metadata))
        logger.info("Resolution: %s " % (self.resolution,))
        logger.info("Bounds: %s " % (self.bounds,))
        self.climate_data = climate_data
        logger.info("Dataset loaded. Use .read() or .read_masks() to access the layers.")
        return self.climate_data

    def close_dataset(self):
        if not self.climate_data.closed:
            self.climate_data.close()
            logger.info("Dataset %s closed. " % self.file_path)

    def get_data(self):
        if not self.climate_data or self.climate_data.closed:
            logger.info("The dataset is closed. Please load it first using .load_dataset()")
            return
        return self.climate_data

    def reproject(self, source_file=None, destination_file=None, resampling=RESAMPLING.nearest, **kwargs):
        if not source_file:
            if not self.file_path:
                raise AttributeError("Please provide a source_file to load the data from.")
            else:
                source_file = self.file_path

        with rasterio.open(source_file) as src:
            affine, width, height = calculate_default_transform(src_crs=src.crs,
                                                                dst_crs=kwargs.get('dst_crs', src.crs),
                                                                width=src.width,
                                                                height=src.height,
                                                                left=src.bounds.left,
                                                                bottom=src.bounds.bottom,
                                                                right=src.bounds.right,
                                                                top=src.bounds.top,
                                                                resolution=kwargs.get('resolution', src.res)
                                                                )
            logger.info("Calculated default transformation:")
            logger.info("Affine:\n{0} \n width={1}, height={2}".format(affine, width, height))

            kwargs = src.meta.copy()
            kwargs.update({'transform': affine,
                           'affine': affine,
                           'width': width,
                           'height': height
                           })

            with rasterio.open(destination_file, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    rasterio.warp.reproject(source=rasterio.band(src, i),
                                            destination=rasterio.band(dst, i),
                                            src_transform=src.affine,
                                            src_crs=src.crs,
                                            dst_transform=affine,
                                            dst_crs=kwargs.get('dst_crs', src.crs),
                                            resampling=resampling
                                            )
            logger.info("Reprojected data in %s " % destination_file)

    def overlay(self, range_map, plot=False):
        """
        To extract data from a raster only where it intersects with a vector feature.
        rasterio provides option to "burn" vector shapes into rasters (rasterize the geometry). Then we create
        a raster mask layer
        """
        if not (isinstance(range_map, GeoSeries)
                or isinstance(range_map, EnvironmentalLayer)
                or isinstance(range_map, IUCNSpecies)):
            raise AttributeError("Please provide a correct rangemap input.")

        # TODO: what about if the overlay is just a shape file?
        # TODO: what if there are multiple geometries in a shapefile? Currently it just returns the last
        # but could make a list and append
        masked_data = None

        if isinstance(range_map, EnvironmentalLayer) or isinstance(range_map, IUCNSpecies):
            for geometry in range_map.get_data()['geometry']:
                if self.climate_data.closed:
                    self.load_dataset(self.file_path)
                with self.climate_data as raster:
                    # get pixel coordinates of the geometry's bounding box
                    ul = raster.index(*geometry.bounds[0:2])
                    lr = raster.index(*geometry.bounds[2:4])

                    # read the subset of the data into a numpy array
                    window = ((lr[0], ul[0] + 1), (ul[1], lr[1] + 1))
                    data = raster.read(1, window=window)
                    # create an affine transform for the subset data
                    t = raster.affine
                    shifted_affine = Affine(t.a, t.b, t.c + ul[1] * t.a, t.d, t.e, t.f + lr[0] * t.e)
                    # rasterize the geometry
                    mask = rasterio.features.rasterize(
                        [(geometry, 0)],
                        out_shape=data.shape,
                        transform=shifted_affine,
                        fill=1,
                        all_touched=True,
                        dtype=np.uint8)

                    # create a masked numpy array
                    masked_data = np.ma.array(data=data, mask=mask.astype(bool))

        self.masked_data = masked_data
        logger.info("Overlayed raster climate data with the given range map.")
        logger.info("Use the .masked_data attribute to access it.")


class LandCoverlayer(EnvironmentalLayer):   # this would e a shapefile, not a raster?
    pass


class LandUseLayer(EnvironmentalLayer):
    pass


class DEMLayer(EnvironmentalLayer):
    pass
