
import logging
from enum import Enum
import rasterio
import pprint
from rasterio.warp import calculate_default_transform, RESAMPLING
from geopandas import GeoSeries
from iSDM.species import IUCNSpecies

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

    def overlay(self, range_map):
        """
        To extract data from a raster only where it intersects with a vector feature.
        rasterio provides option to "burn" vector shapes into rasters (rasterize the geometry). Then we create
        a raster mask layer
        """
        if not (isinstance(range_map, GeoSeries)
                or isinstance(range_map, EnvironmentalLayer)
                or isinstance(range_map, IUCNSpecies)):
            raise AttributeError("Please provide a correct rangemap input.")

        logger.info("Overlayed raster climate data with the given range map.")


class ClimateLayer(EnvironmentalLayer):
    def __init__(self, source=None, file_path=None, **kwargs):
        EnvironmentalLayer.__init__(self, source, file_path, **kwargs)

    def load_dataset(self, file_path=None):
        if file_path:
            self.file_path = file_path

        if not self.file_path:
            raise AttributeError("Please provide a file_path argument to load the data from.")

        logger.info("Loading data from %s " % self.file_path)
        climate_data = rasterio.open(self.file_path)
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


class LandCoverlayer(EnvironmentalLayer):   # this would e a shapefile, not a raster?
    pass


class LandUseLayer(EnvironmentalLayer):
    pass


class DEMLayer(EnvironmentalLayer):
    pass
