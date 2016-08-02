
"""Documentation on environment module."""
import logging
from enum import Enum
import rasterio
import pprint
from rasterio.warp import calculate_default_transform, RESAMPLING
from iSDM.species import IUCNSpecies
import numpy as np
from geopandas import GeoSeries, GeoDataFrame
from rasterio.transform import Affine
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
from shapely.geometry import Point
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch
import shapely
from rasterio import features
from shapely.geometry import Polygon

logger = logging.getLogger('iSDM.environment')
logger.setLevel(logging.DEBUG)


class Source(Enum):
    WORLDCLIM = 1
    GLOBE = 2
    UNKNOWN = 3
    TNC = 4
    ARCGIS = 5


class EnvironmentalLayer(object):
    """
    Some class-level documentaiton on environmental layer class
    """
    def __init__(self, source=None, file_path=None, name_layer=None, **kwargs):
        # you want to be able to agregate at a different resolution
        # and back/forth, right?
        # self.resolution = kwargs['resolution']
        if source:
            if not isinstance(source, Source):
                raise AttributeError("The source can only be one of the following: %s " % list(Source.__members__))
            self.source = source
        else:
            self.source = Source.UNKNOWN
        if file_path:
            self.file_path = file_path
        self.name_layer = name_layer

    def load_data(self, file_path=None):
        pass

    def set_source(self, source):
        if not isinstance(source, Source):
            raise AttributeError("The source can only be one of the following: %s " % list(Source.__members__))
        self.source = source

    def get_source(self):
        return self.source.name

    def save_data(self, full_name=None, dir_name=None, file_name=None):
        pass

    def get_data(self):
        pass


class RasterEnvironmentalLayer(EnvironmentalLayer):
    """
    Some class-level documentation on raster environmental layer
    """
    def __init__(self, source=None, file_path=None, name_layer=None, **kwargs):
        EnvironmentalLayer.__init__(self, source, file_path, name_layer, **kwargs)

    # def load_data(self, file_path=None):
    #     """
    #     Documentation pending on how to load environment raster data
    #     """
    #     if file_path:
    #         self.file_path = file_path

    #     if not self.file_path:
    #         raise AttributeError("Please provide a file_path argument to load the data from.")

    #     logger.info("Loading data from %s " % self.file_path)
    #     raster_reader = rasterio.open(self.file_path, 'r')
    #     self.metadata = raster_reader.meta
    #     self.resolution = raster_reader.res
    #     self.bounds = raster_reader.bounds
    #     pp = pprint.PrettyPrinter(depth=5)

    #     logger.info("Metadata: %s " % pp.pformat(self.metadata))
    #     logger.info("Resolution: %s " % (self.resolution,))
    #     logger.info("Bounds: %s " % (self.bounds,))
    #     self.raster_reader = raster_reader
    #     logger.info("Dataset loaded. Use .read() or .read_masks() to access the layers.")
    #     return self.raster_reader

    def load_data(self, file_path=None):
        """
        Documentation pending on how to load environment aster data
        """
        if file_path:
            self.file_path = file_path
        if not self.file_path:
            raise AttributeError("Please provide a file_path to read raster environment data from.")

        src = rasterio.open(self.file_path)
        logger.info("Loaded raster data from %s " % self.file_path)
        logger.info("Driver name: %s " % src.driver)
        pp = pprint.PrettyPrinter(depth=5)
        self.metadata = src.meta
        logger.info("Metadata: %s " % pp.pformat(self.metadata))
        logger.info("Resolution: x_res={0} y_res={1}.".format(src.width, src.height))
        logger.info("Bounds: %s " % (src.bounds,))
        logger.info("Coordinate reference system: %s " % src.crs)
        logger.info("Affine transformation: %s " % (src.affine.to_gdal(),))
        logger.info("Number of layers: %s " % src.count)
        logger.info("Dataset loaded. Use .read() or .read_masks() to access the layers.")
        self.raster_affine = src.affine
        self.resolution = src.res
        self.bounds = src.bounds
        self.raster_reader = src
        return self.raster_reader

    def pixel_to_world_coordinates(self,
                                   raster_data=None,
                                   no_data_value=0,
                                   filter_no_data_value=True,
                                   band_number=1,
                                   plot=False):
        """
        Map the pixel coordinates to world coordinates. The affine transformation matrix
        is used for this purpose. The convention is to reference the pixel corner. To
        reference the pixel center instead, we translate each pixel by 50%.
        The "no value" pixels (cells) can be filtered out.

        A dataset's pixel coordinate system has its origin at the "upper left" (imagine it displayed on your screen).
        Column index increases to the right, and row index increases downward. The mapping of these coordinates to
        "world" coordinates in the dataset's reference system is done with an affine transformation matrix.

        :returns: a tuple of arrays. The first array contains the latitude values for each
        non-zero cell, the second array contains the longitude values for each non-zero cell.
        """
        if raster_data is None:
            logger.info("No raster data provided, attempting to load default...")
            try:
                raster_data = self.load_data(self.file_path).read(band_number)  # we work on one layer, the first
                logger.info("Succesfully loaded existing raster data from %s." % self.file_path)
            except AttributeError as e:
                logger.error("Could not open raster file. %s " % str(e))

        logger.info("Transforming to world coordinates...")
        # first get the original Affine transformation matrix
        T0 = self.raster_affine
        # convert it to gdal format (it is otherwise flipped)
        T0 = Affine(*reversed(T0.to_gdal()))
        logger.debug("Affine transformation T0:\n %s " % (T0,))
        # shift by 50% to get the pixel center
        T1 = T0 * Affine.translation(0.5, 0.5)
        # apply the shift, filtering out no_data_value values
        logger.debug("Raster data shape: %s " % (raster_data.shape,))
        logger.debug("Affine transformation T1:\n %s " % (T1,))
        raster_data[raster_data == no_data_value] = 0  # reset the nodata values to 0, easier to manipulate
        if filter_no_data_value:
            logger.info("Filtering out no_data pixels.")
            coordinates = (T1 * np.where(raster_data != no_data_value))
        else:
            logger.info("Not filtering any no_data pixels.")
            coordinates = (T1 * np.where(np.ones_like(raster_data)))

        if plot:
            self.plot_world_coordinates(coordinates)

        logger.info("Transformation to world coordinates completed.")
        return coordinates

    @classmethod
    def __geometrize__(cls, data, latitude_col_name='decimallatitude', longitude_col_name='decimallongitude', crs=None, dropna=True):
        try:
            if crs is None:
                crs = {'init': "EPSG:4326"}
            # exclude those points with NaN in plot_world_coordinates
            if dropna:
                geometry = [Point(xy) for xy in zip(
                    data[longitude_col_name].dropna(),
                    data[latitude_col_name].dropna()
                )]
                data_geo = GeoDataFrame(
                    data.dropna(
                        subset=[latitude_col_name, longitude_col_name]),
                    crs=crs,
                    geometry=geometry)
                logger.info("Data geometrized: converted into GeoPandas dataframe.")
                logger.info("Points with NaN coordinates ignored. ")
            else:
                geometry = [Point(xy) for xy in zip(
                    data[longitude_col_name],
                    data[latitude_col_name]
                )]
                data_geo = GeoDataFrame(data, crs=crs, geometry=geometry)
                logger.info("Data geometrized: converted into GeoPandas dataframe.")
        except AttributeError:
            logger.error("No latitude/longitude data to convert into a geometry. Please load the data first.")
        return data_geo

    def polygonize(self, band_number=1):
        raster_data = self.read(band_number)
        mask = raster_data != self.raster_reader.nodata
        T0 = self.raster_reader.affine
        shapes = features.shapes(raster_data, mask=mask, transform=T0)
        df = GeoDataFrame.from_records(shapes, columns=['geometry', 'value'])
        # convert the geometry dictionary from {'coordinates': [[(-73.5, 83.5),
        #   (-73.5, 83.0),
        #   (-68.0, 83.0),
        #   (-68.0, 83.5),
        #   (-73.5, 83.5)]],
        # 'type': 'Polygon'}
        # to a proper shapely polygon format
        df.geometry = df.geometry.apply(lambda row: Polygon(row['coordinates'][0]))
        df.crs = self.raster_reader.crs
        return df   # TODO: maybe return here a VectorEnvironmentalLayer?

    @classmethod
    def plot_world_coordinates(cls,
                               coordinates=None,
                               figsize=(16, 12),
                               projection='merc',
                               facecolor='crimson'):
        if coordinates is None or not isinstance(coordinates, tuple):
            logger.error("Please provide the coordinates to plot, in the correct format.")
            logger.error("Use pixel_to_world_coordinates() for this.")
            return
        # empty dataset
        if coordinates[0].shape[0] == 0:
            logger.error("No data to plot.")
            return
        data = pd.DataFrame([coordinates[0], coordinates[1]]).T
        data.columns = ['decimallatitude', 'decimallongitude']
        data_geometrized = cls.__geometrize__(data)
        mm = Basemap(projection=projection, lat_0=50, lon_0=-100,
                     resolution='l', area_thresh=1000.0,
                     llcrnrlon=data_geometrized.geometry.total_bounds[0],  # lower left corner longitude point
                     llcrnrlat=data_geometrized.geometry.total_bounds[1],  # lower left corner latitude point
                     urcrnrlon=data_geometrized.geometry.total_bounds[2],  # upper right longitude point
                     urcrnrlat=data_geometrized.geometry.total_bounds[3]   # upper right latitude point
                     )

        # prepare longitude/latitude list for basemap
        ax1 = plt.subplots(figsize=figsize)[1]

        plt.title("World coordinates of raster data.")

        mm.drawcoastlines()
        mm.drawcountries()
        mm.drawrivers(color='lightskyblue', linewidth=1.5)
        mm.drawmapboundary(fill_color='lightskyblue')
        mm.fillcontinents(color='cornsilk')
        # draw latitude and longitude
        mm.drawmeridians(np.arange(-180, 180, 10), labels=[False, False, False, True])
        mm.drawparallels(np.arange(-180, 180, 10), labels=[True, True, False, False])

        patches = []
        selection = data_geometrized
        for poly in selection.geometry:
            if poly.geom_type == 'Polygon':
                mpoly = shapely.ops.transform(mm, poly)
                patches.append(PolygonPatch(mpoly))
            elif poly.geom_type == 'MultiPolygon':
                for subpoly in poly:
                    mpoly = shapely.ops.transform(mm, subpoly)
                    patches.append(PolygonPatch(mpoly))
            elif poly.geom_type == "Point":
                patches.append(PolygonPatch(Point(mm(poly.x, poly.y)).buffer(9999)))  # TODO: this buffer thing is tricky around the edges of a map
            else:
                logger.warning("Geometry type %s not supported. Skipping ... " % poly.geom_type)
                continue
        ax1.add_collection(PatchCollection(patches, facecolor=facecolor, match_original=True, zorder=100))
        plt.show()

    def plot(self, figsize=(25, 20), band_number=1):
        if not self.raster_reader or self.raster_reader.closed:
            logger.info("The dataset is closed. Please load it first using .load_data()")
            return
        plt.figure(figsize=figsize)
        plt.title(self.name_layer)
        plt.imshow(self.read(band_number), cmap="flag", interpolation="none")

    def close_dataset(self):
        """
        Documentation pending on what it means to close a dataset
        """
        if not self.raster_reader.closed:
            self.raster_reader.close()
            logger.info("Dataset %s closed. " % self.file_path)

    def get_data(self):
        """
        Documentation on what it means to get the data
        """
        if not self.raster_reader or self.raster_reader.closed:
            logger.info("The dataset is closed. Please load it first using .load_data()")
            return
        return self.raster_reader

    def read(self, band_number=None):
        """
        Documentation pending on reading raster bands
        """
        if not self.raster_reader or self.raster_reader.closed:
            logger.info("The dataset is closed. Please load it first using .load_data()")
            return
        return self.raster_reader.read(band_number)

    def reproject(self, source_file=None, destination_file=None, resampling=RESAMPLING.nearest, **kwargs):
        """
        Documentation pending on how to reproject/resample data
        """
        if not source_file:
            if not self.file_path:
                raise AttributeError("Please provide a source_file to load the data from.")
            else:
                source_file = self.file_path

        with rasterio.open(source_file) as src:
            affine, width, height = calculate_default_transform(src_crs=src.crs,
                                                                dst_crs=kwargs.get('dst_crs', src.crs),
                                                                width=kwargs.get('width', src.width),
                                                                height=kwargs.get('height', src.height),
                                                                left=kwargs.get('left', src.bounds.left),
                                                                bottom=kwargs.get('bottom', src.bounds.bottom),
                                                                right=kwargs.get('right', src.bounds.right),
                                                                top=kwargs.get('top', src.bounds.top),
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
        if not (isinstance(range_map, EnvironmentalLayer) or isinstance(range_map, IUCNSpecies)):
            raise AttributeError("Please provide a correct rangemap input.")

        # TODO: what about if the overlay is just a shape file?
        # TODO: what if there are multiple geometries in a shapefile? Currently it just returns the last
        # but could make a list and append
        masked_data = None

        if isinstance(range_map, EnvironmentalLayer) or isinstance(range_map, IUCNSpecies):
            for geometry in range_map.get_data()['geometry']:
                if self.raster_reader.closed:
                    self.load_data(self.file_path)
                with self.raster_reader as raster:
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
                    mask = features.rasterize(
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

    def sample_pseudo_absences(self,
                               species_raster_data=None,
                               band_number=1,
                               number_of_pseudopoints=1000):
        if not (isinstance(species_raster_data, np.ndarray)) or not (set(np.unique(species_raster_data)) == set({0, 1})):
            logger.error("Please provide the species raster data as a numpy array with pixel values 1 and 0 (presence/absence).")
            return
        try:
            env_raster_data = self.read(band_number)
            logger.info("Succesfully loaded existing raster data from %s." % self.file_path)
        except AttributeError as e:
            logger.error("Could not open raster file. %s " % str(e))

        if species_raster_data.shape != env_raster_data.shape:
            logger.error("Please provide (global) species raster data at the same resolution as the environment")
            logger.error("Environment data has the following shape %s " % (env_raster_data.shape, ))
            return

        logger.info("Sampling %s pseudo-absence points from environmental layer." % number_of_pseudopoints)
        # first set to zero all pixels that have "nodata" values
        env_raster_data[env_raster_data == self.raster_reader.nodata] = 0
        # next get all the overlapping pixels between the species raster and the environment data
        presences_pixels = env_raster_data * species_raster_data
        # what are the unique values left? (these are the distinct "regions" that need to be taken into account)
        # Do NOT take into account the 0-value pixel, which we assigned to all "nodata" pixels
        unique_regions = np.unique(presences_pixels[presences_pixels != 0])
        if len(unique_regions) == 0:
            logger.info("There are no environmental layers to sample pseudo-absences from. ")
            return
        logger.debug("The following unique (pixel) values will be taken into account for sampling pseudo-absences")
        logger.debug(unique_regions)
        # add the pixels of all these regions to a layer array
        regions = []
        for region in unique_regions:
            regions.append(np.where(env_raster_data == region))
        # now "regions" contains a list of tuples, each tuple with separate x/y indexes (arrays thereof) of the pixels
        # make an empty "base" matrix and fill it with the selected regions pixel values
        selected_pixels = np.zeros_like(env_raster_data)
        # pick out only those layers that have been selected and fill in the matrix
        for layer in regions:
            selected_pixels[layer] = env_raster_data[layer]

        # sample from those pixels which are in the selected raster regions, minus those of the species presences
        pixels_to_sample_from = selected_pixels - presences_pixels
        # These are x/y positions of pixels to sample from. Tuple of arrays.
        (x, y) = np.where(pixels_to_sample_from > 0)
        number_pixels_to_sample_from = x.shape[0]  # == y.shape[0] since every pixel has (x,y) position.
        logger.info("There are %s pixels to sample from..." % (number_pixels_to_sample_from))
        sampled_pixels = np.zeros_like(selected_pixels)

        if number_pixels_to_sample_from < number_of_pseudopoints:
            logger.warning("There are less pixels to sample from, than the desired number of pseudo-absences")
            logger.warning("Will select all pixels as psedudo-absences.")
            random_indices = np.arange(0, number_pixels_to_sample_from)
        else:
            # now randomly choose <number_of_pseudopoints> indices to fill in with pseudo absences
            random_indices = np.random.randint(0, number_pixels_to_sample_from, number_of_pseudopoints)
            logger.info("Filling %s random pixel positions..." % (len(random_indices)))

            # fill in those indices with the pixel values of the environment layer
        for position in random_indices:
            sampled_pixels[x[position]][y[position]] = pixels_to_sample_from[x[position], y[position]]

        logger.info("Sampled %s unique pixels as pseudo-absences." % sampled_pixels.nonzero()[0].shape[0])

        return (pixels_to_sample_from, sampled_pixels)


class ClimateLayer(RasterEnvironmentalLayer):
    pass


class DEMLayer(RasterEnvironmentalLayer):
    pass


class VectorEnvironmentalLayer(EnvironmentalLayer):
    """
    Some class-level documentation on vector environmental layers here
    """
    def __init__(self, source=None, file_path=None, name_layer=None, **kwargs):
        EnvironmentalLayer.__init__(self, source, file_path, name_layer, **kwargs)

    def load_data(self, file_path=None):
        """
        Documentation on how to load shapefile environmental layer
        """
        if file_path:
                self.file_path = file_path

        if not self.file_path:
            raise AttributeError("Please provide a file_path argument to load the data from.")

        logger.info("Loading data from %s " % self.file_path)
        self.data_full = GeoDataFrame.from_file(self.file_path)
        self.data_full.columns = [x.lower() for x in self.data_full.columns]
        logger.info("The shapefile contains data on %d environmental regions." % self.data_full.shape[0])
        self.shape_file = self.file_path

    def save_data(self, full_name=None, driver='ESRI Shapefile', overwrite=False):
        """
        Saves the current (geopandas) data as a shapefile.
        The geopandas data needs to have geometry as a column, besides the metadata.
        """

        if not (isinstance(self.data_full, GeoSeries) or isinstance(self.data_full, GeoDataFrame)):
            raise AttributeError("The data is not of a Geometry type (GeoSeries or GeoDataFrame).",
                                 "Please geometrize first!")

        if overwrite:
            full_name = self.shape_file
        elif not overwrite and full_name is None:
            raise AttributeError("Please provide a shape_file location, or set overwrite=True")

        try:
            self.data_full.to_file(full_name, driver="ESRI Shapefile")
            logger.debug("Saved data: %s " % full_name)
        except AttributeError as e:
            logger.error("Could not save data! %s " % str(e))

    def get_data(self):
        return self.data_full

    def set_data(self, data_frame):
        """
        Careful, overwrites the existing raw data!. More documentation
        """
        self.data_full = data_frame

    def rasterize(self, raster_file=None, pixel_size=None, all_touched=True,
                  no_data_value=0,
                  default_value=1,
                  crs=None,
                  cropped=False,
                  *args, **kwargs):
        """
        Documentation pending on how to rasterize geometrical shapes
        """
        if not (pixel_size or raster_file):
            raise AttributeError("Please provide pixel_size and a target raster_file.")

        if not hasattr(self, 'data_full'):
            raise AttributeError("You have not loaded the data.")

        if crs is None:
            crs = {'init': "EPSG:4326"}

        # crop to the boundaries of the shape?
        if cropped:
            # cascaded_union_geometry = shapely.ops.cascaded_union(self.data_full.geometry)
            x_min, y_min, x_max, y_max = self.data_full.geometry.total_bounds
        # else global map
        else:
            x_min, y_min, x_max, y_max = -180, -90, 180, 90

        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        # translate
        transform = Affine.translation(x_min, y_max) * Affine.scale(pixel_size, -pixel_size)
        result = features.rasterize(self.data_full.geometry,
                                    transform=transform,
                                    out_shape=(y_res, x_res),
                                    all_touched=all_touched,
                                    fill=no_data_value,
                                    default_value=default_value
                                    )

        with rasterio.open(raster_file, 'w', driver='GTiff', width=x_res, height=y_res,
                           count=1,
                           dtype=np.uint8,
                           nodata=no_data_value,
                           transform=transform,
                           crs=crs) as out:
            out.write(result.astype(np.uint8), indexes=1)
            out.close()
        logger.info("RASTERIO: Data rasterized into file %s " % raster_file)
        logger.info("RASTERIO: Resolution: x_res={0} y_res={1}".format(x_res, y_res))
        self.raster_file = raster_file
        self.raster_affine = transform
        return result

    def load_raster_data(self, raster_file=None):
        """
        Documentation pending on how to load raster data
        """
        if raster_file:
            self.raster_file = raster_file
        if not self.raster_file:
            raise AttributeError("Please provide a raster_file to read raster data from.")

        src = rasterio.open(self.raster_file)
        logger.info("Loaded raster data from %s " % self.raster_file)
        logger.info("Driver name: %s " % src.driver)
        pp = pprint.PrettyPrinter(depth=5)
        self.metadata = src.meta
        logger.info("Metadata: %s " % pp.pformat(self.metadata))
        logger.info("Resolution: x_res={0} y_res={1}.".format(src.width, src.height))
        logger.info("Bounds: %s " % (src.bounds,))
        logger.info("Coordinate reference system: %s " % src.crs)
        logger.info("Affine transformation: %s " % (src.affine.to_gdal(),))
        logger.info("Number of layers: %s " % src.count)
        logger.info("Dataset loaded. Use .read() or .read_masks() to access the layers.")
        self.raster_affine = src.affine
        self.raster_reader = src
        return self.raster_reader


class ContinentsLayer(VectorEnvironmentalLayer):
    """
    Continents Layer will have a separate treatment. This is mostly because, when
    rasterizing the continents shapefile, there are multiple shapes which should end up with
    a different pixel value for each continents. The typical rasterizing operation rather
    produces a raster with binary values (0s and 1s).

    For overlaying selected pseudo-absence biomes (pixels) with continents, it is best to
    have one individual (raster band) matrix per continent, each filled with 0s and 1s, and loop through
    the list of continents (in total 8 layers, should be fast). Then "overlaying" or clipping
    with biomes would amount to simple arithmetics, like multiplying the matrices (element by element)
    to filter out anything with 0-valued pixels in any matrix.

    """
    def __init__(self, source=None, file_path=None, name_layer=None, **kwargs):
        VectorEnvironmentalLayer.__init__(self, source, file_path, name_layer, **kwargs)

    def rasterize(self, raster_file=None, pixel_size=None, all_touched=False,
                  no_data_value=0,
                  default_value=1,
                  crs=None,
                  cropped=False,
                  *args, **kwargs):
        """
        Documentation pending on how to rasterize continents
        """
        if not (pixel_size or raster_file):
            raise AttributeError("Please provide pixel_size and a target raster_file.")

        if not hasattr(self, 'data_full'):
            raise AttributeError("You have not loaded the data.")

        if crs is None:
            crs = {'init': "EPSG:4326"}

        # crop to the boundaries of the shape?
        if cropped:
            # cascaded_union_geometry = shapely.ops.cascaded_union(self.data_full.geometry)
            x_min, y_min, x_max, y_max = self.data_full.geometry.total_bounds
        # else global map
        else:
            x_min, y_min, x_max, y_max = -180, -90, 180, 90

        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)

        continents_list = self.data_full.continent.unique()
        logger.info("Will rasterize continent-by-continent.")
        # translate
        stacked_layers = []
        for continent_name in continents_list:
            logger.info("Rasterizing continent %s " % continent_name)
            transform = Affine.translation(x_min, y_max) * Affine.scale(pixel_size, -pixel_size)
            result = features.rasterize(self.data_full.geometry[self.data_full.continent == continent_name],
                                        transform=transform,
                                        out_shape=(y_res, x_res),
                                        all_touched=all_touched,
                                        fill=no_data_value,
                                        default_value=default_value
                                        )
            stacked_layers.append(result)

        stacked_layers = np.stack(stacked_layers)

        for i, band in enumerate(stacked_layers, 1):
            with rasterio.open(raster_file, 'w', driver='GTiff', width=x_res, height=y_res,
                               count=stacked_layers.shape[0],
                               dtype=np.uint8,
                               nodata=no_data_value,
                               transform=transform,
                               crs=crs) as out:
                out.write(result.astype(np.uint8), indexes=i)

        out.close()
        logger.info("RASTERIO: Data rasterized into file %s " % raster_file)
        logger.info("RASTERIO: Resolution: x_res={0} y_res={1}".format(x_res, y_res))
        self.raster_file = raster_file
        self.raster_affine = transform
        return stacked_layers

    def load_raster_data(self, raster_file=None):
        """
        Documentation pending on how to load raster data
        """
        if raster_file:
            self.raster_file = raster_file
        if not self.raster_file:
            raise AttributeError("Please provide a raster_file to read raster data from.")

        src = rasterio.open(self.raster_file)
        logger.info("Loaded raster data from %s " % self.raster_file)
        logger.info("Driver name: %s " % src.driver)
        pp = pprint.PrettyPrinter(depth=5)
        self.metadata = src.meta
        logger.info("Metadata: %s " % pp.pformat(self.metadata))
        logger.info("Resolution: x_res={0} y_res={1}.".format(src.width, src.height))
        logger.info("Bounds: %s " % (src.bounds,))
        logger.info("Coordinate reference system: %s " % src.crs)
        logger.info("Affine transformation: %s " % (src.affine.to_gdal(),))
        logger.info("Number of layers: %s " % src.count)
        logger.info("Dataset loaded. Use .read() or .read_masks() to access the layers.")
        self.raster_affine = src.affine
        self.raster_reader = src
        return self.raster_reader


class LandCoverlayer(VectorEnvironmentalLayer):
    pass
