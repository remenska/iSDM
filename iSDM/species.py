"""
Documentation on the species module.
Some more text to document what this module is for.

      .. moduleauthor:: Daniela Remenska <remenska@gmail.com>

"""

from pygbif import species   # http://pygbif.readthedocs.org/en/latest/
from pygbif import occurrences
import copy
import pandas as pd
import logging
import os
import csv
from enum import Enum
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from geopandas import GeoSeries, GeoDataFrame
from osgeo import gdal, ogr
from shapely.geometry import Point
import shapely.ops
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch
from rasterio.transform import Affine
import rasterio
from rasterio import features
from shapely.prepared import prep
import pprint

logger = logging.getLogger('iSDM.species')
logger.setLevel(logging.DEBUG)


class Source(Enum):
    """
    Possible sources of global species data.
    """
    GBIF = 1
    IUCN = 2
    PREDICTS = 3
    MOL = 4


class ObservationsType(Enum):
    PRESENCE_ONLY = 1
    PRESENCE_ABSENCE = 2
    RICHNESS = 3
    ABUNDANCE = 4


class Species(object):
    """
    Species
    A generic Species class used for subclassing different global-scale species data sources.

    :ivar ID: a unique ID for a particular species. For example, for GBIF sources, it is the gbifid metadata field.
    :ivar name_species: initial value: 'Unknown'

    """

    ID = int(0)
    name_species = 'Unknown'

    def __init__(self, **kwargs):

        if 'ID' not in kwargs and 'name_species' not in kwargs:
            raise AttributeError("Cannot initialize species without a 'name_species' or an 'ID' argument supplied")

        if 'name_species' in kwargs:
            self.name_species = kwargs['name_species']
        if 'ID' in kwargs:
            self.ID = kwargs['ID']
        try:
            # Enable C-based speedups available from 1.2.10+
            from shapely import speedups
            speedups.enable()
            logger.debug("Enabled Shapely speedups for performance.")
        except ImportError:
            logger.info("Upgrade Shapely for Performance enhancements")

    def save_data(self, full_name=None, dir_name=None, file_name=None, method="pickle"):
        """
        Serializes the loaded species occurrence filtered dataset (`pandas <http://pandas.pydata.org/pandas-docs/stable/dsintro.html>`_ or `geopandas <http://geopandas.org/user.html>`_ DataFrame) into a binary `pickle <https://en.wikipedia.org/wiki/Pickle_%28Python%29>`_  file.

       :param str full_name: The full path of the file (including the directory and name in one string),
        where the data will be saved.

       :param str dir_name: The directory where the file will be stored. If :attr:`file_name` is not specified, the default one :attr:`name_species` + .pkl is given.

       :param str file_name: The name of the file where the data will be saved. If :attr:`dir_name` is not specified, the current working directory is taken by default.

       :param str method: The type of serialization to use for the data frame. Default is "pickle". Another possibility is "msgpack", as it has shown as 10$ more efficient for the type of data

       :raises: AttributeError: if the data has not been loaded in the object before. See :func:`load_data` and :func:`find_species_occurrences`

       :returns: None

        """
        if full_name is None:
            if file_name is None:
                file_name = str(self.name_species) + (".pkl" if method == "pickle" else ".msg")
            if dir_name is None:
                dir_name = os.getcwd()

            full_name = os.path.join(dir_name, file_name)

        try:
            if method == "msgpack":
                self.data_full.to_msgpack(full_name)
            elif method == "pickle":
                self.data_full.to_pickle(full_name)
            else:
                logger.error("Incorrect method of serializing: %s " % method)
                return

            logger.debug("Saved data: %s " % full_name)
            logger.debug("Type of data: %s " % type(self.data_full))
        except IOError as e:
            logger.error("Could not save data! %s " % str(e))
        except AttributeError as e:
            logger.error("No data to save. Please load it first. %s " % str(e))

    def load_data(self, file_path=None):
        """
        Loads the serialized species pickle file data into a pandas DataFrame.

        :param str file_path: The full path to the file where the data is serialized to.

        :returns: Data loaded into geopandas dataframe.

        :rtype: geopandas.GeoDataFrame
        """
        if file_path is None:
            filename = str(self.name_species) + ".pkl"
            file_path = os.path.join(os.getcwd(), filename)

        logger.info("Loading data from: %s" % file_path)

        try:
            self.data_full = pd.read_pickle(file_path)
            logger.info("Succesfully loaded previously saved data.")
            return self.data_full
        except IOError as e:
            logger.error("Problem loading data! %s " % str(e))

    def find_species_occurrences(self, name_species=None, **kwargs):
        raise NotImplementedError("You need to implement this method!")

    def get_data(self):
        """
        Returns the (pre)loaded species data.

        :rtype: geopandas.GeoDataFrame
        """
        return self.data_full

    def set_data(self, data_frame):
        """
        Careful, overwrites the existing raw data!. More documentation
        """
        self.data_full = data_frame

    def plot_species_occurrence(self, figsize=(16, 12), projection='merc', facecolor='crimson'):
        """
        Documentation pending on plotting
        """
        if not isinstance(self.data_full, GeoDataFrame):
            if not isinstance(self.data_full, pd.DataFrame):
                raise AttributeError("No data to save. Please load it first.")
            else:
                self.geometrize(dropna=True)

        # empty dataset
        if self.data_full.shape[0] == 0:
            logger.error("No data to plot.")
            return
        # dataset with one point (one dimensional) is also problematic for plotting, if no buffer around
        if self.data_full.shape[0] == 1 and isinstance(self.data_full.geometry.iat[0], shapely.geometry.Point):
            logger.error("Only one point in dataset.")
            return

        # Now we have a geometry column (GeoPandas instance). Could be filled with Point/Polygon...
        mm = Basemap(projection=projection, lat_0=50, lon_0=-100,
                     resolution='l', area_thresh=1000.0,
                     llcrnrlon=self.data_full.geometry.total_bounds[0],  # lower left corner longitude point
                     llcrnrlat=self.data_full.geometry.total_bounds[1],  # lower left corner latitude point
                     urcrnrlon=self.data_full.geometry.total_bounds[2],  # upper right longitude point
                     urcrnrlat=self.data_full.geometry.total_bounds[3]   # upper right latitude point
                     )

        # prepare longitude/latitude list for basemap
        ax1 = plt.subplots(figsize=figsize)[1]

        plt.title("%s occurrence records from %s " % (self.name_species, self.source.name))

        mm.drawcoastlines()
        mm.drawcountries()
        mm.drawrivers(color='lightskyblue', linewidth=1.5)
        mm.drawmapboundary(fill_color='lightskyblue')
        mm.fillcontinents(color='cornsilk')
        # draw latitude and longitude
        mm.drawmeridians(np.arange(-180, 180, 10), labels=[False, False, False, True])
        mm.drawparallels(np.arange(-180, 180, 10), labels=[True, True, False, False])

        patches = []
        selection = self.data_full
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


class GBIFSpecies(Species):
    """
    Some class-level documentation on GBIF species. One two three.
    """

    def __init__(self, **kwargs):

        Species.__init__(self, **kwargs)
        self.source = Source.GBIF
        self.observations_type = ObservationsType.PRESENCE_ONLY

    def find_species_occurrences(self, name_species=None, **kwargs):
        """
        Finds and loads species occurrence data into pandas DataFrame.
        Data comes from the GBIF database, based on name or gbif ID
        the occurrences.search(...) returns a list of json structures
        which we load into Pandas DataFrame for easier manipulation.

        """
        if name_species:
            self.name_species = name_species
        if not self.name_species:
            raise AttributeError("You have not provided a name for the species.")

        try:
            species_result = species.name_backbone(name=self.name_species, verbose=False)
            if species_result['matchType'] == 'NONE':
                logger.error("No match for the species %s in GBIF backbone!" % self.name_species)
                return pd.DataFrame()
                # raise ValueError("No match for the species %s " % self.name_species)
                # TODO: maybe just return, no error raising...
            self.ID = species_result['usageKey']
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)

        except AttributeError:   # name not provided, assume at least ID is provided
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)

        full_results = copy.copy(first_res)
        logger.info("Number of occurrences in GBIF backbone: %s " % full_results['count'])
        # http://lists.gbif.org/pipermail/api-users/2015-February/000135.html
        # maximum offset+limit allowed by the API
        if full_results['count'] > 200000:
            logger.warning("There are more than 200000 observations for %s" % self.name_species)
            logger.warning("The GBIF API has a limitation of 200000.")
            logger.warning("You may want to manually request a download from the website.")
            logger.warning("Will continue fetching 200000 results.")
            return

        # results are paginated so we need a loop to fetch them all
        # http://www.gbif.org/developer/occurrence
        # "This API provides services for searching occurrence records that have been indexed by GBIF. In order to retrieve all
        # results for a given search filter you need to issue individual requests for each page, which is limited to a maximum
        # size of 300 records per page. Note that for technical reasons we also have a hard limit for any query of 200,000
        # records. You will get an error if the offset + limit exceeds 200,000. To retrieve all records beyond 200,000 you should
        # use our asynchronous download service instead."
        counter = 1
        limit = 10000
        while first_res['endOfRecords'] is False and (300 * counter) + limit < 200000:
            first_res = occurrences.search(taxonKey=self.ID, offset=300 * counter, limit=limit)
            logger.debug("Page offset: %s, counter %s. Got %s more records ... " % ((300 * counter), counter, len(first_res['results'])))
            full_results['results'].extend(first_res['results'])
            counter += 1

        logger.info("Loading species ... ")
        logger.debug(full_results['count'] == len(full_results['results']))   # match?
        logger.debug("Full results: %s , got: %s " % (full_results['count'], len(full_results['results'])))
        # TODO: do we want a special way of loading? say, suggesting data types in some columns?

        # TODO: should we reformat the dtypes of the columns? at least day/month/year we care?
        # data_cleaned[['day', 'month', 'year']] = data_cleaned[['day', 'month', 'year']].fillna(0.0).astype(int)

        self.data_full = pd.DataFrame(full_results['results'])  # load results in pandas DF
        # self.data_full.columns = map(str.lower, self.data_full.columns)  # convert all column headers to lowercase
        self.data_full.columns = [x.lower() for x in self.data_full.columns]

        if self.data_full.empty:
            logger.info("Could not retrieve any occurrences!")
        else:
            logger.info("Loaded species: %s " % self.data_full['species'].unique())
        return self.data_full

    def load_csv(self, file_path):
        """
        Documentation pending on loading from csv file
        """

        logger.info("Loading data from: %s" % file_path)
        f = open(file_path, 'r')
        try:
            dialect = csv.Sniffer().sniff(f.read(10240))
            self.data_full = pd.read_csv(file_path, sep=dialect.delimiter)
            # self.data_full.columns = map(str.lower, self.data_full.columns)  # convert all column headers to lowercase
            self.data_full.columns = [x.lower() for x in self.data_full.columns]
            logger.info("Succesfully loaded previously saved CSV data.")
            if 'specieskey' in self.data_full and self.data_full['specieskey'].unique().size == 1:
                self.ID = self.data_full['specieskey'].unique()[0]
                logger.info("Updated species ID: %s " % self.ID)

            return self.data_full
        finally:
            f.close()

    # vectorize? better name
    def geometrize(self, dropna=True, longitude_col_name='decimallongitude', latitude_col_name='decimallatitude', crs=None):
        """
        Converts data to geopandas. The latitude/longitude is converted into shapely Point geometry.
        Geopandas/GeoSeries data structures have a geometry column.
        """
        try:
            if crs is None:
                crs = {'init': "EPSG:4326"}
            # exclude those points with NaN in coordinates
            if dropna:
                geometry = [Point(xy) for xy in zip(
                    self.data_full[longitude_col_name].dropna(),
                    self.data_full[latitude_col_name].dropna()
                )]
                self.data_full = GeoDataFrame(
                    self.data_full.dropna(
                        subset=[latitude_col_name, longitude_col_name]),
                    crs=crs,
                    geometry=geometry)
                logger.info("Data geometrized: converted into GeoPandas dataframe.")
                logger.info("Points with NaN coordinates ignored. ")
            else:
                geometry = [Point(xy) for xy in zip(
                    self.data_full[longitude_col_name],
                    self.data_full[latitude_col_name]
                )]
                self.data_full = GeoDataFrame(self.data_full, crs=crs, geometry=geometry)
                logger.info("Data geometrized: converted into GeoPandas dataframe.")
        except AttributeError:
            logger.error("No latitude/longitude data to convert into a geometry. Please load the data first.")

    def polygonize(self,
                   buffer_distance=1,
                   buffer_resolution=16,
                   simplify_tolerance=0.1,
                   preserve_topology=False,
                   with_envelope=False):
        """
        Expand each sample point into its polygon of influence (buffer).
        Merge the polygons that overlap into a cascaded union (multipolygon)
        Return a GeoDataFrame
        """
        if not (isinstance(self.data_full, GeoSeries) or isinstance(self.data_full, GeoDataFrame)):
            self.geometrize(dropna=True)

        data_polygonized = self.data_full.copy(deep=True)
        if with_envelope:
            data_polygonized = data_polygonized.buffer(
                buffer_distance,
                buffer_resolution
            ).simplify(simplify_tolerance, preserve_topology).envelope
            logger.debug("Data polygonized with envelope.")
        else:
            data_polygonized = data_polygonized.buffer(
                buffer_distance,
                buffer_resolution
            ).simplify(simplify_tolerance, preserve_topology)
            logger.debug("Data polygonized without envelope.")

        # TODO: use GeoSeries.cascaded_union directly? does this always return multipolygon? sometimes just one?
        # TODO if not self._usable_area.is_valid only then cascaded_union?
        cascaded_union_multipolygon = shapely.ops.cascaded_union(data_polygonized.geometry)
        logger.debug("Cascaded union of polygons created.")

        # no .tolist for MultiPolygon unfortunatelly
        # TODO: this can sometimes create a single polygon (non iterable)
        df_polygonized = GeoDataFrame(geometry=[pol for pol in cascaded_union_multipolygon])
        return df_polygonized

    def overlay(self, species_range_map):
        """
        Overlays the point records with a species range map. The map can be an instance of IUCNSpecies, or directly
        a GeoSeries datastructure containing geometry data.
        This overlaying effectively crops the point records to the area within the range map.
        """
        if not (isinstance(species_range_map, GeoSeries) or isinstance(species_range_map, IUCNSpecies)):
            raise AttributeError("Please provide a correct species rangemap input.")

        if not isinstance(self.data_full, GeoDataFrame):
            if not isinstance(self.data_full, pd.DataFrame):
                raise AttributeError("No data to save. Please load it first.")
            else:
                self.geometrize(dropna=True)

        try:
            if isinstance(species_range_map, GeoSeries):
                prepped_range_map_union = prep(species_range_map.unary_union)
            elif isinstance(species_range_map, IUCNSpecies):
                prepped_range_map_union = prep(species_range_map.data_full.geometry.unary_union)
            # self.data_full = self.data_full[self.data_full.geometry.intersects(range_map_union)]
            self.data_full = self.data_full[self.data_full.geometry.apply(lambda x: prepped_range_map_union.contains(x))]

            logger.info("Overlayed species occurrence data with the given range map.")
        except ValueError as e:
            logger.error("The rangemap geometries seem to be invalid (possible self-intersections). %s " % str(e))


class IUCNSpecies(Species):
    """
    The input data is in shapefiles(ESRI native format) and contains the known expert range of each species. Ranges are
    depicted as polygons. One large shapefile contains all the distribution maps of that group, i.e., all geometries. We
    will need to rasterize IUCN range polygons to grids with a predefined resolution. "Gridify" the data, and that per
    species rather than all in one region.

    """

    def __init__(self, **kwargs):
        Species.__init__(self, **kwargs)
        self.source = Source.IUCN
        self.observations_type = ObservationsType.PRESENCE_ONLY

    def load_shapefile(self, file_path):
        """
        A GeoDataFrame is a tablular data structure that contains a column called geometry which contains a GeoSeries.
        So data_full will be a geopandas dataframe, you can obtain it by .get_data()
        Supported input formats: Shapefiles, GeoJSON.
        """
        logger.info("Loading data from: %s" % file_path)
        self.data_full = GeoDataFrame.from_file(file_path)   # shapely.geometry type of objects are used
        # self.data_full.columns = map(str.lower, self.data_full.columns)   # convert all column headers to lowercase
        self.data_full.columns = [x.lower() for x in self.data_full.columns]   # python 2

        logger.info("The shapefile contains data on %d species areas." % self.data_full.shape[0])
        self.shape_file = file_path

    def find_species_occurrences(self, name_species=None, **kwargs):
        """
        Documentation pending on filtering species data from a shapefile
        """

        if not hasattr(self, 'data_full'):
            raise AttributeError("You have not loaded the data.")
        if name_species:
            self.name_species = name_species
        if not self.name_species:
            raise AttributeError("You have not provided a name for the species.")

        all_data = self.data_full[self.data_full['binomial'] == self.name_species]

        if all_data.shape[0] == 0:
            raise ValueError("There is no species with the name '%s' in the shapefile"
                             % self.name_species)
        else:
            self.data_full = all_data
            logger.info("Loaded species: %s " % self.data_full['binomial'].unique())

        if self.data_full['id_no'].shape[0] == 1:
            self.ID = int(self.data_full['id_no'].iloc[0])

        return self.data_full

    def save_shapefile(self, full_name=None, driver='ESRI Shapefile', overwrite=False):
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

    def rasterize_old(self, raster_file=None, pixel_size=None, x_res=None, y_res=None, *args, **kwargs):
        # options = ["ALL_TOUCHED=TRUE"]
        # right now it is pixel_size. But we could complicate further with cell width/height.
        # TODO: no need for x_res y_res?

        if not (pixel_size or raster_file):
            raise AttributeError("Please provide pixel_size and a target raster_file.")

        NoData_value = -9999

        # Open the data source and read in the extent
        # TODO: check shapefile exists
        source_ds = ogr.Open(self.shape_file)
        source_layer = source_ds.GetLayer()
        x_min, x_max, y_min, y_max = source_layer.GetExtent()   # boundaries

        # Create the destination data source
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        target_ds = gdal.GetDriverByName('GTiff').Create(raster_file, x_res, y_res, 1, gdal.GDT_Byte)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(NoData_value)
        # Rasterize
        gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[255], *args, **kwargs)

        # WOW: https://trac.osgeo.org/gdal/wiki/PythonGotchas "To save and close GDAL raster datasets or OGR vector
        # datasources, the object needs to be dereferenced, such as setting it to None, a different value, or deleting
        # the object. "
        # From the book "Python Geospatial Analysis Essentials":
        # dstFile.Destroy() This closes the destination file and makes sure everything
        # has been saved to disk.
        band = None
        target_ds = None
        logger.info("Data rasterized into file %s " % raster_file)
        logger.info("Resolution: x_res={0} y_res={1}".format(x_res, y_res))

        self.raster_file = raster_file
        self.x_res = x_res
        self.y_res = y_res

    def rasterize(self, raster_file=None, pixel_size=None, all_touched=False,
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
                                    efault_value=default_value
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

    # def load_raster_data(self, raster_file=None):
    #     """
    #     Documentation pending on how to load raster data
    #     """
    #     if raster_file:
    #         self.raster_file = raster_file
    #     if not self.raster_file:
    #         raise AttributeError("Please rasterize the data first, ",
    #                              "or provide a raster_file to read from.")

    #     with rasterio.open(self.raster_file) as src:
    #         logger.info("Loaded raster data from %s " % self.raster_file)
    #         logger.info("Driver name: %s " % src.driver)
    #         logger.info("Resolution: x_res={0} y_res={1}.".format(src.width, src.height))
    #         logger.info("Coordinate reference system: %s " % src.crs)
    #         logger.info("Affine transformation: %s " % (src.affine.to_gdal(),))
    #         logger.info("Number of layers: %s " % src.count)
    #         self.raster_affine = src.affine
    #         return src.read()

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

    def pixel_to_world_coordinates(self,
                                   raster_data=None,
                                   no_data_value=0,
                                   filter_no_data_value=True,
                                   band_number=1):
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
                raster_data = self.load_raster_data(self.raster_file).read(band_number)
                logger.info("Succesfully loaded existing raster data from %s." % self.raster_file)
            except AttributeError as e:
                logger.error("Could not open raster file. %s " % str(e))

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
            coordinates = (T1 * np.where(raster_data > no_data_value))
        else:
            logger.info("Not filtering any no_data pixels.")
            coordinates = (T1 * np.where(np.ones_like(raster_data)))

        return coordinates

    def random_pseudo_absence_points(self,
                                     buffer_distance=2,
                                     buffer_resolution=16,
                                     simplify_tolerance=1,
                                     preserve_topology=True,
                                     fast=False,
                                     count=100):
        """
        Draw random pseudo-absence points from within a buffer around the geometry. First it simplifies the geometry
        with a buffer around the original geometry. Then calculates the difference between this one, and the original
        geometry, to determine a geometry from which to sample random points. Finally, generates random points one by
        one and tests if they fall in that difference-geometry, until a <count> number of points are generated.
        If the "buffered" geometry is invalid (which could happen), it gradually tries to simplify it by applying a
        bigger value for the simplify_tolerance parameter, until the geometry becomes valid. The reason is that an
        operation like difference/intersection is problematic to apply on an invalid geometry.
        The value is increased by maximum of 100.

        A more efficient approach would be to just generate a <count> number of points from the first step, i.e.,
        from the buffer. Some points will fall within the original shape, and they can be discarded,
        so the number of pseudo-absence points will not actually be equal to <count>.
        If precision is not an issue, we could provide a <count> number that is larger but calculated according
        to the original_area/buffered_convex_hull ratio.
        update: Maybe not even necessary, given that shapely's prep(..) speeds up a factor of 100 to 1000
        """
        # First simplify (necessary) and apply a buffer around the geometry
        simplified_buffer = (self.data_full.geometry.buffer(0)
                             .simplify(simplify_tolerance, preserve_topology)
                             .buffer(buffer_distance, buffer_resolution))

        # try to "fix" an invalid geometry by gradually simplifying
        simplify_tolerance_inc = simplify_tolerance

        # assume that there is one geometry record in the dataframe (single species)
        while not simplified_buffer.iat[0].is_valid and simplify_tolerance_inc < 100:
            logger.info("Buffered geometry is invalid, trying to simplify it using tolerance = %s" % simplify_tolerance_inc)
            simplified_buffer = simplified_buffer.simplify(simplify_tolerance_inc, preserve_topology)
            simplify_tolerance_inc += 1

        # Get the difference between the original geometry and the one above, to determine a geometry
        # from which to draw random points.
        logger.info("Buffered geometry valid. Now creating a buffer difference from which to draw random points")
        buffer_difference = simplified_buffer.geometry.difference(self.data_full.geometry.buffer(0))
        xmin, ymin, xmax, ymax = buffer_difference.total_bounds
        random_count = 0
        pts = []

        # draw pounts until you reach a <count> number of points
        logger.info("Creating random points ... ")
        prepped_buffer_difference = prep(buffer_difference.iat[0])  # immense speedup with prepared geometries!

        while random_count < count:
            random_point = Point(xmin + (xmax - xmin) * np.random.random(), ymin + (ymax - ymin) * np.random.random())
            if prepped_buffer_difference.contains(random_point):
                pts.append(random_point)
                random_count += 1
            else:
                pass

        geo_pts = GeoSeries(pts)
        self.pseudo_absence_points = geo_pts
        return self.pseudo_absence_points

    def drop_extinct_species(self, presence_column_name='presence', discard_bad=False):
        """
        According to the current IUCN Coded Domain Values for Presence:
        Code    Presence
        1       Extant
        2       Probably Extant (discontinued)
        3       Possibly Extant
        4       Possibly Extinct
        5       Extinct (post 1500)
        6       Presence Uncertain

        Species can have both areas(polygons) in which they are extinct(5) AND areas in which they are not.
        We keep such species, and only filter-out species for which all areas are extinct.

        :param str presence_column_name: The column name which contains the presence code values. Default is 'presence'.

        :param bool discard_bad: Whether to keep or discard species with "unknown only" areas (code==0). By default they
        are kept (discard_bad=False).
        There are currently (july 2016) four such problematic species:
        Acipenser baerii, Ambassis urotaenia, Microphysogobio tungtingensis, Rhodeus sericeus

        :returns: None

        """
        logger.info("There are currently %s unique species. \n" % self.data_full.binomial.unique().size)

        # the following creates a pandas series in the form:
        # Aaptosyax grypus                                     [2.0]
        # Abbottina binhi                                      [2.0]
        # Aborichthys elongatus                 [1.0, 1.0, 1.0, 1.0]
        # Aborichthys garoensis                      [1.0, 1.0, 1.0]
        # Aborichthys kempi           [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        # Aborichthys tikaderi                                 [1.0]
        # Abramis brama                                   [1.0, 1.0]
        # Acantharchus pomotis                       [5.0, 1.0, 1.0]
        # Acanthobrama centisquama                        [1.0, 1.0]
        # Acanthobrama lissneri                           [1.0, 1.0]
        # ...
        # which means species are grouped by binomial, and the presence columns are aggregated in a list
        grouped = self.data_full.groupby('binomial')[presence_column_name].apply(lambda x: x.tolist())

        extinct = []
        for row in grouped.iteritems():
            if row[1] == [5.0]:  # all areas are extinct
                extinct.append(row[0])
            if discard_bad:
                if row[0] == [0.0]:  # all areas are with invalid value of 0
                    extinct.append(row[0])

        logger.info("Filtering out the following extinct species: %s \n" % extinct)
        self.data_full = self.data_full[~self.data_full.binomial.isin(extinct)]
        logger.info("There are now %s unique species after dropping out extinct ones." % self.data_full.binomial.unique().size)


class MOLSpecies(Species):
    pass
