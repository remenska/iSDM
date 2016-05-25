"""
Some interesting text.

:synopsis: A useful module indeed.

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


logger = logging.getLogger('iSDM.species')
logger.setLevel(logging.DEBUG)


class Source(Enum):
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
    py:class: Species
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

    def save_data(self, full_name=None, dir_name=None, file_name=None):
        """
        Serializes the loaded species occurrence filtered dataset (`pandas <http://pandas.pydata.org/pandas-docs/stable/dsintro.html>`_ or `geopandas <http://geopandas.org/user.html>`_ DataFrame) into a binary `pickle <https://en.wikipedia.org/wiki/Pickle_%28Python%29>`_  file.

       :param str full_name: The full path of the file (including the directory and name in one string),
        where the data will be saved.
       :param str dir_name: The directory where the file will be stored. If :attr:`file_name` is not specified, the default one :attr:`name_species` + :attr:`ID`.pkl is given.
       :param str file_name: The name of the file where the data will be saved. If :attr:`dir_name` is not specified, the current working directory is taken by default.
       :raises AttributeError: if the data has not been loaded in the object before. See :func:`load_data` and :func:`find_species_occurrences`

        """
        if full_name is None:
            if file_name is None:
                file_name = str(self.name_species) + str(self.ID) + ".pkl"
            if dir_name is None:
                dir_name = os.getcwd()

            full_name = os.path.join(dir_name, file_name)

        try:
            self.data_full.to_pickle(full_name)
            logger.debug("Saved data: %s " % full_name)
            logger.debug("Type of data: %s " % type(self.data_full))
        except IOError as e:
            logger.error("Could not save data! %s " % str(e))
        except AttributeError as e:
            logger.error("No data to save. Please load it first. %s " % str(e))

    def load_data(self, file_path=None):
        """ Loads the serialized species pickle file into a pandas DataFrame.
    """
        if file_path is None:
            filename = str(self.name_species) + str(self.ID) + ".pkl"
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
        """This function does something.

    Args:
       name (str):  The name to use.

    Kwargs:
       state (bool): Current state to be in.RrR

    Returns:
       int.  The return code::

          0 -- Success!
          1 -- No good.
          2 -- Try again.

    Raises:
       AttributeError, KeyError
        """
        return self.data_full

    def set_data(self, data_frame):
        # !!! Careful, overwrites the existing raw data!
        self.data_full = data_frame

    def plot_species_occurrence(self, figsize=(16, 12), projection='merc', facecolor='crimson'):
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
        fig, ax1 = plt.subplots(figsize=figsize)

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
                patches.append(PolygonPatch(Point(mm(poly.x, poly.y)).buffer(9999)))
            else:
                logger.warning("Geometry type %s not supported. Skipping ... " % poly.geom_type)
                continue
        ax1.add_collection(PatchCollection(patches, facecolor=facecolor, match_original=True, zorder=100))
        plt.show()


class GBIFSpecies(Species):

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
                raise ValueError("No match for the species %s " % self.name_species)
            self.ID = species_result['usageKey']
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)

        except AttributeError:   # name not provided, assume at least ID is provided
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)

        full_results = copy.copy(first_res)

        # results are paginated so we need a loop to fetch them all
        counter = 1
        while first_res['endOfRecords'] is False:
            first_res = occurrences.search(taxonKey=self.ID, offset=300 * counter, limit=10000)
            full_results['results'].extend(first_res['results'])
            counter += 1

        logger.info("Loading species ... ")
        logger.info("Number of occurrences: %s " % full_results['count'])
        logger.debug(full_results['count'] == len(full_results['results']))   # match?

        # TODO: do we want a special way of loading? say, suggesting data types in some columns?

        # TODO: should we reformat the dtypes of the columns? at least day/month/year we care?
        # data_cleaned[['day', 'month', 'year']] = data_cleaned[['day', 'month', 'year']].fillna(0.0).astype(int)

        self.data_full = pd.DataFrame(full_results['results'])  # load results in pandas DF
        self.data_full.columns = map(str.lower, self.data_full.columns)  # convert all column headers to lowercase

        if self.data_full.empty:
            logger.info("Could not retrieve any occurrences!")
        else:
            logger.info("Loaded species: %s " % self.data_full['species'].unique())
        return self.data_full

    def load_csv(self, file_path):

        logger.info("Loading data from: %s" % file_path)
        f = open(file_path, 'r')
        try:
            dialect = csv.Sniffer().sniff(f.read(10240))
            self.data_full = pd.read_csv(file_path, sep=dialect.delimiter)
            self.data_full.columns = map(str.lower, self.data_full.columns)  # convert all column headers to lowercase
            logger.info("Succesfully loaded previously saved CSV data.")
            if 'specieskey' in self.data_full and self.data_full['specieskey'].unique().size == 1:
                self.ID = self.data_full['specieskey'].unique()[0]
                logger.info("Updated species ID: %s " % self.ID)

            return self.data_full
        finally:
            f.close()

    # vectorize? better name
    def geometrize(self, dropna=True, longitude_col_name='decimallongitude', latitude_col_name='decimallatitude'):
        """
        Converts data to geopandas. The latitude/longitude is converted into shapely Point geometry.
        Geopandas/GeoSeries data structures have a geometry column.
        """
        try:
            crs = None
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
                logger.info("Points with NaN coordinnates ignored. ")
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

        cascaded_union_multipolygon = shapely.ops.cascaded_union(data_polygonized.geometry)
        logger.debug("Cascaded union of polygons created.")

        # no .tolist for MultiPolygon unfortunatelly
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
                range_map_union = species_range_map.unary_union
            elif isinstance(species_range_map, IUCNSpecies):
                range_map_union = species_range_map.data_full.geometry.unary_union
            self.data_full = self.data_full[self.data_full.geometry.intersects(range_map_union)]
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
        self.data_full.columns = map(str.lower, self.data_full.columns)   # convert all column headers to lowercase
        logger.info("The shapefile contains data on %d species." % self.data_full.shape[0])
        self.shape_file = file_path

    def find_species_occurrences(self, name_species=None, **kwargs):

        if not self.shape_file:
            raise AttributeError("You have not provided a shapefile to load data from.")
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

    def rasterize(self, raster_file=None, pixel_size=None, x_res=None, y_res=None, *args, **kwargs):
        # options = ["ALL_TOUCHED=TRUE"]
        # right now it is pixel_size. But we could complicate further with cell width/height.

        if not (pixel_size or raster_file):
            raise AttributeError("Please provide pixel_size and a target raster_file.")

        NoData_value = -9999

        # Open the data source and read in the extent
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

    def load_raster_data(self, raster_file=None):
        if raster_file:
            self.raster_file = raster_file
        if not self.raster_file:
            raise AttributeError("Please rasterize the data first, ",
                                 "or provide a raster_file to read from.")

        geo = gdal.Open(self.raster_file)
        # sillly gdal python wrappings don't throw exceptions
        if not geo:
            logger.error("Unable to open %s " % raster_file)
            return None

        drv = geo.GetDriver()

        logger.info("Driver name: %s " % drv.GetMetadataItem('DMD_LONGNAME'))
        logger.info("Raster data from %s loaded." % self.raster_file)
        logger.info("Resolution: x_res={0} y_res={1}. GeoTransform: {2}"
                    .format(geo.RasterXSize, geo.RasterYSize, geo.GetGeoTransform()))
        img = geo.ReadAsArray()

        return img


class MOLSpecies(Species):
    pass
