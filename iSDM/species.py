from pygbif import species # http://pygbif.readthedocs.org/en/latest/
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

try:
    import cPickle as pickle
except:
    import pickle

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


#TODO: split into train/test subset (some random manner?) for model evaluation
class Species(object):

    def __init__(self, **kwargs):

        if 'ID' not in kwargs and 'name_species' not in kwargs:
            raise ValueError("Cannot initialize species without a 'species_name' or an 'ID' argument supplied")

        if 'name_species' in kwargs:
            self.name_species=kwargs['name_species']
        if 'ID' in kwargs:
            self.ID=kwargs['ID']

    def save_data(self, full_name=None, dir_name=None, file_name=None): 
        """
        Serializes the loaded GBIF species occurrence filtered dataset (pandas.DataFrame) into a binary pickle file
        """
        #TODO shall we store in HDF5 / PyTables also?
        
        if full_name is None:
            if file_name is None:
                file_name = str(self.ID) + ".pkl"
            if dir_name is None:
                dir_name = os.getcwd()

            full_name = os.path.join(dir_name, file_name)

        f = open(full_name, 'wb')
        try:
            pickle.dump(self.data_full, f)
            logger.debug("Saved data: %s " %full_name)
        except AttributeError as e:
            logger.error("Occurrences data is not loaded! %s " %str(e))
        finally:
            f.close()

    def load_data(self, file_path=None):
        """
        Loads the serialized GBIF species occurrence filtered dataset into a pandas DataFrame
        """

        if file_path is None:
            filename = str(self.ID) + ".pkl"
            file_path = os.path.join(os.getcwd(), filename)

        logger.info("Loading data from: %s" %file_path)
        f = open(file_path, 'rb')
        try:
            self.data_full = pickle.load(f)
            logger.info("Succesfully loaded previously saved data.")
            return self.data_full
        finally:
            f.close()

    def find_species_occurrences(self, **kwargs):
        raise NotImplementedError("You need to implement this method!")

    def get_data(self):
        return self.data_full

    def set_data(self, data_frame):
        # !!! Careful, overwrites the existing raw data!
        self.data_full = data_frame

    def plot_species_occurrence(self, figsize=(16,12), projection='merc'):
        data_clean = self.data_full.dropna(how='any', subset=['decimalLatitude', 'decimalLongitude'])

        # latitude/longitude lists
        data_full_latitude = data_clean.decimalLatitude
        data_full_longitude = data_clean.decimalLongitude

        plt.figure(figsize=figsize)
        plt.title("%s occurrence records from %s " 
            % (self.name_species, self.source.name)
            )

        my_map = Basemap(projection=projection, lat_0=50, lon_0=-100,
                        resolution='l', area_thresh=1000.0, 
                        llcrnrlon=data_full_longitude.min(),# lower left corner longitude point 
                        llcrnrlat=data_full_latitude.min(), # lower left corner latitude point
                        urcrnrlon=data_full_longitude.max(), # upper right longitude point
                        urcrnrlat=data_full_latitude.max() # upper right latitude point
                        )
        # prepare longitude/latitude list for basemap
        df_x, df_y = my_map(data_full_longitude.tolist(), data_full_latitude.tolist())
        my_map.drawcoastlines()
        my_map.drawcountries()
        my_map.drawmapboundary(fill_color='#649eff')
        my_map.fillcontinents(color='#cc9955')
        # draw latitude and longitude
        my_map.drawmeridians(np.arange(0, 360, 30))
        my_map.drawparallels(np.arange(-90, 90, 30))
        my_map.plot(df_x, df_y, 'bo', markersize=5, color="#b01a1a")


class GBIFSpecies(Species):

    def __init__(self, **kwargs):
        Species.__init__(self, **kwargs)
        self.source = Source.GBIF
        self.observations_type = ObservationsType.PRESENCE_ONLY


    def find_species_occurrences(self, **kwargs):
        """
        Finds and loads species occurrence data into pandas DataFrame.
        Data comes from the GBIF database, based on name or gbif ID
        the occurrences.search(...) returns a list of json structures
        which we load into Pandas DataFrame for easier manipulation.

        """

        try:
            species_result = species.name_backbone(name=self.name_species, verbose=False)
            if species_result['matchType']=='NONE':
                raise ValueError("No match for the species %s " % self.name_species)
            self.ID = species_result['usageKey']
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)

        except AttributeError: # name not provided, assume at least ID is provided
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)
        
        #TODO: more efficient way than copying...appending to the same dataframe?
        
        full_results = copy.copy(first_res)

        # results are paginated so we need a loop to fetch them all
        counter = 1
        while first_res['endOfRecords'] is False:
            first_res = occurrences.search(taxonKey=self.ID, offset=300*counter, limit=10000)
            full_results['results'] = copy.copy(full_results['results']) + copy.copy(first_res['results'])
            counter+=1
        
        logger.info("Loading species ... ")
        logger.info("Number of occurrences: %s " % full_results['count'])
        logger.debug(full_results['count'] == len(full_results['results'])) # match?

        #TODO: do we want a special way of loading? say, suggesting data types in some columns?

        #TODO: should we reformat the dtypes of the columns? at least day/month/year we care?
        #data_cleaned[['day', 'month', 'year']] = data_cleaned[['day', 'month', 'year']].fillna(0.0).astype(int)
        
        self.data_full = pd.DataFrame(full_results['results']) # load results in pandas dataframes
        if self.data_full.empty:
            logger.info("Could not retrieve any occurrences!")
        else:   
            logger.info("Loaded species: %s " % self.data_full['species'].unique())
        return self.data_full


    def load_csv(self, file_path):
        logger.info("Loading data from: %s" %file_path)
        f = open(file_path, 'r')
        try:
            dialect = csv.Sniffer().sniff(f.read(10240))
            self.data_full = pd.read_csv(file_path, sep=dialect.delimiter)
            logger.info("Succesfully loaded previously CSV data.")
            if 'specieskey' in self.data_full and self.data_full['specieskey'].unique().size == 1:
                self.ID = self.data_full['specieskey'].unique()[0]
                logger.info("Updated species ID: %s " %self.ID)

            return self.data_full
        finally:
            f.close()


class IUCNSpecies(Species):
    """
    Data are held in shapefiles, the ESRI native format and contains the known range of each species. Ranges are depicted as polygons.
    The maps are available as shapefiles. Not as one layer per species, but one (very) large shapefile 
    that contains all the distribution maps of that group.
    We will need to rasterize IUCN range polygons to grids with a predefined resolution. ("Gridify" the data, and per species rather than all in one region)

    http://www.petrkeil.com/?p=648
    "   It does not matter that much what is (or will be in future) the pattern of species distribution at a single scale. It does not matter that much what predicts species occurrences at the finest grain resolutions. It is naive to say that, for applied and conservation purposes, it is fine grain that matters. All grains 
        (potentially) matter. Appropriate scaling relationships are the baseline that links all of the grains together."

    """
    def __init__(self):
        Species.__init__(self)
        self.source=Source.IUCN


class MOLSpecies(Species):
    pass






