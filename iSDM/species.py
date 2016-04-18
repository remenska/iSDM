from pygbif import species # http://pygbif.readthedocs.org/en/latest/
from pygbif import occurrences
import copy
import pandas as pd
import logging
import os
import csv

try:
    import cPickle as pickle
except:
    import pickle

logger = logging.getLogger('iSDM.species')
logger.setLevel(logging.DEBUG)

class Species(object):

    def __init__(self, **kwargs):

        if 'ID' not in kwargs and 'name_species' not in kwargs:
            raise ValueError("Cannot initialize species without a 'species_name' or an 'ID' argument supplied")

        if 'name_species' in kwargs:
            self.name_species=kwargs['name_species']
        if 'ID' in kwargs:
            self.ID=kwargs['ID']


    def save_data(self, dirname=None, filename=None):
        """
        Serializes the loaded GBIF species occurrence filtered dataset (pandas.DataFrame) into a binary pickle file
        """

        if filename is None:
            filename = str(self.ID) + ".pkl"
        if dirname is None:
            dirname = os.getcwd()

        fullname = os.path.join(dirname, filename)

        f = open(fullname, 'wb')
        try:
            pickle.dump(self.data_full, f)
            logger.debug("Saved data: %s " %fullname)
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


class GBIFSpecies(Species):

    def find_species_occurrences(self, **kwargs):
        """
        Finds and loads species occurrence data into pandas DataFrame.
        Data comes from the GBIF database, based on name or gbif ID

        """

        try:
            species_result = species.name_backbone(name=self.name_species, verbose=False)
            if species_result['matchType']=='NONE':
                raise ValueError("No match for the species %s " % self.name_species)
            self.ID = species_result['usageKey']
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)

        except AttributeError: # name not provided, assume at least ID is provided
            first_res = occurrences.search(taxonKey=self.ID, limit=100000, **kwargs)
        
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
    Data are held in shapefiles, the ESRI native format. Ranges are depicted as polygons.
    The maps are available as shapefiles. Not as one layer per species, but one (very) large shapefile 
    that contains all the distribution maps of that group.
    We will need to rasterize IUCN range polygons to grids with a predefined resolution. 
    """
    

class MOLSpecies(Species):
    pass




