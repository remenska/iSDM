from pygbif import species # http://pygbif.readthedocs.org/en/latest/
from pygbif import occurrences
import copy
import pandas as pd
import logging
logger = logging.getLogger('iSDM')
logger.setLevel(logging.DEBUG)

class Species(object):
	# Let's generalize to list of species

	def __init__(self, **kwargs):

		if 'ID' not in kwargs and 'name_species' not in kwargs:
			raise ValueError("Cannot initialize species without a 'species_name' or an 'ID' argument supplied")

		if 'name_species' in kwargs:
			self.name_species=kwargs['name_species']
		if 'ID' in kwargs:
			self.ID=kwargs['ID']


	def load_species_occurrences(self, **kwargs):
		#self.name_species = name_species
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

		self.df_full = pd.DataFrame(full_results['results']) # load results in pandas dataframes
		if self.df_full.empty:
			logger.info("Could not retrieve any occurrences!")
		else:	
			logger.info("Loaded species: %s " % self.df_full['species'].unique())
		return self.df_full

