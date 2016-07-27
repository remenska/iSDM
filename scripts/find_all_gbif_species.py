#!/usr/bin/env python
"""
script: find_all_gbif_species.py
Description: Scrap all GBIF backend occurrence records, for a (filtered) list of species from the IUCN dataset,
             and store the data locally in separate files, per species.
Input:
 - list of species (binomials, loaded from a file)
 - method of serialization (msgpack or pickle)
 - folder where to save the data in separate files

This script does the following:
 1. Loads a list of non-extinct binomials (non_extinct_binomials.pkl).
    This list contains all the non-extinct species already selected from the IUCN expert range data,
    by filtering out species which have only extinct regions and nothing more.
 2. Loops through the list, and queries (using the GBIF API) the GBIF backend for occurrences for each individual binomial.
 3. Saves (serializes) the result of the query, as occurrences records, in a separate file for each binomial (species).
    Each individual file contains GBIF occurrence records for a particular species. All files are stored in the same folder.
    By default the name of the file corresponds to the name of the species.
"""

from iSDM.species import GBIFSpecies
import logging
import timeit
import pickle
import os

# input
method = "msgpack"  # could also be set to "pickle", which was default but msgpack is slightly better in speed/memory.
non_extinct_binomials_file_path = "./data/fish/selection/non_extinct_binomials.pkl"
save_data_path = "./data/fish/selection"

# logging
logger = logging.getLogger('iSDM.species')
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('./data/fish/selection/' + '/to_' + method + '.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)

# Load the list of non-extinct binomials that is previously filtered from the IUCN expert range data.
non_extinct_binomials = pickle.load(open(non_extinct_binomials_file_path, "rb"))

total_time = 0
save_data_path = os.path.join(save_data_path, method)

logger.info("size = %s species." % len(non_extinct_binomials))
# Loop through the list of non-extinct binomials
for idx, name in enumerate(non_extinct_binomials):
    next_species = GBIFSpecies(name_species=name)
    logger.info("Attempting to query for species: %s " % name)
    logger.info("Species #: %s " % idx)
    if "\'" in name:
        logger.info("Skipping problematic name: %s " % name)
        continue
    # Query the GBIF backend, and if empty just skip
    if next_species.find_species_occurrences().empty:
        logger.info("No data to save on %s " % name)
    # If there are occcurrence records in GBIF...
    else:
        logger.info("Done with GBIF API query, now saving.")
        start_time = timeit.default_timer()
        # 3. Save the data
        next_species.save_data(dir_name=save_data_path, method=method)
        stop_time = timeit.default_timer()
        logger.info("Serializing with %s took %s seconds." % (method, stop_time - start_time))
        total_time += stop_time - start_time

logger.info("Total serializing time %s seconds for %s. " % (total_time, method))
