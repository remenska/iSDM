#!/usr/bin/env python

from iSDM.species import GBIFSpecies
import logging
import pickle

logger = logging.getLogger('iSDM.species')
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler('./find_all_gbif_species3.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.info("Loading binomials list from pickled file...")
non_extinct_binomials = pickle.load(open("./data/fish/selection/non_extinct_binomials.pkl", "rb"))
logger.info("Done with loading binomials list.")
logger.info("size = %s species." % len(non_extinct_binomials))

for name in non_extinct_binomials[275:1000]:
    next_species = GBIFSpecies(name_species=name)
    logger.info("Attempting to query for species: %s " % name)
    if "\'" in name:
        logger.info("Skipping problematic name: %s " % name)
        continue
    next_species.find_species_occurrences()
    logger.info("Done with GBIF API query, now saving.")
    next_species.save_data(dir_name="./data/fish/selection/gbif/")
