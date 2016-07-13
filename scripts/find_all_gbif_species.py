#!/usr/bin/env python

from iSDM.species import GBIFSpecies
import logging
import timeit

method = "msgpack"
size = "large"
logger = logging.getLogger('iSDM.species')
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler('./data/fish/selection/test/again/' + size + '/to_' + method + '.log')
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to thes logger
logger.addHandler(fh)
logger.addHandler(ch)

# logger.info("Loading binomials list from pickled file...")
non_extinct_binomials = pickle.load(open("./data/fish/selection/non_extinct_binomials.pkl", "rb"))
# non_extinct_binomials = ["Salmo trutta", "Anguilla anguilla", "Abramis brama", "Esox lucius", "Gasterosteus aculeatus", "Gobio gobio", "Salmo salar", "Cyprinus carpio", "Barbatula barbatula", "Barbus intermedius", "Scardinius erythrophthalmus", "Platichthys flesus", "Phoxinus phoxinus", "Aborichthys tikaderi", "Aborichthys kempi", "Aborichthys garoensis", "Abbottina binhi"]
# logger.info("Done with loading binomials list.")
# non_extinct_binomials = ["Esox americanus", "Etheostoma blennioides", "Etheostoma caeruleum", "Misgurnus anguillicaudatus",  "Notropis texanus", "Percina nigrofasciata", "Plecoglossus altivelis", "Rhinichthys atratulus", "Rhinichthys cataractae", "Thymallus thymallus"]

# small
total_time = 0
non_extinct_binomials = []
# if size == "small":
#     non_extinct_binomials = ["Cottus princeps", "Mastacembelus erythrotaenia", "Petrocephalus pellegrini", "Polypterus delhezi", "Polypterus retropinnis"]
# elif size == "medium":
#     non_extinct_binomials = ["Trinectes maculatus", "Coregonus lavaretus", "Lepomis microlophus", "Etheostoma olmstedi", "Moxostoma erythrurum"]
# elif size == "large":
#     non_extinct_binomials = ["Salmo trutta", "Anguilla anguilla", "Abramis brama", "Esox lucius", "Gasterosteus aculeatus"]

logger.info("size = %s species." % len(non_extinct_binomials))
for idx, name in enumerate(non_extinct_binomials):
    next_species = GBIFSpecies(name_species=name)
    logger.info("Attempting to query for species: %s " % name)
    logger.info("Species #: %s " % idx)
    if "\'" in name:
        logger.info("Skipping problematic name: %s " % name)
        continue

    if next_species.find_species_occurrences().empty:
        logger.info("No data to save on %s " % name)
    else:
        logger.info("Done with GBIF API query, now saving.")
        start_time = timeit.default_timer()
        if method == "msgpack":
            next_species.save_data(dir_name="./data/fish/selection/test/again/" + size + "/msgpack", method=method)
        elif method == "pickle":
            next_species.save_data(dir_name="./data/fish/selection/test/again/" + size + "/pickle", method=method)
        stop_time = timeit.default_timer()
        logger.info("Serializing with %s took %s seconds." % (method, stop_time - start_time))
        total_time += stop_time - start_time

logger.info("Total serializing time %s seconds for %s. " % (total_time, method))
