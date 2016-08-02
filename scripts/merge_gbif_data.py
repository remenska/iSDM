#!/usr/bin/env python
"""
script: merge_gbif_data.py
Description: Merge all individual files containing occurrence records per species (obtained with the script find_all_gbif_species.py)
             into one large file. Species occurrences can contain various metadata (columns), but the merged file will only
             contain the relevant columns, i.e., those deemed important for further analysis. The rest of the metadata is ignored.
             For example, the latitude and longitude (if available) for each occurrence record are relevant.
Input:
 - method of deserialization (msgpack or pickle, same as the one used for serializing the individual files)
 - folder where the separate files are stored
 - a list of important columns (loaded from a file)

This script does the following:
 1. Loops through all files in the folder where occurrences for each individual species are stored, and deserializes (reads)
    the data in a temporary dataframe.
 2. Appends the deserialized data in one large file, selecting only the important columns (if available).
"""

import logging
import pickle
import os
import pandas as pd
import timeit

# input
method = "pickle"  # this method should be the same as the one used to serialize the GBIF occurrence data in files.
saved_data_path = './data/fish/selection/'
important_columns_file_path = "./data/fish/selection/important_columns.pkl"

# logging
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('./data/fish/selection/' + '/read_' + method + '.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)


important_columns = pickle.load(open(important_columns_file_path, "rb"))
no_occurrences = []
saved_data_path = os.path.join(saved_data_path, method)
total_time = 0

# Loop through each file in the folder
for idx, my_file in enumerate(os.listdir(saved_data_path)):
    logger.info("Opening %s " % my_file)
    logger.info("Species #: %s " % idx)
    start_time = timeit.default_timer()
    # Open and read data into pandas dataframe
    if method == "msgpack":
        df = pd.read_msgpack(os.path.join(saved_data_path, my_file))
    elif method == "pickle":
        df = pd.read_pickle(os.path.join(saved_data_path, my_file))
    stop_time = timeit.default_timer()
    logger.info("Reading with %s took %s seconds." % (method, stop_time - start_time))
    total_time += stop_time - start_time
    # If the dataframe is not empty
    if not df.empty:
        # Make an intersection from the list of important columns and those present in the dataframe
        common_columns = list(important_columns.intersection(set(df.columns.tolist())))
        # Append the data from those columns of the intersection, into one big "merged" file.
        # (beter save the merged file one directory "up" to avoid merging with itself)
        # (memory usage doesn't grow constantly, individual frames are dumped on disk and memory discarded in each loop)
        df[common_columns].to_msgpack(os.path.join(saved_data_path, "../merged.msg"), append=True)
        logger.info("Columns: %s" % len(common_columns))
        logger.info(df.info(memory_usage='deep'))
    elif df.empty:
        # also keep a list of species with no occurrences available.
        no_occurrences.append(my_file)
    del df


logger.info("Total reading time %s seconds for %s. " % (total_time, method))
# also store a list of all species with 0 occurrences in a file.
pickle.dump(no_occurrences, open(os.path.join(saved_data_path, "no_occurrences.pkl"), "wb"))
logger.info("Stored a list of species with zero occurrences in %s " % os.path.join(saved_data_path, "no_occurrences.pkl"))
