import logging
import pickle
import os
import pandas as pd

logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler('./merge_gbif_records.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

path = './data/fish/selection/'
test_dir = os.path.join(path, 'gbif')

df_main = pd.DataFrame()
no_occurrences = []
first = True
for idx, my_file in enumerate(os.listdir(test_dir)):
    logger.info("Opening %s " % my_file)
    logger.info("Species #: %s " % idx)
    df = pd.read_pickle(os.path.join(test_dir, my_file))
    if not df.empty:
        df_main = df_main.append(df, ignore_index=True)
        if first:
            common_columns = df_main.columns
            first = False
        common_columns = df.columns.intersection(common_columns)
        df_main = df_main[common_columns]   # reduce data frame by dropping all but common columns
        logger.info(common_columns.shape)
    elif df.empty:
        no_occurrences.append(my_file)


df_main.to_pickle(os.path.join(test_dir, "main_pickle.pkl"))
logger.info("Stored merged data in %s " % os.path.join(test_dir, "main_pickle.pkl"))
pickle.dump(no_occurrences, open(os.path.join(test_dir, "no_occurrences.pkl"), "wb"))
logger.info("Stored a list of species with zero occurrences in %s " % os.path.join(test_dir, "no_occurrences.pkl"))
