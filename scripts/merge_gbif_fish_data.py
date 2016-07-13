import logging
import pickle
import os
import pandas as pd
import timeit

method = "pickle"
size = "large"
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler('./data/fish/selection/test/again/' + size + '/read_' + method + '.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

path = './data/fish/selection/test/again/' + size
test_dir = os.path.join(path, method)

important_columns = pickle.load(open("./data/fish/selection/test/again/important_columns.pkl", "rb"))
# df_main = pd.DataFrame()
no_occurrences = []
total_time = 0
first = True
for idx, my_file in enumerate(os.listdir(test_dir)):
    logger.info("Opening %s " % my_file)
    logger.info("Species #: %s " % idx)
    start_time = timeit.default_timer()
    if method == "msgpack":
        df = pd.read_msgpack(os.path.join(test_dir, my_file))
    elif method == "pickle":
        df = pd.read_pickle(os.path.join(test_dir, my_file))
    stop_time = timeit.default_timer()
    logger.info("Reading with %s took %s seconds." % (method, stop_time - start_time))
    total_time += stop_time - start_time
    if not df.empty:
        # df_main = df_main.append(df, ignore_index=True)
        # if first:
        #     common_columns_all = set(df.columns.tolist())
        #     first = False
        # common_columns = df.columns.intersection(common_columns)
        # this way memory usage doesn't grow constantly in memory, only individual frames are dumped
        common_columns = list(important_columns.intersection(set(df.columns.tolist())))
        # common_columns_all = common_columns_all.intersection(set(df.columns.tolist()))
        df[common_columns].to_msgpack(os.path.join(path, "merged_msgpack.msg"), append=True)
        # df_main = df_main[common_columns]   # reduce data frame by dropping all but common columns
        logger.info("Columns: %s" % len(common_columns))
        logger.info(df.info(memory_usage='deep'))
    elif df.empty:
        no_occurrences.append(my_file)
    del df


logger.info("Total reading time %s seconds for %s. " % (total_time, method))
# if method == "msgpack":
#     df_main.to_msgpack(os.path.join(path, "merged.msg"))
#     logger.info("Stored merged data in %s " % os.path.join(path, "merged.msg"))
# elif method == "pickle":
#     df_main.to_pickle(os.path.join(path, "merged.pkl"))
#     logger.info("Stored merged data in %s " % os.path.join(path, "merged.pkl"))

# reduce to common columns, if possible?
# df2 = pd.concat([df[common_columns] for df in pd.read_msgpack(os.path.join(path, 'merged_new.msg'))], ignore_index=True)
# df2.to_msgpack(os.path.join(path, "merged_again.msg"))
# df_reduced = pd.DataFrame()
# for df in pd.read_msgpack(os.path.join(path, "merged_new1.msg"), iterator=True):
#    df[common_columns].to_msgpack(os.path.join(path, "merged_reduced.msg"), append=True)

pickle.dump(no_occurrences, open(os.path.join(path, "no_occurrences.pkl"), "wb"))
# pickle.dump(common_columns_all, open(os.path.join(path, "common_columns_all.pkl"), "wb"))

logger.info("Stored a list of species with zero occurrences in %s " % os.path.join(path, "no_occurrences.pkl"))


# FIND ALL WITHOUT DECIMAL LATITUDE
# for df in pd.read_msgpack("/home/daniela/git/iSDM/data/fish/selection/test/again/largest/merged_new2.msg", iterator=True):
#     if "decimallatitude" not in df.columns.tolist():
#         print(df['species'].unique())
