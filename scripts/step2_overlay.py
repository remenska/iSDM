import logging
# import timeit
import pandas as pd
import numpy as np
# from iSDM.environment import RasterEnvironmentalLayer
from iSDM.environment import RealmsLayer, RasterEnvironmentalLayer
from iSDM.environment import Source
from iSDM.environment import ClimateLayer
from iSDM.species import IUCNSpecies
from iSDM.model import Model
from rasterio.transform import Affine

# from iSDM.model import Algorithm
import os
import argparse
import errno
import rasterio

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--realms-location', default=os.path.join(os.getcwd(), "data", "terrestrial_ecoregions"), help='The full path to the folder where the biogeographic realms shapefiles are located.')
parser.add_argument('-l', '--habitat-location', default=os.path.join(os.getcwd(), "data", "GLWD", "downscaled"), help="The folder where the temperature raster files are.")
parser.add_argument('-s', '--species-location', default=os.path.join(os.getcwd(), "data", "fish"), help="The folder where the IUCN species shapefiles are located.")
parser.add_argument('-o', '--output-location', default=os.path.join(os.getcwd(), "data", "fish"), help="Output location (folder) for storing the output of the processing.")
parser.add_argument('--reprocess', action='store_true', help="Reprocess the data, using the already-rasterized individual species rangemaps. Assumes these files are all available.")
parser.set_defaults(reprocess=False)
args = parser.parse_args()

# 0. logging
try:
    os.makedirs(args.output_location)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
fh = logging.FileHandler(os.path.join(args.output_location, "step2_overlay.log"))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

pixel_size = 0.083333333
x_min, y_min, x_max, y_max = -180, -90, 180, 90
x_res = int((x_max - x_min) / pixel_size)
y_res = int((y_max - y_min) / pixel_size)

habitat_model = Model(pixel_size=pixel_size)
base_dataframe = habitat_model.get_base_dataframe()

# Load suitable habitat layer
logger.info("Adding GLWD layer:")
glwd_layer = RasterEnvironmentalLayer(file_path=os.path.join(args.habitat_location, "downscaled_cubic_again.tif"), name_layer="GLWD")
glwd_reader = glwd_layer.load_data()
glwd_data = glwd_reader.read(1)
# hmmm this below is useless, as the data is read from scratch again when adding the layer :-/
# better prepare the layer in the expected form
glwd_data[glwd_data != glwd_reader.nodata] = 1  # unify all pixel values (now ranging [1,..,12])
glwd_data[glwd_data == glwd_reader.nodata] = 0  # cutoff any nodata values (like here 255)
logger.info("GWLD data has %s data pixels. " % np.count_nonzero(glwd_data))
if glwd_data.shape != (y_res, x_res):
    logger.error("The layer is not at the proper resolution! Layer shape:%s " % (glwd_data.shape, ))

habitat_model.add_environmental_layer(glwd_layer)   # not needed to keep the whole layer?
logger.info("Saving base_merged dataframe to csv")
base_merged = habitat_model.get_base_dataframe()
base_merged.to_csv(os.path.join(args.output_location, "base_merged.csv"))
del base_merged
logger.info("STEP 3: LOADING all species rangemaps.")
species = IUCNSpecies(name_species='All')
species.load_shapefile(args.species_location)   # warning, all species data will be loaded, may take a while!!
# 3.1 Get the list of non-extinct binomials, for looping through individual species
species_data = species.get_data()
species.drop_extinct_species()
non_extinct_species = species.get_data()
non_extinct_binomials = non_extinct_species.binomial.unique().tolist()

try:
    os.makedirs(os.path.join(args.output_location, "rasterized"))
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
try:
    os.makedirs(os.path.join(args.output_location, "csv"))
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

# 3.2 LOOP/RASTERIZE/STORE_RASTER/MERGE_WITH_BASE_DATAFRAME
logger.info(">>>>>>>>>>>>>>>>>Looping through species!<<<<<<<<<<<<<<<<")
for idx, name_species in enumerate(non_extinct_binomials):
    species.set_data(species_data[species_data.binomial == name_species])
    species.name_species = name_species
    logger.info("ID=%s Processing species: %s " % (idx, name_species))

    if args.reprocess:
        logger.info("%s Reprocessing individual species data, will use existing raster file for species: %s" % (idx, name_species))
        full_location_raster_file = os.path.join(args.output_location, "rasterized", name_species + ".tif")
        if not os.path.exists(full_location_raster_file):
            logger.error("%s Raster file does NOT exist for species: %s " % (idx, name_species))
            continue
        raster_reader = species.load_raster_data(raster_file=full_location_raster_file)
        rasterized = species.raster_reader.read(1)
    else:
        logger.info("%s Rasterizing species: %s " % (idx, name_species))
        rasterized = species.rasterize(raster_file=os.path.join(args.output_location, "rasterized", name_species + ".tif"), pixel_size=pixel_size)
    # special case with blank map after rasterizing species (could be too small region)
    # we attempt to rasterize it with all_touched=True which means any pixel touching a geometry will be burned
    if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
        logger.warning("%s Rasterizing very small area, will use all_touched=True to avoid blank raster for species %s " % (idx, name_species))
        rasterized = species.rasterize(raster_file=os.path.join(args.output_location, "rasterized", name_species + ".tif"), pixel_size=pixel_size, all_touched=True)
        if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
            logger.error("%s Rasterizing did not succeed for species %s , (raster is empty)    " % (idx, name_species))
            continue
    logger.info("%s Finished rasterizing species: %s " % (idx, name_species))

    # filter out unsuitable habitat (overlay species raster with the glwd)
    expert_range_suitable = glwd_data * rasterized

    # TODO: we will not need this; just for sanity check
    crs = {'init': "EPSG:4326"}
    transform = Affine.translation(x_min, y_max) * Affine.scale(pixel_size, -pixel_size)
    with rasterio.open(os.path.join(args.output_location, "rasterized", name_species + "_suitable.tif"), 'w', driver='GTiff', width=x_res, height=y_res,
                       count=1,
                       dtype=np.uint8,
                       nodata=0,
                       transform=transform,
                       crs=crs) as out:
        out.write(expert_range_suitable.astype(np.uint8), indexes=1)
        out.close()

    # TODO rasterize the corresponding GBIF records and overlay
    # probably will not need this column either, first we need to overlay with rasterized GBIF records to filter out more pixels
    logger.info("%s Pixel-to-world coordinates transformation of suitable habitat for species: %s " % (idx, name_species))
    suitable_habitat_coordinates = species.pixel_to_world_coordinates(raster_data=expert_range_suitable)
    logger.info("%s Finished pixel-to-world coordinates transformation of suitable habitat for species: %s " % (idx, name_species))
    logger.info("%s Constructing a data frame for suitable_habitat_coordinates and merging with base data frame." % idx)
    suitable_habitat_dataframe = pd.DataFrame([suitable_habitat_coordinates[0], suitable_habitat_coordinates[1]]).T
    suitable_habitat_dataframe.columns = ['decimallatitude', 'decimallongitude']
    suitable_habitat_dataframe[species.name_species] = 1   # fill presences with 1's
    suitable_habitat_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
    # merged = base_dataframe.combine_first(presences_dataframe)
    merged = pd.merge(base_dataframe, suitable_habitat_dataframe, how='left', left_index=True, right_index=True)
    del suitable_habitat_dataframe
    logger.info("%s Serializing to storage." % idx)
    merged.to_csv(open(os.path.join(args.output_location, "csv", name_species + ".csv"), "w"))
    logger.info("%s Finished constructing a data frame for presences and merging with base data frame." % idx)
    del merged
logger.info("DONE!")
