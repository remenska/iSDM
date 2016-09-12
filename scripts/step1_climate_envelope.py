"""
script: step1_climate_envelope.py
Description: Creates a species-by-species dataframe to be used as input for step 1 of the 3-step approach to modeling species distributions.
             For this step, all non-extinct binomials from IUCN, are taken into account. Rangemaps are rasterized for each species,
             and 1000 random pseudo-absences are selected according to the following criteria:  all biogeographical realms covered by
             the species rangemaps are taken into account for sampling pseudo-absences.
             All the layers should be at the same resolution (30 arcmin).
             Water temperature layers (min/max/mean) are also added in the dataframe.
             Each row of the dataframe corresponds to one cell in the global
             raster layer. Each column corresponds to an environmental layer; for convenience one column identifies the realm to which each cell
             belongs. And finally, one column contains the individual species presences/absences.
             The index of the dataframe is (decimallatitude, decimallongitude)

Input:
 - Full location of the folder where the biogeographic realms (terrestrial ecoregions) shapefiles are stored.
 - Full location of the folder where the temperature raster layers (files) are location.
 - Full location of the folder where the IUCN species shapefiles are located.
 - Output location (folder) for storing the individual dataframes (csv output of the processing)

This script does the following:
 1. Loads the environment layers (biogeographical realms, temperature). In case of shapefile input, it is rasterized at the same fixed resolution as
    the other layers. Pixels from each layer are converted to world coordinates (middle of pixel location), and each layer is merged into a "base"
    dataframe as one separate column of that datafframe. Each row then represents one cell.
 2. Loads the IUCN species shapefile containing rangemaps of individual binomials in a group (such as, freshwater species, or mammals).
 3. Loops through the individual (non-extinct) binomials, rasterizing all the rangemaps that belong to each individual binomial. Selects an area
    outside the raster pixels, and samples 1000 pseudo-absence pixels (using the criteria as described above). Both the presences and pseudo-absences
    raster layers' pixels are converted to world coordinates (middle of pixel location) and merged with the "base" dataframe that contains only
    latitude/longitude as index.
 4. Serializes the constructed dataframe into a csv format.

"""
import logging
# import timeit
import pandas as pd
import numpy as np
# from iSDM.environment import RasterEnvironmentalLayer
from iSDM.environment import RealmsLayer
from iSDM.environment import Source
from iSDM.environment import ClimateLayer
from iSDM.species import IUCNSpecies
from iSDM.model import Model
# from iSDM.model import Algorithm
import os
import argparse
import errno

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--realms-location', default=os.path.join(os.getcwd(), "data", "terrestrial_ecoregions"), help='The full path to the folder where the biogeographic realms shapefiles are located.')
parser.add_argument('-t', '--temperature-location', default=os.path.join(os.getcwd(), "data", "watertemp"), help="The folder where the temperature raster files are.")
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
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.path.join(args.output_location, "step1_climate_envelope.log"))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

logger.info("Preparing a Model base dataframe")
climate_envelope_model = Model()
base_dataframe = climate_envelope_model.get_base_dataframe()
logger.info("Saving base dataframe to csv")
base_dataframe.to_csv(os.path.join(args.output_location, "base.csv"))
# 1. Temperature layers
logger.info("STEP 1: LOADING Temperature layers.")
water_min_layer = ClimateLayer(file_path=os.path.join(args.temperature_location, "min_wt_2000.tif"), name_layer="MinT")
climate_envelope_model.add_layer(water_min_layer, discard_threshold=0)
water_max_layer = ClimateLayer(file_path=os.path.join(args.temperature_location, "max_wt_2000.tif"), name_layer="MaxT")
climate_envelope_model.add_layer(water_max_layer, discard_threshold=0)
water_mean_layer = ClimateLayer(file_path=os.path.join(args.temperature_location, "mean_wt_2000.tif"), name_layer="MeanT")
climate_envelope_model.add_layer(water_mean_layer, discard_threshold=0)

logger.info("STEP 4: LOADING Biogeographical realms layer")
realms_layer = RealmsLayer(file_path=args.realms_location, source=Source.WWL, name_layer='Realm')
realms_layer.set_classifier('realm')  # which category (column) to use for grouping the polygons
realms_layer.set_raster_file(raster_file=os.path.join(args.realms_location, "realms_raster.tif"))  # destination raster file
realms_layer.set_pixel_size(pixel_size=0.5)
climate_envelope_model.add_layer(realms_layer)

logger.info("Saving base_merged to a csv dataframe...")
base_merged = climate_envelope_model.get_base_dataframe()
base_merged.to_csv(open(os.path.join(args.output_location, "base_merged.csv"), "w"))
# 5. Load entire fish IUCN data. (TODO: maybe we can load one by one, if the data grows beyond RAM)
# download from Google Drive: https://drive.google.com/open?id=0B9cazFzBtPuCSFp3YWE1V2JGdnc
logger.info("STEP 5: LOADING all species rangemaps.")
fish = IUCNSpecies(name_species='All')
fish.load_shapefile(args.species_location)   # warning, 2GB of data will be loaded, may take a while!!
# 5.1 Get the list of non-extinct binomials, for looping through individual species
fish_data = fish.get_data()
fish.drop_extinct_species()
non_extinct_fish = fish.get_data()
non_extinct_binomials = non_extinct_fish.binomial.unique().tolist()

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

# rasterized_species = IUCNSpecies(name_species="Temporary name")

# 5.2 LOOP/RASTERIZE/STORE_RASTER/MERGE_WITH_BASE_DATAFRAME
logger.info(">>>>>>>>>>>>>>>>>Looping through species!<<<<<<<<<<<<<<<<")
for idx, name_species in enumerate(non_extinct_binomials):
    fish.set_data(fish_data[fish_data.binomial == name_species])
    fish.name_species = name_species
    logger.info("ID=%s Processing species: %s " % (idx, name_species))

    if args.reprocess:
        logger.info("%s Reprocessing individual species data, will use existing raster file for species: %s" % (idx, name_species))
        full_location_raster_file = os.path.join(args.output_location, "rasterized", name_species + ".tif")
        if not os.path.exists(full_location_raster_file):
            logger.error("%s Raster file does NOT exist for species: %s " % (idx, name_species))
            continue
        raster_reader = fish.load_raster_data(raster_file=full_location_raster_file)
        rasterized = fish.raster_reader.read(1)
    else:
        logger.info("%s Rasterizing species: %s " % (idx, name_species))
        rasterized = fish.rasterize(raster_file=os.path.join(args.output_location, "rasterized", name_species + ".tif"), pixel_size=0.5)
    # special case with blank map
    if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
        logger.warning("%s Rasterizing very small area, will use all_touched=True to avoid blank raster for species %s " % (idx, name_species))
        rasterized = fish.rasterize(raster_file=os.path.join(args.output_location, "rasterized", name_species + ".tif"), pixel_size=0.5, all_touched=True)
        if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
            logger.error("%s Rasterizing did not succeed for species %s , (raster is empty)    " % (idx, name_species))
            continue
    logger.info("%s Finished rasterizing species: %s " % (idx, name_species))

    logger.info("%s Selecting pseudo-absences for species: %s " % (idx, name_species))
    selected_layers, pseudo_absences = realms_layer.sample_pseudo_absences(species_raster_data=rasterized,
                                                                           number_of_pseudopoints=1000)
    logger.info("%s Finished selecting pseudo-absences for species: %s " % (idx, name_species))

    logger.info("%s Pixel-to-world coordinates transformation of presences for species: %s " % (idx, name_species))
    presence_coordinates = fish.pixel_to_world_coordinates(raster_data=rasterized)
    logger.info("%s Finished pixel-to-world coordinates transformation of presences for species: %s " % (idx, name_species))
    logger.info("%s Constructing a data frame for presences and merging with base data frame." % idx)
    presences_dataframe = pd.DataFrame([presence_coordinates[0], presence_coordinates[1]]).T
    presences_dataframe.columns = ['decimallatitude', 'decimallongitude']
    presences_dataframe[fish.name_species] = 1   # fill presences with 1's
    presences_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
    # merged = base_dataframe.combine_first(presences_dataframe)
    merged = pd.merge(base_dataframe, presences_dataframe, how='left', left_index=True, right_index=True)
    del presences_dataframe
    logger.info("%s Finished constructing a data frame for presences and merging with base data frame." % idx)

    if pseudo_absences is not None:
        logger.info("%s Pixel-to-world coordinates transformation of pseudo-absences for species: %s " % (idx, name_species))
        pseudo_absence_coordinates = fish.pixel_to_world_coordinates(raster_data=pseudo_absences)
        logger.info("%s Finished pixel-to-world coordinates transformation of pseudo-absences for species: %s " % (idx, name_species))
        logger.info("%s Constructing a data frame for pseudo-absences and merging with base data frame." % idx)
        pseudo_absences_dataframe = pd.DataFrame([pseudo_absence_coordinates[0], pseudo_absence_coordinates[1]]).T
        pseudo_absences_dataframe.columns = ['decimallatitude', 'decimallongitude']
        pseudo_absences_dataframe[fish.name_species] = 0   # fill pseudo-absences with 0
        pseudo_absences_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
        # merged = merged.combine_first(pseudo_absences_dataframe)
        merged.update(pseudo_absences_dataframe, overwrite=False)
        del pseudo_absences_dataframe
        logger.info("%s Finished constructing a data frame for pseudo-absences and merging with base data frame." % idx)
    else:
        logger.warning("%s No pseudo absences sampled for species %s " % (idx, name_species))

    logger.info("%s Finished processing species: %s " % (idx, name_species))
    if base_dataframe.shape[0] < merged.shape[0]:
        logger.warning("%s Something is fishy with species %s : merged.shape = %s " % (idx, name_species, merged.shape[0]))
    logger.info("%s Serializing to storage." % idx)
    merged.to_csv(open(os.path.join(args.output_location, "csv", name_species + ".csv"), "w"))
    logger.info("%s Finished serializing to storage." % idx)
    logger.info("%s Shape of dataframe: %s " % (idx, merged.shape,))
    del merged
logger.info("DONE!")
