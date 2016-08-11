#!/usr/bin/env python
"""
script: step1_climate_envelope.py
Description:
Input:
- Biomes raster folder/file
- Continents shapefile

- temperature(/min/max/mean) folder/file
- log file location?
This script does the following:

Output:
- base_dataframe
- individual csv files
"""
import logging
# import timeit
import pandas as pd
import numpy as np
from iSDM.environment import RasterEnvironmentalLayer
from iSDM.environment import ContinentsLayer
from iSDM.environment import Source
from iSDM.environment import ClimateLayer
from iSDM.species import IUCNSpecies
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--biomes-location', default="./data/rebioms/w001001.adf", help='The full location of the folder and biomes raster file.')
parser.add_argument('-c', '--continents-location', default="./data/continents/", help='The full path to the folder where the continents shapefiles are located.')
parser.add_argument('-t', '--temperature-location', default="./data/watertemp/", help="The folder where the temperature raster files are.")
parser.add_argument('-s', '--species-location', default='./data/fish/', help="The folder where the species shapefiles are located.")
parser.add_argument('-o', '--output-location', default="./data/fish/", help="Output location (folder) for storing the output of the processing.")
parser.add_argument('--reprocess', action='store_true', help="Reprocess the data, using the already-rasterized individual species rangemaps. Assumes these files are all available.")
parser.set_defaults(reprocess=False)
args = parser.parse_args()

# 0. logging
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(args.output_location + '/step1_climate_envelope.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

# 1. Biomes layer
logger.info("LOADING Biomes layer.")
biomes_adf = RasterEnvironmentalLayer(file_path=args.biomes_location, name_layer="Biomes")
biomes_adf.load_data()

# 2. Continents layer (vector layer originally)
logger.info("LOADING Continents layer")
continents = ContinentsLayer(file_path=args.continents_location, source=Source.ARCGIS)
continents.load_data()
continents_rasters = continents.rasterize(raster_file=args.continents_location + "/continents_raster.tif", pixel_size=0.5, all_touched=True)
continents_rasters[0] = continents_rasters[0] + continents_rasters[2]   # combine Europe and Asia
continents_rasters[0][continents_rasters[0] > 1] = 1
continents_rasters = np.delete(continents_rasters, 2, 0)
logger.info("Continents rasters shape: %s " % (continents_rasters.shape,))

# 2. 1
# merge continents on a single band, so oceans pixels can be discarded
# immediately in the "base" data frame
continents_flattened = np.zeros_like(continents_rasters[0])
for continent in continents_rasters:
    continents_flattened += continent

# 3. Temperature layers
logger.info("LOADING Temperature layers.")
water_min_layer = ClimateLayer(file_path=args.temperature_location + "/min_wt_2000.tif")
water_min_reader = water_min_layer.load_data()
water_min_data = water_min_reader.read(1)
# cut out anything below 0 Kelvins, absolute zero.
water_min_data[water_min_data < 0] = 0
water_min_coordinates = water_min_layer.pixel_to_world_coordinates(raster_data=water_min_data,
                                                                   filter_no_data_value=True)
mintemp_dataframe = pd.DataFrame([water_min_coordinates[0], water_min_coordinates[1]]).T
mintemp_dataframe.columns = ['decimallatitude', 'decimallongitude']
flattened_watermin_data = water_min_data.reshape(np.product(water_min_data.shape))
mintemp_dataframe['MinT'] = flattened_watermin_data[flattened_watermin_data != 0]
mintemp_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
logger.info("!!! Shape mintemp_dataframe: % s " % (mintemp_dataframe.shape, ))

water_max_layer = ClimateLayer(file_path=args.temperature_location + "/max_wt_2000.tif")
water_max_reader = water_max_layer.load_data()
water_max_data = water_max_reader.read(1)
# cut out anything below 0 Kelvins, absolute zero.
water_max_data[water_max_data < 0] = 0
water_max_coordinates = water_max_layer.pixel_to_world_coordinates(raster_data=water_max_data,
                                                                   filter_no_data_value=True)
maxtemp_dataframe = pd.DataFrame([water_max_coordinates[0], water_max_coordinates[1]]).T
maxtemp_dataframe.columns = ['decimallatitude', 'decimallongitude']
flattened_watermax_data = water_max_data.reshape(np.product(water_max_data.shape))
maxtemp_dataframe['MaxT'] = flattened_watermax_data[flattened_watermax_data != 0]
maxtemp_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
logger.info("!!! Shape maxtemp_dataframe: % s " % (maxtemp_dataframe.shape, ))

water_mean_layer = ClimateLayer(file_path=args.temperature_location + "/mean_wt_2000.tif")
water_mean_reader = water_mean_layer.load_data()
water_mean_data = water_mean_reader.read(1)
# cut out anything below 0 Kelvins, absolute zero.
water_mean_data[water_mean_data < 0] = 0
water_mean_coordinates = water_mean_layer.pixel_to_world_coordinates(raster_data=water_mean_data,
                                                                     filter_no_data_value=True)
meantemp_dataframe = pd.DataFrame([water_mean_coordinates[0], water_mean_coordinates[1]]).T
meantemp_dataframe.columns = ['decimallatitude', 'decimallongitude']
flattened_watermean_data = water_mean_data.reshape(np.product(water_mean_data.shape))
meantemp_dataframe['MeanT'] = flattened_watermean_data[flattened_watermean_data != 0]
meantemp_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
logger.info("!!! Shape meantemp_dataframe: % s " % (meantemp_dataframe.shape, ))

# 4. Construct a base and a merged dataframe
logger.info("Creating a base dataframe.")

# Right now the base dataframe contains all cells (360x720); the advantage of this is that some operations may be a lot
# simpler if all dataframes contain the same number of rows (the index being latitude/longitude)
all_coordinates = biomes_adf.pixel_to_world_coordinates(raster_data=np.zeros_like(water_mean_data), filter_no_data_value=False)
base_dataframe = pd.DataFrame([all_coordinates[0], all_coordinates[1]]).T
# Alternatively, the non-continent (typically ocean) coordinates can be discarded like the commented code below
# all_non_ocean_coordinates = biomes_adf.pixel_to_world_coordinates(raster_data=continents_flattened, filter_no_data_value=True)
# base_dataframe = pd.DataFrame([all_non_ocean_coordinates[0], all_non_ocean_coordinates[1]]).T
# base_dataframe contains only the index, i.e., latitude/longitude of the data
# we will use this dataframe to combine with the individual species presences/absences. Combine with matching on the index.
base_dataframe.columns = ['decimallatitude', 'decimallongitude']
base_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)

base_dataframe.to_csv(args.output_location + "/base.csv")
logger.info("!!! Shape of base: %s " % (base_dataframe.shape, ))
# base_merged contains the columns on min/max/mean temperature and continent ID
# this dataframe can be used to reconstruct the full data about a particular species, as input for modeling
# we don't want to keep duplicate data in each individual species data frame.

base_merged = base_dataframe.combine_first(mintemp_dataframe)
base_merged = base_merged.combine_first(maxtemp_dataframe)
base_merged = base_merged.combine_first(meantemp_dataframe)

# 5. Add continent column, with an ID of the continent
for idx, band in enumerate(continents_rasters):
    continents_coordinates = biomes_adf.pixel_to_world_coordinates(raster_data=band)
    continent_dataframe = pd.DataFrame([continents_coordinates[0], continents_coordinates[1]]).T
    continent_dataframe.columns = ['decimallatitude', 'decimallongitude']
    continent_dataframe['Continent'] = idx
    continent_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
    base_merged = base_merged.combine_first(continent_dataframe)

base_merged.to_csv(open(args.output_location + "/base_merged.csv", "w"))
logger.info("!!! Shape of base_merged: %s " % (base_merged.shape, ))
# release individual frames memory
del maxtemp_dataframe
del mintemp_dataframe
del meantemp_dataframe
del continent_dataframe

# 6. Load entire fish IUCN data. (TODO: maybe we can load one by one, if the data grows beyond RAM)
# download from Google Drive: https://drive.google.com/open?id=0B9cazFzBtPuCSFp3YWE1V2JGdnc
logger.info("LOADING all species rangemaps.")
fish = IUCNSpecies(name_species='All')
fish.load_shapefile(args.species_location)   # warning, 2GB of data will be loaded, may take a while!!
# 4.1 Get the list of non-extinct binomials, for looping through individual species
fish_data = fish.get_data()
fish.drop_extinct_species()
non_extinct_fish = fish.get_data()
non_extinct_binomials = non_extinct_fish.binomial.unique().tolist()

# merged.to_csv(open("./data/fish/full_merged.csv", 'w'))

os.makedirs(args.output_location + "/rasterized/", exist_ok=True)
os.makedirs(args.output_location + "/csv/", exist_ok=True)

# rasterized_species = IUCNSpecies(name_species="Temporary name")

# 4.2 LOOP/RASTERIZE/STORE_RASTER/MERGE_WITH_BASE_DATAFRAME
logger.info(">>>>>>>>>>>>>>>>>Looping through species!<<<<<<<<<<<<<<<<")
for idx, name_species in enumerate(non_extinct_binomials):
    fish.set_data(fish_data[fish_data.binomial == name_species])
    fish.name_species = name_species
    logger.info("ID=%s Processing species: %s " % (idx, name_species))

    if args.reprocess:
        logger.info("Reprocessing individual species data, will use existing raster file for species: %s" % name_species)
        full_location_raster_file = args.output_location + "/rasterized/" + name_species + ".tif"
        if not os.path.exists(full_location_raster_file):
            logger.error("Raster file does NOT exist for species: %s " % name_species)
            continue
        raster_reader = fish.load_raster_data(raster_file=full_location_raster_file)
        rasterized = fish.raster_reader.read(1)
    else:
        logger.info("Rasterizing species: %s " % name_species)
        rasterized = fish.rasterize(raster_file=args.output_location + "/rasterized/" + name_species + ".tif", pixel_size=0.5)
    # special case with blank map
    if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
        logger.warning("Rasterizing very small area, will use all_touched=True to avoid blank raster for species %s " % name_species)
        rasterized = fish.rasterize(raster_file=args.output_location + "/rasterized/" + name_species + ".tif", pixel_size=0.5, all_touched=True)
        if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
            logger.error("Rasterizing did not succeed for species %s , (raster is empty)    " % name_species)
            continue
    logger.info("Finished rasterizing species: %s " % name_species)

    logger.info("Selecting pseudo-absences for species: %s " % name_species)
    selected_layers, pseudo_absences = biomes_adf.sample_pseudo_absences(species_raster_data=rasterized,
                                                                         continents_raster_data=continents_rasters,
                                                                         number_of_pseudopoints=1000)
    logger.info("Finished selecting pseudo-absences for species: %s " % name_species)

    logger.info("Pixel-to-world coordinates transformation of presences for species: %s " % name_species)
    presence_coordinates = fish.pixel_to_world_coordinates(raster_data=rasterized)
    logger.info("Finished pixel-to-world coordinates transformation of presences for species: %s " % name_species)
    logger.info("Constructing a data frame for presences and merging with base data frame.")
    presences_dataframe = pd.DataFrame([presence_coordinates[0], presence_coordinates[1]]).T
    presences_dataframe.columns = ['decimallatitude', 'decimallongitude']
    presences_dataframe[fish.name_species] = 1   # fill presences with 1's
    presences_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
    merged = base_dataframe.combine_first(presences_dataframe)
    del presences_dataframe
    logger.info("Finished constructing a data frame for presences and merging with base data frame.")

    if pseudo_absences is not None:
        logger.info("Pixel-to-world coordinates transformation of pseudo-absences for species: %s " % name_species)
        pseudo_absence_coordinates = biomes_adf.pixel_to_world_coordinates(raster_data=pseudo_absences)
        logger.info("Finished pixel-to-world coordinates transformation of pseudo-absences for species: %s " % name_species)
        logger.info("Constructing a data frame for presences and merging with base data frame.")
        pseudo_absences_dataframe = pd.DataFrame([pseudo_absence_coordinates[0], pseudo_absence_coordinates[1]]).T
        pseudo_absences_dataframe.columns = ['decimallatitude', 'decimallongitude']
        pseudo_absences_dataframe[fish.name_species] = 0   # fill pseudo-absences with 0
        pseudo_absences_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
        merged = merged.combine_first(pseudo_absences_dataframe)
        del pseudo_absences_dataframe
        logger.info("Finished constructing a data frame for pseudo-absences and merging with base data frame.")
    else:
        logger.warning("No pseudo absences sampled for species %s " % name_species)

    logger.info("Finished processing species: %s " % name_species)
    if base_dataframe.shape[0] < merged.shape[0]:
        logger.warning("Something is fishy with species %s : merged.shape = %s " % (name_species, merged.shape[0]))
    logger.info("Serializing to storage.")
    merged.to_csv(open(args.output_location + "/csv/" + name_species + ".csv", "w"))
    logger.info("Finished serializing to storage.")
    logger.info("Shape of dataframe: %s " % (merged.shape,))
    del merged
logger.info("DONE!")
