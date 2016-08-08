#!/usr/bin/env python
"""
script: step1_climate_envelope.py
Description:
Input:
This script does the following:
"""
import logging
# import timeit
import pandas as pd
import numpy as np
from iSDM.environment import RasterEnvironmentalLayer
# import argparse

# 0. logging
logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler('./data/fish/' + '/step1_climate_envelope.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

# 1. Biomes layer
logger.info("LOADING Biomes layer.")
biomes_adf = RasterEnvironmentalLayer(file_path="./data/rebioms/w001001.adf", name_layer="Biomes")
biomes_adf.load_data()

# 2. Continents layer (vector layer originally)
logger.info("LOADING Continents layer")
from iSDM.environment import ContinentsLayer
from iSDM.environment import Source
continents = ContinentsLayer(file_path="./data/continents/continent.shp", source=Source.ARCGIS)
continents.load_data()
continents_rasters = continents.rasterize(raster_file="./data/continents/continents_raster.tif", pixel_size=0.5)
continents_rasters[0] = continents_rasters[0] + continents_rasters[2]   # combine Europe and Asia
continents_rasters[0][continents_rasters[0] > 1] = 1
continents_rasters = np.delete(continents_rasters, 2, 0)
logger.info("Continents rasters shape: %s " % (continents_rasters.shape,))

# 3. Temperature layers
logger.info("LOADING Temperature layers.")
from iSDM.environment import ClimateLayer
water_min_layer = ClimateLayer(file_path="./data/watertemp/min_wt_2000.tif")
water_min_reader = water_min_layer.load_data()
water_min_coordinates = water_min_layer.pixel_to_world_coordinates(filter_no_data_value=False)
mintemp_dataframe = pd.DataFrame([water_min_coordinates[0], water_min_coordinates[1]]).T
mintemp_dataframe.columns = ['decimallatitude', 'decimallongitude']
water_min_matrix = water_min_reader.read(1)
mintemp_dataframe['MinT'] = water_min_matrix.reshape(np.product(water_min_matrix.shape))
mintemp_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)

water_max_layer = ClimateLayer(file_path="./data/watertemp/max_wt_2000.tif")
water_max_reader = water_max_layer.load_data()
water_max_coordinates = water_max_layer.pixel_to_world_coordinates(filter_no_data_value=False)
maxtemp_dataframe = pd.DataFrame([water_max_coordinates[0], water_max_coordinates[1]]).T
maxtemp_dataframe.columns = ['decimallatitude', 'decimallongitude']
water_max_matrix = water_max_reader.read(1)
maxtemp_dataframe['MaxT'] = water_max_matrix.reshape(np.product(water_max_matrix.shape))
maxtemp_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)

water_mean_layer = ClimateLayer(file_path="./data/watertemp/mean_wt_2000.tif")
water_mean_reader = water_mean_layer.load_data()
water_mean_coordinates = water_mean_layer.pixel_to_world_coordinates(filter_no_data_value=False)
meantemp_dataframe = pd.DataFrame([water_mean_coordinates[0], water_mean_coordinates[1]]).T
meantemp_dataframe.columns = ['decimallatitude', 'decimallongitude']
water_mean_matrix = water_mean_reader.read(1)
meantemp_dataframe['MeanT'] = water_mean_matrix.reshape(np.product(water_mean_matrix.shape))
meantemp_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)

# 4. Base data frame to merge all data into
logger.info("Creating a base dataframe.")
all_coordinates = biomes_adf.pixel_to_world_coordinates(raster_data=np.zeros_like(water_min_matrix), filter_no_data_value=False)
base_dataframe = pd.DataFrame([all_coordinates[0], all_coordinates[1]]).T
base_dataframe.columns = ['decimallatitude', 'decimallongitude']
base_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
base_merged = base_dataframe.combine_first(mintemp_dataframe)
base_merged = base_merged.combine_first(maxtemp_dataframe)
base_merged = base_merged.combine_first(meantemp_dataframe)
base_merged.to_msgpack(open("./data/fish/dataframes/full_merged.msg", "wb"))

# release individual frames memory
del maxtemp_dataframe
del mintemp_dataframe
del meantemp_dataframe
# 5. Load entire fish IUCN data. (TODO: maybe we can load one by one, if the data grows beyond RAM)
# download from Google Drive: https://drive.google.com/open?id=0B9cazFzBtPuCSFp3YWE1V2JGdnc
logger.info("LOADING all species rangemaps.")
from iSDM.species import IUCNSpecies
fish = IUCNSpecies(name_species='All')
fish.load_shapefile('./data/fish/FW_FISH.shp')   # warning, 2GB of data will be loaded, may take a while!!
# 4.1 Get the list of non-extinct binomials, for looping through individual species
fish_data = fish.get_data()
fish.drop_extinct_species()
non_extinct_fish = fish.get_data()
non_extinct_binomials = non_extinct_fish.binomial.unique().tolist()

# merged.to_csv(open("./data/fish/full_merged.csv", 'w'))

# 4.2 LOOP/RASTERIZE/STORE_RASTER/MERGE_WITH_BASE_DATAFRAME
logger.info(">>>>>>>>>>>>>>>>>Looping through species!<<<<<<<<<<<<<<<<")
for name_species in non_extinct_binomials:
    fish.set_data(fish_data[fish_data.binomial == name_species])
    fish.name_species = name_species
    logger.info("Processing species: %s " % name_species)
    logger.info("Rasterizing species: %s " % name_species)
    rasterized = fish.rasterize(raster_file="./data/fish/rasterized/" + name_species + ".tif", pixel_size=0.5)
    # special case with blank map
    if not (isinstance(rasterized, np.ndarray)) or not (set(np.unique(rasterized)) == set({0, 1})):
        logger.warning("Rasterizing very small area, will use all_touched=True to avoid blank raster for species %s " % name_species)
        rasterized = fish.rasterize(raster_file="./data/fish/rasterized/" + name_species + ".tif", pixel_size=0.5, all_touched=True)
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
    logger.info("Constructing a data frame for presences and merging with main data frame.")
    presences_dataframe = pd.DataFrame([presence_coordinates[0], presence_coordinates[1]]).T
    presences_dataframe.columns = ['decimallatitude', 'decimallongitude']
    presences_dataframe[fish.name_species] = 1   # fill presences with 1's
    presences_dataframe.set_index(['decimallatitude', 'decimallongitude'], inplace=True, drop=True)
    merged = base_merged.combine_first(presences_dataframe)    # del presences_dataframe
    del presences_dataframe
    logger.info("Finished constructing a data frame for presences and merging with base data frame.")

    if pseudo_absences is not None:
        logger.info("Pixel-to-world coordinates transformation of pseudo-absences for species: %s " % name_species)
        pseudo_absence_coordinates = biomes_adf.pixel_to_world_coordinates(raster_data=pseudo_absences)
        logger.info("Finished pixel-to-world coordinates transformation of pseudo-absences for species: %s " % name_species)
        logger.info("Constructing a data frame for presences and merging with main data frame.")
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
    logger.info("Serializing to storage.")
    merged[name_species].to_msgpack(open("./data/fish/dataframes/" + name_species + "_merged.msg", "wb"))
    logger.info("Finished serializing to storage.")
    del merged
# merged.to_csv(open("./data/fish/full_merged.csv", 'a'))
logger.info("DONE!")
