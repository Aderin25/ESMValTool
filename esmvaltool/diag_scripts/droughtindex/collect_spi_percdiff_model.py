#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Collects data produced by diag_save_spi.r to plot/process them further.

###############################################################################
droughtindex/collect_spi.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Collects data produced by diag_save_spi.r to plot/process them further.

Configuration options
---------------------
    TBD

###############################################################################

"""
#The imported modules and libraries below allow this script to run accordingly

import os #allows the script to interface with the underlying operating system that pyton is using on Linux
import glob #find all the pathnames matching a specified pattern according to the rules used by the Unix shell
import iris #is a dataset in the python library used for analysing and visualizing Earth science data, usually based on matplot, cartopy, numpy and dask.
from iris.util import rolling_window
from iris.analysis import Aggregator
from iris.time import PartialDateTime
import cf_units as unit
import datetime
import numpy as np #package for scientific computing(N-dimensional array)
import cartopy.crs as cart #to draw maps for visualization and analysis 
import matplotlib.pyplot as plt #plotting library that provides an object oriented API for embedding plots into application
import matplotlib.dates as mda #MDA means measured data analysis-for the time series
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_spi_func import *


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    ###########################################################################
    # Read recipe data
    ###########################################################################

    # Make an aggregator from the user function.
    SPELL_COUNT = Aggregator('spell_count',count_spells,
                            units_func=lambda units: 1)

    # Define the parameters of the test.
    threshold_spi = -2.0
    number_drought_charac = 4
    first_run = 1

    # Get filenames of input files produced by diag_spi.r
    # "cfg[n.INPUT_FILES]" is produced by the ESMValTool and contains
    # information on the SPI input files produced by diag_spi.r
    input_filenames = (cfg[n.INPUT_FILES])[0] + "/*.nc"

    # Write out search pattern for input file names
    # print("input_filenames")
    # print(input_filenames)

    # For loop: "glob.iglob" findes all files which match the
    # pattern of "input_filenames".
    # It writes the resulting exact file name onto spi_file
    # and runs the following indented lines for all possibilities
    # for spi_file.
    for iii, spi_file in enumerate(glob.iglob(input_filenames)):
        # Loads the file into a special structure (IRIS cube),
        # which allows us to access the data and additional information
        # with python.
        cube = iris.load(spi_file)[0]
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        time = cube.coord('time')
        
        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time, i.e. data you need
        # to plot the average over the time series of SPI on a map
        coords = ('time')
        cube2 = cube.collapsed(coords, iris.analysis.MEAN) # 3D to 2D (Mean over time)

        if first_run == 1:
            #shape_all = cube2.data.shape + (number_drought_charac,) + (len(input_filenames),)
            shape_all = cube2.data.shape + (number_drought_charac,) + (len(os.listdir((cfg[n.INPUT_FILES])[0])),)
            print("shape_all")
            print(shape_all)
            all_drought_hist = np.zeros(shape_all)
            all_drought_rcp85 = np.zeros(shape_all)
            print("iii")
            print(iii)
            first_run = 0
        
        # Test if time series goes until 2100/12
        timecheck = time.units.date2num(datetime.datetime(2100, 11, 30, 0, 0, 0))
        lasttime = cube.coord('time').points[-1]

        if lasttime > timecheck:
            # extract time series from 1950-2000 historical model data
            start = datetime.datetime(1950, 1, 15, 0, 0, 0)
            end = datetime.datetime(2000, 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]

            # make a new cube to increase the size of the data array
            # by an additional dimension to get two (instead of one) values back from the aggregator SPELL_COUNT
            new_shape = tscube.shape + (number_drought_charac,)
            new_data = iris.util.broadcast_to_shape(tscube.data, new_shape, [0,1,2])
            new_cube = iris.cube.Cube(new_data)

            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('time').points, long_name='time'), 0)
            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('latitude').points, long_name='latitude'), 1)
            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('longitude').points, long_name='longitude'), 2)
            new_cube.add_dim_coord(iris.coords.DimCoord(np.arange(0, number_drought_charac, 1), long_name='z'), 3)

            # calculate the number of drought events and their average duration
            drought_periods = new_cube.collapsed('time', SPELL_COUNT,
                                          threshold=threshold_spi)
            drought_periods.rename('Drought characteristics')
            # length of time series
            time_length = len(new_cube.coord('time').points) / 12.0
            # Convert number of droughtevents to frequency (per year)
            drought_periods.data[:,:,0] = drought_periods.data[:,:,0] / time_length
            all_drought_hist[:,:,:,iii] = drought_periods.data
        
            # plot the number of drought events 
            # drought_numbers_level = np.arange(0 , 60 , 6) # set color levels
            # use cube2 to get metadata
            # cube2.data = drought_periods.data[:,:,0]

            drought_numbers_level = np.arange(0 , 0.7 , 0.05) # set color levels
            # Put the data on cube2 because it contains metadata plot_map_spi needs
            cube2.data = drought_periods.data[:,:,0]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Historic_Number_of_Events_per_year', name = 'Historic_Number of Events per year')
        
            # plot the average duration of drought events 
            drought_numbers_level = np.arange(0 , 6 , 1) # set color levels
            # Put the data on cube2 because it contains metadata plot_map_spi needs
            cube2.data = drought_periods.data[:,:,1]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Historic_Duration_of_Events', name = 'Historic_Duration of Events(month)')

            # plot the average severity index of drought events 
            drought_numbers_level = np.arange(0 , 9 , 1) # set color levels
            # Put the data on cube2 because it contains metadata plot_map_spi needs
            cube2.data = drought_periods.data[:,:,2]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Historic_Severity_index_of_Events', name = 'Historic_Severity Index of Events')

            # plot the average spi of drought events 
            drought_numbers_level = np.arange(-3.5, -1.8, 0.2) # set color levels
            # Put the data on cube2 because it contains metadata plot_map_spi needs
            cube2.data = drought_periods.data[:,:,3]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Historic_Average_spi_of_Events', name = 'Historic_Average spi of Events')

            
            # extract time series from 2050-2100 rcp model data
            start = datetime.datetime(2050, 1, 15, 0, 0, 0)
            end = datetime.datetime(2100, 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]

            # make a new cube to increase the size of the data array
            # by an additional dimension to get two (instead of one) values back from the aggregator SPELL_COUNT
            new_shape = tscube.shape + (number_drought_charac,)
            new_data = iris.util.broadcast_to_shape(tscube.data, new_shape, [0,1,2])
            new_cube = iris.cube.Cube(new_data)

            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('time').points, long_name='time'), 0)
            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('latitude').points, long_name='latitude'), 1)
            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('longitude').points, long_name='longitude'), 2)
            new_cube.add_dim_coord(iris.coords.DimCoord(np.arange(0, number_drought_charac, 1), long_name='z'), 3)

            # calculate the number of drought events and their average duration
            drought_periods = new_cube.collapsed('time', SPELL_COUNT,
                                          threshold=threshold_spi)
            drought_periods.rename('Drought characteristics')
            # length of time series
            time_length = len(new_cube.coord('time').points) / 12.0
            drought_periods.data[:,:,0] = drought_periods.data[:,:,0] / time_length
            all_drought_rcp85[:,:,:,iii] = drought_periods.data
        
            # plot the number of drought events 
            drought_numbers_level = np.arange(0 , 0.7 , 0.05) # set color levels
            # use cube2 to get metadata
            cube2.data = drought_periods.data[:,:,0]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'RCP85_Number_of_Events_per_year', name = 'RCP85_Number of Events per year')
        
            # plot the average duration of drought events 
            drought_numbers_level = np.arange(0 , 6 , 1) # set color levels
            # use cube2 to get metadata
            cube2.data = drought_periods.data[:,:,1]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'RCP85_Duration_of_Events', name = 'RCP85_Duration of Events(month)')


            # plot the average severity index of drought events 
            drought_numbers_level = np.arange(0 , 9 , 1) # set color levels
            # use cube2 to get metadata
            cube2.data = drought_periods.data[:,:,2]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'RCP85_Severity_index_of_Events', name = 'RCP85_Severity Index of Events')

            # plot the average spi of drought events 
            drought_numbers_level = np.arange(-3.5, -1.8, 0.2) # set color levels
            # use cube2 to get metadata
            cube2.data = drought_periods.data[:,:,3]
            plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'RCP85_Average_spi_of_Events', name = 'RCP85_Average spi of Events')

    # Calculating multi model mean and plot it

    print("all_drought_hist")
    print(all_drought_hist)
    all_drought_hist_mean = np.mean(all_drought_hist, axis = -1) # 4D for all models to 3D multi model mean
    all_drought_rcp85_mean = np.mean(all_drought_rcp85, axis = -1) # 4D for all models to 3D multi model mean
    print("all_drought_rcp85_mean")
    print(all_drought_rcp85_mean)
    
    # Historic MultiModelMean
    data_dict = {}
    data_dict['data'] = all_drought_hist_mean[:,:,0]
    data_dict['datasetname'] = 'MultiModelMean'
    data_dict['latitude'] = lats
    data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Historic'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'Historic_Number_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0 , 0.7 , 0.05) # set color levels

    
    plot_map_spi_multi(cfg, data_dict)

    data_dict['data'] = all_drought_hist_mean[:,:,1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'Historic_Duration_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0 , 6 , 1) # set color levels
    plot_map_spi_multi(cfg, data_dict)

    data_dict['data'] = all_drought_hist_mean[:,:,2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'Historic_Severity_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0 , 9 , 1) # set color levels
    plot_map_spi_multi(cfg, data_dict)

    data_dict['data'] = all_drought_hist_mean[:,:,3]
    data_dict['drought_char'] = 'Average SPI of Events'
    data_dict['filename'] = 'Historic_Average_spi_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-3.5, -1.8, 0.2) # set color levels
    plot_map_spi_multi(cfg, data_dict)

    # RCP85 MultiModelMean
    data_dict['data'] = all_drought_rcp85_mean[:,:,0]
    data_dict['model_kind'] = 'RCP85'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'RCP85_Number_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0 , 0.7 , 0.05) # set color levels
    plot_map_spi_multi(cfg, data_dict)

    data_dict['data'] = all_drought_rcp85_mean[:,:,1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'RCP85_Duration_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0 , 6 , 1) # set color levels
    plot_map_spi_multi(cfg, data_dict)

    data_dict['data'] = all_drought_rcp85_mean[:,:,2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'RCP85_Severity_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0 , 9 , 1) # set color levels
    plot_map_spi_multi(cfg, data_dict)

    data_dict['data'] = all_drought_rcp85_mean[:,:,3]
    data_dict['drought_char'] = 'Average SPI of Events'
    data_dict['filename'] = 'RCP85_Average_spi_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-3.5, -1.8, 0.2) # set color levels
    plot_map_spi_multi(cfg, data_dict)
  


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
