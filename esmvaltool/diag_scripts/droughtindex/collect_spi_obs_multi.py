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
    for spi_file in glob.iglob(input_filenames):
        print("hello for 1")
        # Loads the file into a special structure (IRIS cube),
        # which allows us to access the data and additional information
        # with python.
        cube = iris.load(spi_file)[0]

        print("hello for 2")
        
        
        # Prints the information about the data in the iris cube
        # print("cube")
        # print(cube)

        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time, i.e. data you need
        # to plot the average over the time series of SPI on a map
        coords = ('time')
        cube2 = cube.collapsed(coords, iris.analysis.MEAN) #does it mean iris.analysi.mean compressed the dat from 3d TO 2D?
        # print("cube2")
        # print(cube2)#Prints information about the average cube #does it mean for cube2,the coordinate which is time was separately compressed to 2D and same for cube 3,the coordinate which is longitude and latitude experienced 3d to 1d?

        # Calls a python program which makes the plot
        # of 2D data on a map (see above).
        
        # Set levels to display
        spi_levels = [-4.0, -3.0, -2.0, -1.0,
                      0, 1.0, 2.0,  3.0,  4.0]
        plot_map_spi(cfg, cube2, spi_levels, name = 'SPI')

        # Here is an average over both, latitude and longitude
        # to get a mean global time series of SPI
        coords = ('longitude', 'latitude')
        cube3 = cube.collapsed(coords, iris.analysis.MEAN)
        # print("cube3")
        # print(cube3) #Prints information about the average cube

        plot_time_series_spi(cfg, cube3)
       
        # Pick an area and average the time series for this area
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        # Bremen
        index_lat = get_latlon_index(lats, 52, 53)
        index_lon = get_latlon_index(lons, 7, 9)
        add_to_filename = 'Bremen'
        # print("index_lat")
        # print(index_lat)
        # print("lats[index_lat[0]:index_lat[-1]]")
        # print(lats[index_lat[0]:index_lat[-1] + 1])
        # print("index_lon")
        # print(index_lon)
        # print("lons[index_lon[0]:index_lon[-1]]")
        # print(lons[index_lon[0]:index_lon[-1] + 1])
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        cube_grid_areas = iris.analysis.cartography.area_weights(cube[:, index_lat[0]:index_lat[-1] + 1, index_lon[0]:index_lon[-1] + 1])
        cube4 = (cube[:, index_lat[0]:index_lat[-1] + 1, index_lon[0]:index_lon[-1] + 1]).collapsed(coords, iris.analysis.MEAN, weights=cube_grid_areas)
  
        # print("cube4")
        # print(cube4) #Prints information about the average cube
 
        # extract time series from 1901-2001 from observations and (historical) model data
        start = datetime.datetime(1901, 1, 15, 0, 0, 0)
        end = datetime.datetime(2001, 12, 16, 0, 0, 0)
        time = cube4.coord('time')
        stime = time.nearest_neighbour_index(time.units.date2num(start))
        etime = time.nearest_neighbour_index(time.units.date2num(end))
        tscube4 = cube4[stime:etime]
        # Calls a python program which makes the plot
        # of time series (1D) of SPI (see above).
        plot_time_series_spi(cfg, tscube4, add_to_filename)

        # extract time series from 1901-2001 from observations and (historical) model data
        start = datetime.datetime(1901, 1, 15, 0, 0, 0)
        end = datetime.datetime(2001, 12, 16, 0, 0, 0)
        time = cube4.coord('time')
        stime = time.nearest_neighbour_index(time.units.date2num(start))
        etime = time.nearest_neighbour_index(time.units.date2num(end))
        tscube4 = cube4[stime:etime]
        # Calls a python program which makes the plot
        # of time series (1D) of SPI (see above).
        plot_time_series_spi(cfg, tscube4, add_to_filename)

        index_lat = get_latlon_index(lats, 7, 9)
        index_lon = get_latlon_index(lons, 8, 10)
        add_to_filename = 'Nigeria'
        cube_grid_areas = iris.analysis.cartography.area_weights(cube[:, index_lat[0]:index_lat[-1] + 1, index_lon[0]:index_lon[-1] + 1])
        cube5 = (cube[:, index_lat[0]:index_lat[-1] + 1, index_lon[0]:index_lon[-1] + 1]).collapsed(coords, iris.analysis.MEAN, weights=cube_grid_areas)
        #print("cube5")
        #print(cube5) #Prints information about the average cube

        start = datetime.datetime(1901, 1, 15, 0, 0, 0)
        end = datetime.datetime(2001, 12, 16, 0, 0, 0)
        time = cube5.coord('time')
        stime = time.nearest_neighbour_index(time.units.date2num(start))
        etime = time.nearest_neighbour_index(time.units.date2num(end))
        tscube5 = cube5[stime:etime]
        # Calls a python program which makes the plot
        # of time series (1D) of SPI (see above).

        plot_time_series_spi(cfg, tscube5, add_to_filename)

        start = datetime.datetime(1901, 1, 15, 0, 0, 0)
        end = datetime.datetime(2001, 12, 16, 0, 0, 0)
        time = cube5.coord('time')
        stime = time.nearest_neighbour_index(time.units.date2num(start))
        etime = time.nearest_neighbour_index(time.units.date2num(end))
        tscube5 = cube5[stime:etime]
        # Calls a python program which makes the plot
        # of time series (1D) of SPI (see above).

        plot_time_series_spi(cfg, tscube5, add_to_filename)

        # Make an aggregator from the user function.
        SPELL_COUNT = Aggregator('spell_count',
                                 count_spells,
                                 units_func=lambda units: 1)

        # Define the parameters of the test.
        threshold_spi = -2.0
        # spell_months = 1

        #print("cube.coord('time').points")
        #print(cube.coord('time').points)
        #print(cube.coord('time'))
        #print("type(cube.coord('time'))")
        #print(type(cube.coord('time')))

        # extract time series from 1901-2001 from observations and (historical) model data
        start = datetime.datetime(1901, 1, 15, 0, 0, 0)
        end = datetime.datetime(2001, 12, 16, 0, 0, 0)
        time = cube.coord('time')
        stime = time.nearest_neighbour_index(time.units.date2num(start))
        etime = time.nearest_neighbour_index(time.units.date2num(end))
        tscube = cube[stime:etime, :, :]

        # make a new cube to increase the size of the data array
        # by an additional dimension to get three (instead of one) values back from the aggregator SPELL_COUNT
        new_shape = tscube.shape + (4,)
        new_data = iris.util.broadcast_to_shape(tscube.data, new_shape, [0,1,2])
        new_cube = iris.cube.Cube(new_data)

        new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('time').points, long_name='time'), 0)
        new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('latitude').points, long_name='latitude'), 1)
        new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('longitude').points, long_name='longitude'), 2)
        new_cube.add_dim_coord(iris.coords.DimCoord([0, 1, 2, 3], long_name='z'), 3)

        # calculate the number of drought events and their average duration
        drought_periods = new_cube.collapsed('time', SPELL_COUNT,
                                          threshold=threshold_spi)
        drought_periods.rename('Number of 1 month drought')

        #print("drought_periods")
        #print(drought_periods)
        
        # plot the number of drought events 
        # drought_numbers_level = np.arange(0 , 60 , 6) # set color levels
        # use cube2 to get metadata
        # cube2.data = drought_periods.data[:,:,0]

        drought_numbers_level = np.arange(0 , 0.7 , 0.05) # set color levels
        # length of time series
        time_length = len(new_cube.coord('time').points) / 12.0
        cube2.data = drought_periods.data[:,:,0] / time_length
        plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Number_of_Events', name = 'Number of Events per year')
   
        # drought_periods.rename('Average duration of droughts')
        
        # plot the average duration of drought events 
        drought_numbers_level = np.arange(0 , 6 , 1) # set color levels
        # use cube2 to get metadata
        cube2.data = drought_periods.data[:,:,1]
        plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Duration_of_Events', name = 'Duration of Events(month)')

        # plot the average severity index of drought events 
        drought_numbers_level = np.arange(0 , 9 , 1) # set color levels
        # use cube2 to get metadata
        cube2.data = drought_periods.data[:,:,2]
        plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Severity_index_of_Events', name = 'Severity Index of Events')

        # plot the average spi of drought events 
        drought_numbers_level = np.arange(-3.5, -1.8, 0.2) # set color levels
        # use cube2 to get metadata
        cube2.data = drought_periods.data[:,:,3]
        plot_map_spi(cfg, cube2, drought_numbers_level, add_to_filename = 'Average_spi_of_Events', name = 'Average spi of Events')


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
