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
# The imported modules and libraries below allow this script to run accordingly

import os
import glob
import datetime
import iris
# from iris.util import rolling_window
from iris.analysis import Aggregator
# from iris.time import PartialDateTime
# import cf_units as unit
import numpy as np
import cartopy.crs as cart
# import matplotlib.pyplot as plt
# import matplotlib.dates as mda
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_drought_func import (
    count_spells, plot_map_spei_multi, plot_map_spi_multi,
    plot_map_spi, plot_map_spei)


def main(cfg):
    """Run the diagnostic.
    Parameters
    ----------
    cfg : dict
    """
    ###########################################################################
    # Read recipe data
    ###########################################################################

    # Make an aggregator from the user function.
    spell_no = Aggregator('spell_count', count_spells,
                          units_func=lambda units: 1)

    # Define the parameters of the test.
    threshold_drought = -2.0
    number_drought_charac = 4
    first_run = 1
    end_time = 2050

    # Get filenames of input files produced by diag_spi.r
    # "cfg[n.INPUT_FILES]" is produced by the ESMValTool and contains
    # information on the SPI input files produced by diag_spi.r
    input_filenames = (cfg[n.INPUT_FILES])[0] + "/*.nc"

    # Write out search pattern for input file names
    # print("input_filenames")
    # print(input_filenames)

    # For loop: "glob.iglob" findes all files which match the
    # pattern of "input_filenames".
    # It writes the resulting exact file name onto drought_file
    # and runs the following indented lines for all possibilities
    # for drought_file.
    for iii, drought_file in enumerate(glob.iglob(input_filenames)):
        # Loads the file into a special structure (IRIS cube),
        # which allows us to access the data and additional information
        # with python.
        cube = iris.load(drought_file)[0]
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        time = cube.coord('time')
        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time, i.e. data you need
        # to plot the average over the time series of SPI on a map
        coords = ('time')
        cube2 = cube.collapsed(coords, iris.analysis.MEAN)  # 3D to 2D

        if first_run == 1:
            shape_all = (cube2.data.shape + (number_drought_charac,)
                         +(len(os.listdir((cfg[n.INPUT_FILES])[0])),))
            print("shape_all")
            print(shape_all)
            all_drought_hist = np.zeros(shape_all)
            all_drought_rcp85 = np.zeros(shape_all)
            print("iii")
            print(iii)
            first_run = 0
        # Test if time series goes until 2050/12
        timecheck = time.units.date2num(datetime.datetime(end_time, 11, 30,
                                                          0, 0, 0))
        lasttime = cube.coord('time').points[-1]

        if lasttime > timecheck:
            # extract time series from 1950-2000 historical model data
            start = datetime.datetime(1950, 1, 15, 0, 0, 0)
            end = datetime.datetime(2000, 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]

            # make a new cube to increase the size of the data array
            # aggregator spell number
            new_shape = tscube.shape + (number_drought_charac,)
            new_data = iris.util.broadcast_to_shape(tscube.data, new_shape,
                                                    [0, 1, 2])
            new_cube = iris.cube.Cube(new_data)

            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('time').points,
                                                        long_name='time'), 0)
            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('latitude').points,
                                                        long_name='latitude'), 1)
            new_cube.add_dim_coord(iris.coords.DimCoord(tscube.coord('longitude').points,
                                                        long_name='longitude'), 2)
            new_cube.add_dim_coord(iris.coords.DimCoord(np.arange(0,
                                                                  number_drought_charac, 1),
                                                        long_name='z'), 3)
            # calculate the number of drought events and average duration
            drought_show = new_cube.collapsed('time', spell_no,
                                              threshold=threshold_drought)
            drought_show.rename('Drought characteristics')
            # length of time series
            time_len = len(new_cube.coord('time').points) / 12.0
            # Convert number of droughtevents to frequency (per year)
            drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                           0] / time_len
            all_drought_hist[:, :, :, iii] = drought_show.data
            # plot the number of drought events
            # drought_numbers_level = np.arange(0 , 60 , 6)
            # set color levels
            # use cube2 to get metadata
            # cube2.data = drought_show.data[:, :, 0]
            drought_numbers_level = np.arange(0, 0.6, 0.05)
            # Put the data on cube2 as it contains metadata plot_map_spi needs
            cube2.data = drought_show.data[:, :, 0]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_No_of_Events_per_year',
                         name='Historic_Number of Events per year')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_No_of_Events_per_year',
                         name='Historic_Number of Events per year')
            # plot the average duration of drought events
            drought_numbers_level = np.arange(0, 6, 1)  # set color levels
            # Put the data on cube2,it contains metadata plot_map_spi needs
            cube2.data = drought_show.data[:, :, 1]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Dur_of_Events',
                         name='Historic_Duration of Events(month)')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Dur_of_Events',
                         name='Historic_Duration of Events(month)')

            # plot the average severity index of drought events
            drought_numbers_level = np.arange(0, 9, 1)
            # set color levels
            cube2.data = drought_show.data[:, :, 2]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Sev_index_of_Events',
                         name='Historic_Severity Index of Events')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Sev_index_of_Events',
                         name='Historic_Severity Index of Events')

            # plot the average spi of drought events
            drought_numbers_level = np.arange(-3.0, -1.8, 0.2)
            # set color levels
            cube2.data = drought_show.data[:, :, 3]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Avr_spi_of_Events',
                         name='Historic_Average spi of Events')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Avr_spei_of_Events',
                         name='Historic_Average spei of Events')
            # extract time series from 2050-2100 rcp model data
            start = datetime.datetime(2000, 1, 15, 0, 0, 0)
            end = datetime.datetime(end_time, 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]

            # make a new cube to increase the size of the data array
            # get two (instead of one) values back from the
            # aggregator spell_no
            new_shape = tscube.shape + (number_drought_charac,)
            new_data = iris.util.broadcast_to_shape(
                tscube.data, new_shape, [0, 1, 2])
            new_cube = iris.cube.Cube(new_data)

            new_cube.add_dim_coord(iris.coords.DimCoord(
                tscube.coord('time').points, long_name='time'), 0)
            new_cube.add_dim_coord(iris.coords.DimCoord(
                tscube.coord('latitude').points, long_name='latitude'), 1)
            new_cube.add_dim_coord(iris.coords.DimCoord(
                tscube.coord('longitude').points, long_name='longitude'), 2)
            new_cube.add_dim_coord(iris.coords.DimCoord(
                np.arange(0, number_drought_charac, 1), long_name='z'), 3)

            # calculate the number of drought events and their average duration
            drought_show = new_cube.collapsed('time', spell_no,
                                              threshold=threshold_drought)
            drought_show.rename('Drought characteristics')
            # length of time series
            time_len = len(new_cube.coord('time').points) / 12.0
            drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                           0] / time_len
            all_drought_rcp85[:, :, :, iii] = drought_show.data

            # plot the number of drought events
            drought_numbers_level = np.arange(0, 0.6, 0.05)
            # set color levels
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 0]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_No_of_Events_per_year',
                         name='RCP85_Number of Events per year')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_No_of_Events_per_year',
                         name='RCP85_Number of Events per year')

            # plot the average duration of drought events
            drought_numbers_level = np.arange(0, 6, 1)  # set color levels
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 1]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_Dur_of_Events',
                         name='RCP85_Duration of Events(month)')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_Dur_of_Events',
                         name='RCP85_Duration of Events(month)')

            # plot the average severity index of drought events
            drought_numbers_level = np.arange(0, 9, 1)  # set color levels
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 2]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_Sev_index_of_Events',
                         name='RCP85_Severity Index of Events')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_Sev_index_of_Events',
                         name='RCP85_Severity Index of Events')

            # plot the average spi of drought events
            drought_numbers_level = np.arange(-3.0, -1.8, 0.2)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 3]
            plot_map_spi(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_Avr_spi_of_Events',
                         name='RCP85_Average spi of Events')
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='RCP85_Avr_spei_of_Events',
                         name='RCP85_Average spei of Events')
    # Calculating multi model mean and plot it
    print("all_drought_hist")
    print(all_drought_hist)
    all_drought_hist_mean = np.nanmean(all_drought_hist, axis=-1)
    # to 3D multi model mean
    all_drought_rcp85_mean = np.nanmean(all_drought_rcp85, axis=-1)
    # to 3D multi model mean
    perc_diff = ((all_drought_rcp85_mean - all_drought_hist_mean)
                 / (all_drought_rcp85_mean + all_drought_hist_mean) * 200)
    print("all_drought_rcp85_mean")
    print(all_drought_rcp85_mean)

    # Historic 
    data_dict = {}
    data_dict['data'] = all_drought_hist_mean[:, :, 0]
    data_dict['datasetname'] = 'MultiModelMean'
    data_dict['latitude'] = lats
    data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Historic'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'Historic_No_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0, 0.4, 0.05)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'Historic_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 6, 1)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplott')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'Historic_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 3]
    data_dict['drought_char'] = 'Average SPI of Events'
    data_dict['filename'] = 'Historic_Average_spi_of_Events'
    data_dict['drought_char'] = 'Average SPEI of Events'
    data_dict['filename'] = 'Historic_Average_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # RCP85 
    data_dict['data'] = all_drought_rcp85_mean[:, :, 0]
    data_dict['model_kind'] = 'RCP85'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'RCP85_No_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0, 0.4, 0.05)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_rcp85_mean[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'RCP85_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 6, 1)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_rcp85_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'RCP85_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_rcp85_mean[:, :, 3]
    data_dict['drought_char'] = 'Average SPI of Events'
    data_dict['filename'] = 'RCP85_Avr_spi_of_Events'
    data_dict['drought_char'] = 'Average SPEI of Events'
    data_dict['filename'] = 'RCP85_Avr_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spi_multi(cfg, data_dict, colormap='gnuplot')
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # RCP85 Percentage difference
    data_dict['data'] = perc_diff[:, :, 0]
    data_dict['datasetname'] = 'Percentage'
    # data_dict['latitude'] = lats
    # data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Difference'
    data_dict['drought_char'] = 'Number of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_No_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-100, 110, 10)
    plot_map_spi_multi(cfg, data_dict, colormap='rainbow')
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-100, 110, 10)
    plot_map_spi_multi(cfg, data_dict, colormap='rainbow')
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Sev_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spi_multi(cfg, data_dict, colormap='rainbow')
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 3]
    data_dict['drought_char'] = 'Average SPI of Events [%]'
    data_dict['drought_char'] = 'Average SPEI of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Avr_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spi_multi(cfg, data_dict, colormap='rainbow')
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
