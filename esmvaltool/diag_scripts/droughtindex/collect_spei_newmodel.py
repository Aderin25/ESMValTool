#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Collects data produced by diag_save_sspei.r to plot/process them further.

###############################################################################
droughtindex/collect_spei.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Collects data produced by diag_save_spei.r to plot/process them further.

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
# import cartopy.crs as cart
# import matplotlib.pyplot as plt
# import matplotlib.dates as mda
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.droughtindex.collect_spei_newfunc import (
    count_spells, plot_map_spei_multi, plot_map_spei,
    plot_histogram_spei_multi)


def main(cfg):
    """Run the diagnostic.
    Parameters
    ----------
    cfg : dict
    """
    ######################################################################
    # Read recipe data
    ######################################################################

    # Make an aggregator from the user function.
    spell_no = Aggregator('spell_count', count_spells,
                          units_func=lambda units: 1)

    # Define the parameters of the test.
    threshold_spei = -2.0
    number_drought_charac = 4
    first_run = 1
    start_year_hist = 1950
    end_year_hist = 2000
    start_year_future = 2050
    end_year_future = 2100

    # Get filenames of input files produced by diag_spei.r
    input_filenames = (cfg[n.INPUT_FILES])[0] + "/*.nc"

    for iii, spei_file in enumerate(glob.iglob(input_filenames)):
        # Loads the file into a special structure (IRIS cube),
        # which allows us to access the data and additional information
        # with python.
        cube = iris.load(spei_file)[0]
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        time = cube.coord('time')
        # The data are 3D (time x latitude x longitude)
        # To plot them, we need to reduce them to 2D or 1D
        # First here is an average over time.
        coords = ('time')
        cube2 = cube.collapsed(coords, iris.analysis.MEAN)  # 3D to 2D
        
        
        # Get model name from cube    
        try:
            modelname = cube.metadata.attributes['model_id']
        except KeyError:
            try:
                modelname = cube.metadata.attributes['source_id']
            except KeyError:
                modelname = 'Unknown'

        if first_run == 1:
            files = os.listdir((cfg[n.INPUT_FILES])[0])
            ncfiles = list(filter(lambda f: f.endswith('.nc'), files))
            shape_all = cube2.data.shape + (number_drought_charac,) + \
                (len(ncfiles),)
            all_drought_hist = np.full(shape_all, np.nan)
            all_drought_ssp585 = np.full(shape_all, np.nan)
            modelnames = (modelname,)
            first_run = 0
        else:
            modelnames = modelnames + (modelname,)
        # Test if time series goes until 2100/12
        timecheck = time.units.date2num(datetime.datetime(end_year_future, 11,
                                                          30, 0, 0, 0))
        lasttime = cube.coord('time').points[-1]

        if lasttime > timecheck:
            # extract time series from 1950-2000 historical model data
            start = datetime.datetime(start_year_hist, 1, 15, 0, 0, 0)
            end = datetime.datetime(end_year_hist, 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]
            print("tscube.data.shape hist")
            print(tscube.data.shape)
            print(locals().keys())
            #print(locals().values())
            if 'all_spei_hist' in locals().keys():
                all_spei_hist[:, :, :, iii] = tscube.data
            else:
                shape_histo = tscube.data.shape + (len(ncfiles),)
                all_spei_hist = np.full(shape_histo, np.nan)
                all_spei_hist[:, :, :, iii] = tscube.data

            # make a new cube to increase the size of the data array
            # aggregator spell number
            new_shape = tscube.shape + (number_drought_charac,)
            new_data = iris.util.broadcast_to_shape(tscube.data, new_shape,
                                                    [0, 1, 2])
            new_cube = iris.cube.Cube(new_data)

            new_cube.add_dim_coord(iris.coords.DimCoord(
                tscube.coord('time').points, long_name='time'), 0)
            new_cube.add_dim_coord(iris.coords.DimCoord(
                tscube.coord('latitude').points, long_name='latitude'), 1)
            new_cube.add_dim_coord(iris.coords.DimCoord(
                tscube.coord('longitude').points, long_name='longitude'), 2)
            new_cube.add_dim_coord(iris.coords.DimCoord(
                np.arange(0, number_drought_charac, 1), long_name='z'), 3)

            # calculate the number of drought events and average duration
            drought_show = new_cube.collapsed('time', spell_no,
                                              threshold=threshold_spei)
            drought_show.rename('Drought characteristics')
            time_len = len(new_cube.coord('time').points) / 12.0
            # Convert number of droughtevents to frequency (per year)
            drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                           0] / time_len
            all_drought_hist[:, :, :, iii] = drought_show.data
            drought_numbers_level = np.arange(0, 0.4, 0.05)
            cube2.data = drought_show.data[:, :, 0]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_No_of_Events_per_year',
                         name='Historic_Number of Events per year')

            # plot the average duration of drought events
            drought_numbers_level = np.arange(0, 6, 1)
            cube2.data = drought_show.data[:, :, 1]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Dur_of_Events',
                         name='Historic_Duration of Events(month)')

            # plot the average severity index of drought events
            drought_numbers_level = np.arange(0, 9, 1)
            cube2.data = drought_show.data[:, :, 2]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Sev_index_of_Events',
                         name='Historic_Severity Index of Events')

            # plot the average spei of drought events
            drought_numbers_level = np.arange(-2.8, -1.8, 0.2)
            cube2.data = drought_show.data[:, :, 3]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='Historic_Avr_spei_of_Events',
                         name='Historic_Average spei of Events')
            # extract time series from 2050-2100 ssp model data
            start = datetime.datetime(start_year_future, 1, 15, 0, 0, 0)
            end = datetime.datetime(end_year_future, 12, 16, 0, 0, 0)
            stime = time.nearest_neighbour_index(time.units.date2num(start))
            etime = time.nearest_neighbour_index(time.units.date2num(end))
            tscube = cube[stime:etime, :, :]
            
            print("tscube.data.shape hist")
            print(tscube.data.shape)
            if 'all_spei_ssp585' in locals().keys():
                all_spei_ssp585[:, :, :, iii] = tscube.data
            else:
                shape_histo = tscube.data.shape + (len(ncfiles),)
                all_spei_ssp585 = np.full(shape_histo, np.nan)
                all_spei_ssp585[:, :, :, iii] = tscube.data
            

            # make a new cube to increase the size of the data array
            # get two (instead of one) values from the aggregator spell_no
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
                                              threshold=threshold_spei)
            drought_show.rename('Drought characteristics')
            # length of time series
            time_len = len(new_cube.coord('time').points) / 12.0
            drought_show.data[:, :, 0] = drought_show.data[:, :,
                                                           0] / time_len
            all_drought_ssp585[:, :, :, iii] = drought_show.data

            # plot the number of drought events
            drought_numbers_level = np.arange(0, 0.4, 0.05)
            # set color levels
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 0]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='SSP585_No_of_Events_per_year',
                         name='SSP585_Number of Events per year')

            # plot the average duration of drought events
            drought_numbers_level = np.arange(0, 6, 1)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 1]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='SSP585_Dur_of_Events',
                         name='SSP585_Duration of Events(month)')

            # plot the average severity index of drought events
            drought_numbers_level = np.arange(0, 9, 1)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 2]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='SSP585_Sev_index_of_Events',
                         name='SSP585_Severity Index of Events')

            # plot the average spei of drought events
            drought_numbers_level = np.arange(-2.8, -1.8, 0.2)
            # use cube2 to get metadata
            cube2.data = drought_show.data[:, :, 3]
            plot_map_spei(cfg, cube2, drought_numbers_level,
                         add_to_filename='SSP585_Avr_spei_of_Events',
                         name='SSP585_Average spei of Events')
    # Calculating multi model mean and plot it
    print("all_drought_hist")
    print(all_drought_hist)
    all_drought_hist_mean = np.nanmean(all_drought_hist, axis=-1)
    # Variance
    all_drought_hist_std = np.nanstd(all_drought_hist, axis=-1)
    # to 3D multi model mean
    all_drought_ssp585_mean = np.nanmean(all_drought_ssp585, axis=-1)
    # Variance
    all_drought_ssp585_std = np.nanstd(all_drought_ssp585, axis=-1)
    # to 3D multi model mean
    perc_diff = ((all_drought_ssp585_mean - all_drought_hist_mean)
                 / (all_drought_ssp585_mean + all_drought_hist_mean) * 200)
    perc_diff_std = ((all_drought_ssp585_std - all_drought_hist_std)
                 / (all_drought_ssp585_std + all_drought_hist_std) * 200)
    print("all_drought_ssp585_mean")
    print(all_drought_ssp585_mean)
    
    # Histogram
    histo_dict = {}
    histo_dict['data'] = all_spei_hist
    histo_dict['model_kind'] = 'Historic'
    histo_dict['start'] = start_year_hist
    histo_dict['end'] = end_year_hist
    histo_dict['modelnames'] = modelnames
    histo_dict['filename'] = 'Historic_histo'
    histo_dict['bins'] = [-3.0, -2.0, -1.0, 0,
              1, 2, 3]
    plot_histogram_spei_multi(cfg, histo_dict)
    
    
    histo_dict['data'] = all_spei_ssp585
    histo_dict['model_kind'] = 'SSP585'
    histo_dict['start'] = start_year_future
    histo_dict['end'] = end_year_future
    histo_dict['filename'] = 'SSP585_histo'
    plot_histogram_spei_multi(cfg, histo_dict)

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
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'Historic_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 6, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'Historic_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_hist_mean[:, :, 3]
    data_dict['drought_char'] = 'Average SPEI of Events'
    data_dict['filename'] = 'Historic_Average_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # SSP585
    data_dict['data'] = all_drought_ssp585_mean[:, :, 0]
    data_dict['model_kind'] = 'SSP585'
    data_dict['drought_char'] = 'Number of Events per year'
    data_dict['filename'] = 'SSP585_No_of_Events_per_year'
    data_dict['drought_numbers_level'] = np.arange(0, 0.4, 0.05)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_ssp585_mean[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [month]'
    data_dict['filename'] = 'SSP585_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 6, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_ssp585_mean[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events'
    data_dict['filename'] = 'SSP585_Sev_index_of_Events'
    data_dict['drought_numbers_level'] = np.arange(0, 9, 1)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    data_dict['data'] = all_drought_ssp585_mean[:, :, 3]
    data_dict['drought_char'] = 'Average SPEI of Events'
    data_dict['filename'] = 'SSP585_Avr_spei_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-2.8, -1.8, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='gnuplot')

    # SSP585 Percentage difference
    data_dict['data'] = perc_diff[:, :, 0]
    data_dict['datasetname'] = 'Percentage'
    # data_dict['latitude'] = lats
    # data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Difference'
    data_dict['drought_char'] = 'Number of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_No_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-200, 210, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 1]
    data_dict['drought_char'] = 'Duration of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Dur_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-100, 110, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 2]
    data_dict['drought_char'] = 'Severity Index of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Sev_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff[:, :, 3]
    data_dict['drought_char'] = 'Average SPEI of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Avr_of_Events'
    data_dict['drought_numbers_level'] = np.arange(-50, 60, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')
    
    # Standard deviation
    # Historic
    data_dict['data'] = all_drought_hist_std[:, :, 0]
    data_dict['datasetname'] = 'MultiModelStd'
    data_dict['model_kind'] = 'Historic'
    data_dict['drought_char'] = 'Standard deviation of the Number of Events per year'
    data_dict['filename'] = 'Historic_No_of_Events_per_year_std'
    data_dict['drought_numbers_level'] = np.arange(0, 0.2, 0.025)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    data_dict['data'] = all_drought_hist_std[:, :, 1]
    data_dict['drought_char'] = 'Standard deviation of the Duration of Events [month]'
    data_dict['filename'] = 'Historic_Dur_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(0, 3, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    data_dict['data'] = all_drought_hist_std[:, :, 2]
    data_dict['drought_char'] = 'Standard deviation of the Severity Index of Events'
    data_dict['filename'] = 'Historic_Sev_index_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(0, 4, 0.25)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    data_dict['data'] = all_drought_hist_std[:, :, 3]
    data_dict['drought_char'] = 'Standard deviation of the Average SPEI of Events'
    data_dict['filename'] = 'Historic_Average_spei_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(0, 0.5, 0.05)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    # SSP585
    data_dict['data'] = all_drought_ssp585_std[:, :, 0]
    data_dict['model_kind'] = 'SSP585'
    data_dict['drought_char'] = 'Standard deviation of the Number of Events per year'
    data_dict['filename'] = 'SSP585_No_of_Events_per_year_std'
    data_dict['drought_numbers_level'] = np.arange(0, 0.2, 0.025)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    data_dict['data'] = all_drought_ssp585_std[:, :, 1]
    data_dict['drought_char'] = 'Standard deviation of the Duration of Events [month]'
    data_dict['filename'] = 'SSP585_Dur_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(0, 3, 0.2)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    data_dict['data'] = all_drought_ssp585_std[:, :, 2]
    data_dict['drought_char'] = 'Standard deviation of the Severity Index of Events'
    data_dict['filename'] = 'SSP585_Sev_index_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(0, 4, 0.25)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    data_dict['data'] = all_drought_ssp585_std[:, :, 3]
    data_dict['drought_char'] = 'Standard deviation of the Average SPEI of Events'
    data_dict['filename'] = 'SSP585_Avr_spei_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(0, 0.5, 0.05)
    plot_map_spei_multi(cfg, data_dict, colormap='cool')

    # SSP585 Percentage difference
    data_dict['data'] = perc_diff_std[:, :, 0]
    data_dict['datasetname'] = 'Percentage'
    # data_dict['latitude'] = lats
    # data_dict['longitude'] = lons
    data_dict['model_kind'] = 'Difference'
    data_dict['drought_char'] = 'Standard deviation of the Number of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_No_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(-200, 210, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff_std[:, :, 1]
    data_dict['drought_char'] = 'Standard deviation of the Duration of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Dur_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(-200, 210, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff_std[:, :, 2]
    data_dict['drought_char'] = 'Standard deviation of the Severity Index of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Sev_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(-200, 210, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')

    data_dict['data'] = perc_diff_std[:, :, 3]
    data_dict['drought_char'] = 'Standard deviation of the Average SPEI of Events [%]'
    data_dict['filename'] = 'Percentage_difference_of_Avr_of_Events_std'
    data_dict['drought_numbers_level'] = np.arange(-200, 210, 10)
    plot_map_spei_multi(cfg, data_dict, colormap='rainbow')


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
