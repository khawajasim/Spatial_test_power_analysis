#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:51:36 2021

@author: khawaja
"""
import numpy
import pandas
from csep.core.poisson_evaluations import spatial_test
from csep.core.regions import QuadtreeGrid2D
from csep.core.forecasts import GriddedForecast
from csep.core.catalogs import CSEPCatalog
from csep.utils.time_utils import decimal_year_to_utc_epoch
import os





n_EQs = 2**numpy.arange(0,10)#numpy.round(numpy.arange(0.06,0.26, 0.001)*1024).astype(int)
#n_EQ50 = 2**numpy.arange(0,10)
#n_EQ25 = 2**numpy.arange(0,10)
#n_EQ10 = 2**numpy.arange(0,10)
#n_EQ5 = 2**numpy.arange(0,10)
#n_EQ1 = 2**numpy.arange(0,10)
#n_L11 = 2**numpy.arange(6,16)

#Store above number of earthquakes for each Grid in a Dictionary. 
grid = {"EQ100L11":n_EQs,
        "EQ50L11":n_EQs,
        "EQ25L11":n_EQs,
        "EQ10L11":n_EQs,
        "EQ5L11":n_EQs,
        "EQ1L11":n_EQs }


grid_fn = ['EQ100L11_power', 'EQ50L11_power', 'EQ25L11_power', 'EQ10L11_power', 'EQ5L11_power', 'EQ1L11_power' ]
zoom = ['EQ100L11', 'EQ50L11', 'EQ25L11', 'EQ10L11', 'EQ5L11', 'EQ1L11']
mbins = numpy.array([5.95])
n_cat = 100

os.mkdir('power_stest_grids_GEAR_multi_res')

for z in range(len(grid_fn)):

    zoom_level = zoom[z]   
    print('Grid Zoom-level :',zoom_level)
#    /home/khawaja/GFZ/Testing forecasts
    qk = numpy.genfromtxt('grids/qk_zoom='+zoom_level+'.txt', dtype='str')
    r = QuadtreeGrid2D.from_quadkeys(qk, magnitudes=mbins)
#    r = QuadtreeGrid2D.from_regular_grid(zoom_level, magnitudes=mbins)
    r.get_cell_area()
    power_stest_zoom = [0,0]

    #get into the folder of every N earthquakes.
    N_eqs = grid[str(zoom_level)]
    for N in N_eqs:
        path = 'catalogs_gear_forecast'+'/N='+str(N)+'/'
        print(path)
        stest_fail_N = 0
        #Counter over the number of catalogs for every N
        for c in range(n_cat):
            data_file = path + 'test_catalog_'+str(c)+'.csv'          
            #Read Catalogs
            dfcat = pandas.read_csv(data_file)

            column_name_mapper = {'index': 'id'}
            # maps the column names to the dtype expected by the catalog class
            dfcat = dfcat.reset_index().rename(columns=column_name_mapper)
            # create the origin_times from decimal years
            dfcat['origin_time'] = dfcat.apply(lambda row: decimal_year_to_utc_epoch(row.year), axis=1) #----
    
            catalog = CSEPCatalog.from_dataframe(dfcat)
#            print(catalog)

            #Generate uniform forecast
            fcst = (r.cell_area/ sum(r.cell_area)) * catalog.event_count
            fcst = fcst.reshape(-1,1)
            forecast = GriddedForecast(data = fcst, region = r, magnitudes = mbins, name = zoom_level)
#            print('sum forecast =', forecast.sum())


            catalog.region = forecast.region
            stest = spatial_test(forecast, catalog, verbose=False, seed = 123456)
            if stest.observed_statistic < numpy.percentile(stest.test_distribution, 2.5):
                    stest_fail_N = stest_fail_N+1

        power_value_N = stest_fail_N / n_cat
        print('------Power : ', power_value_N)
        power_stest_zoom = numpy.row_stack((power_stest_zoom, [N, power_value_N]))
    numpy.savetxt('power_stest_grids_GEAR_multi_res/'+grid_fn[z]+'.csv', power_stest_zoom, delimiter=',')
 
#Generate spatial grid 
