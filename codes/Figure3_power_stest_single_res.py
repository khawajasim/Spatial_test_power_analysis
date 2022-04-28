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

#from csep.utils.plots import plot_poisson_consistency_test




n_EQs = 2**numpy.arange(6,16)#numpy.round(numpy.arange(0.06,0.26, 0.001)*1024).astype(int)


#Store above number of earthquakes for each Grid in a Dictionary. 
grid = {"L5":n_EQs,
        "L6":n_EQs,
        "L7":n_EQs,
        "L8":n_EQs,
        "L9":n_EQs,
        "L10":n_EQs,
        "L11":n_EQs
        }

#grid = {
##        "L9":n_EQs,
##        "L10":n_EQs,
#        "L11":n_EQs
#        }

grid_fn = ['L5_power', 'L6_power', 'L7_power', 'L8_power', 'L9_power', 'L10_power', 'L11_power' ]
#grid_fn = ['L11_power' ] #'L9_power', 'L10_power', 

zoom = [5, 6, 7, 8, 9, 10, 11]
#zoom = [11] #9, 10, 
mbins = numpy.array([5.95])
n_cat = 2 #100
#power_stest_zoom = [] #[0,0,0]

for z in range(len(grid_fn)):

    zoom_level = zoom[z]   
    print('Grid Zoom-level :',zoom_level)
    r = QuadtreeGrid2D.from_single_resolution(zoom_level, magnitudes=mbins)
    r.get_cell_area()
    
    power_stest_zoom = []
    #get into the folder of every N earthquakes.
    N_eqs = grid['L'+str(zoom_level)]
    for N in N_eqs:
        path = '../Data/catalogs_gear_forecast'+'/N='+str(N)+'/'
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
            catalog.filter_spatial
#            print(catalog)
            catalog.filter('latitude < ' + str(85.05))
            catalog.filter('latitude > ' + str(-85.05))
            print('Events after Filter 85.05 :', catalog.event_count)
            #Generate uniform forecast
            fcst = (r.cell_area/ sum(r.cell_area)) * catalog.event_count
            fcst = fcst.reshape(-1,1)
            forecast = GriddedForecast(data = fcst, region = r, magnitudes = mbins, name = 'L5')
#            print('sum forecast =', forecast.sum())


            catalog.region = forecast.region
            stest = spatial_test(forecast, catalog, verbose=False, seed = 123456)
            if stest.observed_statistic < numpy.percentile(stest.test_distribution, 5):
                    print('S-test Failed')
                    stest_fail_N = stest_fail_N+1
            
        power_value_N = stest_fail_N / n_cat
        print('------Power : ', power_value_N)
        power_stest_zoom.append(numpy.array([zoom[z], N, power_value_N]))
#        power_stest_zoom = numpy.row_stack((power_stest_zoom, [zoom[z], N, power_value_N]))
    numpy.savetxt('../Data/power_stest/'+grid_fn[z]+'.csv', numpy.array(power_stest_zoom), delimiter=',')
 
#Generate spatial grid 
