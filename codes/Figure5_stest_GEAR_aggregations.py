#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:51:36 2021

@author: khawaja
"""
import numpy
import pandas
from csep.core.poisson_evaluations import spatial_test
from csep.core.regions import QuadtreeGrid2D, CartesianGrid2D
from csep.core.forecasts import GriddedForecast
from csep.core.catalogs import CSEPCatalog
from csep.utils.time_utils import decimal_year_to_utc_epoch
from csep.utils.plots import plot_poisson_consistency_test
import matplotlib.pyplot as plt
plt.ioff()

n_L01 = 2**numpy.arange(0,16) #16
#
##Store above number of earthquakes for each Grid in a Dictionary. 
grid = {"L01":n_L01
        }


#folder = ['L5_power', 'L6_power', 'L7_power', 'L8_power', 'L9_power', 'L10_power', 'L11_power' ]
#folder = 'catalogs_gear_forecast'


mbins = numpy.array([5.95])
n_cat = 100

#    Read 0.1x0.1 forecast...

#forecast_path = '/home/khawaja/GFZ/Testing forecasts/global testing experiment/forecasts/GEAR1/'
forecast_path = 'GEAR1/'
# names = ['Uniform', 'WHEEL', 'KJSS', 'SHIFT2F_GSRM', 'TEAMr', 'GEAR1'  ] 
#names = ['5','9','6','7','EQ100L11','EQ50L11','EQ25L11','EQ10L11','EQ5L11','EQ1L11','8','9','10'  ]
names = ['EQ100L11','EQ50L11','EQ25L11','EQ10L11','EQ5L11', 'EQ1L11']
#names = ['EQ1L11']
for name in names:
    print('-----------', name)
    print(forecast_path+'GEAR1_rate_zoom='+name+'.csv')
    data = numpy.loadtxt(forecast_path+'GEAR1_rate_zoom='+name+'.csv', delimiter=',')
#    coords = data[:,:2]
#    qk_path = '/home/khawaja/GFZ/Testing forecasts/global testing experiment/grids/qk_zoom='+name+'.txt'
    qk_path = 'grids/qk_zoom='+name+'.txt'
    print('Quadkeys :', qk_path)
    qk = numpy.genfromtxt(qk_path, dtype = 'str')
    r =QuadtreeGrid2D.from_quadkeys(qk, magnitudes=mbins)        
#    fcst = data[:,4]
    
    fcst = data.reshape(-1,1)
    forecast = GriddedForecast(data = fcst, region = r, magnitudes = mbins, name = name)
    print('Total forecast: ',forecast.sum())
    
    power_stest_zoom = [0.1 , 0, 0]
    
 #get into the folder of every N earthquakes.
    N_eqs = grid['L01']
    for N in N_eqs:
        path = 'catalogs_gear_forecast/N='+str(N)+'/'
        print(path)
        stest_fail_N = 0
        #Counter over the number of catalogs for every N
        for c in range(n_cat):
            data_file = path + 'test_catalog_'+str(c)+'.csv'          
            #Read Catalogs
#            print(data_file)
            dfcat = pandas.read_csv(data_file)

            column_name_mapper = {'index': 'id'}
            # maps the column names to the dtype expected by the catalog class
            dfcat = dfcat.reset_index().rename(columns=column_name_mapper)
            # create the origin_times from decimal years
            dfcat['origin_time'] = dfcat.apply(lambda row: decimal_year_to_utc_epoch(row.year), axis=1) #----
    
            catalog = CSEPCatalog.from_dataframe(dfcat)
#            print(catalog)

            #Generate uniform forecast
#            fcst = (r.cell_area/ sum(r.cell_area)) * catalog.event_count
#            print('sum forecast =', forecast.sum())


            catalog.region = forecast.region
            stest = spatial_test(forecast, catalog, verbose=False, seed = 123456)
#            plot_poisson_consistency_test([stest], 
#                                plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=True, one_sided_lower=True)
#            
#            plt.savefig(name+'_stest_N='+str(N)+'_'+str(c)+'.png', dpi = 250)
#            if stest.observed_statistic < numpy.percentile(stest.test_distribution, 5): #for one_sided_lower=True
            if stest.observed_statistic < numpy.percentile(stest.test_distribution, 2.5): #for one_sided_lower=False (default)
                    stest_fail_N = stest_fail_N+1
                    print('Stest Fail---Catalog # :',c, stest.quantile)
            else:
                    print('Stest Pass---Catalog # :',c, stest.quantile)

        power_value_N = stest_fail_N / n_cat
        print('------Power : ', power_value_N)
#        numpy.savetxt(path+'power_value.txt', [power_value_N])      
        power_stest_zoom = numpy.row_stack((power_stest_zoom, [0.1, N, power_value_N]))
    numpy.savetxt('power_stest_GEAR_zoom='+name+'.csv', power_stest_zoom, delimiter=',')
 
#Generate spatial grid 