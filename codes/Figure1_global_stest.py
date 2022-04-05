#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:11:12 2022

@author: khawaja
"""
import numpy
import pickle
import pandas
from csep.core.regions import CartesianGrid2D
from csep.core.forecasts import GriddedForecast
from csep.core.catalogs import CSEPCatalog
from csep.utils.time_utils import decimal_year_to_utc_epoch
from csep.core.poisson_evaluations import spatial_test
from csep.utils.plots import plot_poisson_consistency_test

mbins = numpy.array([5.95])
stest_all = []
folder = '../Data/forecasts/' #, #'forecasts/wheel/wheel_' #'forecasts/uniform/uniform_' 
names = ['KJSS', 'SHIFT2F_GSRM', 'WHEEL', 'TEAM', 'GEAR1', 'Uniform']     # 'Uniform' 
#names = ['uniform']


mbins = numpy.array([5.95])
#forecast_name = 'forecasts/helmstetter/helmstetter_' #, #'forecasts/wheel/wheel_' #'forecasts/uniform/uniform_' 
catalog_name = '../Data/cat_test.csv'

#------Read Tono Forecast and Convert it into CSEP format
#----Read Catalog format
dfcat = pandas.read_csv(catalog_name)

column_name_mapper = {
    'lon': 'longitude',
    'lat': 'latitude',
    'mag': 'magnitude'
    }  #'index': 'id'

# maps the column names to the dtype expected by the catalog class
dfcat = dfcat.reset_index().rename(columns=column_name_mapper)
# create the origin_times from decimal years
dfcat['origin_time'] = dfcat.apply(lambda row: decimal_year_to_utc_epoch(row.year), axis=1) #----
# add depth
#df['depth'] = 5
# create catalog from dataframe
catalog = CSEPCatalog.from_dataframe(dfcat)
print(catalog)


for name in names:
    data = numpy.loadtxt(folder+name+'/'+name+'_rate_zoom=01.csv', delimiter=',')
    coords = data[:,:2]
    forecast = data[:,4]
    print('Forecast Name:', name)
    print('Forecast / year:',sum(forecast))
    
    r = CartesianGrid2D.from_origins(coords, magnitudes=mbins)
    
    
    forecast = forecast.reshape(-1,1)
    
    forecast_gridded = GriddedForecast(data = forecast, region = r, magnitudes = mbins, name=name)
    print(f"expected event count before scaling: {forecast_gridded.event_count}")
    forecast_gridded.scale(6)
    #print(f"expected event count after scaling: {forecast_gridded.event_count}")
    
    #Binding region (spatial + magnitude)
    catalog.region = forecast_gridded.region
    
    stest_01 = spatial_test(forecast_gridded, catalog, verbose=False, seed = 123456)
#    
#    print(stest_01.observed_statistic)
#    print(stest_01.quantile)
    stest_all.append(stest_01)
    #Out[25]: -4556.393591179189
#    ntest = number_test(forecast_gridded, catalog)
#    ntest_all.append(ntest)
    
   
    del forecast_gridded, forecast, r, data, coords


plot_poisson_consistency_test(stest_all, 
                                plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=True, one_sided_lower=True)


#with open(folder+'ntest_01_all.dat', 'wb') as f:
#        pickle.dump(ntest_all, f)

#with open(folder+'stest_01_all.dat', 'wb') as f:
#        pickle.dump(stest_all, f)
    
with open(folder+'/stest_01_all.dat', 'rb') as f:
         reordered_results = pickle.load(f)
#         