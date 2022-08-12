#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:16:30 2021

@author: khawaja
"""
import pickle
import numpy
import pandas
from csep.core.poisson_evaluations import spatial_test, number_test, conditional_likelihood_test, magnitude_test
from csep.utils.plots import plot_poisson_consistency_test
from csep.core.regions import CartesianGrid2D 
from csep.core.catalogs import CSEPCatalog
from csep.utils.time_utils import decimal_year_to_utc_epoch
from csep.core.forecasts import GriddedForecast
import datetime
import matplotlib.pyplot as plt
from csep.core.regions import QuadtreeGrid2D


def main():
    mbins = numpy.array([5.95])
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
    
    
    #--------------- Just Spatial counts ....
    folder = '../Data/forecasts/' #, #'forecasts/wheel/wheel_' #'forecasts/uniform/uniform_' 
    names = ['WHEEL', 'Uniform','GEAR1', 'KJSS', 'SHIFT2F_GSRM', 'TEAM' ]     # ,  
    result_folder = '../Data/power_stest/'
    
    for name in names:
        stest_all = []
    
        print('------     ',name)
        data = numpy.loadtxt(folder+name+'_rate_zoom=01.csv', delimiter=',')
        coords = data[:,:2]
        forecast = data[:,4]
        r = CartesianGrid2D.from_origins(coords, magnitudes=mbins)
        
        
        forecast = forecast.reshape(-1,1)
        
        forecast_gridded = GriddedForecast(data = forecast, region = r, magnitudes = mbins, name=name+' - 0.1x0.1' )
        print(f"expected event count before scaling: {forecast_gridded.event_count}")
        forecast_gridded.scale(6)
        print(f"expected event count after scaling: {forecast_gridded.event_count}")
        
        #Binding region (spatial + magnitude)
        catalog.region = forecast_gridded.region
        
        stest_01 = spatial_test(forecast_gridded, catalog, verbose=False, seed = 123456)
        
        print(stest_01.observed_statistic)
        print(stest_01.quantile)
        stest_all.append(stest_01)
    #Out[25]: -4556.393591179189
        
        #plot_poisson_consistency_test([stest_01], 
    #                                plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=False, one_sided_lower=True)
    
    
        del forecast_gridded, forecast, r, data, coords
    
    #-----------Wheel Zoom = 11
        zoom = ['11','10', '9', '8','7','6','5','N100L11', 'N50L11', 'N25L11', 'N10L11', 'N5L11', 'N1L11' ]  #
        for z in zoom:
            print('Zoom level :',z)
            qk = numpy.genfromtxt('../Data/grids/qk_zoom='+z+'.txt', dtype = 'str')
            r = QuadtreeGrid2D.from_quadkeys(qk, magnitudes=mbins)
    
            forecast = numpy.loadtxt(folder+name+'_rate_zoom='+z+'.csv', delimiter=',')
            #Reshape forecast as Nx1 array
            forecast = forecast.reshape(-1,1)     
    
            forecast_gridded = GriddedForecast(data = forecast, region = r, magnitudes = mbins, name = name+ ' Grid: '+z)
            print(f"expected event count before scaling: {forecast_gridded.event_count}")
            forecast_gridded.scale(6)
            print(f"expected event count after scaling: {forecast_gridded.event_count}")
            catalog.region = forecast_gridded.region
            
            stest = spatial_test(forecast_gridded, catalog, verbose=False, seed = 123456)
            
            
            print(stest.observed_statistic)
            print(stest.quantile)
            stest_all.append(stest)
            
       
            del forecast_gridded, forecast, r
        
        plot_poisson_consistency_test(stest_all, 
                                    plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=True, one_sided_lower=True)
        plt.savefig(result_folder+'stest_'+name+'_normalize.png', dpi = 250)
    #    plot_poisson_consistency_test(stest_all, 
    #                                plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=False, one_sided_lower=True)
    #    plt.savefig(folder+name+'/stest_'+name+'.png', dpi = 250)
        with open(result_folder+'stest_'+name+'.dat', 'wb') as f:
            pickle.dump(stest_all, f)
    
if __name__ == "__main__":
    main()   
    
    
            
    #----------------------
    #----Load saved results and improve plots
    #import pickle    
    #from csep.utils.plots import plot_poisson_consistency_test
    #import matplotlib.pyplot as plt
    #import numpy
    # 
    #names = ['WHEEL', 'Uniform', 'GEAR1', 'KJSS', 'SHIFT2F_GSRM', 'TEAM' ]     # 'Uniform' wheel  'uniform', 
    #
    ##name = 'WHEEL'
    #for name in names:
    #    with open(result_folder+'/stest_'+name+'.dat', 'rb') as f:
    #         reordered_results = pickle.load(f)
    #    ax =  plot_poisson_consistency_test(reordered_results, 
    #                                plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=True, one_sided_lower=True)
    #    
    #    ticks = numpy.array(plt.yticks()[0])
    #    
    ##    labels_left = ['6480000' , str(4**11), str(4**10), str(4**9), str(4**8), str(4**7), str(4**6), str(4**5), '922', '1780', '3502', '8089', '14782','39811']
    #    labels_left = ['0.1$^\circ$x0.1$^\circ$ (6480000)', 'L11 (4194304)','L10 (1048576)', 'L9 (262144)', 'L8 (65536)','L7 (16384)','L6 (4096)','L5 (1024)',
    #                   'N100L11 (922)', 'N50L11 (1780)', 'N25L11 (3502)', 'N10L11 (8089)', 'N5L11 (14782)', 'N1L11 (39811)' ]
    #  
    #    ticks_left = numpy.flip(ticks)
    #    plt.yticks(ticks_left,labels_left )
    #    
    #    plt.ylabel('Test grid and number of cells in each grid')
    #    plt.title(name+' : S-Test')
    #    plt.tight_layout()
    #    plt.savefig(folder+name+'/stest_'+name+'_normalized.png', dpi = 250)
    
    
