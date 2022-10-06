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
#from csep.utils.plots import plot_poisson_consistency_test
import matplotlib.pyplot as plt
#plt.ioff()
import utils

def main():
    n_EQs = 2**numpy.arange(0,12) #16
    mbins = numpy.array([5.95])
    n_cat = 10
    
    #    Read 0.1x0.1 forecast...
    
    #forecast_path = '/home/khawaja/GFZ/Testing forecasts/global testing experiment/forecasts/GEAR1/'
    forecast_path = '../Data/forecasts/'
    names = ['N100L11','N50L11','N25L11','N10L11','N5L11', 'N1L11']

    
    run_simulations_again = input('Run Simulations again (Y/N) :')
    
    if run_simulations_again =='Y' or run_simulations_again =='y':
        results_folder = '../Data/results/'
    else: 
        results_folder = '../Data/generated_results/'
    
    if run_simulations_again =='Y' or run_simulations_again == 'y':  
        
        
        for name in names:
            print('-----------', name)
            print(forecast_path+'GEAR1_rate_zoom='+name+'.csv')
            data = numpy.loadtxt(forecast_path+'GEAR1_rate_zoom='+name+'.csv', delimiter=',')
            qk_path = '../Data/grids/qk_zoom='+name+'.txt'
            print('Grid path :', qk_path)
            qk = numpy.genfromtxt(qk_path, dtype = 'str')
            r =QuadtreeGrid2D.from_quadkeys(qk, magnitudes=mbins)        
    #        fcst = data[:,4]
            
            fcst = data.reshape(-1,1)
            forecast = GriddedForecast(data = fcst, region = r, magnitudes = mbins, name = name)
            print('Total forecast: ',forecast.sum())
            
            power_stest_zoom = []
        
     #get into the folder of every N earthquakes.
    #    N_eqs = grid['L01']
            for N in n_EQs:
                path = '../Data/catalogs_gear_forecast/N='+str(N)+'/'
                print(path)
                stest_fail_N = 0
                #Counter over the number of catalogs for every N
                for c in range(n_cat):
                    data_file = path + 'test_catalog_'+str(c)+'.csv'          
                    #Read Catalogs
    #                print(data_file)
                    dfcat = pandas.read_csv(data_file)
                    
                    column_name_mapper = {'index': 'id'}
                    # maps the column names to the dtype expected by the catalog class
                    dfcat = dfcat.reset_index().rename(columns=column_name_mapper)
                    # create the origin_times from decimal years
                    dfcat['origin_time'] = dfcat.apply(lambda row: decimal_year_to_utc_epoch(row.year), axis=1) #----
                    
                    catalog = CSEPCatalog.from_dataframe(dfcat)
                    catalog.region = forecast.region
                    stest = spatial_test(forecast, catalog, verbose=False, seed = 123456)
    #                plot_poisson_consistency_test([stest], 
    #                                    plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=True, one_sided_lower=True)
    #                
    #                plt.savefig(name+'_stest_N='+str(N)+'_'+str(c)+'.png', dpi = 250)
    #                if stest.observed_statistic < numpy.percentile(stest.test_distribution, 5): #for one_sided_lower=True
                    if stest.observed_statistic < numpy.percentile(stest.test_distribution, 2.5): #for one_sided_lower=False (default)
                            stest_fail_N = stest_fail_N+1
                            print('Stest Fail---Catalog # :',c, stest.quantile)
                    else:
                            print('Stest Pass---Catalog # :',c, stest.quantile)
    
                power_value_N = stest_fail_N / n_cat
                print('------Power : ', power_value_N)
                power_stest_zoom.append(numpy.array([N, power_value_N])) #[zoom[z], 
            numpy.savetxt(results_folder+'GEAR_aggregation_zoom='+name+'.csv', power_stest_zoom, delimiter=',')
     
    #Generate spatial grid 
#    print('--Plotting Power for Single-resolution grids')
    power = []
    for name in names:
    #    print(fn)
        pp = numpy.loadtxt(results_folder+'GEAR_aggregation_zoom='+name+'.csv', delimiter =',')
        power.append(pp[1:,1])
        
    power = (1 - numpy.array(power)) #Calculating Pass-ratio
    col_eqs = pp[1:,0].astype(int)
    row_name = ['N100L11 (922)', 'N50L11 (1780)', 'N25L11 (3502)', 
                  'N10L11 (8089)', 'N5L11 (14782)', 'N1L11 (39811)']
    
    fig, ax = plt.subplots()
    ax.set_xlabel('Number of earthquakes in test catalogs',fontsize=16)
    ax.set_ylabel('Number of cells in multi-resolution grids', fontsize = 16)
#    ax.set_title('S-test Pass ratio for GEAR1 aggregations', fontsize=22)
    
    im, cbar = utils.heatmap(power, row_name, col_eqs, ax=ax,
                       cmap="YlGn", cbarlabel="Power")
    texts = utils.annotate_heatmap(im, valfmt="{x:.3f}") #,textcolors=["blue", "blue"]
    
    fig.set_size_inches(32, 18)
    ax.figure.tight_layout()
    fig.savefig('../Data/Figures/Figure5_GEAR1_aggregations_pass_ratio.png',  bbox_inches='tight', dpi = 400)
    
if __name__ == "__main__":
    main()