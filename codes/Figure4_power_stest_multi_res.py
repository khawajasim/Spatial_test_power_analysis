#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:51:36 2021

@author: khawaja
"""
import numpy
import pandas
import matplotlib.pyplot as plt
from csep.core.poisson_evaluations import spatial_test
from csep.core.regions import QuadtreeGrid2D
from csep.core.forecasts import GriddedForecast
from csep.core.catalogs import CSEPCatalog
from csep.utils.time_utils import decimal_year_to_utc_epoch
import utils
#import os

def main():
    n_EQs = 2**numpy.arange(0,10)
    
    #Store above number of earthquakes for each Grid in a Dictionary. 
    grid = {"N100L11":n_EQs,
            "N50L11":n_EQs,
            "N25L11":n_EQs,
            "N10L11":n_EQs,
            "N5L11":n_EQs,
            "N1L11":n_EQs }
    
    
    grid_fn = ['N100L11_power', 'N50L11_power', 'N25L11_power', 'N10L11_power', 'N5L11_power', 'N1L11_power' ]
    zoom = ['N100L11', 'N50L11', 'N25L11', 'N10L11', 'N5L11', 'N1L11']
    mbins = numpy.array([5.95])
    
    run_simulations_again = input('Run Simulations again (Y/N) :')
    
    if run_simulations_again =='Y' or run_simulations_again == 'y':  
    
        n_cat = 100
    
        #os.mkdir('power_stest_grids_GEAR_multi_res')
    
        for z in range(len(grid_fn)):
    
            zoom_level = zoom[z]   
            print('Grid Zoom-level :',zoom_level)
            #    /home/khawaja/GFZ/Testing forecasts
            qk = numpy.genfromtxt('../Data/grids/qk_zoom='+zoom_level+'.txt', dtype='str')
            r = QuadtreeGrid2D.from_quadkeys(qk, magnitudes=mbins)
            #    r = QuadtreeGrid2D.from_regular_grid(zoom_level, magnitudes=mbins)
            r.get_cell_area()
            power_stest_zoom = []
            
            #get into the folder of every N earthquakes.
            N_eqs = grid[str(zoom_level)]
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
    #                print(catalog)
                    
                    #Generate uniform forecast
                    fcst = (r.cell_area/ sum(r.cell_area)) * catalog.event_count
                    fcst = fcst.reshape(-1,1)
                    forecast = GriddedForecast(data = fcst, region = r, magnitudes = mbins, name = zoom_level)
    #                print('sum forecast =', forecast.sum())
                    
                    
                    catalog.region = forecast.region
                    stest = spatial_test(forecast, catalog, verbose=False, seed = 123456)
                    if stest.observed_statistic < numpy.percentile(stest.test_distribution, 5): #2.5
                            stest_fail_N = stest_fail_N+1
    
                power_value_N = stest_fail_N / n_cat
                print('------Power : ', power_value_N)
    #            power_stest_zoom = numpy.row_stack((power_stest_zoom, [N, power_value_N]))
                power_stest_zoom.append(numpy.array([N, power_value_N])) #zoom[z],
            numpy.savetxt('../Data/power_stest/'+grid_fn[z]+'.csv', numpy.array(power_stest_zoom), delimiter=',')
    #        numpy.savetxt('power_stest_grids_GEAR_multi_res/'+grid_fn[z]+'.csv', power_stest_zoom, delimiter=',')
     
    #Generate spatial grid 
    
    print('--Plotting Power for Single-resolution grids')
    power = []
    for fn in grid_fn:
    #    print(fn)
        pp = numpy.loadtxt('../Data/power_stest/'+fn+'.csv', delimiter =',')
        power.append(pp[:,1])
        
    power = numpy.array(power)
    col_eqs = pp[:,0].astype(int)
    row_name = ['N100L11 (922)', 'N50L11 (1780)', 'N25L11 (3502)', 
                  'N10L11 (8089)', 'N5L11 (14782)', 'N1L11 (39811)']
    
    fig, ax = plt.subplots()
    ax.set_xlabel('Number of earthquakes in test catalogs',fontsize=16)
    ax.set_ylabel('Number of cells in multi-resolution grids', fontsize = 16)
#    ax.set_title('Statistical power of Spatial-test for multi-resolution grids', fontsize=22)
    
    im, cbar = utils.heatmap(power, row_name, col_eqs, ax=ax,
                       cmap="YlGn", cbarlabel="Power")
    texts = utils.annotate_heatmap(im, valfmt="{x:.2f}")
    
    fig.set_size_inches(32, 18)
    ax.figure.tight_layout()
    fig.savefig('../Figures/Figure4_power_multi_resolution_grids.png',  bbox_inches='tight', dpi=400)
    
if __name__ == "__main__":
    main()