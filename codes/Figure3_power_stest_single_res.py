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

def main():
    n_EQs = 2**numpy.arange(6,16)
    
    
    #Store above number of earthquakes for each Grid in a Dictionary. 
    grid = {"L5":n_EQs,
            "L6":n_EQs,
            "L7":n_EQs,
            "L8":n_EQs,
            "L9":n_EQs,
            "L10":n_EQs,
            "L11":n_EQs
            }
    
    grid_fn = ['L5_power', 'L6_power', 'L7_power', 'L8_power', 'L9_power', 'L10_power', 'L11_power' ]
    #grid_fn = ['L11_power' ] #'L9_power', 'L10_power', 
    
    zoom = [5, 6, 7, 8, 9, 10, 11]
    #zoom = [11] #9, 10, 
    mbins = numpy.array([5.95])
    
    
    run_simulations_again = input('Run Simulations again (Y/N) :')
    
    if run_simulations_again =='Y' or run_simulations_again == 'y':    
        n_cat = 2 #00   
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
    #                print('Events after Filter 85.05 :', catalog.event_count)
                    #Generate uniform forecast
                    fcst = (r.cell_area/ sum(r.cell_area)) * catalog.event_count
                    fcst = fcst.reshape(-1,1)
                    forecast = GriddedForecast(data = fcst, region = r, magnitudes = mbins, name = 'L5')
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
        
    print('--Plotting Power for Single-resolution grids')
    power = []
    for fn in grid_fn:
    #    print(fn)
        pp = numpy.loadtxt('../Data/power_stest/'+fn+'.csv', delimiter =',')
        power.append(pp[:,2])
        
    power = numpy.array(power)
    col_eqs = pp[:,1].astype(int)
    row_name = numpy.array(['L5 (1024)', 'L6 (4096)', 'L7 (16384)', 'L8 (65536)', 
                                'L9 (262144)', 'L10 (1048576)', 'L11 (4194304)'])
        
    fig, ax = plt.subplots()
    ax.set_xlabel('Number of earthquakes in test catalogs',fontsize=16)
    ax.set_ylabel('Number of cells in single-resolution grids', fontsize = 16)
#    ax.set_title('Statistical power of Spatial-test for single-resolution grids', fontsize=22)
    
    im, cbar = utils.heatmap(power, row_name, col_eqs, ax=ax,
                       cmap="YlGn", cbarlabel="Power")
    texts = utils.annotate_heatmap(im, valfmt="{x:.2f}")
    
    fig.set_size_inches(32, 18)
    ax.figure.tight_layout()
    fig.savefig('../Figures/Figure3_power_single_resolution_grids.png',  bbox_inches='tight', dpi =400)
    
    
if __name__ == "__main__":
    main()