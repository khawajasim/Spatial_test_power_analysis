#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:16:30 2021

@author: khawaja
"""
import pickle
import numpy
import pandas
# from csep.core.poisson_evaluations import spatial_test, number_test, conditional_likelihood_test, magnitude_test
from csep.utils.plots import plot_poisson_consistency_test
from csep.core.regions import CartesianGrid2D 
from csep.core.catalogs import CSEPCatalog
from csep.utils.time_utils import decimal_year_to_utc_epoch
from csep.core.forecasts import GriddedForecast
import datetime
import matplotlib.pyplot as plt
from csep.core.regions import QuadtreeGrid2D
from utils import binomial_spatial_test

def main():
    
    run_evaluations_again = input('Run Evaluations again (Y/N) :')
    
    if run_evaluations_again =='Y' or run_evaluations_again =='y':
        results_folder = '../Data/results/'
    else: 
        results_folder = '../Data/generated_results/'
        
    
    plot_normalize= True
    mbins = numpy.array([5.95])
    #forecast_name = 'forecasts/helmstetter/helmstetter_' #, #'forecasts/wheel/wheel_' #'forecasts/uniform/uniform_' 
    catalog_name = '../Data/cat_test.csv'
    
    #------Read Tono Forecast and Convert it into CSEP format
    
    
    
    #--------------- Just Spatial counts ....
    folder = '../Data/forecasts/' #, #'forecasts/wheel/wheel_' #'forecasts/uniform/uniform_' 
    names = ['GEAR1','KJSS','SHIFT_GSRM2f','TEAM','WHEEL', 'Uniform']     # ,  
    # result_folder = '../Data/power_stest/'
      
    if run_evaluations_again =='Y' or run_evaluations_again =='y':  
        print('Step 1: Use Toolkit Directly for forecast Evaluation')
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
        catalog = CSEPCatalog.from_dataframe(dfcat)
        print(catalog)
        for name in names:
            binary_stest_all = []
           
            print('------     ',name)
            data = numpy.loadtxt(folder+name+'_rate_zoom=01.csv', delimiter=',')
            coords = data[:,:2]
            forecast = data[:,4]
            r = CartesianGrid2D.from_origins(coords, magnitudes=mbins)
            
            
            forecast = forecast.reshape(-1,1)
            
            forecast_gridded = GriddedForecast(data = forecast, region = r, magnitudes = mbins, name='0.1x0.1 ('+str(len(coords))+')' )
            print(f"expected event count before scaling: {forecast_gridded.event_count}")
            forecast_gridded.scale(6)
            print(f"expected event count after scaling: {forecast_gridded.event_count}")
            
            #Binding region (spatial + magnitude)
            catalog.region = forecast_gridded.region
            
            stest_01 = binomial_spatial_test(forecast_gridded, catalog, verbose=False, seed = 123456)
            
            print(stest_01.observed_statistic)
            print(stest_01.quantile)
            binary_stest_all.append(stest_01)
        #Out[25]: -4556.393591179189
            
            #plot_poisson_consistency_test([stest_01], 
        #                                plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, normalize=False, one_sided_lower=True)
        
        
            del forecast_gridded, forecast, r, data, coords
        
        #-----------Wheel Zoom = 11
            zoom = ['11','10', '9', '8','7','6','5','4','3','2','1','N100L11', 'N50L11', 'N25L11', 'N10L11', 'N5L11', 'N1L11' ] #
            # zoom = ['5','4']
            for z in zoom:
                print('Zoom level :',z)
                qk = numpy.genfromtxt('../Data/grids/qk_zoom='+z+'.txt', dtype = 'str')
                r = QuadtreeGrid2D.from_quadkeys(qk, magnitudes=mbins)
        
                forecast = numpy.loadtxt(folder+name+'_rate_zoom='+z+'.csv', delimiter=',')
                #Reshape forecast as Nx1 array
                forecast = forecast.reshape(-1,1)     
        
                forecast_gridded = GriddedForecast(data = forecast, region = r, magnitudes = mbins, name = z+' ('+str(len(qk))+')' )
                print(f"expected event count before scaling: {forecast_gridded.event_count}")
                forecast_gridded.scale(6)
                print(f"expected event count after scaling: {forecast_gridded.event_count}")
                catalog.region = forecast_gridded.region
                
                stest = binomial_spatial_test(forecast_gridded, catalog, verbose=False, seed = 123456)
                
                
                print(stest.observed_statistic)
                print(stest.quantile)
                binary_stest_all.append(stest)
                
           
                del forecast_gridded, forecast, r
            
            ax = plot_poisson_consistency_test(binary_stest_all, 
                                        normalize=plot_normalize, one_sided_lower=True)
            ax.set_title(name+': Binary S-test')
            
            ticks = numpy.array(plt.yticks()[0])
            labels_left = ['0.1$^\circ$x0.1$^\circ$ (6480000)', 'L11 (4194304)','L10 (1048576)', 'L9 (262144)', 'L8 (65536)',
                        'L7 (16384)','L6 (4096)','L5 (1024)', 'L4 (256)', 'L3 (64)', 'L2 (16)', 'L1 (4)',
                        'N100L11 (922)', 'N50L11 (1780)', 'N25L11 (3502)',
                        'N10L11 (8089)', 'N5L11 (14782)', 'N1L11 (39811)' ]
            labels_left = numpy.flip(labels_left)
            ax.set_yticklabels(labels_left )
            ax.set_ylabel('Test grid and number of cells in each grid')
            
            if plot_normalize == True:
                ax.set_xlabel('Simulated log-likelihood - Observed log-likelihood')
                ax.figure.tight_layout()
                ax.figure.savefig('../Data/Figures/binary_stest_'+name+'_normalize.png', dpi = 250)
                ax.figure.savefig('../Data/Figures/binary_stest_'+name+'_normalize.svg')
            else:
                ax.set_xlabel('Log-likelihood')
                ax.figure.tight_layout()
                ax.figure.savefig('../Data/Figures/binary_stest_'+name+'.png', dpi = 250)
                ax.figure.savefig('../Data/Figures/binary_stest_'+name+'.svg')
                
            # with open(folder+name+'/binary_stest_'+name+'.dat', 'wb') as f:
            with open(results_folder+'binary_stest_'+name+'.dat', 'wb') as f:
                pickle.dump(binary_stest_all, f)


    else:
        for name in names:
            with open(results_folder+'binary_stest_'+name+'.dat', 'rb') as f:  #binary_stest_WHEEL.dat
                reordered_results = pickle.load(f)

            ax = plot_poisson_consistency_test(reordered_results, 
                                     normalize=plot_normalize, one_sided_lower=True)
            ax.set_title(name+': Binary S-test')
            ticks = numpy.array(plt.yticks()[0])

            labels_left = ['0.1$^\circ$x0.1$^\circ$ (6480000)', 'L11 (4194304)','L10 (1048576)', 'L9 (262144)', 'L8 (65536)',
                     'L7 (16384)','L6 (4096)','L5 (1024)', 'L4 (256)', 'L3 (64)', 'L2 (16)', 'L1 (4)',
                     'N100L11 (922)', 'N50L11 (1780)', 'N25L11 (3502)',
                     'N10L11 (8089)', 'N5L11 (14782)', 'N1L11 (39811)' ]
            labels_left = numpy.flip(labels_left)
            ax.set_yticklabels(labels_left )
            ax.set_ylabel('Test grid and number of cells in each grid')
             # plot_args={'xlabel': 'Log-Likelihood', 'title': 'S-Test'}, 
            if plot_normalize == True:
                 
                 ax.set_xlabel('Simulated log-likelihood - Observed log-likelihood')
                 ax.figure.tight_layout()
                 ax.figure.savefig('../Data/Figures/binary_stest_'+name+'_normalize.png', dpi = 250)
                 ax.figure.savefig('../Data/Figures/binary_stest_'+name+'_normalize.svg')
        
if __name__ == "__main__":
    main() 


        
# #----------------------
# ----Load saved results and improve plots
# import pickle    
# from csep.utils.plots import plot_poisson_consistency_test
# import matplotlib.pyplot as plt
# import numpy

# folder = '../Data/forecasts/' #, #'forecasts/wheel/wheel_' #'forecasts/uniform/uniform_' 
# names = ['WHEEL', 'Uniform', 'GEAR1', 'KJSS', 'SHIFT2F_GSRM', 'TEAMr' ]     # 'Uniform' wheel  'uniform', 

# #name = 'WHEEL'
# for name in names:
#     # folder+name+'/binary_stest_modified'+name+'.dat'
#     with open(folder+name+'/binary_stest_'+name+'.dat', 'rb') as f:
#         reordered_results = pickle.load(f)


# ----Load saved results and improve plots
# import pickle    
# from csep.utils.plots import plot_poisson_consistency_test
# import matplotlib.pyplot as plt
# import numpy
#-----Select the grids we want to plot

# plot_normalize = True


# ax = plot_poisson_consistency_test(reordered_results, 
#                                     normalize=plot_normalize, one_sided_lower=True)
# ax.set_title(name+': Binary S-test')
        
# ticks = numpy.array(plt.yticks()[0])
# labels_left = ['0.1$^\circ$x0.1$^\circ$ (6480000)', 'L11 (4194304)','L10 (1048576)', 'L9 (262144)', 'L8 (65536)',
#                     'L7 (16384)','L6 (4096)','L5 (1024)', 'L4 (256)', 'L3 (64)', 'L2 (16)', 'L1 (4)',
#                     'N100L11 (922)', 'N50L11 (1780)', 'N25L11 (3502)',
#                     'N10L11 (8089)', 'N5L11 (14782)', 'N1L11 (39811)' ]
# ticks_left = numpy.flip(ticks)
# ax.set_yticklabels(ticks_left,labels_left )
# ax.set_ylabel('Test grid and number of cells in each grid')

# if plot_normalize == True:
#     ax.set_xlabel('Simulated log-likelihood - Observed log-likelihood')
#     ax.figure.savefig(folder+'binary_stest_'+name+'allgrids_normalize.png', dpi = 250)
# else:
#     ax.set_xlabel('Log-likelihood')
#     ax.figure.savefig(folder+'binary_stest_'+name+'allgrids.png', dpi = 250)



# ---------

#   ----------------------
#     # ----Load saved results and improve plots
#     import pickle    
#     from csep.utils.plots import plot_poisson_consistency_test
#     import matplotlib.pyplot as plt
#     import numpy
#     folder = '../Data/forecasts/'
#     plot_normalize =True
#     names = ['Uniform', 'GEAR1', 'KJSS', 'SHIFT_GSRM2f', 'TEAM', 'WHEEL' ]     # 'Uniform' wheel  'uniform', 
    
#     #name = 'WHEEL'
