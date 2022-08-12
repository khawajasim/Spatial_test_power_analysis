#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 17:59:53 2022

@author: khawaja
"""
import Figure1_global_stest
import Figure2_power_analysis_showcase
import Figure3_power_stest_single_res
import Figure4_power_stest_multi_res
import Figure5_stest_GEAR_aggregations
import Figure6_forecast_stest_poisson

#Call all plots
print('Generating Figure 1')
#Figure1_global_stest.main()

print('Generating Figure 2')
Figure2_power_analysis_showcase.main()

print('Generating Figure 3')
Figure3_power_stest_single_res.main()

print('Generating Figure 4')
Figure4_power_stest_multi_res.main()

print('Generating Figure 5')
Figure5_stest_GEAR_aggregations

print('Generating Figure 6')
Figure6_forecast_stest_poisson.main()
