#import os
#import sys
import scipy.stats
import numpy as np
from scipy.special import erf  #See what this is doing?
from scipy.special import erfinv
import matplotlib.pyplot as plt
import matplotlib

def simulation(model, N):
    """
    Generate simiulated gridded catalogs
    Inputs: 
            Model : probabalistic model
            N : Number of earthquake in the simulated catalog
    Output:
            nobs: Simulated gridded catalog
    """
    nobs = np.zeros(len(model))
    nn = np.random.choice(len(model), N, p=model/np.sum(model))
    ii, ni = np.unique(nn, return_counts=True)
    for k, i in enumerate(ii):
        nobs[i] = ni[k]
    return nobs

def uniformmodel(N, binedges):  #Explore. !!!!!
    dx = binedges[1:] - binedges[:-1]
    res = N * dx / np.sum(dx)
    return res

def gaussweight(x0, x1, mu, sig):
    """
    A gussian distribution is defined by its mu and sig.
    It takes the the EDGES (x0 and x1) and returns the Percentage of data contained in these two bins.
    
    OR: Assigns Gaussian weights to the Bins Edges
    Inputs: 
        x0: Starting Edge of bin
        x1: Ending Edge of bin
        mu: Gaussian Parameter for Mean ?
        Sig: Gaussian parameter for Sigma ? 
    """
    
    y0 = (x0 - mu) / (np.sqrt(2.0) * sig)
    y1 = (x1 - mu) / (np.sqrt(2.0) * sig)
    res = 0.5 * (erf(y1) - erf(y0))
    return res
    
def gaussmodel(N, binedges, mu, sig):
    """
    Generates a model for a spatial grid using Gaussian distribution.
    Calculates percentage of data contained each spatial bin for a hypothetical grid.
    
    Inputs: 
        N: Number of Earthquakes to generate model
        binedges: Bin edges of (hypothetical) Spatial grid
        mu: Gaussian Parameter for Mean ?
        Sig: Gaussian parameter for Sigma ?
    
    """
    #Percentalge of data contained in whole hypothetical grid. Should be =1.
    norm = gaussweight(binedges[0], binedges[-1], mu, sig)  
    
    res = []
    for i in range(len(binedges)-1):
        res = np.append(res, N * gaussweight(binedges[i], binedges[i+1], mu, sig) / norm)
    return res


def gaussian_binedges(Z, Nc2, mu, sig):
    """
    Generates Bin Edges based on Gaussian distribution. 
    The philosopy is that each  bin should receive equal percentage of data as determined by
    Gaussian Distribution (mu & Sig). Increase the bin size on Edges, and decrease in center.
    
    Inputs:
        Z: Determines the the last Edge of bin on positive and nagative side.
        Nc2: Total number of Bins on Positive and Negative side. That is why we pass 2*Nc2 in when used.
        mu: Mean of Gaussian Distribution
        sig: Std of Gaussian distribution
    """
    norm = 2 * gaussweight(0, Z, 0, sig) #*sig
    yi = np.linspace(0, norm, Nc2+1)
    xi = erfinv(yi) * np.sqrt(2.0) * sig
    if np.isinf(xi[-1]):
        xi[-1] = Z #*sig
    res = mu + xi
    res = np.append(res, mu - xi)
    res = np.sort(np.unique(res))
    return res

def LLvalue(model, nobs):
    LL = np.sum(np.log(scipy.stats.poisson.pmf(nobs, model)))
    return LL

def calculate_power(model0, model1, N, Nsim):
    """
    Testing model0.
    Model1 is only given for generating Observed Catalog
    """
    power = 0
    for nsim in range(Nsim):
        #Generate observed gridded catalog using model1 with N earthquakes
        nobs = simulation(model1, N)
        
        #Compute Observed LL-value using observed gridded catalog generated above with model0
        LLobs = LLvalue(model0, nobs)
        LL = []
        for nsim2 in range(Nsim):
            #Use model0 to generate test simulations [inside of S-test] 
            nsim = simulation(model0, N)
            LL = np.append(LL, LLvalue(model0, nsim))
        #Lower limit of Simulated LL interval 
        LL0 = np.percentile(LL, 2.5)    #2.5 or 5. Depends, If its two-sided Or One-sided
        
        #Compute power: (S-test failing correctly) / Nsim
        if LLobs < LL0:
            power += 1.0/(1.0*Nsim)
    return power



def main():
    # -- Parameters
    mu = 0.0
    sigma_1 = 1
    sigma_2=0.5
    Nsim = 1000
    
    #1 # Array, in case if running code of multiple sigma.
    # Use Z as three times of sigma to Cover 99.7% of the data it has to offer.
    Nc2_all = Nc2_all = np.arange(1,51,1)
    Neqs_fix = 10 #Use this number of earthquakes to plot trend (Frame c)
    Ncell_fix = 20 #Use this number of cells to plot trend (Frame d)
    Z = 3 #*sigma #10.0
    print('Standard Deviation :', sigma_1)
    gridtype = ['Uniform grid', 'Density grid']
    
    run_sim = input("Run simulation again ?  (Y / N):  ")
#    Nsim = input("Choose the number of simulations (By Default 1000): ")
#    Nsim = int(Nsim)
    sigmas = [sigma_1, sigma_2]
    
    #Results folder name, where to look for results.
    if run_sim =='Y' or run_sim =='y':
        results_folder = '../Data/results/power_analysis_showcase/'
    else: 
        results_folder = '../Data/generated_results/power_analysis_showcase/'
        
        
    if run_sim == 'Y' or run_sim == 'y':   
        for sigma in sigmas:
            for Nc2 in Nc2_all:
                print('-----Number of Cells = ', 2*Nc2)
                uniformbinedges = np.linspace(-Z, Z, 2 * Nc2 + 1) 
                print('Uniform grid size :', len(uniformbinedges)-1) 
                densitybinedges = gaussian_binedges(Z, Nc2, mu, sigma_1) 
                
                outname = results_folder+'powertest-Z%.0f-Ncell%.0f-Nsim%d-sigma%.1f.out' % (Z, 2*Nc2, Nsim,sigma)
                f = open(outname, 'w')
                f.write('# N    power_uniformgrid   power_densitygrid\n')
                f.close()
                #Below conditional Nc2 assignment is to reduce time for unnecessary calculation.
                #Doing only those calculations, that we want to show in figure. 
                if Nc2 ==  (Ncell_fix/2):  #10: #10*2= Ncells
                    Neqs = np.arange(1,101,1) 
                else:
                    Neqs = [Neqs_fix]
        #        print('Total Earthquakes :', Neqs)
                for N in  Neqs: 
        #            print('Number of Cells :', Nc2)
        #            print('---EQS :,', N)
                    for ng, grid in enumerate(gridtype):
                        if grid == 'Uniform grid':
                            
                            model0 = uniformmodel(N, uniformbinedges)
                            model1 = gaussmodel(N, uniformbinedges, mu, sigma)
                            poweruniform = calculate_power(model0, model1, N, Nsim)
                        elif grid == 'Density grid':
                            model0 = uniformmodel(N, densitybinedges)
                            model1 = gaussmodel(N, densitybinedges, mu, sigma)
                            powerdensity = calculate_power(model0, model1, N, Nsim)
                    print('Ncell=%d  N=%d  power(uniformgrid)=%f  power(densitygrid)=%f' % (2*Nc2,  N, poweruniform, powerdensity))
                    f = open(outname, 'a')
                    f.write('%d\t%f\t%f\n' % (N, poweruniform, powerdensity))
                    f.close()
                print('\n\t OUTPUT: %s\n' % (outname))
    
    
    
    # Plot Figure
    figname = 'fig-power.png'
    plt.rc('font', size=14)
    plt.rc('font', family='sans-serif')
    plt.rc('legend', fontsize=14)
    plt.rc('axes', labelsize=16)
    fig = plt.figure(1, figsize=(35, 35))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    
    Nc2 = int(Ncell_fix/2) # 10
    Neqs = Neqs_fix #Earthquakes
    uniformbinedges = np.linspace(-Z, Z, 2 * Nc2 + 1) #*sigma
    densitybinedges = gaussian_binedges(Z, Nc2, mu, sigma_1)
    
    #Frame a and b - Demonstrate single and multi-resolution grids 
    for ng, grid in enumerate(gridtype):
        print('ng :', ng)
        print('grid :', grid)
        if grid == 'Uniform grid':
            binedges = uniformbinedges
            titlename = '%s $\mathregular{(N_{cell}=%.0f)}$' % (grid, 2*Nc2)
            ax = plt.subplot2grid((2,2), (0,0), colspan=1, rowspan=1)
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        elif grid == 'Density grid':
            binedges = densitybinedges
            titlename = '%s $\mathregular{(N_{cell}=%.0f)}$' % (grid, 2*Nc2)     #r'$e^{-t}$'
            ax = plt.subplot2grid((2,2), (0,1), colspan=1, rowspan=1)
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            
        xb = 0.5 * (binedges[:-1] + binedges[1:])
        model0 = uniformmodel(Neqs, binedges)
        model1 = gaussmodel(Neqs, binedges, mu, sigma_1)
        nobs = simulation(model1, Neqs)     #Using Gaussian model to generate Observed gridded catalog
        plt.plot(xb, model0, c='k', label='Uniform model')
        plt.plot(xb, model1, '--', c='b', label='Gaussian model ($\sigma=%.1f$)'%(sigma_1))
        plt.scatter(xb, nobs, alpha=0.5, label='Test Catalog ($\mathregular{N_{eq}=%d}$)' % (Neqs))
        plt.xlim(-Z, Z) #*sigma
        plt.xlabel('Bin mid-point position')
        plt.ylabel('Events per bin')
#        plt.title(titlename)
        if grid == 'Uniform grid':
            plt.legend(fontsize=14, loc = 'upper right')
            ax.text(min(xb), max(max(model1), max(nobs)) -0.05, '(a)', fontsize=16)
        if grid == 'Density grid':
            ax.text(min(xb)-0.4, max(max(model0), max(nobs))-0.05, '(b)', fontsize=16)
        
        #----- Set Xticks specific to Gaussian distribution Frame b
        tik = [""]
        tik_all = len(xb)*tik
        tik_all[0] = str(np.round(xb[0],2))
        tik_all[-1] = str(np.round(xb[-1],2))
        tik_all[0+int(len(xb)/6)] = str(np.round(xb[0+int(len(xb)/6)],2))
        tik_all[-1-int(len(xb)/6)] = str(np.round(xb[-1-int(len(xb)/6)],2))
        tik_all[0+int(len(xb)/3)] = str(np.round(xb[0+int(len(xb)/3)],2))
        tik_all[-1-int(len(xb)/3)] = str(np.round(xb[-1-int(len(xb)/3)],2))
        plt.xticks(xb, tik_all)
    
    
    #Frame c - Ncells vs Power --- Using Fixed Neqs.
    ax = plt.subplot2grid((2,2), (1,0), colspan=1, rowspan=1)
    Ncell_all = 2*np.arange(1,51,1)
    power_cell = np.array([0,0,0])
    power_cell_2 = np.array([0,0,0])
    for Ncell in Ncell_all:
    #    print('---Number of Cells :', Ncell)
        outname = results_folder+'powertest-Z%.0f-Ncell%.0f-Nsim%d-sigma%.1f.out' % (Z, Ncell, Nsim,sigma_1) #Z=3
        data = np.loadtxt(outname, skiprows=1)
        #The below IF is for the reduced computations we did to reduce computations only to the values we need for figure.
        if data.ndim == 2:
            power_cell = np.row_stack((power_cell,data[data[:,0]==Neqs]))
        else:
            power_cell = np.row_stack((power_cell,data))
    
        #Data 2---
        outname_2 = results_folder+'powertest-Z%.0f-Ncell%.0f-Nsim%d-sigma%.1f.out' % (Z, Ncell, Nsim,sigma_2) #Z=3
        data_2 = np.loadtxt(outname_2, skiprows=1)
        #The below IF is for the reduced computations we did to reduce computations only to the values we need for figure.
        if data_2.ndim == 2:
            power_cell_2 = np.row_stack((power_cell_2,data_2[data_2[:,0]==Neqs]))
        else:
            power_cell_2 = np.row_stack((power_cell_2,data_2))
    
    
    
    power_cell = power_cell[1:]
    power_cell_2 = power_cell_2[1:]
    
    plt.plot(Ncell_all, power_cell[:,1], c='k', label='Uniform Grid')
    plt.plot(Ncell_all, power_cell[:,2], '--', c='b', label='Density Grid ($\sigma=%.1f$)'% (sigma_1))
    plt.plot(Ncell_all, power_cell_2[:,2], ':', c='b', label='Density Grid ($\sigma=%.1f$)'% (sigma_2))
    ax.set_xscale('log', basex=2)
    plt.legend(fontsize=14,loc='center right')
    plt.xlabel('$\mathrm{N_{cell}}$', fontsize = 14) #
    plt.ylabel('Power', fontsize = 14)  #, 
#    plt.title('$\mathregular{N_{cell}}$ vs Power ($\mathregular{N_{eq}=%.0f}$)' % (Neqs)) #, fontsize = 14
    plt.xticks([ 2,  4,  8, 16, 32, 64],[2,  4,  8, 16, 32, 64])
    ax.text(2,0.95, '(c)', fontsize=16)
    ax.set_ylim(0,1.02)
    
    
    #Frame d ---- Neqs vs Power -- Using Fixed Ncells
    ax = plt.subplot2grid((2,2), (1,1), colspan=1, rowspan=1)
    Ncell = 2*Nc2
    outname = results_folder+'powertest-Z%.0f-Ncell%.0f-Nsim%d-sigma%.1f.out' % (Z, Ncell,Nsim, sigma_1) #Z=3
    data = np.loadtxt(outname, skiprows=1)
    
    outname_2 = results_folder+'powertest-Z%.0f-Ncell%.0f-Nsim%d-sigma%.1f.out' % (Z, Ncell, Nsim, sigma_2) #Z=3
    data_2 = np.loadtxt(outname_2, skiprows=1)
    
    plt.plot(data[:,0], data[:,1], c='k', label='Unifrom Grid')
    plt.plot(data[:,0], data[:,2], '--', c='b', label='Density Grid ($\sigma=%.1f$)'% (sigma_1))
    plt.plot(data_2[:,0], data_2[:,2], ':', c='b', label='Density Grid ($\sigma=%.1f$)'% (sigma_2))
    ax.set_xscale('log', basex=2)
    #plt.legend()
    plt.xlabel('$\mathrm{N_{eq}}$', fontsize = 16) #, 
    plt.ylabel('Power', fontsize = 14) #, 
    #plt.title('$\mathregular{N_{eq}}$ vs Power ($\mathregular{N_{cell}=%.0f}$)' % (Ncell)) #, fontsize = 14
    plt.xticks([ 1,  2,  4,  8, 16, 32, 64],[ 1,  2,  4,  8, 16, 32, 64])
    ax.text(1,0.95, '(d)', fontsize=16)
    ax.set_ylim(0,1.02)
    #plt.tight_layout()
    plt.savefig('../Data/Figures/Figure2_stest_power_showcase_sigma_'+str(sigma_1)+'_'+str(sigma_2)+'.png', dpi = 400)
    plt.savefig('../Data/Figures/Figure2_stest_power_showcase_sigma_'+str(sigma_1)+'_'+str(sigma_2)+'.svg')
    
if __name__ == "__main__":
    main()