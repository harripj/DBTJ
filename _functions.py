#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 12:25:20 2018

@author: pjh523
"""
import numpy as np
import multiprocessing as mp

# Boltzmann constant
kB = 1.38e-23
# electronic charge
e = 1.602e-19

def remove_line_from_axes(ax, label):
    '''removes a labelled line from an axes.
    Axes will need redrawing'''
    for line in ax.get_lines():
        if line.get_label() == label:
            ax.lines.remove(line)

def sort_ylim(ax, padding=False, pad=0.025):
    '''Calculates and sets y axes limits given the data on the plot.
    
    Parameters
    ----------
    ax: MPL axes
        Axes instance to adjust.
    padding: bool
        If True add whitespace to edge of axes around data. Default is False.
    pad: float
        Fraction of axes to add as whitespace.

    Returns
    -------
    ymin: float
        New axes lower y limit.
    ymax: float
        New axes upper y limit.
    '''
    for index, line in enumerate(ax.get_lines()):
        ymin_temp = min(line.get_ydata())
        ymax_temp = max(line.get_ydata())
        
        if index == 0:
            ymin = ymin_temp
            ymax = ymax_temp
        else:
            if ymin_temp < ymin:
                ymin = ymin_temp
            if ymax_temp > ymax:
                ymax = ymax_temp
    
    if padding:
        # add padding to axes limits
        yrange = ymax - ymin
        ymin = ymin - pad*yrange
        ymax = ymax + pad*yrange
    
    ax.set_ylim(ymin, ymax)
    
    return ymin, ymax

#pm is plus minus ±, takes values of +1 for plus and -1 for mius
def gamma(R, dE, T):
    '''The Fermi's golden rule tunneling constant.
    
    Parameters
    ----------
    R: float
        The resistance in Ohms across the tunneling junction.
    dE: float
        The change in energy expected from a change in occupancy across the 
        tunneling junction.
    T: float
        Temperature in  K.
    
    Returns
    -------
    gamma: float
        The Fermi's golden rule constant across the tunneling junction.
    '''
    
    
    
    # if dE is positive then gamma will be essentially zero: avoid /0 errors
    if dE > 0:
        # set tunneling rate small- typical gamma when dE<0 is 10**9 larger
        return 1e-1
    else:
        return (-dE/(1-np.exp(dE/(kB*T)))) / (R*e**2)
    
def dE(pm1, pm2, n, V, C, Csum, Q0):
    '''The change in energy expected from a change in occupancy across the
    tunneling junction.
    
    Parameters
    ----------
    pm1, pm2: ints
        +1 for plus, -1 for minus.
    n: int
        Number of electrons on central electrode.
    V: float
        Bias across the central electrode, eg. source-sink, tip-substrate.
    C: float
        Capacitance in Farads of the tunneling junction in question.
    Csum: float
        Total capacitance in Farads across both tunneling junctions.
        Csum = C1 + C2.
    Q0: float
        Value of fractional charge on central electrode in units of electron
        charge. Should be in the range -0.5: 0.5.
        
    Returns
    -------
    dE: float
        The change in energy of moving a charge across the tunneling junction.
    '''
    return (e/Csum) * (e/2 + np.sign(pm1)*(n*e - Q0) + np.sign(pm2)*C*V)

def calculate_tunneling_on_off_ratio(n, bias, R1, R2, C1, C2, Csum, Q0, T):
    '''Calculates the ratio of tunneling on and off of the center electrode.'''
    # energy change across J1 from             
    dE1_minus = dE(-1, -1, n, bias, C2, Csum, Q0)
    dE1_plus = dE(+1, +1, n-1, bias, C2, Csum, Q0)
            
    dE2_minus = dE(-1, +1, n, bias, C1, Csum, Q0)
    dE2_plus = dE(+1, -1, n-1, bias, C1, Csum, Q0)
            
    gamma1_minus = gamma(R1, dE1_minus, T)            
    gamma1_plus = gamma(R1, dE1_plus, T)
    
    gamma2_minus = gamma(R2, dE2_minus, T)
    gamma2_plus = gamma(R2, dE2_plus, T)
            
    # in steady state, sigma(n)*tunnel_off = sigma(n-1)*tunnel_on
    tunnel_on = gamma1_plus + gamma2_plus
    tunnel_off = gamma1_minus + gamma2_minus
    
    return tunnel_on/ tunnel_off

def dbtj_mp(biases, R1, R2, C1, C2, Q0, T=4.2, n_e=(-10, 10)):
    '''
    Double Barrier Tunneling Junction Model
    Hanna, Tinkham (1991), Phys. Rev. B 44 (11), p. 5919-5922.
    
    Parameters
    ----------
    biases: array-like
        Bias range over which to model tunneling current behaviour.
    R1, R2, C1, C2: floats
        Values for resistances and capacitances of tunneling junctions 1 and 2.
        Resistances in Ohm, Capacitances in Farads.   
    Q0: float
        Value of fractional charge on central electrode in units of electron
        charge. Should be in the range -0.5: 0.5.
    T: float
        Temperature of the system in K. Default is 4.2.
    n_e: tuple of ints
        The range of occupancies of electrons on the central electrode over
        which to sum the current contributions, eg. (-5, 5) sums the current
        contributions through the system in which the central electrode has
        -5 to +5 (+Q0) fractional charge state.
        
    Returns
    -------
    I: numpy array
        The modelled DBTJ current at each bias in biases.
    '''
    #sum of capacitances in the system
    Csum = C1 + C2
    n_min, n_max = n_e
    
    #Q0 in same units as electron charge
    Q0 *= e
    if Q0>e/2:
        Q0 = e/2
        print('Q0 has been set to 0.5e.')
    elif Q0<-e/2:
        Q0 = -e/2
        print('Q0 has been set to -0.5e.')
        
    #current array
    I = np.zeros_like(biases)
    # I2 = I.copy()
    
    #loop over bias range of interest
    for index, bias in enumerate(biases):
        print(index)
        
        #loop no electrons on center electrode from -10 to 10 in integer steps
        sigma = np.zeros((n_max-n_min, 2))
        
        pool = mp.Pool(processes=mp.cpu_count()-1)
        tunnel_ratio = pool.starmap(calculate_tunneling_on_off_ratio,
                                  ((n, bias, R1, R2, C1, C2, Csum, Q0, T) \
                                   for n in range(n_min, n_max, 1)))
        
        #loop over n to calculate sigma n for this given voltage
        for index_n, n in enumerate(range(n_min, n_max, 1)):
            if index_n == 0:
                # sigma starts at 1 -> will be normalised
                sigma[index_n] = np.array((n, 1))
                continue
            
            sigma[index_n] = np.array((n, sigma[index_n-1, 1] \
                                          * tunnel_ratio[index_n]))
        
        # normalise sigma after looping over all n values
        sigma[:,1] /= np.sum(sigma[:,1])
        
        for n, sigma_val in sigma:
#            Useful to compare tunnel_on vs tunnel_off currents-
#            they end up shifted by one inde from each other
            
#            dE1_plus = dE(+1, +1, n, bias, C2, Csum ,Q0)
#            dE1_minus = dE(-1, -1, n, bias, C2, Csum, Q0)
#            tunnel_on_J1 = gamma(R1, dE1_plus, T)
#            tunnel_off_J1 = gamma(R1, dE1_minus, T)
#            
#            I2[index_sigma] += e * sigma_val * (tunnel_on_J1 - tunnel_off_J1)
            
            dE2_plus = dE(+1, -1, n, bias, C1, Csum, Q0)
            dE2_minus = dE(-1, +1, n, bias, C1, Csum, Q0)
            tunnel_on_J2 = gamma(R2, dE2_plus, T)
            tunnel_off_J2 = gamma(R2, dE2_minus, T)
            
            I[index] += e * sigma_val * (tunnel_on_J2 - tunnel_off_J2)
    #print((I==I2).all(), (I==-I2).all())
    return I# , I2

def dbtj(biases, R1, R2, C1, C2, Q0, T=4.2, n_e=(-10, 10)):
    '''
    Double Barrier Tunneling Junction Model
    Hanna, Tinkham (1991), Phys. Rev. B 44 (11), p. 5919-5922.
    
    Parameters
    ----------
    biases: array-like
        Bias range over which to model tunneling current behaviour.
    R1, R2, C1, C2: floats
        Values for resistances and capacitances of tunneling junctions 1 and 2.
        Resistances in Ohm, Capacitances in Farads.   
    Q0: float
        Value of fractional charge on central electrode in units of electron
        charge. Should be in the range -0.5: 0.5.
    T: float
        Temperature of the system in K. Default is 4.2.
    n_e: tuple of ints
        The range of occupancies of electrons on the central electrode over
        which to sum the current contributions, eg. (-5, 5) sums the current
        contributions through the system in which the central electrode has
        -5 to +5 (+Q0) fractional charge state.
        
    Returns
    -------
    I: numpy array
        The modelled DBTJ current at each bias in biases.
    '''
    #sum of capacitances in the system
    Csum = C1 + C2
    n_min, n_max = n_e
    
    #Q0 in same units as electron charge
    Q0 *= e
    if Q0>e/2:
        Q0 = e/2
        print('Q0 has been set to 0.5e.')
    elif Q0<-e/2:
        Q0 = -e/2
        print('Q0 has been set to -0.5e.')
        
    #current array
    I = np.zeros_like(biases)
    # I2 = I.copy()
    
    #loop over bias range of interest
    for index, bias in enumerate(biases):
        #loop no electrons on center electrode from -10 to 10 in integer steps
        sigma = np.zeros((n_max-n_min, 2))
        
        #loop over n to calculate sigma n for this given voltage
        for index_n, n in enumerate(range(n_min, n_max, 1)):
            if index_n == 0:
                # sigma starts at 1 -> will be normalised
                sigma[index_n] = np.array((n, 1))
                continue
            
            tunneling_ratio = calculate_tunneling_on_off_ratio(n, bias, R1, \
                                                       R2, C1, C2, Csum, Q0, T)
            
            sigma[index_n] = np.array((n, sigma[index_n-1, 1]*tunneling_ratio))
        
        # normalise sigma after looping over all n values
        sigma[:,1] /= np.sum(sigma[:,1])
        
        for n, sigma_val in sigma:
#            Useful to compare tunnel_on vs tunnel_off currents-
#            they end up shifted by one inde from each other
            
#            dE1_plus = dE(+1, +1, n, bias, C2, Csum ,Q0)
#            dE1_minus = dE(-1, -1, n, bias, C2, Csum, Q0)
#            tunnel_on_J1 = gamma(R1, dE1_plus, T)
#            tunnel_off_J1 = gamma(R1, dE1_minus, T)
#            
#            I2[index_sigma] += e * sigma_val * (tunnel_on_J1 - tunnel_off_J1)
            
            dE2_plus = dE(+1, -1, n, bias, C1, Csum, Q0)
            dE2_minus = dE(-1, +1, n, bias, C1, Csum, Q0)
            tunnel_on_J2 = gamma(R2, dE2_plus, T)
            tunnel_off_J2 = gamma(R2, dE2_minus, T)
            
            I[index] += e * sigma_val * (tunnel_on_J2 - tunnel_off_J2)
    #print((I==I2).all(), (I==-I2).all())
    return I# , I2