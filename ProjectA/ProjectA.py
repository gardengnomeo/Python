# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 12:24:07 2024

@author: lucas Simpson
"""

import math
import matplotlib.pyplot as plt


def main():

    # Use these initial conditions for the application
    
    #Initializing variables.
    initial_B1_fast, initial_B1_slow = 200, 200
    initial_B2_fast, initial_B2_slow = 1, 1
    initial_P_fast, initial_P_slow = 0.01, 0.01
    sB2fast = 0.2  # B2 sensitivity parameter for fast compartment
    sB2slow = 0.1  # B2 sensitivity parameter for slow compartment
    sB1fast = 0.01  # B1 effectiveness parameter for fast compartment
    sB1slow = 0.005  # B1 effectiveness parameter for slow compartment
    B1p = 80  # Plasma B1 concentration
    B1t = 20  # Target B1 concentration
    delay_fast = 5  # Time delay for B2 action in fast compartment
    delay_slow = 20 # Time delay for B2 action in slow compartment
    beta = 1 # use 3 for normal and 1 for disease
    B1thres = 150
    
    dt = 0.1
    total_time = 60
    
    simulate(initial_B1_fast, initial_B1_slow,initial_B2_fast, 
                  initial_B2_slow,initial_P_fast, initial_P_slow,
                  sB2fast, sB2slow, sB1fast, sB1slow, B1p, B1t, 
                  delay_fast,delay_slow,beta, B1thres,dt,total_time)
    
#defining all of the functions that will be used.
        
def molecule1_dynam_fast(sB1fast, B1fast, B1p, B2fast, Pfast):
    """
    Update rate of change of bio-molecule 1 value in fast compartment. 

    Parameters
    ----------
    sB1fast : float
        Sensitivity parameter for bio-molecule 1 in the fast compartment.
    B1fast : float
        Bio-molecule 1 concentration in the fast compartment.
    B1p : int
        Plasma bio-molecule 1 concentration.
    B2fast : float
        Bio-molecule 1 concentration in the fast compartment.
    Pfast : int
        Bio-molecule 2 secretion rate in the fast compartment.

    Returns
    -------
    float
        Rate of change of B1 fast value.
    """
    #How B1 levels change in fast compartment using eqauation 1.
    return (-1 * sB1fast * (B1fast - B1p) - (B2fast * B1fast) + Pfast)

def molecule1_dynam_slow(sB1slow, B1slow, B1p, B2slow, Pslow):
    """
    Update rate of change of bio-molecule 1 value in slow compartment. 


    Parameters
    ----------
    sB1slow : float
        Sensitivity parameter for bio-molecule 1 in the slow compartment.
    B1slow : float
        Bio-molecule 1 concentration in the slow compartment.
    B1p : int
        Plasma bio-molecule 1 concentration.
    B2slow : float
        Bio-molecule 1 concentration in the slow compartment.
    Pslow : int  
        Bio-molecule 2 secretion rate in the slow compartment.

    Returns
    -------
    float
        Rate of change of B1 slow value.
    """
    #How B1 levels change in slow compartment using eqauation 2.
    return -1 * sB1slow * (B1slow - B1p) - (B2slow * B1slow) + Pslow
    
def molecule2_dynam_fast(delay_fast,sB2fast,B1fast, B2fast, B1t):
    """
    Update rate of change of bio-molecule 2 value in fast compartment. 

    Parameters
    ----------
    delay_fast : int
        Time delay for the fast compartment.
    sB2fast : float
        Sensitivity parameter for bio-molecule 2 in the fast compartment.
    B1fast : float
        Concentration of the bio-molecule 1 in the fast compartment.
     B2fast : float
        Concentration of bio-molecule 2 interacting with the bio-
molecule 1 in the fast compartment.
    B1t : int
        Target bio-molecule 1 concentration.

    Returns
    -------
    float
        Rate of change of B2 fast value.
    """
    #How B2 levels change in fast compartment in responce to B1 using eqauation 3.
    return (1/delay_fast)* (sB2fast *(B1fast - B1t) - (B2fast))

def molecule2_dynam_slow(delay_slow, sB2slow, B1slow, B2slow, B1t):
    """
    Update rate of change of bio-molecule 2 value in slow compartment.

    Parameters
    ----------
    delay_slow : int
        Time delay parameter for the slow compartment.
    sB2slow : flaot
        Sensitivity parameter for bio-molecule 2 in the slow compartment.
    B1slow : float
        Concentration of the bio-molecule 1 in the slow compartment.
    B2slow : float
        Concentration of bio-molecule 2 interacting with the biomolecule 1 
        in the slow compartment.
    B1t : int
        Target concentration for the bio-molecule 1.

    Returns
    -------
    float
        Rate of change of B2 slow value.
    """
    #How B2 levels change in slow compartment in responce to B1 using eqauation 4.
    return (1/delay_slow)* (sB2slow *(B1slow - B1t) - (B2slow))

def molecule2_sec_fast(beta, B1fast, B1thres, Pfast, t):
    """
    Update the rate of change in the secretion rate of bio-molecule 2
    affecting the concentration of the bio-molecule 1 in the fast compartment.

    Parameters
    ----------
    beta : int
        Sensitivity parameter for bio-molecule 2 secretion.
    B1fast : float
        Bio-molecule 1 concentration in the fast compartment.
    B1thres : int
        Threshold of bio-molecule-1 concentration for bio-molecule 2 release.
    Pfast : flaot
        Bio-molecule 2 secretion rate in the fast compartment.
    t : int
        Simulation time variable.

    Returns
    -------
    float
        Rate of change of Pfast.
    """
    #B2 secretion rate in fast compartment using eqauation 5.
    return beta * max(0, B1fast - B1thres) * Pfast * math.sin(0.1 * t)

def molecule2_sec_slow(beta, B1slow, B1thres, Pslow, t):
    """
    Update the rate of change in the secretion rate of bio-molecule 2
    affecting the concentration of the bio-molecule 1 in the slow compartment.

    Parameters
    ----------
    beta : int
        Sensitivity parameter for bio-molecule 2 secretion.
    B1slow : float
        DESCRIPTION.
    B1thres : int
        Threshold of bio-molecule-1 concentration for bio-molecule 2 release.
    Pslow : float
        Bio-molecule 2 secretion rate in the slow compartment.
    t : int
        Simulation time variable.

    Returns
    -------
    float
        Rate of change of Pslow.
    """
    #B2 secretion rate in slow compartment using eqauation 6.
    return beta * max(0, B1slow - B1thres) * Pslow * math.cos(0.1 * t)

def update_mol1_fast(fun, sB1fast, B1fast, B1p, B2fast, Pfast, dt):
    """
    Update concentration of the bio-molecule 1 in the fast compartment.

    Parameters
    ----------
    fun : fun
        Updates rate of change of bio-molecule 1 value in fast compartment.
    sB1fast : float
        Sensitivity parameter for bio-molecule 1 in the fast compartment.
    B1fast : float
        Bio-molecule 1 concentration in the fast compartment.
    B1p : int
        Plasma bio-molecule 1 concentration.
    B2fast : float
        Concentration of bio-molecule 2 interacting with the biomolecule 1 
        in the fast compartment.
    Pfast : float
        Bio-molecule 2 secretion rate in the fast compartment.
    dt : float
        Time step.

    Returns
    -------
    float
        Concentration of the bio-molecule 1 in the fast compartment.
    """
    #Updataes the concentrations of B1 fast (molecule1_dynam_fast)
    return (B1fast + fun(sB1fast, B1fast, B1p, B2fast, 
                                          Pfast) * dt)

def update_mol1_slow(fun, sB1slow, B1slow, B1p, B2slow, Pslow, dt):
    """
    Update concentration of the bio-molecule 1 in the slow compartment.

    Parameters
    ----------
    fun : fun
    Updates rate of change of bio-molecule 1 value in slow compartment.
    sB1slow : float
        Sensitivity parameter for bio-molecule 1 in the slow compartment.
    B1slow : float
        Bio-molecule 1 concentration in the slow compartment.
    B1p : int
        Plasma bio-molecule 1 concentration.
    B2slow : float
        Concentration of bio-molecule 2 interacting with the biomolecule 1 
        in the slow compartment..
    Pslow : float
        Bio-molecule 2 secretion rate in the slow compartment.
    dt : float
        Time step.

    Returns
    -------
    float
        Concentration of the bio-molecule 1 in the slow compartment.
    """
    #Updataes the concentrations of B1 slow (molecule1_dynam_slow)
    return (B1slow + fun(sB1slow, B1slow, B1p, B2slow, 
                                           Pslow) * dt)

def update_mol2_fast(fun, delay_fast, sB2fast, B1fast, B2fast, B1t, dt):
    """
    Update concentration of the bio-molecule 2 in the fast compartment.


    Parameters
    ----------
    fun : fun
        Updates rate of change of bio-molecule 2 value in fast compartment.
    delay_fast : int
        Time delay for the fast compartment.
    sB2fast : float
        Sensitivity parameter for bio-molecule 2 in the fast compartment.
    B1fast : float
        Bio-molecule 1 concentration in the fast compartment.
    B2fast : float
        Concentration of bio-molecule 2 interacting with the biomolecule 1 
        in the fast compartment.
    B1t : int
        Target concentration for the bio-molecule 1.
    dt : float
        Time Step.

    Returns
    -------
    float
        Updated concentration of the bio-molecule 2 in the fast compartment.
    """
    #Updataes the concentrations of B2 fast (molecule2_dynam_fast)
    return (B2fast + fun(delay_fast,sB2fast,B1fast,
                                           B2fast, B1t) * dt)

def update_mol2_slow(fun, delay_slow, sB2slow, B1slow, B2slow, B1t, dt):
    """
    Update concentration of the bio-molecule 2 in the slow compartment.


    Parameters
    ----------
    fun : fun
        Update rate of change of bio-molecule 2 value in slow compartment.
    delay_slow : int
        Time delay for the slow compartment.
    sB2slow : float
        Sensitivity parameter for bio-molecule 2 in the slow compartment.
    B1slow : flaot
        Bio-molecule 1 concentration in the slow compartment.
    B2slow : float
        Bio-molecule 2 concentration in the slow compartment.
    B1t : int
        Target concentration for the bio-molecule 1.
    dt : float
        Tiime step.

    Returns
    -------
    float
        Updated concentration of the bio-molecule 2 in the slow compartment.
    """
    #Updtaes the concentrations of B2 slow (molecule2_dynam_slow)
    return (B2slow + fun(delay_slow, sB2slow, B1slow,
                                           B2slow, B1t) * dt)

def update_mol2_sec_fast(fun, beta, B1fast, B1thres, Pfast, t, dt):
    """
    Update rate of change in the secretion rate of bio-molecule 2 in the 
    fast compartment.

    Parameters
    ----------
    fun : fun
        Update the rate of change in the secretion rate of bio-molecule 2
        affecting the concentration of the bio-molecule 1 in the fast 
        compartment.
    beta : int
        Sensitivity parameter for bio-molecule 2 secretion.
    B1fast : float
        Bio-molecule 1 concentration in the fast compartment.
    B1thres : int
        Threshold of bio-molecule-1 concentration for bio-molecule 2 release.
    Pfast : float
        Bio-molecule 2 secretion rate in the fast compartment.
    t : int
        Simulation time variable.
    dt : float
        Time step.

    Returns
    -------
    float
        Updated rate of change of the secretation rate of biomolecule 2 in 
        the fast compartment.
    """
    #Updataes the secretion rate of B2 fast (molecule2_sec_fast)
    return (Pfast + fun(beta, B1fast, B1thres, Pfast, t) * dt)

def update_mol2_sec_slow(fun, beta, B1slow, B1thres, Pslow, t, dt):
    """
    Update rate of change in the secretion rate of bio-molecule 2 in the 
    slow compartment.

    Parameters
    ----------
    fun : fun
        Update rate of change in the secretion rate of bio-molecule 2 in the 
        fast compartment.
    beta : int
        Sensitivity parameter for bio-molecule 2 secretion.
    B1slow : float
        Bio-molecule 1 concentration in the fast compartment.
    B1thres : int
        Threshold of bio-molecule-1 concentration for bio-molecule 2 release.
    Pslow : float
        Bio-molecule 2 secretion rate in the fast compartment.
    t : int
        Simulation time variable.
    dt : float
        Time step.

    Returns
    -------
    float
        Updated rate of change of the secretation rate of biomolecule 2 in 
        the fast compartment.
    """
    #Updataes the secretion rate of B1 fast (molecule2_sec_slow)
    return (Pslow + fun(beta, B1slow, B1thres, Pslow, t) * dt)

'''Simulates the dynamics of the system with B1 and B2 in the
fast and slow compartments. Visualizes the data'''
def simulate(initial_B1_fast, initial_B1_slow,initial_B2_fast, 
             initial_B2_slow,initial_P_fast, initial_P_slow,
             sB2fast, sB2slow, sB1fast, sB1slow, B1p, B1t, 
             delay_fast,delay_slow,beta, B1thres,dt,total_time):
    
    
    fig, ax = plt.subplots(1,2,figsize=(15, 7))
    
    #### ADD YOUR CODE BELOW ####
    
    '''initializing the variables'''
    B1fast = initial_B1_fast
    B1slow = initial_B1_slow
    B2fast = initial_B2_fast
    B2slow = initial_B2_slow
    Pfast = initial_P_fast
    Pslow = initial_P_slow
    t = 0
    
    #Finding the values for all values of t from 0 to total_time. 
    while t < (total_time + dt):
        #plots the point for the t value.
        plot_result(ax, B1fast, B1slow, B2fast, B2slow, t)
        
        #Updating the variables with their new values for each value of t.
        B1fast_new = (update_mol1_fast(molecule1_dynam_fast, sB1fast, B1fast, 
                                       B1p, B2fast, Pfast, dt))
        B1slow_new = update_mol1_slow(molecule1_dynam_slow, sB1slow, B1slow,
                                      B1p, B2slow, Pslow, dt)
        B2fast_new = update_mol2_fast(molecule2_dynam_fast, delay_fast,
                                      sB2fast, B1fast, B2fast, B1t, dt)
        B2slow_new = update_mol2_slow(molecule2_dynam_slow, delay_slow,
                                      sB2slow, B1slow, B2slow, B1t, dt)
        Pfast_new = update_mol2_sec_fast(molecule2_sec_fast, beta, B1fast,
                                         B1thres, Pfast, t, dt)
        Pslow_new = update_mol2_sec_slow(molecule2_sec_slow, beta, B1slow,
                                         B1thres, Pslow, t, dt)
        #Setting each variable to the new value. 
        B1fast = B1fast_new
        B1slow = B1slow_new
        B2fast = B2fast_new
        B2slow = B2slow_new
        Pfast = Pfast_new
        Pslow = Pslow_new
        
        #Increasing the value of t by dt.
        t += dt
    
    # Do not modify this return statement
    # This should ramain the last statement in
    # the simulate function
    return ax

### DO NOT CHANGE THIS FUNCTION ####  
### Use this function  in your simulate function 
### to plot the concentrations of B1 and B2 molecules 

def plot_result(ax,B1fast,B1slow,B2fast,B2slow,t):
    
    """
    Plot the concentrations of B1 and B2 over time.

    Parameters
    ----------
    ax : list of matplotlib.axes.Axes
        List containing two Axes objects for plotting B1 and B2 concentrations.
    B1fast : float
        Concentration of B1 in the fast compartment at time t.
    B1slow : float
        Concentration of B1 in the slow compartment at time t.
    B2fast : float
        Concentration of B2 in the fast compartment at time t.
    B2slow : float
        Concentration of B2 in the slow compartment at time t.
    t : float
        Time value for which concentrations are plotted.

    Returns
    -------
    None

    Notes
    -----
    This function plots the concentrations of B1 and B2 in the fast and slow compartments at a particular time point.
    """
    
    
    ax[0].plot(t, B1fast, 'bo')
    ax[0].plot(t, B1slow, 'go')
    ax[1].plot(t, B2fast, 'ro')
    ax[1].plot(t, B2slow, 'mo')
    ax[0].set_xlabel('Time')
    ax[0].set_ylabel('B1 Concentration')
    
    ax[0].set_title(f'Complex B1 Response Model at Time: {t:.1f}')
    
    ax[1].set_xlabel('Time')
    ax[1].set_ylabel('B2 Concentration')
    
    ax[1].set_title(f'Complex B2 Response Model at Time: {t:.1f}')
    
    # Uncomment to see an animation of B1 and B2 
    # concentrations at each time point as opposed
    # to see all at the end of the simulation
    #plt.pause(0.0001)  # Pause
    
    ax[0].legend(['B1 (Fast)','B1 (Slow)'],loc=1)
    ax[1].legend(['B2 (Fast)','B2 (Slow)'],loc=1)
   



if __name__ == '__main__':
    
    main()
  
