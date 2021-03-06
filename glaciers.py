#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:44:02 2018

Python2 code for Glacier Mass Balance

@author: nadine
"""
#%% import modules
import numpy as np
import matplotlib.pyplot as plt
import os

#%% ANALYTIC SOLUTION FOR ICE DISCHARGE (Q)

gamma = 0.01 # m/yr per m elevation
ELA = 1500. # meters

xmax = 20000. # meters
dx = 150. # meters
x = np.arange(0,xmax,dx)

zbmax = 2000. # meters
zbmin = 500. # meters
zb = np.linspace(zbmax,zbmin,num=len(x))

slope = np.arctan(zbmax/xmax)
b = gamma*(zb-ELA)/(365.25*24*60*60)


Q = gamma*x*(zbmax-ELA) -(0.5*gamma*slope*(x**2))
Q = np.maximum(0,Q)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1,sharex=True)
plt.tight_layout()


ax1.plot(x,zb,'gray')
ax1.fill_between(x,0,zb,facecolor='dimgray')
ax1.set_ylabel('elevation [m]')
ax1.set_title('Analytical Solution for Ice Discharge, Q')


ax2.plot(x,b,'blue')
ax2.set_ylabel('b [m/yr]')
    
ax3.plot(x,Q,'red')
ax3.grid(color='lightgray',linestyle='--')
ax3.set_ylabel('Q [m^2/sec]')
#ax3.set_ylim(0,300)
ax3.set_xlabel('distance [m]')


#%% FINITE DIFFERENCE SOLUTION

# define variables
xmax = 20000. # meters
dx = 150. # meters
x = np.arange(0,xmax,dx)
xstar = 2000

zbmax = 1500. # meters
zbmax2 = 500
slope = .05 # meters/meters
zbmin = 0. # meters
#zb = np.linspace(zbmax,zbmin,num=len(x)) # bedrock topography [meters]
zb = zbmax - (slope*x)
zb = zb+(zbmax2*(np.exp(-x/xstar)))

#dz = zb[0]-zb[1]

rho_ice = 917. # density of ice in kg/m^3
rho_water = 1000. # density of water in kg/m^3
g = 9.81 # gravity in meters/sec^2

Ao = 2.1*(10**-16) * (1/(365.25*24*60*60))  # 1/yr Pa^3 --> convert to seconds
A = Ao
#Ea = 61. * (10**3)   # activate energy of ice in Joule/mole
#R = 8.31            # universal gas constant in Joule/mole-K
#Tk =  273. # make this vary through time!  --> array length of time          # absolute temperature in Kelvin
#A = Ao*np.exp(-Ea/(R*Tk))


ELA = 1500. # meters
gamma = .01 # slope of mass balance gradient [meters/year per meter of elevation]
b = gamma * (zb - ELA)/(365.25*24*60*60) # glacier mass balance rule
#dbdz = np.diff(b) / np.diff(zb) # lapse rate? change in b with elevation


years = 1000. # number of years to run for
tmax = years * 365.25*24*60*60 # years in seconds
dt = 365.25*24*60*60 * .002 # seconds - 86400 seconds in one day
time = np.arange(0,tmax+dt,dt) # time array in seconds


# initialize arrays for loop
zi = np.zeros(shape=(len(x),len(time)),dtype=float)
H = np.zeros(shape=(len(x),len(time)),dtype=float)
Hedge = np.zeros(len(x)-1)
Q = np.zeros(shape=(len(x)-1,len(time)),dtype=float)
V = np.zeros(len(time),dtype=float)

dQdx = np.zeros(len(x)-2)
dHdt = np.zeros(len(x)-2)
dzidx = np.zeros(len(x)-1)


#%% run loop

for i in range(len(time)-1):
    V[i] = np.sum(H[:,i] * (dx**2)) # calculate ice volume at this time
    dzidx = np.diff(zi[:,i]) / dx # calculate dzidx for all x at this time
    Hedge = H[0:-1,i] # assigning upstream value of H to Hedge at this time
    #Hedge = np.maximum(0,Hedge[:])
    Q[:,i] = A * ((rho_ice*g*np.abs(dzidx))**3)*(-np.sign(dzidx)) * ((Hedge**5)/5) # calculate Q for all x at this time
    dQdx = np.diff(Q[:,i]) / dx # calculate dQdx for all x at this time
    b = gamma * (zi[:,i] - ELA)/(365.25*24*60*60) # calculate b based on ice elevation for this time
    dHdt = b[1:-1] - dQdx # calculate dHdt based on b and Q for all x at this time
    H[1:-1,i+1] = H[1:-1,i] + (dHdt * dt) # update H, ice thickness, with H from last time plus dH this time for all x
    H[1:-1,i+1] = np.maximum(0,H[1:-1,i+1]) # take maximum of 0 and H to ensure no negative H
    zi[:,i+1] = zb + H[:,i+1] # update ice elev for next time = bedrock elev + ice thickness
    
#%% plot ice volume through time 
    
# calculate characteristic time scale
chartimescale = 0.63212 * V[-10]

# plot ice voume with characteristic time scale
plt.plot((time[:-1]/(365.25*24*60*60)),V[:-1])
#plt.plot(chartimescale,'.r')
plt.axhline(chartimescale,0,tmax,linestyle='--',color='k')
plt.plot(187,chartimescale,'.r')
plt.xlabel('time [years]')
plt.ylabel('ice volume [m^3]')
plt.title('Ice Volume')
plt.grid(linestyle='--',color='lightgrey')

#equiltime = 187*3.5

#%% plot many frames
    
plots = 2500
    
for i in range(len(time)-1):
    if i % plots == 0:
        
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        plt.tight_layout()

        ax1.plot(x,zi[:,i],'blue')
        ax1.fill_between(x,zb,zi[:,i],facecolor='powderblue')
        ax1.fill_between(x,0,zb,facecolor='dimgray')
        ax1.set_ylabel('elevation [m]')
        ax1.set_xlabel('distance [m]')
        ax1.set_xlim(0,xmax)
        ax1.set_ylim(0,2200)
        ax1.grid(color='lightgray',linestyle='--')
        ax1.text(10000,1500,str(np.around((time[i]/(365.25*24*60*60)),decimals=1))+' years')

        ax2.plot(x[1:],Q[:,i],'red')
        ax2.grid(color='lightgray',linestyle='--')
        ax2.set_xlim(0,xmax)
        ax2.set_ylim(0,.0005)
        ax2.set_ylabel('Q [m^2/sec]')
        ax2.set_xlabel('distance [m]')

        plt.savefig('tmp'+str(i/plots)+'.png',bbox_inches="tight",dpi=200)
        plt.close()

    
#%% make a movie with ffmpeg
## -r = fps
os.system("rm movie.mp4") # remove a previous movie.mp4 file so don't get overwrite problems
os.system("ffmpeg -r 15 -pattern_type sequence -i tmp'%d'.png -vcodec mpeg4 movie.mp4") 
os.system("rm tmp*.png")



#%% plot a couple frames
        
#i = 0,1,2,3,4,5
#
#for i in range(len(i)):
#    fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True)
#    plt.tight_layout()
#
#    ax1.plot(x,zi[:,i],'blue')
#    ax1.fill_between(x,zb,zi[:,i],facecolor='powderblue')
#    ax1.fill_between(x,0,zb,facecolor='dimgray')
#    ax1.set_ylabel('elevation [m]')
#    #ax1.text(600,400,str(np.around((i/dt),decimals=1))+' years')
#    
#    ax2.plot(x[1:],Q[:,i],'red')
#    ax2.grid(color='lightgray',linestyle='--')
#    ax2.set_ylabel('Q [m^2/sec]')
#    ax2.set_xlabel('distance [m]')

