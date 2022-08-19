# Copyright (C) 2015-2020: The University of Edinburgh
#                 Authors: Craig Warren and Antonis Giannopoulos
#
# This file is part of gprMax.
#
# gprMax is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gprMax is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gprMax.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import os
import sys

import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from gprMax.exceptions import CmdInputError
from .outputfiles_merge import get_output_data
from .outputfiles_merge import get_output_all_data

from tqdm import tqdm 
import time 


# Normalized scan data
def norm(_inp):
    _max = np.max(_inp)
    _min = np.min(_inp)
    _zeroLev = 0
    if _max - _min > 0 :
        _zeroLev = (_inp[0] - _min)/(_max - _min)
    
    norm = []
    for i in range(len(_inp)):
        tmp_Val = 0
        if _max - _min > 0 :
            tmp_Val = (_inp[i] - _min)/(_max - _min) - _zeroLev
        norm.append(tmp_Val)
    return norm
def cm_to_inch(value):
    return value/2.54
def mpl_plot(filename, outputdata, dt, dx, er, gmin, gmax, start, rxnumber, rxcomponent, width, height, isnorm = False, rawdata = False):  
    #Converting data
    #tcut = 4.15e-9
    #c = 299792458.0 # light speed (m/s)
    c = 3e+8 # light speed (m/s)
    d = [] # depth
    p = [] # position as follows x-direction [m]
    t = []
    dx = float(0.01)
    start = 0.5	
    gain = np.linspace(float(gmin),float(gmax),len(outputdata)) # gain range value [gmin, gmax]
    #dx = float(dx)
    #start = float(start)
    #print("Viet %f",start)
    for k in range(len(outputdata[0])): # calculate position (x-direction)
        tmp = start + k*(float)(dx)
        p.append(tmp)
    # Calculate Depth range [m]
    for i in range(len(outputdata)): #number of row
        if er == -1:
            t.append(i*dt*1e+9)
        else:
            depthval = (c*i*dt)/(2*math.sqrt((float)(er))); # Depth in metre [m]
            d.append(depthval)
    #Time data
    time = np.array(t, dtype = np.float64)
    time = time.reshape(len(time), 1) 
    #Simulated text file - Raw data
    if rawdata:
        a_file = open(filename + "_raw.txt", "w")
        savedata = np.array(outputdata,dtype = np.float64)
        bscan_d = np.hstack((time, outputdata))
        #print('' + str(len(bscan_d)) + ',' + str(bscan_d.shape[1]))
        np.savetxt(a_file, bscan_d, delimiter=",")
        a_file.close()
    
    #Transpose matrix outputdata
    caldata = np.transpose(outputdata)
    if gmin is not gmax:       
        for i in tqdm(range(len(caldata)),  
                desc="Applying gain function…",  
                ascii=False, ncols=75): #number of column
            #time.sleep(0.01) 
            rowdata = caldata[i]
            for j in range(len(rowdata)): #loop for row
                scaleval = gain[j] # assuming scale factor equal row index
                dataval = rowdata[j]         
                rowdata[j] = dataval* scaleval # already converted as follows scale factor
                
            #print("scale={}, before={}, after={}".format(scaleval,dataval,rowdata[j] ))         
            caldata[i] = rowdata
        # Back to transformation
    if isnorm:
        for i in range(len(caldata)):
            rowdata = caldata[i]
            caldata[i] = norm(rowdata)
    outputdata = np.transpose(caldata)
    
    #Simulated text file - Normalized data
    if rawdata:
        a_file = open(filename + "_norm.txt", "w")
        savedata = np.array(outputdata,dtype = np.float64)
        bscan_d = np.hstack((time, outputdata))
        #print('' + str(len(bscan_d)) + ',' + str(bscan_d.shape[1]))
        np.savetxt(a_file, bscan_d, delimiter=",")
        a_file.close()

    z = outputdata
    
    #Save matrix to text file
    #a_file = open("test.txt", "w")
    #for row in outputdata:
    #    np.savetxt(a_file, row)
    #a_file.close()
    
    # Plot contour
    
    #plt.plots(num='rx' + str(rxnumber), figsize=(20, 10), facecolor='w', edgecolor='w')

    # Calibrate axis
    if er == -1:
        tmin = np.min(t)
        tmax = np.max(t)
    else:
        dmax = np.max(d)
        dmin = np.min(d)
    pmax = np.max(p)
    pmin = np.min(p)

    # plot filled contour map with 100 levels

    cmapname = 'gist_gray' #'binary'

    # MODELLING B-SCAN
    #(path, filename) = os.path.split(filename)
    #plt.plot(num='Bscan',figsize=(float(width), float(height)), facecolor='w', edgecolor='w')
    plt.figure(figsize=(cm_to_inch(float(width)), cm_to_inch(float(height))))
    levels1 = MaxNLocator(nbins=200).tick_values(-np.max(abs(z)), np.max(abs(z)))
    if er == -1:
        cs1 = plt.contourf(p, t, z, 200, cmap=cmapname,vmin=-np.max(abs(z)),vmax=np.max(abs(z)),levels=levels1)
    else:
        cs1 = plt.contourf(p, d, z, 200, cmap=cmapname,vmin=-np.max(abs(z)),vmax=np.max(abs(z)),levels=levels1)
    
    plt.title('Simulated B-scan (3D)',style='normal', fontsize=18)
    
    plt.xlim(pmin,pmax)
    if er == -1:
        plt.ylim(tmax,tmin)
    else:
        plt.ylim(dmax,dmin)
    plt.tick_params(labelsize=14)
    if er == -1:
        plt.ylabel('Time [s]',style='normal', fontsize=14)
    else:
        plt.ylabel('Depth [m]',style='normal', fontsize=14)
    
    plt.xlabel('Position [m]',style='normal', fontsize=14)

    # add default colorbar for the map
    cb = plt.colorbar()
    plt.savefig(filename + "_Bscan.png", dpi=300, bbox_inches='tight')
    #cb.setlabel('Normalized Amplitude [-]',style='normal', fontsize=14, size='large', weight='bold')
    #cb.tick_params(labelsize='large')
    return plt


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Plots a B-scan image.', usage='cd gprMax; python -m tools.plot_Bscan_gain outputfile output')
    parser.add_argument('outputfile', help='name of output file including path')
    parser.add_argument('rx_component', help='name of output component to be plotted', choices=['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', 'Ix', 'Iy', 'Iz'])
    parser.add_argument('-er',  default=-1, help='Dielectric constant of the medium')
    parser.add_argument('-norm', action='store_true', help='Is Normalized', default=False)
    parser.add_argument('-gmin',  default=1, help='Minimum Gain Factor')
    parser.add_argument('-gmax', default=1, help='Maximum Gain Factor')
    parser.add_argument('-width', default=36, help='Width of figure')
    parser.add_argument('-height', default=12, help='Heigh of figure')
    parser.add_argument('-rawdata', action='store_true', help='Write rawdata', default=False)
    args = parser.parse_args()

    f = h5py.File(args.outputfile, 'r')
    nrx = f.attrs['nrx']
    f.close()

    # Check there are any receivers
    if nrx == 0:
        raise CmdInputError('No receivers found in {}'.format(outputfile))
    
    for rx in range(1, nrx + 1):
        outputdata, dt, dx, start = get_output_all_data(args.outputfile, rx, args.rx_component)
        #print("dt=%f, dx=%f, start=%f",dt,dx,start)
        print("dt={},dx={}, len[output]={}, {}".format(dt,dx,len(outputdata),len(outputdata[0])))
        #Convert to depth
        plthandle = mpl_plot(args.outputfile, outputdata, dt, dx, args.er, args.gmin,args.gmax, start, rx, args.rx_component, args.width, args.height,isnorm = args.norm, rawdata = args.rawdata)
    plthandle.show()
