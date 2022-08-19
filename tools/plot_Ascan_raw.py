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
from gprMax.receivers import Rx

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
#mpl_plot(args.tracenumber,scanfile, args.timescan,outputfile, outputdata, dt,args.ts, rx, args.rx_component)
def mpl_plot(filename, outputs=Rx.defaultoutputs, rawdata = False):
   
    #Read A-scan file
        # Open output file and read some attributes
    f = h5py.File(filename, 'r')
    nrx = f.attrs['nrx']
    dt = f.attrs['dt']
    iterations = f.attrs['Iterations']
    time = np.linspace(0, (iterations - 1) * dt, num=iterations)   
    # Check there are any receivers
    if nrx == 0:
        raise CmdInputError('No receivers found in {}'.format(filename))
    # New plot for each receiver
    for rx in range(1, nrx + 1):
        path = '/rxs/rx' + str(rx) + '/'
        availableoutputs = list(f[path].keys())
        # Check for polarity of output and if requested output is in file
        if outputs[0][-1] == '-':
            polarity = -1
            outputtext = '-' + outputs[0][0:-1]
            output = outputs[0][0:-1]
        else:
            polarity = 1
            outputtext = outputs[0]
            output = outputs[0]

        if output not in availableoutputs:
            raise CmdInputError('{} output requested to plot, but the available output for receiver 1 is {}'.format(output, ', '.join(availableoutputs)))

        outputdata = f[path + output][:] * polarity
        #Simulated text file
        if rawdata:
            a_file = open(filename + "_raw.csv", "w")
            a_file.writelines("Time(ns),Ez(V/m)\n")
            for i in range(len(outputdata)):
                #if time[i] >= 0: 
                rowdata = str(time[i]) + "," + str(outputdata[i]) + "\n"
                a_file.writelines(rowdata)
            a_file.close()
        outputdata = norm(outputdata)
        time = time*1e+9
        #Simulated text file
        if rawdata:
            a_file = open(filename + "_norm.csv", "w")
            a_file.writelines("Time(ns),Normalized Amplitude\n")
            for i in range(len(outputdata)):
                #if time[i] >= 0: 
                rowdata = str(time[i]) + "," + str(outputdata[i]) + "\n"
                a_file.writelines(rowdata)
            a_file.close()
        #a_file = open("measured.txt", "w")
        #np.savetxt(a_file, zs)
        #a_file.close()
        #a_file = open("tmeasured.txt", "w")
        #np.savetxt(a_file, te1)
        #a_file.close()
        # Plot time history of output component
       
        # plotting the points  
        plt.plot(time, outputdata, color='red', linewidth = 2,label='Simulated') 
        # setting x and y axis range 
        plt.ylim(-1,1) 
        plt.xlim(0,np.max(time)) 
        
        plt.xticks(style='normal', fontsize=18)
        plt.yticks(style='normal', fontsize=18)
        # naming the x axis 
        plt.xlabel('Time [ns]', style='normal', fontsize=24) 
        # naming the y axis 
        plt.ylabel('Normalized Amplitude [-]', style='normal', fontsize=24) 
        plt.legend(frameon=True, loc='upper center', ncol=2, fontsize=18)
        # giving a title to my graph 
        #plt.title('Comparison of A-scan', style='normal', fontsize=30)
            # Close Stream File
    f.close()

    # plot filled contour map with 100 levels
    # Save a PDF/PNG of the figure
    # savefile = os.path.splitext(filename)[0]
    # fig.savefig(path + os.sep + savefile + '.pdf', dpi=None, format='pdf', bbox_inches='tight', pad_inches=0.1)
    # fig.savefig(path + os.sep + savefile + '.png', dpi=150, format='png', bbox_inches='tight', pad_inches=0.1)

    return plt


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Plots a comparison of simulated and measured A-scan image.', usage='cd gprMax; python -m tools.plot_Ascan outputfile output')
    parser.add_argument('outputfile', help='name of output file including path A-scan')
    parser.add_argument('--outputs', help='outputs to be plotted', default=Rx.defaultoutputs, choices=['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz', 'Ix', 'Iy', 'Iz', 'Ex-', 'Ey-', 'Ez-', 'Hx-', 'Hy-', 'Hz-', 'Ix-', 'Iy-', 'Iz-'], nargs='+')
    parser.add_argument('-rawdata', action='store_true', help='Write rawdata', default=False)
    args = parser.parse_args()

    
    # Open output file and read number of outputs (receivers)
    plthandle = mpl_plot(args.outputfile,args.outputs,rawdata = args.rawdata)
       

    plthandle.show()
