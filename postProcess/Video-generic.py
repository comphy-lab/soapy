# Author: Vatsal Sanjay
# vatsalsy@comphy-lab.org
# CoMPhy Lab
# Physics of Fluids Department
# Last updated: Mar 8, 2025

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import multiprocessing as mp
from functools import partial
import argparse  # Add at top with other imports
import sys


matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def gettingFacets(filename,includeCoat='true'):
    exe = ["./getFacet", filename, includeCoat]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2, z2)))
                    segs.append(((-r1, z1),(-r2, z2)))
                    skip = True
    return segs

def gettingfield(filename, zmin, zmax, rmax, nr):
    exe = ["./getData", filename, str(zmin), str(0), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    # print(temp2) #debugging
    Rtemp, Ztemp, Ttemp  = [],[],[]

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            Ttemp.append(float(temp3[2]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    T = np.asarray(Ttemp)
    nz = int(len(Z)/nr)

    # print("nr is %d %d" % (nr, len(R))) # debugging
    print("nz is %d" % nz)

    R.resize((nz, nr))
    Z.resize((nz, nr))
    T.resize((nz, nr))

    return R, Z, T, nz
# ----------------------------------------------------------------------------------------------------------------------

def process_timestep(ti, caseToProcess, folder, tsnap, GridsPerR, rmin, rmax, zmin, zmax, lw):
    t = tsnap * ti
    place = f"{caseToProcess}/intermediate/snapshot-{t:.4f}"
    name = f"{folder}/{int(t*1000):08d}.png"

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return

    if os.path.exists(name):
        print(f"{name} Image present!")
        return

    segs = gettingFacets(place)

    if not segs:
        print(f"Problem in the available file {place}")
        return

    nr = int(GridsPerR * rmax)
    R, Z, T, nz = gettingfield(place, zmin, zmax, rmax, nr)
    zminp, zmaxp, rminp, rmaxp = Z.min(), Z.max(), R.min(), R.max()

    # Plotting
    AxesLabel, TickLabel = 50, 20
    fig, ax = plt.subplots()
    fig.set_size_inches(19.20, 10.80)

    ax.plot([0, 0], [zmin, zmax], '-.', color='grey', linewidth=lw)
    ax.plot([rmin, rmin], [zmin, zmax], '-', color='black', linewidth=lw)
    ax.plot([rmin, rmax], [zmin, zmin], '-', color='black', linewidth=lw)
    ax.plot([rmin, rmax], [zmax, zmax], '-', color='black', linewidth=lw)
    ax.plot([rmax, rmax], [zmin, zmax], '-', color='black', linewidth=lw)

    line_segments = LineCollection(segs, linewidths=4, colors='blue', linestyle='solid')
    ax.add_collection(line_segments)

    cntrl2 = ax.imshow(T, interpolation='Bilinear', cmap='hot_r', origin='lower', extent=[rminp, rmaxp, zminp, zmaxp], vmax=1.0, vmin=0.0)

    ax.set_aspect('equal')
    ax.set_xlim(rmin, rmax)
    ax.set_ylim(zmin, zmax)
    ax.set_title(f'$t/\\tau_\\gamma$ = {t:4.3f}', fontsize=TickLabel)

    ax.axis('off')

    plt.savefig(name, bbox_inches="tight")
    plt.close()

def main():
    # Get number of CPUs from command line argument, or use all available
    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=8, help='Number of CPUs to use')
    parser.add_argument('--nGFS', type=int, default=4000, help='Number of restart files to process')
    parser.add_argument('--GridsPerR', type=int, default=128, help='Number of grids per R')
    parser.add_argument('--ZMAX', type=float, default=4.0, help='Maximum Z value')
    parser.add_argument('--RMAX', type=float, default=4.0, help='Maximum R value')
    parser.add_argument('--ZMIN', type=float, default=0.0, help='Minimum Z value')
    parser.add_argument('--tsnap', type=float, default=0.01, help='Time snap')
    parser.add_argument('--caseToProcess', type=str, default='../simulationCases/soapBubble-full', help='Case to process')  
    parser.add_argument('--folderToSave', type=str, default='Video', help='Folder to save')
    args = parser.parse_args()

    num_processes = args.CPUs
    nGFS = args.nGFS
    tsnap = args.tsnap
    ZMAX = args.ZMAX
    RMAX = args.RMAX
    ZMIN = args.ZMIN
    
    rmin, rmax, zmin, zmax = [-RMAX, RMAX, ZMIN, ZMAX]
    GridsPerR = args.GridsPerR


    lw = 2
    folder = args.folderToSave
    caseToProcess = args.caseToProcess

    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        # Create partial function with fixed arguments
        process_func = partial(process_timestep, caseToProcess=caseToProcess,
                             folder=folder, tsnap=tsnap,
                             GridsPerR=GridsPerR, rmin=rmin, rmax=rmax, 
                             zmin=zmin, zmax=zmax, lw=lw)
        # Map the process_func to all timesteps
        pool.map(process_func, range(nGFS))

if __name__ == "__main__":
    main()
