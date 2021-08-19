#!/usr/bin/env python3
"""
"""
# Stack frequency plot (subplots)

# Load files from each of the folders using pickle and use those to make RTI plots.
#3new1. Keep only those schedules which occur at least 5 times  throughout the day
#3new2. Keep only those schedules which occur at least 3 times continuously throughout the day.

import numpy as n
import matplotlib.pyplot as plt
import glob
import h5py
import scipy.constants as c
import chirp_config as cc
import chirp_det as cd
import sys
import os
import time
import shutil
import datetime
from datetime import timezone
#from datetime import datetime
from numpy import unravel_index
import pickle
import matplotlib as mpl
import math


import cartopy.crs as ccrs
import ipdb

freqlist = [60, 80, 100, 120, 160, 180]
#rootdir = '/home/dev/Downloads/chirp_juha2b'
rootdir = '/media/dev/Seagate Backup Plus Drive/lfm_files'
# for subdir, dirs, files in os.walk(rootdir):
#dirs = sorted(os.listdir(rootdir))

output_dir1 = "/home/dev/Downloads/chirp_juha2b/Plots20"
output_dir2 = "/home/dev/Downloads/chirp_juha2b/Plots20/AllRTI"

output_dir21 = "/home/dev/Downloads/chirp_juha2b/Plots23"
dirs = sorted(os.listdir(output_dir1))


def k_largest_index_argsort(S, k):
    idx = n.argsort(S.ravel())[:-k-1:-1]
    return n.column_stack(n.unravel_index(idx, S.shape))


# def save_var(img_fname="img_1b.png"):
def save_var(DataDict):
    print('check')
    
    path1 = rootdir + '/' + dirs1 + '/' + dirs1[5:10] + 'b.data'
    path2 = output_dir1 + '/' + dirs1 + '/' + dirs1[5:10] + 'k.data'

    img_fname1 = "%s/%s/RTIk-%s.png" % (output_dir1, dirs1, dirs1[0:10])
    img_fname2 = "%s/RTI-%sk.png" % (output_dir2, dirs1[0:10])

  
    with open(path2, 'rb') as f:
        #T03, dB3, dB3a, dB3b, dB3c, range_gates, range_gates2, range_gates3, freqs = pickle.load(f)
        DataDict = pickle.load(f)
    
    # Check if there is change in schedule -- get the modulus, round it and take the differences of all values -- if any change, one or more of the differences will not be zero

    T03 = DataDict['Time']
    print(len(T03))
    T03old = len(T03)
    T03old1 = T03
    path3 = os.path.join(output_dir21, cd.unix2dirname(T03[0]))  + 'a'
    path4 = os.path.join(output_dir21, cd.unix2dirname(T03[0]))  + 'b'
    if not os.path.exists(path3):
            os.makedirs(path3)
    if not os.path.exists(path4):
            os.makedirs(path4)
    
    for jj,ii in enumerate(T03):
            img_fname2a  = glob.glob("%s/lfm*-%1.2f.png"%(path,T03[jj]))
            if img_fname2a[0] not in path3:
                shutil.copy(img_fname2a[0], path3)   
    
    x1 = [x % 720 for x in T03]
    x3 = n.array([int(round(x)) for x in x1])
    x3old = x3
    (x3u,C) = n.unique(x3,return_counts = True)
    x3c = x3u[C > 5]
    
    ## Check if x3c have three continous values for its elements
    x3c1 = {}
    for i, j in enumerate(x3c):
        AN = n.where(x3==x3c[i])[0]
        AN1 = n.where(n.diff(AN)==1)
        AN2 = AN[AN1[0]]
        if len(AN2) > 3:
            x3c1[i] = x3c[i]
    
    # Get the elements from x3c1                
    x3c2 = []
    for j in x3c1.keys():
           x3c2.append(x3c1[j])
    
    # Rename x3c2 as x3c as we would want to work with x3c        
    x3c = n.array(x3c2)

    
    x3n = []
    for i, j in enumerate(x3c):
            x3n.append(n.where(x3 == x3c[i]))    
    
    if len(x3n) == 0: 
        print('No useful data')     
        sys.exit()

    x3n = n.concatenate(x3n,axis=None).astype(int)    
    x3in = n.array(range(0,len(x3)))
    x3diff = n.setdiff1d(x3in, x3n)
        
    range_gates3 = DataDict['range_gates3']
    range_gates2 = DataDict['range_gates2']
    freqlist = DataDict['freqlist']
    freqs = DataDict['freqs']
    
    T03 = n.delete(T03, x3diff)

    for jj,ii in enumerate(T03):
        img_fname2b  = glob.glob("%s/lfm*-%1.2f.png"%(path,T03[jj]))
        if img_fname2b[0] not in path4:
            shutil.copy(img_fname2b[0], path4)  
    
    for k in [j for j in DataDict['DBall'].keys()]:
        DataDict['DBall'][k]  = n.delete(DataDict['DBall'][k],x3diff,1)
        #DataDict['DBall'][k]  = DataDict['DBall'][k][:,x3n]
        
    print(len(T03))
    range_gates3 = n.delete(range_gates3,x3diff,1)
    #range_gates3 = range_gates3[:,x3n]
    x1 = [x % 720 for x in T03]
    x3 = [round(y, 2) for y in x1]
    x4 = n.diff([x3])[0]
    x3new = x3
    
    dtz = datetime.datetime.utcfromtimestamp(T03[-1])
    dtz2 = dtz.replace(hour=23,minute=59, second=59)
    T03z = dtz2.replace(tzinfo=timezone.utc).timestamp()

    T03a = n.arange(T03[0], T03z, 720)
    T03b = T03a

    if len(n.where(abs(x4) > 1)[0]) > 0:
        x5 = n.array([n.where(abs(x4) > 1)])[0][0]+1

    sch_ch1 = n.where(abs(x4)>1)[0]
    SM = n.split(sch_ch1,n.argwhere(n.diff(sch_ch1)>1)[:,0]+1)   

    # If T03a doesn't span the whole day, build it. 
    if len(T03a) < 120:
        dtest = int(datetime.datetime.utcfromtimestamp(T03[0]).strftime("%d"))
        jj1 = 0
        while True:
            # subtract 720 until it gets to beginning of the day and in next while loop, add 720 until it gets to the end of the day
            print('jj1=%d' %(jj1))
            if abs((int(datetime.datetime.utcfromtimestamp(T03[0] - 720*jj1).strftime("%d")) - dtest)) > 0:
                T03b = n.insert(T03b, 0, T03[0] - 720*(jj1-1))
                print('jj1= %1.2f' % (jj1))
                break
            jj1 += 1

        T03a = n.arange(T03b[0], T03z, 720)
        
    # And, apply the corrections for schedule change if it happens.     
    if (len(n.where(abs(x4) > 1)[0]) > 0):    
        # The change has happened at x5 and hence, I am looking for equivalent element in T03a corresponding to (x5-1) [which is SM] position of T03. That's because
        # upto that position, T03 is unchanged and I want to keep T03a unchanged upto that equivalent position which will be given by nn1!
        TT3 = []
        NNJ = []
        for j in range(0,len(SM)):
            NN = abs(T03a - T03[SM[j][0]])
            nn1 = n.where(NN == NN.min())[0][0]
            NNJ.append(nn1)
            if j == 0:
            # Going upto nn1 requires it to write [0:(nn1+1)] as the last element is not taken in Python. So, this portion is unchanged.
                T03a1 = T03a[0:(nn1+1)]
                TT3.append(T03a1)
            if j > 0:  
            # For more than one schedule change. We begin where we had left. That is one element ahead of the last element of 'previous TT3' and add the change to all elements 
            # for this batch. The length of the batch will be determined by len(T03aa1) : which is [NNJ[j-1] + 1]th position  upto new nn1 
            # And repeat the process. 
                T03aa1 = T03a[NNJ[j-1]+1:(nn1+1)] + x4[x5[j-1]-1]
                Element = TT3[j-1][-1]+720 + x4[x5[j-1]-1]
                T03a1 = n.arange(Element, Element + 720*(len(T03aa1)-1)+20, 720)
                TT3.append(T03a1)

        # Apply the last change here : Start it from the (final+1 = TT3a[-1] + 720) position of TT3 nn1 position. And, add the "change = x4[x5-1]" at (nn1 + 1)th position
        TT3a = n.hstack(TT3)
        #T03a2 = n.arange(T03a[nn1+1] + x4[x5[-1]-1], T03z, 720)
        T03a2 = n.arange(TT3a[-1] + 720 + x4[x5[-1]-1], T03z, 720)
        T03a = n.concatenate((TT3a, T03a2))
        # If for some reasons (it can happen as the 'space' might have shrunk due to schedule change to accomodate more than 120 elements between the start and the end of the day), 
        # the length of T03a exceeds 120 , keep only upto 120. 
        if len(T03a) > 120:
            T03a = T03a[0:120]
    
    # Construct full dB3 and ranges_gatesnew with NaNs
    dB3test = n.full([3999, 120], None)
    dB3test[:] = n.NaN

    DataDict['DBallnew'] = {}
    for k in [j for j in DataDict['DBall'].keys()]:
        DataDict['DBallnew'][k] = n.full([3999, 120], None)

    range_gatestest = n.full([3999, 120], None)
    range_gatestest[:] = n.NaN
    range_gatesnew = n.full([3999, 120], None)
    
    CT = 0
    CT1 = 0
    for i, x in enumerate(T03a):
        #print(i)
        DIFF = abs(T03 - x)
        MIN = min(abs(T03-x))
        if MIN < 2:

            CT += 1
            ij = n.where(DIFF == n.amin(DIFF))[0][0]
            
            for k in [j for j in DataDict['DBall'].keys()]:
                DataDict['DBallnew'][k][:, i] = DataDict['DBall'][k][:, ij]

            range_gatesnew[:, i] = range_gates3[:, ij]
            # it (ij) comes out as a tuple which contains an array. So, the first index [0] gets
            # the array out of the tuple. And the second index [0] gets the index out of the array.
        else:
            CT1 += 1
            for k in [j for j in DataDict['DBall'].keys()]:
                DataDict['DBallnew'][k][:, i] = dB3test[:, i]

            range_gatesnew[:, i] = range_gates2

    FileName = os.path.join(path3,"Var.txt")
    file = open(FileName, "w")
    file.write("%s = %s\n" %("Date", str(cd.unix2dirname(T03[0]))))
    file.write("%s = %s\n" %("CT", str(CT)))
    file.write("%s = %s\n" %("CT1", str(CT1)))
    file.write("%s = %s\n" %("T03_old", str(T03old)))
    file.write("%s = %s\n" %("T03_new", str(len(T03))))
    file.write("%s = %s\n" %("T03_old1", str(T03old1)))
    file.write("%s = %s\n" %("T03_new1", str(T03)))
    file.write("%s = %s\n" %("x3_old", str(x3old)))
    file.write("%s = %s\n" %("x3_new", str(x3new)))
    file.close()
    
    #fig = plt.figure(figsize=(1.5*10, 1.5*3))
    fig = plt.figure(figsize=(1.5*6, 1.5*12))
  
    fig.suptitle("RTI Plots : %s UTC" % datetime.datetime.utcfromtimestamp(T03[1]).strftime('%Y-%m-%d'),weight='bold',fontsize=14)
    new_times = [datetime.datetime.utcfromtimestamp(x) for x in T03a]
    new_times = n.array(new_times)
    new_times1 = [datetime.datetime.fromtimestamp(x) for x in T03a]  # local-time
    new_times1 = n.array(new_times1)


    for j, k in enumerate(reversed([j for j in DataDict['DBall'].keys()])):
        #ax1 = fig.add_subplot(4, 1, k.astype(n.int64)-2)
        ax1 = fig.add_subplot(6, 1, j+1)
        for ja in range(0, 118):
            plt.pcolormesh(new_times[ja:ja+2], n.column_stack((range_gatesnew[:, ja], range_gatesnew[:, ja])),DataDict['DBallnew'][k].astype(n.float)[:-1, ja:ja+1], vmin=-3, vmax=30.0,cmap="inferno")
        
        cb = plt.colorbar()
        cb.set_label("SNR (dB)")
        #plt.title("RTI plot for %1.2f MHz" %(k))
        #plt.title("RTI plot for %1.2f MHz" %(freqs[freqlist[k.astype(n.int64)]]/1e6))
        #plt.xlabel("Time (UTC)")
        #plt.ylabel("One-way range offset (km)")
        plt.ylabel("%1.2f MHz\n Range (km)" %(k), weight='bold',fontsize=12)
        plt.ylim([0, 4000])
        plt.xlim(new_times[0], new_times[-1])
        plt.tight_layout(rect=[0, 0.07, 1, 0.99],pad=1.0)
        #plt.tight_layout(rect=[0.1, 0.4, 0.85, 0.99],pad=1.0)
        #plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        #[left, bottom, right, top] in normalized (0, 1) figure coordinates
        #plt.tight_layout(pad=1.0)
        #plt.savefig(img_fname1, bbox_inches='tight')
        #plt.savefig(img_fname2, bbox_inches='tight')

    #cb = plt.colorbar()
    plt.xlabel("Time (UTC)", weight='bold', fontsize=12)
    plt.savefig(img_fname1, bbox_inches='tight')
    plt.savefig(img_fname2, bbox_inches='tight')
    #plt.show()
    # plt.savefig(img_fname1)
    #ipdb.set_trace()
    plt.close()
    plt.clf()
    # ho.close()

if __name__ == "__main__":
    if len(sys.argv) == 2:
        conf = cc.chirp_config(sys.argv[1])
    else:
        conf = cc.chirp_config()

    if conf.realtime:
        while True:
            fl = glob.glob("%s/*/lfm*.h5" % (conf.output_dir))
            fl.sort()
            for f in fl:
                plot_ionogram(conf, f)
            time.sleep(10)
            save_var()
    else:
        for j in range(0, len(dirs)):
            dirs1 = dirs[j]
            dtt1 = datetime.datetime.strptime('2021-05-09','%Y-%m-%d').date()
            dtt2 = datetime.datetime.strptime(dirs1[0:10],'%Y-%m-%d').date()

            #if dirs1[0:10] == '2021-05-02':
            #dir1 = output_dir1 + '/' + dirs[j]
            #for x in os.listdir(dir1):
            #    if x.endswith("h.data"):

            #if dirs1[0:4] == '2021': 
            if dtt2 > dtt1 :
                    path = os.path.join(rootdir, dirs1)
                    print(dirs1)
                    os.chdir(path)
                    fl = glob.glob("%s/lfm*.h5" % (path))
                    fl.sort()

                    DataDict = {}

                # fl = glob.glob("%s/*/lfm*.h5" % (conf.output_dir))
                # fl = glob.glob("%s/lfm*.h5" % (conf.output_dir))
                
                    if len(fl) > 1:
                        save_var(DataDict) 
