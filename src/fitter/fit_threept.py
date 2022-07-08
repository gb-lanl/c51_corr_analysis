import os, sys, argparse
import numpy as np
np.set_printoptions(linewidth=180)
import tables as h5
import warnings
warnings.simplefilter('ignore', h5.NaturalNameWarning)
from glob import glob
import importlib

import os
import numpy as np
import tables as h5
from glob import glob
import sys
import random

def time_reverse(corr,phase=1,time_axis=1):
    '''
    Performe time-reversal of correlation function accounting for BC in time
    phase     = +1 or -1
    time_axis = time_axis in array
    '''
    if len(corr.shape) > 1:
        if time_axis in [0,1]:
            #cr = phase * np.roll(corr[:,::-1],1,axis=time_axis)
            cr = phase * np.roll(np.flip(corr,axis=time_axis),1,axis=time_axis)
            if time_axis == 0:
                cr[0] = phase * cr[0]
            elif time_axis == 1:
                cr[:,0] = phase * cr[:,0]
            else:
                pass
        else:
            print('time_axis != 0,1 not supported at the moment')
            sys.exit()
    else:
        cr = phase * np.roll(np.flip(corr,axis=0),1)
        cr[0] = phase * cr[0]
    return cr

def src_split(src):
    x0 = src.split('x')[1].split('y')[0]
    y0 = src.split('y')[1].split('z')[0]
    z0 = src.split('z')[1].split('t')[0]
    return 'x%s_y%s_z%s_t%s' %(x0,y0,z0,t0)

def xyzt(src):
    x = src.split('x')[1].split('y')[0]
    y = src.split('y')[1].split('z')[0]
    z = src.split('z')[1].split('t')[0]
    t = src.split('t')[1]
    return x,y,z,t

def mom_avg(h5_data,state,mom_lst,weights=False):
    '''
    perform a momentum average of a state from an open h5 file
    data file is assumed to be of shape [Nt,Nz,Ny,Nx,[re,im]]
    data_mom = h5_data[state][:,pz,py,px]
    '''
    d_lst = []
    w = []
    for mom in mom_lst:
        px,py,pz = mom['momentum']
        w.append(mom['weight'])
        #print(state)
        d_lst.append(h5_data[state][:,pz,py,px])
    d_lst = np.array(d_lst)
    w = np.array(w)
    if weights:
        for wi,we in enumerate(w):
            d_lst[wi] = we*d_lst[wi]
        d_avg = np.sum(d_lst,axis=0) / np.sum(w)
    else:
        d_avg = np.mean(d_lst,axis=0)
    return d_avg

mom_lst = []
for px in range(-2,3):
    for py in range(-2,3):
        for pz in range(-2,3):
            if px**2 + py**2 + pz**2 <= 5:
                mom_lst.append('px'+str(px)+'_py'+str(py)+'_pz'+str(pz))

particles = ['proton','proton_SP']

files = glob('./data/C13/*.h5')
srcs = []
for f in files:
    srcs.append(f.split('_')[-1].split('.')[0])
random.shuffle(srcs)
print(srcs)
if not os.path.exists('spec_src_avg'):
    os.makedirs('spec_src_avg')
fout = h5.open_file('spec_src_avg/spec_tslice_src_avg_px0py0pz0_Nsnk1_src_avg.h5','w')
for particle in particles:
    for mom in mom_lst:
        print(particle,mom)
        data = []
        for src in srcs:
            x0,y0,z0,t0 = xyzt(src)
            fin = h5.open_file('spec/spec_px0py0pz0_Nsnk1_'+src+'.h5','r')
            for spin in ['spin_up','spin_dn']:
                h5_path = '/pt/spec/'+particle+'/'+spin+'/'+src_split(src)+'/'+mom
                src_data = fin.get_node(h5_path).read()
                data.append(src_data)
            fin.close()
            #src_data = np.roll(src_data,-int(t0))
        data = np.array(data)
        data_avg = data.mean(axis=0)
        if '_np' in particle:
            print('    time_reversing')
            print('        spec needs time flip')
            data_avg = time_reverse(data_avg,phase=-1,time_axis=0)
        data_slice = data_avg[0:5]
        h5_avg_path = '/pt/spec/'+particle+'/spin_avg'
        fout.create_group(h5_avg_path,mom,createparents=True)
        fout.create_array(h5_avg_path+'/'+mom,'src_avg',data_slice)
        fout.flush()
fout.close()


# PROTON
# within 2pt/ dir, collect correlator dsets for proton, proton_sp, pion, pion_sp
sp   = ['sh','pt']
spin = ['spin_up','spin_dn']
par  = ['proton','proton_sp']
get_data = False
spec_dir = './data/C13'
for s in spin:
    if src not in f5.get_node(spec_dir+'2pt/'+'/proton/'+s) or src not in f5.get_node(spec_dir+'/proton_np/'+s):
        get_data = True
    if args.o and (src in f5.get_node(spec_dir+'/proton/'+s) or src in f5.get_node(spec_dir+'/proton_np/'+s)):
        get_data = True
if get_data:
    fin = h5.open_file(ftmp,'r')
    for corr in par:
        for s in spin:
            tmp_dir = spec_dir+'/'+corr+'/'+s
            p_lst = utils.p_lst(params['BARYONS_PSQ_MAX'])
            for mom in p_lst:
                if mom not in f5.get_node(tmp_dir):
                    f5.create_group(tmp_dir,mom)
                mom_dir = tmp_dir + '/'+mom
                pt = fin.get_node('/pt/'+corr+'/'+s+'/'+src_split+'/'+mom).read()
                sh = fin.get_node('/sh/'+corr+'/'+s+'/'+src_split+'/'+mom).read()
                nt = len(pt)
                data = np.zeros([nt,2,1],dtype=dtype)
                data[:,0,0] = sh
                data[:,1,0] = pt
                if not np.any(np.isnan(data)):
                    if args.v:
                        print(no,corr,s,src,mom)
                    if src not in f5.get_node(mom_dir):
                        f5.create_array(mom_dir,src,data)
                    elif src in f5.get_node(mom_dir) and args.o:
                        f5.get_node(mom_dir+'/'+src)[:] = data
                    elif src in f5.get_node(mom_dir) and not args.o:
                        print('  skipping proton: overwrite = False',no,src)
                else:
                    print('  NAN',no,src)
    fin.close()
    f5.close()
else:
    f5.close()
Footer
