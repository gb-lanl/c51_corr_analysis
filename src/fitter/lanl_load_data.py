import tables as h5
#import h5py as h5
import numpy as np
import gvar as gv
import sys
import os
from pathlib import Path
import copy
'''
TIME REVERSE
'''

# Read group structure
#-------------------------
# data = {}
# data_file_list = 
# cnf_data = np.load(data_file_list[0], allow_pickle=True).tolist()
# for g1 in cnf_data.keys():
#   data[g1] = {}
#   for g2 in cnf_data[g1].keys():
#     data[g1][g2] = {}
#     for g3 in cnf_data[g1][g2].keys():
#       data[g1][g2][g3] = {}
#       for g4 in cnf_data[g1][g2][g3].keys():
#         data[g1][g2][g3][g4] = [None]*N_cnf

def read_dataset(inputfiles, grep=None, keys=None, h5group=None, binsize=1, tcol=0, Gcol=1,):
    if not hasattr(inputfiles, 'keys'):
        # inputfiles is a filename or list of filenames (or files)
        if h5group == [] or h5group is None:   # not needed after next gvar update
            h5group = '/'
        dset = gv.dataset.Dataset(
            inputfiles, binsize=binsize, grep=grep,
            h5group=h5group, keys=keys,
            )
    return dset 


def time_reverse(corr, reverse=True, phase=1, time_axis=1):
    ''' assumes time index is second of array
        assumes phase = +- 1
    '''
    if reverse:
        if len(corr.shape) > 1:
            cr = phase * np.roll(corr[:, ::-1], 1, axis=time_axis)
            cr[:, 0] = phase * cr[:, 0]
        else:
            cr = phase * np.roll(corr[::-1], 1)
            cr[0] = phase * cr[0]
    else:
        cr = phase * corr
    return cr

def load_h5(f5_file, corr_dict, return_gv=True, rw=None, bl=1, uncorr_corrs=False, uncorr_all=False, verbose=True):
    corrs = []

    # check if f5_file is list
    if not isinstance(f5_file, list):
        f5_files = [f5_file]
    else:
        f5_files = f5_file
    # check for re-weighting
    if rw:
        rw_file, rw_path = rw
        if not isinstance(rw_file, list):
            rw_files = [rw_file]
        else:
            rw_files = rw_files
        if len(rw_files) != len(f5_files):
            sys.exit(
                'You must supply the same number of re-weighting files as data files')
        with h5.open_file(rw_files[0], 'r') as rw5:
            reweight = rw5.get_node('/'+rw_path).read()
        for f_i in range(1, len(f5_files)):
            with h5.open_file(rw_files[f_i], 'r') as rw5:
                reweight = np.concatenate(
                    (reweight, rw5.get_node('/'+rw_path).read()), axis=0)
        # normalize rw factors
        reweight = reweight / reweight.sum()

    # collect correlators
    for corr in corr_dict:
        print(corr)
        weights   = corr_dict[corr]['weights']
        t_reverse = corr_dict[corr]['t_reverse']
        # dsets = corr_dict[corr]['dsets']
        # check if data is in an array or single correlators
        if 'corr_array' not in corr_dict[corr]:
            corr_array = True
        else:
            corr_array = corr_dict[corr]['corr_array']
        if corr_array:
            # get first data
            with h5.open_file(f5_files[0], 'r') as f5:
                dsets = corr_dict[corr]['dsets']
                data = np.zeros_like(f5.get_node('/'+dsets[0]).read())
                # print(data)
                for i_d, dset in enumerate(dsets):
                    if 'phase' in corr_dict[corr]:
                        phase = corr_dict[corr]['phase'][i_d]
                    else:
                        phase = 1
                    d_tmp = f5.get_node('/'+dset).read()
                    #print(d_tmp)
                    data = d_tmp
                    print(data.shape)
                    
                    data += weights[i_d] * \
                        time_reverse(
                            d_tmp, reverse=t_reverse[i_d], phase=phase)
            # if we have more than 1 data file
            if len(f5_files) > 1:
                for f_i in range(1, len(f5_files)):
                    with h5.open_file(f5_files[f_i], 'r') as f5:
                        tmp = np.zeros_like(f5.get_node('/'+dsets[0]).read())
                        for i_d, dset in enumerate(dsets):
                            if 'phase' in corr_dict[corr]:
                                phase = corr_dict[corr]['phase'][i_d]
                            else:
                                phase = 1
                            d_tmp = f5.get_node('/'+dset).read()
                            tmp += weights[i_d] * time_reverse(
                                d_tmp, reverse=t_reverse[i_d], phase=phase)
                        # NOTE - we assume the cfg axis == 0
                        data = np.concatenate((data, tmp), axis=0)

            # if fold
            if corr_dict[corr]['fold']:
                data = 0.5*(data + time_reverse(data))
            # populate into [Ncfg, Nt] arrays
            for i,snk in enumerate(corr_dict[corr]['snks']):
                for j,src in enumerate(corr_dict[corr]['srcs']):
                    #print(i,j,src,snk)
                    if 'normalize' in corr_dict[corr] and corr_dict[corr]['normalize']:
                        corrs[corr+'_'+snk+src] = data[:, :, i, j] / \
                            data.mean(axis=0)[0, i, j]
                    else:
                        corrs[corr+'_'+snk+src] = data[:,:,i,j]

        else:  # load individual corrs
            if corr_dict[corr]['type'] == 'mres':
                with h5.open_file(f5_files[0], 'r') as f5:
                    data_MP = f5.get_node('/'+corr_dict[corr]['dset_MP'][0]).read()
                    data_PP = f5.get_node('/'+corr_dict[corr]['dset_PP'][0]).read()
                    # stack the data so it can be treated like other dsets
                    data = np.stack((data_MP,data_PP),axis=-1)
                    if len(f5_files) > 1:
                        for f_i in range(1, len(f5_files)):
                            with h5.open_file(f5_files[f_i], 'r') as f5:
                                data_MP = f5.get_node('/'+corr_dict[corr]['dset_MP'][0]).read()
                                data_PP = f5.get_node('/'+corr_dict[corr]['dset_PP'][0]).read()
                                tmp     = np.stack((data_MP,data_PP),axis=-1)
                                data    = np.concatenate((data, tmp), axis=0)
                # if fold
                if corr_dict[corr]['fold']:
                    data = 0.5*(data + time_reverse(data))
                # add MP and PP data
                corrs[corr+'_MP'] = data[...,0]
                corrs[corr+'_PP'] = data[...,1]

            # load lanl 2pt data, which is structured differently from callat data
            # need to stack SS, PS dsets together    
            if corr_dict[corr]['type'] == 'cosh':
            # if 'pion' in corr_dict[corr]:
                with h5.open_file(f5_files[0], 'r') as f5:
                    data_SS_ = f5.get_node('/'+corr_dict[corr]['dsets'][0]).read()
                    data_PS_ = f5.get_node('/'+corr_dict[corr]['dsets'][1]).read()
                    data = np.stack((data_SS_,data_PS_),axis=0)
                    # print(data[1])
                    # data = np.stack((data_SS,data_PS),axis=1)
                for corr in corrs:

                    corrs[corr+'_SS'] = data[0]
                    corrs[corr+'_PS'] = data[1]

            elif corr_dict[corr]['stack']:
            # if 'proton' in corr_dict[corr]:
                with h5.open_file(f5_files[0], 'r') as f5:
                    data_SS = f5.get_node('/'+corr_dict[corr]['dsets'][0]).read()
                    data_PS = f5.get_node('/'+corr_dict[corr]['dsets'][1]).read()
                    data = np.stack((data_SS,data_PS),axis=0)
                    # print(data[1])
                    # data = np.stack((data_SS,data_PS),axis=1)
                for corr in corrs:
                    corrs[corr+'_SS'] = data[0]
                    corrs[corr+'_PS'] = data[1]
                
            # spin-averaging: combining spin u to spin u and spin dn to spin dn corr fcns 
            # (V4: + , A3: -) reduces stochastic uncertainty of the data 
            # https://arxiv.org/abs/2104.05226                       
            elif corr_dict[corr]['axial']:
                with h5.open_file(f5_files[0], 'r') as f5:
                    data_D = f5.get_node('/'+corr_dict[corr]['dsets'][0]).read()
                    data_U = f5.get_node('/'+corr_dict[corr]['dsets'][1]).read()
                    # print(data_D)
                    data=np.stack((data_U,data_D))
                    print(data)
                for corr in corrs:
                    corrs[corr] = data

                    #     print(np.subtract(data_U[i][1] ,data_D[i][1]))
                        # print(np.array(data_a3,dtype=np.float64))
                    # import IPython 
                    # IPython.embed()
                    # print(data_U,data_D)

            elif corr_dict[corr]['vector']:
                with h5.open_file(f5_files[0], 'r') as f5:
                    data_D = f5.get_node('/'+corr_dict[corr]['dsets'][0]).read()
                    data_U = f5.get_node('/'+corr_dict[corr]['dsets'][1]).read()
                    # print(data_D)A3
                    data=np.stack((data_U,data_D))
                for corr in corrs:
                    corrs[corr] = data

                

            # elif corr_dict['V4']['q_bilinear'] == True:
            #     with h5.open_file(f5_files[0], 'r') as f5:
            #         data_D_ = f5.get_node('/'+corr_dict['V4']['dsets'][0]).read()
            #         data_U_ = f5.get_node('/'+corr_dict['V4']['dsets'][1]).read()
            #         for i in range(len(data_D_)):
            #             data_v4 = np.add(data_U_[i][0] ,data_D_[i][0])
            #             corrs['V4'] = data_v4[:]

            else:
                for i, snk in enumerate(corr_dict[corr]['snks']):
                    for j, src in enumerate(corr_dict[corr]['srcs']):
                        with h5.open_file(f5_files[0], 'r') as f5:
                            d_set = dsets[0] % {'SNK': snk, 'SRC': src}
                            data = np.zeros_like(f5.get_node('/'+d_set).read())
                            for i_d, dset in enumerate(dsets):
                                d_set = dset % {'SNK': snk, 'SRC': src}
                                if 'phase' in corr_dict[corr]:
                                    phase = corr_dict[corr]['phase'][i_d]
                                else:
                                    phase = 1
                                d_tmp = f5.get_node('/'+d_set).read()
                                if corr_dict[corr]['t_reverse']:
                                    data += weights[i_d] * time_reverse(d_tmp, 
                                            reverse=t_reverse[i_d], phase=phase)
                                else:
                                    data = d_tmp

                        # if we have more than 1 data file
                        if len(f5_files) > 1:
                            for f_i in range(1, len(f5_files)):
                                with h5.open_file(f5_files[f_i], 'r') as f5:
                                    d_set = dsets[0] % {'SNK': snk, 'SRC': src}
                                    tmp = np.zeros_like(
                                        f5.get_node('/'+d_set).read())
                                    for i_d, dset in enumerate(dsets):
                                        d_set = dset % {'SNK': snk, 'SRC': src}
                                        if 'phase' in corr_dict[corr]:
                                            phase = corr_dict[corr]['phase'][i_d]
                                        else:
                                            phase = 1
                                        d_tmp = f5.get_node('/'+d_set).read()
                                        tmp += weights[i_d] * time_reverse(
                                            d_tmp, reverse=t_reverse[i_d], phase=phase)
                                    # NOTE - we assume the cfg axis == 0
                                    data = np.concatenate((data, tmp), axis=0)
                        # if fold
                        if corr_dict[corr]['fold']:
                            data = 0.5*(data + time_reverse(data))
                        # normalize?
                        if 'normalize' in corr_dict[corr] and corr_dict[corr]['normalize']:
                            #print('normalizing %s %s %s' %(corr,snk,src))
                            corrs[corr+'_'+snk+src] = data / data.mean(axis=0)[0]
                        else:
                            #print('not normalizing %s %s %s' %(corr,snk,src))
                            corrs[corr+'_'+snk+src] = data

    # re-weight?
    if rw:
        for k in corrs:
            corrs[k] = corrs[k] * reweight[:, None]

    # block/bin data
    if bl != 1:
        print('blocking data in units of saved configs: block length = %d' % bl)
        corrs_bl = {}
        for corr in corrs:
            corrs_bl[corr] = block_data(corrs[corr], bl)
        corrs = corrs_bl

    # return correlators
    if verbose:
        for corr in corrs:
            print(corr, corrs[corr].shape)
    if return_gv:
        if uncorr_corrs or uncorr_all:
            corrs_gv = {}
            if uncorr_all:
                for k in corrs:
                    corrs_gv[k] = gv.dataset.avg_data(corrs[k])
            else:
                for corr in corr_dict:
                    corrs_corr = {k: v for k, v in corrs.items() if corr in k}
                    tmp_gv = gv.dataset.avg_data(corrs_corr)
                    for k in tmp_gv:
                        corrs_gv[k] = tmp_gv[k]
        else:
            corrs_gv = gv.dataset.avg_data(corrs)
        return corrs_gv
    else:
        print(corrs)
        return corrs


def block_data(data, bl):
    ''' data shape is [Ncfg, others]
        bl = block length in configs
    '''
    ncfg, nt_gf = data.shape
    if ncfg % bl == 0:
        nb = ncfg // bl
    else:
        nb = ncfg // bl + 1
    corr_bl = np.zeros([nb, nt_gf], dtype=data.dtype)
    for b in range(nb-1):
        corr_bl[b] = data[b*bl:(b+1)*bl].mean(axis=0)
    corr_bl[nb-1] = data[(nb-1)*bl:].mean(axis=0)

    return corr_bl

