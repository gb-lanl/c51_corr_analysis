import gvar as gv
import numpy as np
import os 
import h5py
import time
import re
import glob 
directory = './data/C13/'
N_cnf = len([name for name in os.listdir(directory) if os.path.isfile(name)])

dirs = os.listdir( directory )
print(dirs)
data_file = directory+dirs[0]
print(data_file)
fit_states = ['proton', 'pion', 'nucleon']
#bs_seed = 'a12m220XL

def get_corr(file_h5,fit_type, species, species_momentum, cnfg_num):
    data = {}
    with h5py.File(data_file,"r") as f:
        path = "/"+ fit_type + "/" + species + "/src5.0_snk5.0" + "/" + species_momentum + "/C13.b_"+ cnfg_num
    
        cfgs = np.array(f[path][()])
        print(cfgs)
        data[species] = cfgs.real # .real for the _hp ensembles 
    return data
get_corr(data_file,'2pt', 'proton', 'proton', '5682')
corr_lst = {
    'proton':{
        'dsets':['2pt/proton/src5.0_snk5.0/proton/C13.b_5682/AMA'], #zero momentum
        'weights'  :[1],
        't_reverse':[False],
        'fold'     :True,
        'snks'     :['S', 'P'],
        'srcs'     :['S'],
        'xlim'     :[0,48.5],
        'ylim'     :[0.12,0.169],
        'colors'   :{'SS':'#70bf41','SP':'k'},
        'type'     :'cosh',
        'ztype'    :'z_snk z_src',
        'z_ylim'   :[0.055,0.26],
        # fit params
        'n_state'  :2,
        'T'        :96,
        't_range'  :np.arange(8,48),
        't_sweep'  :range(3,28),
        'n_sweep'  :range(1,6),
        'eff_ylim' :[0.133,0.1349]
    },
    'pion':{
        'dsets':['2pt/pion/src5.0_snk5.0/pion/C13.b_5682/AMA'], #zero momentum
        'weights'  :[1],
        't_reverse':[False],
        'fold'     :True,
        'snks'     :['S', 'P'],
        'srcs'     :['S'],
        'xlim'     :[0,48.5],
        'ylim'     :[0.12,0.169],
        'colors'   :{'SS':'#70bf41','SP':'k'},
        'type'     :'cosh',
        'ztype'    :'z_snk z_src',
        'z_ylim'   :[0.055,0.26],
        # fit params
        'n_state'  :3,
        'T'        :96,
        't_range'  :np.arange(5,48),
        't_sweep'  :range(2,28),
        'n_sweep'  :range(1,6),
        'eff_ylim' :[0.133,0.1349]
    },
}

priors = gv.BufferDict()
x      = dict()

priors['proton_E_0']  = gv.gvar(0.5, .06)
priors['proton_zS_0'] = gv.gvar(2.0e-5, 1.e-5)
priors['proton_zP_0'] = gv.gvar(2.5e-3, 1.e-3)

priors['pion_E_0']  = gv.gvar(0.14, .006)
priors['pion_zS_0'] = gv.gvar(5e-3, 5e-4)
priors['pion_zP_0'] = gv.gvar(0.125,  0.015)

for corr in corr_lst:#[k for k in corr_lst if 'mres' not in k]:
    for n in range(1,10):
        # use 2 mpi splitting for each dE

        # E_n = E_0 + dE_10 + dE_21 +...
        # use log prior to force ordering of dE_n
        priors['log(%s_dE_%d)' %(corr,n)] = gv.gvar(np.log(2*priors['pion_E_0'].mean), 0.7)

        # for z_P, no suppression with n, but for S, smaller overlaps
        priors['%s_zP_%d' %(corr,n)] = gv.gvar(priors['%s_zP_0' %(corr)].mean, 2*priors['%s_zP_0' %(corr)].sdev)
        zS_0 = priors['%s_zS_0' %(corr)]
        if n <= 2:
            priors['%s_zS_%d' %(corr,n)] = gv.gvar(zS_0.mean, 2*zS_0.sdev)
        else:
            priors['%s_zS_%d' %(corr,n)] = gv.gvar(zS_0.mean/2, zS_0.sdev)
    # x-params
    for snk in corr_lst[corr]['snks']:
        sp = corr_lst[corr]['srcs'][0]+snk
        state = corr+'_'+sp
        x[state] = dict()
        x[state]['state'] = corr
        for k in ['type', 'T', 'n_state', 't_range', 'eff_ylim', 'ztype']:
            if k in corr_lst[corr]:
                x[state][k] = corr_lst[corr][k]
        if 't0' in corr_lst[corr]:
            x[state]['t0'] = corr_lst[corr]['t0']
        else:
            x[state]['t0'] = 0
        if 'mres' not in corr:
            x[state]['color'] = corr_lst[corr]['colors'][sp]
            x[state]['snk']   = snk
            x[state]['src']   = corr_lst[corr]['srcs'][0]
        else:
            x[state]['color'] = corr_lst[corr]['colors']

    print(x)

