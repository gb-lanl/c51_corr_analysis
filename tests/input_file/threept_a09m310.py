import gvar as gv
import numpy as np

data_file = './data/hyper_spec_axial.h5'

#fit_states = ['mres-L','mres-S', 'pion', 'kaon', 'proton', 'omega']
fit_states = ['pion','proton', 'gA', 'gV']
bs_seed = 'a09m310'

corr_lst = {
    # PION
    'pion':{
        'dsets':['a09m310/pion'],
        'weights'  :[1],
        't_reverse':[False],
        'fold'     :True,
        'snks'     :['S', 'P'],
        'srcs'     :['S'],
        'xlim'     :[0,48.5],
        'ylim'     :[0.12,0.169],
        'colors'   :{'SS':'#70bf41','PS':'k'},
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
    # PROTON
    'proton':{
        'dsets':['a09m310/proton'],
        'weights'  :[1],
        't_reverse':[False],
        'fold'     :True,
        'snks'     :['S', 'P'],
        'srcs'     :['S'],
        'xlim'     :[0,48.5],
        'ylim'     :[0.12,0.169],
        'colors'   :{'SS':'#70bf41','PS':'k'},
        'type'     :'exp',
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
    # PROTON_A3
    'gA':{
        'dsets':['a09m310/proton_A3'],
        'weights'  :[1],
        't_reverse':[False],
        'fold'     :True,
        'snks'     :['S', 'P'],
        'srcs'     :['S'],
        'xlim'     :[0,48.5],
        'ylim'     :[0.12,0.169],
        'colors'   :{'SS':'#70bf41','PS':'k'},
        'type'     :'exp',
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
    # PROTON_V4
    'gV':{
        'dsets':['a09m310/proton_V4'],
        'weights'  :[1],
        't_reverse':[False],
        'fold'     :True,
        'snks'     :['S', 'P'],
        'srcs'     :['S'],
        'xlim'     :[0,48.5],
        'ylim'     :[0.12,0.169],
        'colors'   :{'SS':'#70bf41','PS':'k'},
        'type'     :'exp',
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
#proton priors
priors['proton_E_0']  = gv.gvar(0.5, .06)
priors['proton_zS_0'] = gv.gvar(2.0e-5, 1.e-5)
priors['proton_zP_0'] = gv.gvar(2.5e-3, 1.e-3)
#gA, gV priors
priors['gA_nm'] = np.array([gv.gvar(1.2, 0.2), gv.gvar(0,1), gv.gvar(0,1), gv.gvar(0,1)])
priors['gV_nm'] = np.array([gv.gvar(1.0, 0.2),gv.gvar(0,1),gv.gvar(0,1),gv.gvar(0,1)]) 
priors['dA_PS'] = np.array([gv.gvar(-4.4e-6,4.4e-6),gv.gvar(-4.4e-6,4.4e-6),gv.gvar(-4.4e-6,4.4e-6),gv.gvar(-4.4e-6,4.4e-6)]) 
priors['dA_SS'] = np.array([gv.gvar(-6.4e-7,6e-6),gv.gvar(-6.4e-7,6e-6),gv.gvar(-6.4e-7,6e-6),gv.gvar(-6.4e-7,6e-6)]) 
priors['dV_PS'] = np.array([gv.gvar(6.3e-6,6.3e-6),gv.gvar(6.3e-6,6.3e-6),gv.gvar(6.3e-6,6.3e-6),gv.gvar(6.3e-6,6.3e-6)])
priors['dV_SS'] = np.array([gv.gvar(1.2e-6,1.2e-6),gv.gvar(1.2e-6,1.2e-6),gv.gvar(1.2e-6,1.2e-6),gv.gvar(1.2e-6,1.2e-6)])
#pion priors
priors['pion_E_0']  = gv.gvar(0.14, .006)
priors['pion_zS_0'] = gv.gvar(5e-3, 5e-4)
priors['pion_zS_1'] = gv.gvar(5e-3, 5e-4)
priors['pion_zP_0'] = gv.gvar(0.125,  0.015)
priors['pion_zP_1'] = gv.gvar(0.125,  0.015)

# for i in range(5):
#     for j in range(5):
#         if i+j >= 1:
#             if j < i:  
#                 priors['proton_A3_'+str(j)+str(i)] = gv.gvar(0, 1)
#                 priors['proton_V4_'+str(j)+str(i)] = gv.gvar(0, 1)

#             elif j == i:
#                 priors['proton_A3_'+str(j)+str(i)] = gv.gvar(0, 1)
#                 priors['proton_V4_'+str(j)+str(i)] = gv.gvar(1, 0.2)


for corr in corr_lst:#[k for k in corr_lst if 'mres' not in k]:
    if corr in ['proton','pion']:
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
        sp = snk+corr_lst[corr]['srcs'][0]
        state = corr+'_'+sp
        x[state] = dict()
        x[state]['state'] = corr
        if corr in ['gA']:
            x['d']    = corr_lst[corr]['dA_'+snk]
            x['g_nm'] = 'gA_nm'
        elif corr in ['gV']:
            x['d']    = 'dV_'+snk 
            x['g_nm'] = 'gV_nm'
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
