import gvar as gv
import numpy as np
import fitter.fastfit as ffit

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
        'tsep' : 0,
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
        'type'     :'exp_3pt',
        'ztype'    :'z_snk z_src',
        'z_ylim'   :[0.055,0.26],
        # fit params
        'tsep'     : 8,
        'n_state'  :3,
        'tsep'     :12,
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
        'tsep'     : 12,
        'T'        : 96,  
        'n_state'  :3,
        't_snk'        :96,
        't_range'  :np.arange(5,48),
        't_sweep'  :range(2,28),
        'n_sweep'  :range(1,6),
        'eff_ylim' :[0.133,0.1349]
    },

}

# we construct priors as Lepage does, plays nicely with corrfitter

prior = gv.BufferDict()

N=1
gA = gv.gvar('1.2(2)')
prior['log(a)'] = gv.log(gv.gvar(N * ['0.3(3)']))
prior['log(dE)'] =  gv.log(gv.gvar(N * ['0.5(5)']))
prior['log(dE)'][0] = gv.log(gA)

# Ds -- oscillating part
prior['log(ao)'] = gv.log(gv.gvar(N * ['0.1(1)']))
prior['log(dEo)'] = gv.log(gv.gvar(N * ['0.5(5)']))
prior['log(dEo)'][0] = gv.log(gA + gv.gvar('0.3(3)'))

# V
nV = int((N * (N + 1)) / 2)
prior['Vnn'] = gv.gvar(nV * ['0.0(5)'])
prior['Voo'] = gv.gvar(nV * ['0.0(5)'])
prior['Vno'] = gv.gvar(N * [N * ['0.0(5)']])

    # for corr in corr_lst:#[k for k in corr_lst if 'mres' not in k]:
    #     for snk in corr_lst[corr]['snks']:
    #         sp = snk+corr_lst[corr]['srcs'][0]
    #         state = corr+'_'+sp
    #     if state in ['pion_PS','pion_SS','proton_SS','proton_PS']:
    #         amps = ['a', 'b', 'ao', 'bo']
    #         tag = state
    #         n_decay = 1
    #         n_oscillate = 0
    #         # prior = {}
    #         # V
    #         nV = int((N * (N + 1)) / 2)
    #         prior['Vnn'] = gv.gvar(nV * ['0.0(5)'])
    #         prior['Voo'] = gv.gvar(nV * ['0.0(5)'])
    #         prior['Vno'] = gv.gvar(N * [N * ['0.0(5)']])
    #         # ffit_ = ffit.FastFit(data)
    #         # Decaying energies and amplitudes
    #         n = range(n_decay)
    #         prior['dE'] = [gv.gvar('1.0(1.0)')] +\
    #             [gv.gvar('0.6(0.6)') for _ in range(1, n_decay)]
    #         if 'a' in amps:
    #             prior['a'] = [gv.gvar('0.1(1.0)') for _ in n]
    #         if 'b' in amps:
    #             prior['b'] = [gv.gvar('0.1(1.0)') for _ in n]

    #         # Oscillating eneriges and amplitudes
    #         if n_oscillate > 0:
    #             no = range(0, n_oscillate)
    #             prior['dEo'] = [gv.gvar('1.65(50)')] +\
    #                             [gv.gvar('0.6(0.6)') for _ in range(1, n_oscillate)]
    #             if 'ao' in amps:
    #                 prior['ao'] = [gv.gvar('0.1(1.0)') for _ in no]
    #             if 'bo' in amps:
    #                 prior['bo'] = [gv.gvar('0.1(1.0)') for _ in no]

    #         # Extract guesses for the ground-state energy and amplitude
    #         # if ffit is not None:
    #         #     dE_guess = gv.mean(ffit.E)
    #         #     amp_guess = gv.mean(ffit.ampl)
    #         #     prior['dE'][0] = gv.gvar(dE_guess, 0.5 * dE_guess)
    #         #     if 'a' in amps:
    #         #         prior['a'][0] = gv.gvar(amp_guess, 2.0 * amp_guess)
    #         #     elif 'b' in amps:
    #         #         prior['b'][0] = gv.gvar(amp_guess, 2.0 * amp_guess)
    #         #     else:
    #         #         msg = "Error: Unrecognized amplitude structure?"
    #         #         raise ValueError(msg)

    #         # Convert to arrays
    # keys = list(prior.keys())
    # tag = 'proton_PS'
    # if tag is None:
    #     # Just convert to arrays
    #     for key in keys:
    #         prior[key] = np.asarray(prior[key])
    # else:
    #     # Prepend keys with 'tag:' and then convert
    #     for key in keys:
    #         new_key = "{0}:{1}".format(tag, key)
    #         prior[new_key] = np.asarray(prior.pop(key))

#proton
# priors['proton_E_0']  = gv.gvar(0.5, .06)
# priors['proton_zS_0'] = gv.gvar(2.0e-5, 1.e-5)
# priors['proton_zP_0'] = gv.gvar(2.5e-3, 1.e-3)

# #pion
# priors['pion_E_0']  = gv.gvar(0.14, .006)
# priors['pion_zS_0'] = gv.gvar(5e-3, 5e-4)
# priors['pion_zS_1'] = gv.gvar(5e-3, 5e-4)
# priors['pion_zP_0'] = gv.gvar(0.125,  0.015)
# priors['pion_zP_1'] = gv.gvar(0.125,  0.015)

# #gA, gV priors
# #physical values
# priors['gA_00'] = gv.gvar(1.2, 0.2)
# priors['gV_00'] = gv.gvar(1.0, 0.2)

# # set gA_nm, gV_nm priors 
# for i in range(5):
#     for j in range(5):
#         if i+j >= 1:
#             if j < i:  
#                 priors['gA_'+str(j)+str(i)] = gv.gvar(0, 1)
#                 priors['gV_'+str(j)+str(i)] = gv.gvar(0, 1)
#             elif j == i:
#                 priors['gA_'+str(j)+str(i)] = gv.gvar(0, 1)
#                 priors['gV_'+str(j)+str(i)] = gv.gvar(1, 0.2)



# for corr in corr_lst:#[k for k in corr_lst if 'mres' not in k]:
#     if corr in ['proton','pion']:
#         for n in range(1,10):
#             # use 2 mpi splitting for each dE

#             # E_n = E_0 + dE_10 + dE_21 +...
#             # use log prior to force ordering of dE_n
#             priors['log(%s_dE_%d)' %(corr,n)] = np.array(gv.gvar(np.log(2*priors['pion_E_0'].mean), 0.7))

#             # for z_P, no suppression with n, but for S, smaller overlaps
#             priors['%s_zP_%d' %(corr,n)] = gv.gvar(priors['%s_zP_0' %(corr)].mean, 2*priors['%s_zP_0' %(corr)].sdev)
#             zS_0 = priors['%s_zS_0' %(corr)]
#             if n <= 2:
#                 priors['%s_zS_%d' %(corr,n)] = gv.gvar(zS_0.mean, 2*zS_0.sdev)
#             else:
#                 priors['%s_zS_%d' %(corr,n)] = gv.gvar(zS_0.mean/2, zS_0.sdev)
x = dict()
for corr in corr_lst:
    for snk in corr_lst[corr]['snks']:
        sp = snk+corr_lst[corr]['srcs'][0]
        state = corr+'_'+sp
        x[state] = dict()
        x[state]['state'] = corr
        x[state]['z'] = corr+'_z'+snk+'_0'
        # if corr == 'gA':
        #     x[state]['tsep']    = corr_lst[corr]['tsep']
        #     x[state]['g_nm'] = 'gA_nm'        
        # elif x[state]['state'] == 'gV':
        #     x[state]['tsep']    = corr_lst[corr]['tsep']
        #     x[state]['d']    = 'dV_'+sp 
        #     x[state]['g_nm'] = 'gV_nm'
        for k in ['type', 'tsep', 'n_state', 't_range', 'eff_ylim', 'ztype']:
            if k in corr_lst[corr]:
                x[state][k] = corr_lst[corr][k]
        if 't0' in corr_lst[corr]:
            x[state]['t0'] = corr_lst[corr]['t0']
        else:
            x[state]['t0'] = 0
        if 't_snk' in corr_lst[corr]:
            x[state]['t_snk'] = corr_lst[corr]['t_snk']
        else:
            x[state]['t_snk'] = 80
        
        if 'mres' not in corr:
            x[state]['color'] = corr_lst[corr]['colors'][sp]
            x[state]['snk']   = snk
            x[state]['src']   = corr_lst[corr]['srcs'][0]
        else:
            x[state]['color'] = corr_lst[corr]['colors']
