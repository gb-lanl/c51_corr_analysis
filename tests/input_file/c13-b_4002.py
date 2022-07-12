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
cnf_abbr = [files.split(".ama.h5",0) for files in dirs]
print(cnf_abbr)
# data_file_list = os.path.realpath(dirs)
data_file_list = list()
for dirpath,_,filenames in os.walk(directory):
    for f in filenames:
        data_file_list.append(os.path.abspath(os.path.join(dirpath, f)))
        #print(data_file_list)
# Read group structure
#-------------------------
data = {}
with h5py.File(data_file_list[0],'r') as cnf_data:
    for g1 in cnf_data.keys():
        data[g1] = {}
        for g2 in cnf_data[g1].keys():
            data[g1][g2] = {}
            for g3 in cnf_data[g1][g2].keys():
                data[g1][g2][g3] = {}
                for g4 in cnf_data[g1][g2][g3].keys():
                    data[g1][g2][g3][g4] = {}
                    for g5 in cnf_data[g1][g2][g3][g4].keys():
                        data[g1][g2][g3][g4][g5] = {}
                        for g6 in cnf_data[g1][g2][g3][g4][g5].keys():
                            data[g1][g2][g3][g4][g5][g6] = [None]*N_cnf
                            # print(data)
                        #print(data)


def visitor_func(name, node):
    dsets = list()
    timestr = time.strftime("%Y%m%d-%H")
    with open("visit_h5_"+timestr+".txt", "a") as f:
        if isinstance(node, h5py.Group):
            print(node.name, 'is a Group',file=f)
        elif isinstance(node, h5py.Dataset):
            if (node.dtype == 'object') :
                # f[node].read_direct(data)
                print (node.name, 'is an object Dataset',file=f)
            else:
                
                print(node.name, 'is a Dataset',file=f)
                # dsets.append(node.name)
                # print(dsets)
                # n1 = np.array(f[node.name][:])
                # print(n1)

        else:
            print(node.name, 'is an unknown type',file=f)
   

print ('checking hdf5 file')
# with h5py.File(data_file_list[0], 'r') as h5f:
#     h5f.visititems(visitor_func) 
# for i_cnf in range(1,len(dirs)):
#         # cnf_data_= np.fromfile(data_file_list[i_cnf])
#     with h5py.File(data_file_list[i_cnf],'r') as h5f:
#         h5f.visititems(visitor_func)

# for i_cnf in range(1,len(dirs)):
#         # cnf_data_= np.fromfile(data_file_list[i_cnf])
#     with h5.File(data_file_list[i_cnf],'r') as cnf_data_:

# getting generated .txt file from above for data path inspection
# folder_path = os.getcwd()
# file_type = r'\*txt'
# files = glob.glob(folder_path + file_type)
# max_file = max(files, key=os.path.getctime)
# print(max_file)


# for i_cnf in range(1,len(dirs)):
#         # cnf_data_= np.fromfile(data_file_list[i_cnf])
#     with h5py.File(data_file_list[i_cnf],'r') as h5f:
#         print(threept_10_paths[i_cnf])
# # with h5py.File(data_file_list[0], 'r') as h5f:
# #     print(threept_10_paths[0])
# #     n1 = np.array(h5f[threept_10_paths[0]][:]) 
# #     print(n1)

fit_states = ['proton', 'pion', 'nucleon']
#bs_seed = 'a12m220XL

corr_lst = {
    'proton':{
        'q': 'U',
        'state': 'axial',
        'dsets':['2pt/proton/src5.0_snk5.0/proton/C13.b_4002/AMA'],
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
        'n_state'  :2,
        'T'        :96,
        't_range'  :np.arange(8,48),
        't_sweep'  :range(3,28),
        'n_sweep'  :range(1,6),
        'eff_ylim' :[0.133,0.1349]
    },
}



# filename = 'visit_h5_20220707-1230.txt'
# #filename = max_file
# string_fnd_1 = 'is a Dataset'
# with open(filename, "r") as f:
#     for line in f:
#         if string_fnd_1 in line:
#             # cleaning the bad chars in line
#             line = line.strip()
#             # line = line.strip(" \\n")
#             line = re.sub(r"is a Dataset\s*", "", line)
#             line = re.sub(r"(-[0-9]+\.)", r" \1", line)

#             values = [str(value) for value in line.split()]
#             for key in values:
#                 if '2pt' in key:
#                     if 'proton' in key:
#                         twopt_paths.append(key)
#                     # print(twopt_paths)
#                 elif '3pt_tsep8' in key:
#                     if '_g0' in key:
#                         if '_U' in key:
#                             threept_8_paths['U']['scalar'].append(key)
#                         elif '_D' in key:
#                             threept_8_paths['D']['scalar'].append(key)
#                     elif '_g1'in key:
#                         if '_U' in key:
#                             threept_8_paths['U']['axial'].append(key)
#                         elif '_D' in key:
#                             threept_8_paths['D']['axial'].append(key)
#                     elif '_g3' in key:
#                         if '_U' in key:
#                             threept_8_paths['U']['tensor'].append(key)
#                         elif '_D' in key:
#                             threept_8_paths['D']['tensor'].append(key)
#                     elif '_g8' in key:
#                         if '_U' in key:
#                             threept_8_paths['U']['vector'].append(key)
#                         elif '_D' in key:
#                             threept_8_paths['D']['vector'].append(key)
                       




#             #         threept_8_paths.append(key)
#             #         # print(threept_8_paths)
#             #     elif '3pt_tsep10' in key:
#             #         threept_10_paths.append(key)
#             #     elif '3pt_tsep12' in key: 
#             #         threept_12_paths.append(key)
#             #     elif '3pt_tsep14' in key:
#             #         threept_14_paths.append(key)

#             # # array.append(values)
#             print(threept_8_paths)

# with h5py.File(data_file_list[0],'r') as h5f:
#     for i in range(len(threept_10_paths)):
#         n= np.array(h5f[threept_10_paths[i]][:]) 
#         print(n)

def is_numeric_real( val):
  """
  Check if a value is a real numeric type.
  Also checks against numpy native types.

  :param val: value to check

  :returns: `True` if value is a real numerical type, otherwise `False`

  :rtype: `boolean`
  """
  return isinstance( val, (int, float, np.int64, np.float32))

def correlator_simple_3pt( Esrc_m, Esnk_n, Vins_nm, Asrc_im, Asnk_jn, tins, tsnk):
    '''
    :param Esrc_m:
    A shape `(M,)` array of the state energies propagating between source and insertion

    :type Esrc_m: :class:`numpy.ndarray` [ `float` ]

    :param Esnk_n:
        A shape `(N,)` array of the state energies propagating between insertion and sink

    :type Esnk_n: :class:`numpy.ndarray` [ `float` ]

    :param Vins_nm: A shape `(N,M)` array of overlap parameters to fit at the source

    :type Vins_nm: :class:`numpy.ndarray` [ `float` ]

    :param Asrc_im: A shape `(N,)` array of overlap parameters to fit at the source

    :type Asrc_im: :class:`numpy.ndarray` [ `float` ]

    :param Asnk_jn: A shape `(N,)` array of overlap parameters to fit at the sink

    :type Asnk_jn: :class:`numpy.ndarray` [ `float` ]

    :param tins: A shape `(Tins,)` array of the sink Euclidean times

    :type tins: :class:`numpy.ndarray` [ `int` ]

    :param tsnk: A shape `(Tsnk,)` array of the sink Euclidean times

    :type tsnk: :class:`numpy.ndarray` [ `int` ]

    :returns: The shape `(Tins,Tsnk)` correlation function at times `t`

    :rtype: :class:`numpy.ndarray` [ `float` ]
'''
    if is_numeric_real( tins):
        if is_numeric_real( tsnk):
        ## capture fixed time behavior
            expEmt = np.linalg.matrix_power( np.diag( gv.exp( -Esrc_m)), tins)
            expEnt = np.linalg.matrix_power( np.diag( gv.exp( -Esnk_n)), (tsnk-tins))
            return Asnk_jn.dot( expEnt).dot( Vins).dot( expEmt).dot( Asrc_im.T)
        ## loop over tsnk second
        return np.array([ correlator_simple_3pt(
        Esrc_m, Esnk_n, Vins_nm, Asrc_im, Asnk_jn, tins, ti) for ti in tsnk ])
    ## loop over tins first
    return np.array([ correlator_simple_3pt(
        Esrc_m, Esnk_n, Vins_nm, Asrc_im, Asnk_jn, ti, tsnk) for ti in tins ])

def make_fit_params(fp,states,gv_data):
    x = copy.deepcopy(fp.x)
    y = {k: v[x[k]['t_range']]
        for (k, v) in gv_data.items() if k.split('_')[0] in states}
    for k in y:
        if 'exp_r' in x[k]['type']:
            sp = k.split('_')[-1]
            y[k] = y[k] / gv_data[x[k]['denom'][0]+'_'+sp][x[k]['t_range']]
            y[k] = y[k] / gv_data[x[k]['denom'][1]+'_'+sp][x[k]['t_range']]
    if any(['mres' in k for k in y]):
        mres_lst = [k.split('_')[0] for k in y if 'mres' in k]
        mres_lst = list(set(mres_lst))
        for k in mres_lst:
            y[k] = y[k+'_MP'] / y[k+'_PP']
    
    n_states = dict()
    for state in states:
        for k in x:
            if state in k:
                if state in k and 'mres' not in k:
                    n_states[state] = x[k]['n_state']
    priors = dict()
    for k in fp.priors:
        for state in states:
            if 'mres' not in k:
                k_n = int(k.split('_')[-1].split(')')[0])
                if state == k.split('(')[-1].split('_')[0] and k_n < n_states[state]:
                    priors[k] = gv.gvar(fp.priors[k].mean, fp.priors[k].sdev)
            else:
                mres = k.split('_')[0]
                if mres in states:
                    priors[k] = gv.gvar(fp.priors[k].mean, fp.priors[k].sdev)
    return x,y,n_states,priors

if __name__ == "__main__":
  print( "performing unit tests")
  An0 = np.array([0.1, 0.2])
  Bn0 = np.array([0.3, 0.5, 0.7])
  Ea0 = np.array([.50, .75])
  Eb0 = np.array([.63, .82, .91])
  Vba = np.array([[.13, .23], [.17, .29], [.19, .31]])
  ta = np.array( range( 3))
  tb = np.array( range( 5)) +5
  Ct_test0 = np.array([[
      0.3 *0.13 *0.1 *np.exp( -.63*(t1-t0)) *np.exp( -.50*t0)
    + 0.5 *0.17 *0.1 *np.exp( -.82*(t1-t0)) *np.exp( -.50*t0)
    + 0.7 *0.19 *0.1 *np.exp( -.91*(t1-t0)) *np.exp( -.50*t0)
    + 0.3 *0.23 *0.2 *np.exp( -.63*(t1-t0)) *np.exp( -.75*t0)
    + 0.5 *0.29 *0.2 *np.exp( -.82*(t1-t0)) *np.exp( -.75*t0)
    + 0.7 *0.31 *0.2 *np.exp( -.91*(t1-t0)) *np.exp( -.75*t0)
    for t1 in tb ] for t0 in ta ])
  Ct_test1 = correlator_simple_3pt( Ea0, Eb0, Vba, An0, Bn0, ta, tb)
  if not( np.all( np.isclose( Ct_test0, Ct_test1))):
    raise ValueError( "test 0 failed")
  print( "unit tests passed")



# Read group structure
#-------------------------
# data = {}
# with h5.File(data_file_list[0],'r') as cnf_data:
#     for g1 in cnf_data.keys():
#         data[g1] = {}
#         for g2 in cnf_data[g1].keys():
#             data[g1][g2] = {}
#             for g3 in cnf_data[g1][g2].keys():
#                 data[g1][g2][g3] = {}
#                 for g4 in cnf_data[g1][g2][g3].keys():
#                     data[g1][g2][g3][g4] = {}
#                     for g5 in cnf_data[g1][g2][g3][g4].keys():
#                         data[g1][g2][g3][g4][g5] = {}
#                         for g6 in cnf_data[g1][g2][g3][g4][g5].keys():
#                             data[g1][g2][g3][g4][g5][g6] = [None]*N_cnf
#                             print(data)
#                         #print(data)
#data_file = './data/C13/C13-b_4002.ama.h5'
#data_file = os.path.abspath(dirs[0]) 

#data_file = 'data/C13/C13-b_4002.ama.h5'

fit_states = ['proton', 'proton_SP', 
              'pion',   'pion_SP', 
              'ext_current', 'ext_current_SP', 
              'local_axial_SP']
#bs_seed = 'a12m220XL'

corr_lst = {
    # PION
    'proton':{
        'dsets':['2pt/proton/src5.0_snk5.0/proton/C13.b_4002/AMA'],
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
        'n_state'  :2,
        'T'        :96,
        't_range'  :np.arange(8,48),
        't_sweep'  :range(3,28),
        'n_sweep'  :range(1,6),
        'eff_ylim' :[0.133,0.1349]
    },
}