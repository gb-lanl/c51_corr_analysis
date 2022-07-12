import numpy as np
import gvar as gv
import sys
import copy
import tables as h5
import h5py
import os 
import time
import re
# import fitter.corr_functions as cf
# import fitter.fit_twopt

directory = './data/C13/'
N_cnf = len([name for name in os.listdir(directory) if os.path.isfile(name)])

dirs = os.listdir( directory )

cnf_abbr = [files.split(".ama.h5",0) for files in dirs]

# data_file_list = os.path.realpath(dirs)
data_file_list = list()
for dirpath,_,filenames in os.walk(directory):
    for f in filenames:
        data_file_list.append(os.path.abspath(os.path.join(dirpath, f)))
        #print(data_file_list)
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
   

# print ('checking hdf5 file')
# with h5py.File(data_file_list[0], 'r') as h5f:
#     h5f.visititems(visitor_func) 
def three_pt_corr(Esrc, Esnk, Vins, Asrc, Asnk, tins, tsnk):
    if isinstance(tins, (int,float,np.int64,np.float32)):
        if isinstance(tsnk, (int,float,np.int64,np.float32)):
            expEmt = np.linalg.matrix_power( np.diag( gv.exp( -Esrc)), tins)
            expEnt = np.linalg.matrix_power( np.diag( gv.exp( -Esnk)), (tsnk-tins))
            Asnk.dot( expEnt).dot( Vins).dot( expEmt).dot( Asrc.T)
            ## loop over tsnk second
        return np.array([ three_pt_corr(
            Esrc, Esnk, Vins, Asrc, Asnk, tins, ti) for ti in tsnk ])
    ## loop over tins first
    return np.array([ three_pt_corr(
        Esrc, Esnk, Vins, Asrc, Asnk, ti, tsnk) for ti in tins ])
filename = '/home/gbradley/c51_corr_analysis/tests/visit_h5_20220709-17.txt'
threept_8_paths = []
threept_8_paths.append('/3pt_tsep8/NUCL_D_MIXED_NONREL_l0_g0/src5.0_snk5.0/qz+0_qy+0_qx+0/C13.b_5682/AMA')
with h5py.File(data_file_list[0],'r') as h5f:
    for i in range(len(threept_8_paths)):
        n= np.array(h5f[threept_8_paths[i]][:][:]) 
        # print(n[0][1])
        for j in range(len(n)):
            #source @ t=0, tins @ t=t, tsnk @ t=Tau
            # for k in range(j+1):
            Esrc_re = n[j][0] #real part of 3pt corr fcn
            Esnk_im = n[j][1] #imag part of 3pt corr fcn 

def three_pt_fitfcn(pt3_nstates,pt3_tsep_A3, pt3_tsep_V4, pt3_tau_A3, pt3_tau_V4, p,include_ga=True,include_gv=True):
    '''
    Construct the basic three point correlation function to then be summed over n_states in three_pt_summation 
    '''
    E = {} #initialize state energies to be filled 
    for i in range(pt3_nstates):      
            E['E'+str(i)] = p['E0']

    for i in range(1, pt3_nstates): #define Ei    
        for j in range(1, i+1):
                E['E'+str(i)] += p['dE'+str(j)]
    
    result = {}

    if include_ga:
        result['pt3_A3'] = p['A3_00']*p['z0']*p['z0']*np.exp(-E['E0']*pt3_tsep_A3) 

    if include_gv:
        result['pt3_V4'] = p['V4_00']*p['z0']*p['z0']*np.exp(-E['E0']*pt3_tsep_V4) 

    for i in range(pt3_nstates):    
        for j in range(pt3_nstates):
            if i+j >= 1:
                if j == i:
                    if include_ga :
                        result['pt3_A3'] += p['A3_'+str(j)+str(i)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E['E'+str(j)]*pt3_tsep_A3)

                    if include_gv :
                        result['pt3_V4'] += p['V4_'+str(j)+str(i)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E['E'+str(j)]*pt3_tsep_V4)

                else:
                    mi = np.minimum(j, i)
                    ma = np.maximum(j, i)
                    if include_ga:
                        result['pt3_A3'] += p['A3_'+str(mi)+str(ma)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E['E'+str(j)]*pt3_tsep_A3)*np.exp((E['E'+str(j)]-E['E'+str(i)])*pt3_tau_A3)

                    if include_gv:
                        result['pt3_V4'] += p['V4_'+str(mi)+str(ma)]*p['z'+str(j)]*p['z'+str(i)]*np.exp(-E['E'+str(j)]*pt3_tsep_V4)*np.exp((E['E'+str(j)]-E['E'+str(i)])*pt3_tau_V4)

    return result

    def summation(A3_t, V4_t, p,sum_nstates,sum_tau_cut):
        '''
        sum over n_staes 
        '''
        E_list = {}
        for i in range(sum_nstates): #initialize       
            E_list['E'+str(i)] = p['E0']

        for i in range(1, sum_nstates): #define Ei      
            for j in range(1, i+1):
                    E_list['E'+str(i)] += p['dE'+str(j)]

        if sum_nstates == 1:
            cut = sum_tau_cut

            result = {}

            if include_ga :
                result['sum_A3'] = p['z0'] * p['A3_00'] * p['z0'] * np.exp(-E_list['E0'] * A3_t) * (A3_t - 2*cut + 1)

            if include_gv :
                result['sum_V4'] = p['z0'] * p['V4_00'] * p['z0'] * np.exp(-E_list['E0'] * V4_t) * (V4_t - 2*cut + 1)

            return result

        E_list['E_sum'] = E_list['E'+ str(sum_nstates - 2)] + p['E_sum'] 

        E = np.zeros([sum_nstates-1], dtype=gv.GVar)
        z = np.zeros([sum_nstates-1], dtype=gv.GVar)
        for i in range(sum_nstates-1):
            E[i] = E_list['E'+str(i)]
            z[i] = p['z'+str(i)]

        E_sum = E_list['E_sum']
        z_sum = p['z_sum']

        D = np.zeros([sum_nstates-1, sum_nstates-1], dtype=gv.GVar)
        for i in range(sum_nstates-1):
            for j in range(sum_nstates-1):
                D[i][j] = E[i] - E[j]
        
        A3 = np.zeros([sum_nstates-1, sum_nstates-1], dtype=gv.GVar)
        V4 = np.zeros([sum_nstates-1, sum_nstates-1], dtype=gv.GVar)
        sumA3 = np.zeros([sum_nstates], dtype=gv.GVar)
        sumV4 = np.zeros([sum_nstates], dtype=gv.GVar)
        for i in range(sum_nstates-1):
            for j in range(sum_nstates-1):
                mi = np.minimum(j, i)
                ma = np.maximum(j, i)
                if include_ga :
                    A3[i][j] = p['A3_'+str(mi)+str(ma)]
                if include_gv :
                    V4[i][j] = p['V4_'+str(mi)+str(ma)]

        for i in range(sum_nstates):
            if include_ga :
                sumA3[i] = p['sum_A3_'+str(i)]
            if include_gv :
                sumV4[i] = p['sum_V4_'+str(i)]

        cut = sum_tau_cut

        result = {}

        if include_ga :
            result['sum_A3'] = z[0] * A3[0][0] * z[0] * np.exp(-E[0] * A3_t) * (A3_t - 2*cut + 1)
        if include_gv :
            result['sum_V4'] = z[0] * V4[0][0] * z[0] * np.exp(-E[0] * V4_t) * (V4_t - 2*cut + 1)

        for i in range(sum_nstates-1):
            for j in range(sum_nstates-1):
                if i+j >= 1:
                    if j == i: 
                        if include_ga :
                            result['sum_A3'] += z[j] * A3[j][i] * z[i] * np.exp(-E[j] * A3_t) * (A3_t - 2*cut + 1)
                        if include_gv :
                            result['sum_V4'] += z[j] * V4[j][i] * z[i] * np.exp(-E[j] * V4_t) * (V4_t - 2*cut + 1) 

                    else:
                        if include_ga :
                            result['sum_A3'] += z[j] * A3[j][i] * z[i] * np.exp(-E[j] * A3_t) * ((np.exp(cut * D[j][i]) * (1 - np.exp((A3_t - 2*cut + 1) * D[j][i]) )) / (1 - np.exp(D[j][i]) ))
                        if include_gv :
                            result['sum_V4'] += z[j] * V4[j][i] * z[i] * np.exp(-E[j] * V4_t) * ((np.exp(cut * D[j][i]) * (1 - np.exp((V4_t - 2*cut + 1) * D[j][i]) )) / (1 - np.exp(D[j][i]) ))
                        

        for i in range(sum_nstates-1):
            if include_ga:
                result['sum_A3'] += z_sum * sumA3[i] * z[i] * np.exp(-E_sum * A3_t) * (((np.exp(cut * (E_sum - E[i])) ) * (1 - np.exp((A3_t - 2*cut + 1) * (E_sum - E[i])) )) / (1 - np.exp(E_sum - E[i]) ))
                result['sum_A3'] += z[i] * sumA3[i] * z_sum * np.exp(-E[i] * A3_t) * (((np.exp(cut * (E[i] - E_sum)) ) * (1 - np.exp((A3_t - 2*cut + 1) * (E[i] - E_sum)) )) / (1 - np.exp(E[i] - E_sum) ))
            if include_gv:
                result['sum_V4'] += z_sum * sumV4[i] * z[i] * np.exp(-E_sum * V4_t) * (((np.exp(cut * (E_sum - E[i])) ) * (1 - np.exp((V4_t - 2*cut + 1) * (E_sum - E[i])) )) / (1 - np.exp(E_sum - E[i]) ))
                result['sum_V4'] += z[i] * sumV4[i] * z_sum * np.exp(-E[i] * V4_t) * (((np.exp(cut * (E[i] - E_sum)) ) * (1 - np.exp((V4_t - 2*cut + 1) * (E[i] - E_sum)) )) / (1 - np.exp(E[i] - E_sum) ))
            
        if include_ga :
           result['sum_A3'] += z_sum * sumA3[sum_nstates-1] * z_sum * np.exp(-E_sum * A3_t) * (A3_t - 2*cut + 1)
        if include_gv :
            result['sum_V4'] += z_sum * sumV4[self.sum_nstates-1] * z_sum * np.exp(-E_sum * V4_t) * (V4_t - 2*cut + 1)

        return result 
    
    def ratio_fcn(x,p):
        '''
        Import 2pt correlation fcn in order to construct the ratio 
        <N | A_mu/V_mu | N>|_tTau
        _________________________
                C^2pt(Tau)
        Goal is to extract charges A3, V4 (axial, axial vector charge) from simultaneous fits to c^3pt and C^2pt   
        '''

        result = {}

        if include_2pt:
            pt2_t = x['pt2']

            # result['pt2'] = pt2_fit_function(pt2_t, p)['pt2']
            result['2pt'] = fit_twopt.fit 

        if include_3pt:
            pt3_tsep_A3 = []
            pt3_tau_A3 = []
            pt3_tsep_V4 = []
            pt3_tau_V4 = []

            if include_ga:
                pt3_tsep_A3 = x['pt3_A3'][0]
                pt3_tau_A3 = x['pt3_A3'][1]

            if include_gv:
                pt3_tsep_V4 = x['pt3_V4'][0]
                pt3_tau_V4 = x['pt3_V4'][1]
                
            if include_ga:
                result['pt3_A3'] = pt3_fit_function(pt3_tsep_A3, pt3_tsep_V4, pt3_tau_A3, pt3_tau_V4, p)['pt3_A3']
                
            if include_gv:
                result['pt3_V4'] = pt3_fit_function(pt3_tsep_A3, pt3_tsep_V4, pt3_tau_A3, pt3_tau_V4, p)['pt3_V4']





    E = np.array([np.sum([np.exp(log_dE[j]) for j in range(n+1)]) for n in range(n_states)]) + E0

        # output = 0
        # for n in range(n_states):
        #     for m in range(n_states):
        #         if n == m:
        #             #if m > n: g_nm[n, m] = g_nm[m, n]
        #             output += ((t-1)*wf[n]*g_nm[n, m] + d[n]) * np.exp(-E[n] * t)
        #         else:
        #             pass
        #             E_n = E[n]
        #             E_m = E[m]
        #             dE_nm = E_n - E_m
        #             dE_mn = -dE_nm

        #             output += (wf[n]*g_nm[n, m]) * ((np.exp(dE_nm/2 - E_n*t) - np.exp(dE_mn/2 - E_m*t)) /
        #                                             (np.exp(dE_mn/2) - np.exp(dE_nm/2)))
        # return output

            # # print(Esrc,Esnk)
            # Vins = np.array([[.13, .23], [.17, .29], [.19, .31]])
            # Asrc = np.array([0.1, 0.2, 0.1, 0.2,0.1, 0.2,0.1, 0.2,0.1])
            # Asnk = np.array([0.3, 0.5, 0.7,0.3, 0.5, 0.7,0.3, 0.5, 0.7])
            # tins = np.array(range(3))
            # tsnk = np.array(range(5)) + 5

            # # print(Esrc,Esnk)

            # Ct_test1 = three_pt_corr(Esrc, Esnk, Vins, Asrc, Asnk, tins, tsnk)
            # print(Ct_test1)


        # self.n_states = n_states

        # # keys (strings) used to find the wf_overlap and energy in a parameter dictionary
        # self.param_keys = param_keys




        # wf = p[self.param_keys['wf']]
        # E0 = p[self.param_keys['E0']]
        # log_dE = np.append(-np.inf, p[self.param_keys['log(dE)']])
        # d = p[self.param_keys['d']]
        # g_nm = p[self.param_keys['g_nm']]

        # E = np.array([np.sum([np.exp(log_dE[j]) for j in range(n+1)]) for n in range(self.n_states)]) + E0

        # output = 0
        # for n in range(self.n_states):
        #     for m in range(self.n_states):
        #         if n == m:
        #             #if m > n: g_nm[n, m] = g_nm[m, n]
        #             output += ((t-1)*wf[n]*g_nm[n, m] + d[n]) * np.exp(-E[n] * t)
        #         else:
        #             pass
        #             E_n = E[n]
        #             E_m = E[m]
        #             dE_nm = E_n - E_m
        #             dE_mn = -dE_nm

        #             output += (wf[n]*g_nm[n, m]) * ((np.exp(dE_nm/2 - E_n*t) - np.exp(dE_mn/2 - E_m*t)) /
        #                                             (np.exp(dE_mn/2) - np.exp(dE_nm/2)))
        # return output

            # return correlation fcn at times t with shape (tins, tsnk)
            # Esnk = n[j][j+1]
            # Esnk = n[1][1]
#filename = max_file
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
#                     elif '_g3' in key:UCL_U_MIXED_NONREL_l0_g9
#                         if '_U' in key:
#                             threept_8_paths['U']['tensor'].append(key)
#                         elif '_D' in key:
#                             threept_8_paths['D']['tensor'].append(key)
#                     elif '_g8' in key:
#                         if '_U' in key:
#                             threept_8_paths['U']['vector'].append(key)
#                         elif '_D' in key:
#                             threept_8_paths['D']['vector'].append(key)
                       




#                     threept_8_paths.append(key)
#                     # print(threept_8_paths)
#             #     elif '3pt_tsep10' in key:
#             #         threept_10_paths.append(key)
#             #     elif '3pt_tsep12' in key: 
#             #         threept_12_paths.append(key)
#             #     elif '3pt_tsep14' in key:
#             #         threept_14_paths.append(key)

#             # # array.append(values)
#             with h5py.File(data_file_list[0],'r') as h5f:
#                 for i in range(len(threept_10_paths)):
#                     n= np.array(h5f[threept_10_paths[i]][:]) 
#                     print(n)

#             print(threept_10_paths)



# class load_3pt_data():

#     def __init__(self):
#         self.data_dir = directory # LANL or CalLat data directory; import as module
#         # self.data_path = data_path # path to .h5 file or list of files 
#         # self.twopt_range = twopt_range # 2pt data range 
#         # self.threept_range = threept_range # 3 pt data range 
#         self.incl_ga = True
#         self.incl_gv = True

    


#     # def store_keys(self):
#     #     '''
#     #     locate and store keys for 3 pt correlators:
#     #     eg. positive/negative parity with 
#     #     combinations of DD_dn_dn UU_up_up etc 
#     #     '''
#     #     if self.data_dir is 'lanl':
#     #         t_sep = fp['tsep']
#     #         for i in range(1,15):
#     #             NUCL_U_MIXED_NONREL_l0_g{i}
#     #             if key in ['ext_current','ext_current_SP','local_axial_SP']:
        
                

#     #         if self.incl_ga:
#     #             val['pt3_A3'] = p['A3_00']*p['z0']*p['z0']*np.exp(-E_list['E0']*pt3_tsep_A3) 
#     #         else:
#     #             warnings.warn('Are you sure you do not want to include gA in 3pt analysis?')

#     #         if self.incl_gv:
#     #             val['pt3_V4'] = p['V4_00']*p['z0']*p['z0']*np.exp(-E_list['E0']*pt3_tsep_V4) 
#     #         else:
#     #             warnings.warn('Are you sure you do not want to include gV in 3pt analysis?')



#     #     three_point_ml0p00951 = myfile['gf1p0_w3p5_n45_M51p1_L56_a1p5']['formfac']['ml0p00951']

#     #     proton_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_DD_dn_dn')
#     #     proton_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_DD_up_up')
#     #     proton_np_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_DD_dn_dn')
#     #     proton_np_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_np_DD_up_up')

#     #     proton_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_UU_dn_dn')
#     #     proton_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_UU_up_up')
#     #     proton_np_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_UU_dn_dn')
#     #     proton_np_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_np_UU_up_up')

#     def read_data_avg(self):
#         '''
#         Define file naming scheme which uses chroma conventions for proper datset traversal and retrieval.   
#         '''
#         files = os.listdir(self.data_dir)
#         # file = h5.File(self.data_path,'r')
#         correlator_type = ['2pt','3pt']
#         states_2pt      = ['proton','proton_sp','pion','pion_sp']
#         local_curr      = ['ext_current','ext_current_SP','local_axial_SP']
#         cnf_abbr        = [i.split(".ama.h5",0) for i in files]
#         states_3pt      = ['NUCL_D','NUCL_U']
#         cfgs_src_snk    = ['src5.0_snk5.0'] 
#         tsep            = ['tsep8', 'tsep10', 'tsep12','tsep14'] #tsep values for 3pt data, each is a hdf5 group
#         spin_proj_type  = ['MIXED']
#         rel             = 'NONREL'

#         # chroma dictionary: gamma matrices in QDP++ in the DeGrand-Rossi Basis
#         # gamma(N) = g_0^d g_1^c g_2^b g_3^a ; N = a b c d (binary) 
#         g    = {
#             'scalar' : 'g0', #gamma(1)
#             'axial'  : 'g1', #gamma(2)
#             'tensor' : 'g3', #gamma(0)gamma(1)
#             'vector' : 'g8'  #gamma(3)
#                 }

#         # separation of quarks of the bilinear operator: 
#         # 0 == local operator 
#         # 1 == quark bilinear 
#         l = ['l0']#, 'l1'] 

#         spin = {
#             'spin_up': 'U',
#             'spin_dn': 'D'   }
#         #list of combinations of momenta to be summed over (x,y,z) 
#         momentum = [] 
#         data =  {}
#         with h5py.File(data_file_list[0],'r') as cnf_data:
#             for g1 in cnf_data.keys():
#                 data[g1] = {}
#                 print(g1)
#                 for g2 in cnf_data[g1].keys():
#                     data[g1][g2] = {}
#                     for g3 in cnf_data[g1][g2].keys():
#                         data[g1][g2][g3] = {}
#                         # print(g3)
#                         for g4 in cnf_data[g1][g2][g3].keys():
#                             momentum.append(g4)
#                             # print(momentum)

#         #momentum.append(str.split(proton_))
#         first_data = True #if gathering first "raw" data prior to preparation
#         # to handle multiple files in list in input file 
#         for i in range(1,len(dirs)):
#             with h5.File(data_file_list[0],'r') as h5f:
#                 for cfgs in cfgs_src_snk:
#                     for abbs in cnf_abbr:
#                         two_point_pion  = h5f.File(['2pt']['pion'][cfgs]['pion'][abbs]['AMA']) # different data file need different path here
#                         two_point_pion_  = np.array([two_point_pion[:]])

#                 for corr_type in '2pt':
#                     for particle in states_2pt:
#                         for cfgs in cfgs_src_snk:
#                             for mom in momentum:
#                                 for abbs in cnf_abbr:
                                    
#             #         spin_data = dict()
#             #         have_spin = False
#             #         for s in spin:
#             #             spec = np.array([],dtype=dtype)
#             #             if srcs == None:

        
                                    
#                                     for state in states_3pt:
#                                         for spin_type in spin_proj_type:
#                                             for g in g:

#                                                 three_pt_   = h5f[correlator_type+'_'+tsep][state+'_'+spin_type+'_'+rel+'_'+l+'_'+g][cfgs_src_snk][momentum][cnf_abbr]['AMA']
#         # two_point_proton = np.array([two_point_['proton']['px0_py0_pz0']['spin_dn'][:],
#         # two_point_ml0p00951['proton']['px0_py0_pz0']['spin_up'][:]])
#         # two_point_proton_np = np.array([two_point_ml0p00951['proton_np']['px0_py0_pz0']['spin_dn'][:],
#         # two_point_ml0p00951['proton_np']['px0_py0_pz0']['spin_up'][:]])
        
#         # two_point_data = (two_point_proton[0] + two_point_proton[1] + two_point_proton_sp[0] + two_point_proton_sp[1])/4
#         # two_point_data = np.squeeze(two_point_data)
#         # twopt_proton = 
#         # twopt_proton_sp = 

#         # twopt_pion = 
#         # twopt_pion_sp = 

#         C
# def main():
#     load_3pt_data().read_data_avg()


