import numpy as np
import gvar as gv
import sys
import copy
import tables as h5
import os 


class load_3pt_data():
    def __init__(self):
        self.data_dir = data_dir # LANL or CalLat data directory; import as module
        self.data_path = data_path # path to .h5 file or list of files 
        self.twopt_range = twopt_range # 2pt data range 
        self.threept_range = threept_range # 3 pt data range 

        self.incl_ga = 
        self.incl_gv = 

    


    def store_keys(self):
        '''
        locate and store keys for 3 pt correlators:
        eg. positive/negative parity with 
        combinations of DD_dn_dn UU_up_up etc 
        '''
        if self.data_dir is 'lanl':
            t_sep = fp['tsep']
            for i in range(1,15):
                NUCL_U_MIXED_NONREL_l0_g{i}
                if key in ['ext_current','ext_current_SP','local_axial_SP']:
        
                

            if self.incl_ga:
                val['pt3_A3'] = p['A3_00']*p['z0']*p['z0']*np.exp(-E_list['E0']*pt3_tsep_A3) 
            else:
                warnings.warn('Are you sure you do not want to include gA in 3pt analysis?')

            if self.incl_gv:
                val['pt3_V4'] = p['V4_00']*p['z0']*p['z0']*np.exp(-E_list['E0']*pt3_tsep_V4) 
            else:
                warnings.warn('Are you sure you do not want to include gV in 3pt analysis?')



        three_point_ml0p00951 = myfile['gf1p0_w3p5_n45_M51p1_L56_a1p5']['formfac']['ml0p00951']

        proton_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_DD_dn_dn')
        proton_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_DD_up_up')
        proton_np_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_DD_dn_dn')
        proton_np_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_np_DD_up_up')

        proton_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_UU_dn_dn')
        proton_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_UU_up_up')
        proton_np_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_UU_dn_dn')
        proton_np_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_np_UU_up_up')

    def read_data_avg(self):
        file = h5.File(self.data_path,'r')
        correlator_type = ['2pt','3pt']
        spin = ['spin_up', 'spin_dn']
        par = ['proton','proton_sp','pion','pion_sp']
        momentum = {}   #list of combinations of momenta to be summed over (x,y,z) 
        #momentum.append(str.split(proton_))
        first_data = True #if gathering first "raw" data prior to preparation
        cfgs_src = ['src5.0_snk5.0']
        for corr_type in correlator_type:
            for corr in par:
                for mom in momentum:


        
        two_point_ = myfile[correlator_type][par][cfgs_src][momentum] # different data file need different path here

        two_point_proton = np.array([two_point_['proton']['px0_py0_pz0']['spin_dn'][:],
        two_point_ml0p00951['proton']['px0_py0_pz0']['spin_up'][:]])
        two_point_proton_np = np.array([two_point_ml0p00951['proton_np']['px0_py0_pz0']['spin_dn'][:],
        two_point_ml0p00951['proton_np']['px0_py0_pz0']['spin_up'][:]])
        
        two_point_data = (two_point_proton[0] + two_point_proton[1] + two_point_proton_sp[0] + two_point_proton_sp[1])/4
        two_point_data = np.squeeze(two_point_data)
        twopt_proton = 
        twopt_proton_sp = 

        twopt_pion = 
        twopt_pion_sp = 

        three_point_data = myfile[]

