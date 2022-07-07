import numpy as np
import gvar as gv
import sys
import copy
import tables as h5
import os 


class load_3pt_data():
    def __init__(self):
        self.data_origin = data_origin # LANL or CalLat data
        self.incl_ga = 
        self.incl_gv = 

    


    def store_keys(self):
        '''
        locate and store keys for 3 pt correlators:
        eg. positive/negative parity with 
        combinations of DD_dn_dn UU_up_up etc 
        '''
        if self.data_origin is 'lanl':
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


        
        elif self.data_origin is 'callat':


        three_point_ml0p00951 = myfile['gf1p0_w3p5_n45_M51p1_L56_a1p5']['formfac']['ml0p00951']

        proton_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_DD_dn_dn')
        proton_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_DD_up_up')
        proton_np_DD_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_DD_dn_dn')
        proton_np_DD_up_up = self.find_key(three_point_ml0p00951, 'proton_np_DD_up_up')

        proton_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_UU_dn_dn')
        proton_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_UU_up_up')
        proton_np_UU_dn_dn = self.find_key(three_point_ml0p00951, 'proton_np_UU_dn_dn')
        proton_np_UU_up_up = self.find_key(three_point_ml0p00951, 'proton_np_UU_up_up')



