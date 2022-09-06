import gvar as gv
import numpy as np
import os
# import fitter.fastfit as ffit

file_params = {}
file_params['data_dir'] = '/home/gbradley/c51_corr_analysis/tests/data/E7/' #all configurations
file_params['data_file'] = './data/E7/E7-0_1000.ama.h5'

# def ens_base():
#     ens,stream = os.getcwd().split('/')[-1].split('_')
#     return ens,stream 
# print(ens_base())

def parse_cfg_arg(cfg_arg,params):
    allowed_cfgs = range(params['cfg_i'],params['cfg_f']+1,params['cfg_d'])
    if not cfg_arg:
        ci = params['cfg_i']
        cf = params['cfg_f']
        dc = params['cfg_d']
    else:
        if len(cfg_arg) == 3:
            cfgs = range(cfg_arg[0], cfg_arg[1]+1, cfg_arg[2])
        elif len(cfg_arg) == 1:
            cfgs = range(cfg_arg[0], cfg_arg[0]+1, 1)
        if not all([cfg in allowed_cfgs for cfg in cfgs]):
            print('you selected configs not allowed for %s:' %params['ENS_S'])
            print('       allowed cfgs = range(%d, %d, %d)' %(params['cfg_i'], params['cfg_f'], params['cfg_d']))
            sys.exit('  your choice: cfgs = range(%d, %d, %d)' %(cfgs[0],cfgs[-1],cfgs[1]-cfgs[0]))
        elif len(cfg_arg) == 1:
            ci = int(cfg_arg[0])
            cf = int(cfg_arg[0])
            dc = 1
        elif len(cfg_arg) == 3:
            ci = int(cfg_arg[0])
            cf = int(cfg_arg[1])
            dc = int(cfg_arg[2])
        else:
            print('unrecognized use of cfg arg')
            print('cfg_i [cfg_f cfg_d]')
            sys.exit()
    return range(ci,cf+1,dc)



def ensemble(params):
    cfg     = params['ENS_BASE']+'-'+params['STREAM']+'.ama'

params = dict()
params['cfg_i'] = 600  # initial config number
params['cfg_f'] = 2232 # final config number
params['cfg_d'] = 4    # config step value 

params['seed'] = dict()
params['seed']['0'] = '0'
params['seed']['a'] = 'a'
params['seed']['b'] = 'b'
params['seed']['c'] = 'c'

params['ENS_ABBR'] =  'a071m170'
params['ENS']      =  'E7'
params['t_seps']  = [13,15,17,19,21]
params['flavs']   = ['U','D']
params['spins']   = ['up_up','dn_dn']
params['snk_mom'] = ['0 0 0']
params['SS_PS']   = 'SS'
params['particles'] = ['proton','proton_SP','pion','pion_SP']
params['curr_4d'] = ['A3','V4','A1','A2','A4','V1','V2','V3','P','S']
params['curr_0p'] = ['A3','V4','A1','A2','A4','V1','V2','V3','S','T34','T12','T13','T14','T23','T24','CHROMO_MAG']

names = dict()
names['flow']             = 'cfg_flow_%(ENS_LONG)s%(STREAM)s_%(CFG)s_wflow%(FLOW_TIME)s'
names['src']              = 'src_%(ENS_S)s_%(CFG)s_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s_%(SRC)s'
names['prop']             = 'prop_%(ENS_S)s_%(CFG)s_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s'
names['prop']            += '_M5%(M5)s_L5%(L5)s_a%(alpha5)s_mq%(MQ)s_%(SRC)s'
''' the xml generation may incluce multiple quark masses, so no mq info '''
names['prop_xml']         = 'prop_%(ENS_S)s_%(CFG)s_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s'
names['prop_xml']        += '_M5%(M5)s_L5%(L5)s_a%(alpha5)s_%(SRC)s'
names['spec']             = 'spec_%(ENS_S)s_%(CFG)s_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s'
names['spec']            += '_M5%(M5)s_L5%(L5)s_a%(alpha5)s_mq%(MV_L)s_%(SRC)s'
names['spec_4D']            = names['spec'].replace('spec_','spec_4D_')
names['spec_4D_tslice']     = names['spec'].replace('spec_','spec_4D_tslice_')
names['spec_4D_tslice_avg'] = names['spec'].replace('spec_','spec_4D_tslice_avg_')
names['hyperspec']        = 'hyperspec_%(ENS_S)s_%(CFG)s_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s'
names['hyperspec']       += '_M5%(M5)s_L5%(L5)s_a%(alpha5)s_ml%(MV_L)s_ms%(MV_S)s_%(SRC)s'
names['h_spec']           = names['hyperspec']
names['hisq_spec']        = 'hisq_spec_%(ENS_S)s_ml%(ML)s_ms%(MS)s_%(CFG)s_%(SRC)s'
names['seqsrc']           = 'seqsrc_%(ENS_S)s_%(CFG)s_%(PARTICLE)s_%(FLAV_SPIN)s'
names['seqsrc']          += '_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s_M5%(M5)s_L5%(L5)s_a%(alpha5)s_mq%(MQ)s'
names['seqsrc']          += '_%(MOM)s_%(SRC)s_%(SS_PS)s'
names['coherent_seqsrc']  = 'seqsrc_%(ENS_S)s_%(CFG)s_%(PARTICLE)s_%(FLAV_SPIN)s'
names['coherent_seqsrc'] += '_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s_M5%(M5)s_L5%(L5)s_a%(alpha5)s_mq%(MQ)s'
names['coherent_seqsrc'] += '_%(MOM)s_dt%(T_SEP)s_Nsnk%(N_SEQ)s_%(SS_PS)s'
names['seqprop']          = 'seqprop_%(ENS_S)s_%(CFG)s_%(PARTICLE)s_%(FLAV_SPIN)s'
names['seqprop']         += '_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s_M5%(M5)s_L5%(L5)s_a%(alpha5)s_mq%(MQ)s'
names['seqprop']         += '_%(MOM)s_dt%(T_SEP)s_Srcs%(SRC_SET)s_%(SS_PS)s'
names['formfac']          = 'formfac_%(ENS_S)s_%(CFG)s'
names['formfac']         += '_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s_M5%(M5)s_L5%(L5)s_a%(alpha5)s_mq%(MQ)s'
#names['formfac']         += '_%(MOM)s_dt%(T_SEP)s_Nsnk%(N_SEQ)s_%(SRC)s_%(SS_PS)s'
names['formfac']         += '_%(MOM)s_dt%(T_SEP)s_Srcs%(SRC_SET)s_%(SRC)s_%(SS_PS)s'
names['formfac_4D']                = names['formfac'].replace('formfac','formfac_4D')
names['formfac_4D_tslice']         = names['formfac'].replace('formfac','formfac_4D_tslice')
names['formfac_4D_tslice_src_avg'] = names['formfac'].replace('formfac','formfac_4D_tslice_src_avg')

names['mixed_corr']       = 'dwf_hisq_spec_%(ENS_S)s_wflow%(FLOW_TIME)s_M5%(M5)s_L5%(L5)s'
names['mixed_corr']      += '_a%(alpha5)s_cfg_%(CFG)s_src%(SRC)s_%(SMR)s_ml%(MQ_L)s_ms%(MQ_S)s.corr'
names['pipi_scat']        = 'pipi_%(ENS_S)s_%(CFG)s_gf%(FLOW_TIME)s_w%(WF_S)s_n%(WF_N)s'
names['pipi_scat']       += '_M5%(M5)s_L5%(L5)s_a%(alpha5)s_%(MV_LS)s_%(SRC)s'