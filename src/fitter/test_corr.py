import os
import h5py
import sys
import gvar as gv
import pytest
import argparse
import importlib
import corrfitter
import collections
import numpy as np
import lsqfit

import fitter.corr_functions as cf 
import fitter.load_data as ld
import fitter.plotting as plot
import fitter.fastfit as prelim #LePage's preliminary corrfitter to generate p0
import fitter.priors as priors 
sys.path.insert(0, '/home/gbradley/nucleon_elastic_FF')
from nucleon_elastic_ff.data.h5io import get_dsets 


SVDCUT = 0.002
Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
def main():
    parser = argparse.ArgumentParser(
        description='Perform analysis of two-point correlation function')
    parser.add_argument('fit_params',    help='input file to specify fit')
    parser.add_argument('--fit',         default=True, action='store_true',
                        help=            'do fit? [%(default)s]')
    parser.add_argument('--prior_override',default=True, action='store_true',
                        help=            'override generated priors with priors from input file? [%(default)s]')
    parser.add_argument('-b', '--block', default=1, type=int,
                        help=            'specify bin/blocking size in terms of saved configs')
    parser.add_argument('--uncorr_corrs', default=False, action='store_true',
                        help=            'uncorrelate different correlation functions? [%(default)s]')
    parser.add_argument('--uncorr_all',  default=False, action='store_true',
                        help=            'uncorrelate all snk,src for each correlation function? [%(default)s]')
    parser.add_argument('--states',      nargs='+',
                        help=            'specify states to fit?')    
    parser.add_argument('--svdcut',      type=float, help='add svdcut to fit')
    parser.add_argument('--svd_test',    default=True, action='store_false',
                        help=            'perform gvar svd_diagnosis? [%(default)s]')
    parser.add_argument('--svd_nbs',     type=int, default=50, help='number of BS samples for estimating SVD cut [%(default)s]')                    

    args = parser.parse_args()
    # if args.save_figs and not os.path.exists('figures'):
    #     os.makedirs('figures')
    # print(args)
    # add path to the input file and load it
    sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
    fp = importlib.import_module(
        args.fit_params.split('/')[-1].split('.py')[0])

    # h5fname = '/home/gbradley/c51_corr_analysis/tests/data/C13/C13-b_4002.ama.h5'
    
    # data = {}
    # corrs_gv = {}
    # with h5py.File(h5fname, 'r') as h5f:
    #     dsets = get_dsets(h5f)
    #     # print(dsets)
    #     for key in dsets.keys():
        
    #         # print(key)
    #         data[key] = h5f[key][:]
    #     ordered_datatags = sorted(data.keys(),key=str)
    #     # print(ordered_datatags)
    #     sizes = [data[tag].shape[0] for tag in ordered_datatags]
    # #     # print(sizes)
    # #         # corrs_gv[key] = gv.dataset.avg_data(data[key])
    # #         # print(corrs_gv.keys())
    # # # print(data['2pt/ext_current/src5.0_snk5.0/ext_axial_A1_A1/C13.b_4002/AMA'])
    # datatag = '2pt/proton/src5.0_snk5.0/proton/C13.b_4002/AMA'
    # datatag_SP = '2pt/proton_SP/src5.0_snk5.0/proton/C13.b_4002/AMA' 
    # corr_ss = data[datatag]
    # corr_ps = data[datatag_SP]
    
    # for j in range(len(corr_ss)):
    #     #source @ t=0, tins @ t=t, tsnk @ t=Tau
    #     # for k in range(j+1):
    #     Esrc_re[i] = np.array(corr_ss[j][0]) #real part of 3pt corr fcn
    #     Esrc_im[i] = (corr_ss[j][1]) #imag part of 3pt corr fcn 
    #     print(Esrc_re)
    # for i in range(len(corr_ps)):
    #     #source @ t=0, tins @ t=t, tsnk @ t=Tau
    #     # for k in range(j+1):
    #     Esnk_re[i] = np.array(corr_ps[i][0]) #real part of 3pt corr fcn
    #     Esnk_im[i] = (corr_ps[i][1]) #imag part of 3pt corr fcn 
        
    # ydata = {}
    # ydata['ss'] = Esrc_re
    # ydata['ps'] = Esnk_re
    # print(ydata)
   
    # # # nt = corr_ss.shape
    # # # print(nt)
    # # # corr_ss = corr_ss[..., np.newaxis]
    # # _, nt = corr_ss.shape
    # # print(nt)
    # # t = np.arange(nt)
    # # # cr = corr_ss[:, :nt // 2 +1]
    # # # print(cr)
    # # corr_ss_ = corr_ss.real
    # # c = corr_ss.imag
    # fit_out = test_NPoint(tag,Esrc_re,prior=priors)
    # # fit_ = test_NPoint_snk(datatag_SP,data,prior=priors)
    # # print(fit_)
    # print(fit_out)
    
    # print(gv.dataset.avg_data(out))
    # # cov = correct_covariance(data, **kwargs)
    # # dset = corrfitter.read_dataset(h5fname,h5group=datatag)
    # # print(dset.size())
    # def time_reverse(corr, reverse=True, phase=1, time_axis=1):
    #     ''' assumes time index is second of array
    #         assumes phase = +- 1
    #     '''
    #     if reverse:
    #         if len(corr.shape) > 1:
    #             cr = phase * np.roll(corr[:, ::-1], 1, axis=time_axis)
    #             cr[:, 0] = phase * cr[:, 0]
    #         else:
    #             cr = phase * np.roll(corr[::-1], 1)
    #             cr[0] = phase * cr[0]
    #     else:
    #         cr = phase * corr
    #     return cr
    # def fold(data):
    #     data = 0.5*(data + time_reverse(data))
    #     return data

    # # # fold(corr_ss)
    # # for key, value in dset.items():
    # #     print(key,value)
    # # # mass = prelim.FastFit(data[datatag])

   
    
    # # `2pt/pion/src5.0_snk5.0/pion/C13.b_4002/AMA`
    # # `2pt/pion_SP/src5.0_snk5.0/pion/C13.b_4002/AMA`
    # # `3pt_tsep10/NUCL_D_MIXED_NONREL_l0_g0/src5.0_snk5.0/qz+0_qy+0_qx+0/C13.b_4002/AMA`


    # def parse_baryon_tag(datatag):
    #     datatag_split = datatag.split('/')
    #     corr_type     = datatag_split[0]
    #     tsep          = int(corr_type.split('_tsep')[1])
    #     buffer        =  datatag_split[1]
    #     channel       = buffer.split('_')[0]
    #     quark_ins       = buffer.split('_')[1]
    #     spin_proj       = buffer.split('_')[2]
    #     quark_sep       = buffer.split('_')[3]
    #     gamma           = buffer.split('_')[4] #gamma matrix of quark bilinear operator in the CHROMA convention , value accessed via dict
    #     src_snk_sep     = datatag_split[2]
    #     mom         = datatag_split[3]
    #     mom0       = mom.split('_')[0]
    #     mom1       = mom.split('_')[1]
    #     mom2       = mom.split('_')[2]
    #     momentum        = (mom0,mom1,mom2)
    #     config   = datatag_split[4]

    #     data_dict = dict()
    #     data_dict['corr_type']   = corr_type
    #     data_dict['tsep']        = tsep
    #     data_dict['buffer']      = buffer
    #     data_dict['channel']     = channel
    #     data_dict['quark_ins']   = quark_ins
    #     data_dict['spin_proj']   = spin_proj
    #     data_dict['quark_sep']   = quark_sep
    #     data_dict['gamma']       = gamma
    #     data_dict['src_snk_sep'] = src_snk_sep
    #     data_dict['mom']         = momentum
    #     data_dict['config']      = config
    #     return data_dict

    # re-weight correlators?
    try:
        reweight = fp.reweight
    except:
        reweight = False
    # block data
    bl = args.block
    if 'block' in dir(fp) and fp.block != 1:
        bl = fp.block
    if args.block != 1: # allow cl override
        bl = args.block
    if reweight:
        rw_files = fp.rw_files
        rw_path = fp.rw_path
        gv_data = ld.load_h5(fp.data_file, fp.corr_lst, rw=[rw_files, rw_path], bl=bl,
                             uncorr_corrs=args.uncorr_corrs, uncorr_all=args.uncorr_all)
        data_cfg = ld.load_h5(fp.data_file, fp.corr_lst, rw=[rw_files, rw_path], bl=bl,
                             uncorr_corrs=args.uncorr_corrs, uncorr_all=args.uncorr_all, return_gv=False, verbose=False)
    else:
        gv_data = ld.load_h5(fp.data_file, fp.corr_lst, bl=bl,
                             uncorr_corrs=args.uncorr_corrs, uncorr_all=args.uncorr_all)
        data_cfg = ld.load_h5(fp.data_file, fp.corr_lst, bl=bl,
                             uncorr_corrs=args.uncorr_corrs, uncorr_all=args.uncorr_all, return_gv=False, verbose=False)
    if args.states:
        states = args.states
    else:
        states = fp.fit_states

    """The number of configurations present.
        Use for time history plot."""

    nconfigs = [val.shape[0] for val in data_cfg.values()]
    nconfigs = np.unique(nconfigs).item()
    # print(nconfigs)

    ds = {key: gv.dataset.avg_data(data_cfg) for key, val in data_cfg.items()}
    # print(ds)
    x,y,n_states,priors_  = ld.make_fit_params(fp=fp,states=states,gv_data=gv_data)
    fit_funcs = cf.FitCorr()

     # if args.prior_override:
    #     priors = 

    fit_lst, p0, x_fit, y_fit = fit_funcs.get_fit(priors=priors_, states=states,x=x,y=y)
    if x_fit.keys() != y_fit.keys():
        raise ValueError("keys: fit_states should be shared in both fit dicts")
   
    axial_num_gv = {}
    vector_num_gv = {}
    corr_gv = {}
    axial_num_gv['SS'] = gv_data['gA_SS']
    axial_num_gv['PS'] = gv_data['gA_PS']
    vector_num_gv['SS'] = gv_data['gV_SS']
    vector_num_gv['PS'] = gv_data['gV_PS']
    corr_gv['SS'] = gv_data['proton_SS']
    corr_gv['PS'] = gv_data['proton_PS']
    plot.plot_effective_g00(axial_num_gv, corr_gv, 1, 14,observable='gA')
    plot.plot_effective_g00(vector_num_gv, corr_gv, 1, 14,observable='gV')
    plot.plot_effective_mass(corr_gv)
    # plot.plot_correlator_summary(corr_gv)
    mass = prelim.FastFit(gv_data['proton_PS'])
    # print(mass.E)
    # print(y)
    # print(y.items())
    data_chop = dict()
    if args.svd_test:
        for d in y:
            if d in x_fit:
                if d in x_fit and 'mres' not in d:
                    data_chop[d] = data_cfg[d][:,x_fit[d]['t_range']]
                    # print(data_chop)
   
        svd_test = gv.dataset.svd_diagnosis(data_chop, nbstrap=args.svd_nbs)
        svdcut = svd_test.svdcut
        print(svdcut)
    #     has_svd = True
    #     if args.svdcut is not None:
    #         print('    s.svdcut = %.2e' %svd_test.svdcut)
    #         print(' args.svdcut = %.2e' %args.svdcut)
    #         use_svd = input('   use specified svdcut instead of that from svd_diagnosis? [y/n]\n')
    #         if use_svd in ['y','Y','yes']:
    #                 svdcut = args.svdcut
    # if has_svd:
    #         fit = lsqfit.nonlinear_fit(data=(x_fit, y_fit), prior=priors, p0=None, fcn=fit_funcs.fit_function,
    #                                    svdcut=svdcut)
    # else:
    #     fit = lsqfit.nonlinear_fit(
    #         data=(x_fit, y_fit), prior=priors, p0=None, fcn=fit_funcs.fit_function)
    # if args.verbose_fit:
    #     print(fit.format(maxline=True))
    # else:
    #     print(fit)
    # y = {}
    print(x_fit.items())
    ydict = {tag: val for tag, val in x_fit.items() if isinstance(tag, int)}
    print(ydict)
    
    # c_t = gv_data['proton_SS'][0::1][:-1]
    # c_tpdt = gv_data['proton_SS'][0::1][1:]
    # print(np.array((1/2)*np.log(c_t/c_tpdt)))
    # # cf.effective_mass_local(gv_data)
    # # ds = {key: val for key, val in gv_data.items()}
    # cosh_m = (data[2:] + data[:-2]) / (2.0 * data[1:-1])
    # meff = np.zeros(len(cosh_m), dtype=gv._gvarcore.GVar)
    # # The domain of arccosh is [1, Infinity).
    # # Set entries outside of the domain to nans.
    # domain = (cosh_m > 1)
    # meff[domain] = np.arccosh(cosh_m[domain])
    # meff[~domain] = gv.gvar(np.nan)
    # print(meff)
    # ydict = {tag: val for tag, val in gv_data.items() if isinstance(tag, int)}
    # print(ydict)

    # y[t_snk] = np.ones(len(t), dtype=object)
    # denom_src = np.ones(len(t), dtype=object)
    # denom_snk = np.ones(len(t), dtype=object)
    # print(corr_gv)
    Tags = ['proton_SS','proton_PS']
    # ydict = {}
    # for t in Tags:
    #     ydict = x_fit[t]['t_range']
    # t_snks = np.sort(ydict)
    # t_snks
    # Tags_ff = ['gA_SS','gA_PS']
    # tag_3pt = {'gA_SS','gA_PS'}
    # print(data[tag].shape)
    # print(y_fit.keys())
    # Tags = y_fit.keys()
    # print(Tags)
    nt = np.unique([len(arr) for arr in y_fit.values()]).item()
    print(nt)
   
    # c2 = {}
    # for tag in Tags:
    #     c2[tag] = cf.C_2pt(
    #         tag, data_chop[tag], noise_threshy=0.03, nt=nt,skip_fastfit=True)
    # # c2_ = {tag: c2[tag].avg() for tag in c2}
    # # print(c2_)
    # c2_snk = c2['proton_PS']
    # print(c2_snk.mass)
    # c2_src = c2['proton_SS']
    # nts = [c2_src.times.nt,
    #            c2_snk.times.nt,
    #            c3.times.nt]
    # for nt in nts:
    #         if not np.all(nt == nts[0]):
    #             raise ValueError('tdata does not match across correlators')
    # fitter = C_2pt_Analysis(c2_src)
    # # fit = fitter.run_fit()
    
    # fit = fitter.run_fit()
    # print(priors.MesonPrior())
    
    # these should go in input file 
    

    ydict = {
        8:  [0,9],
        10: [0,11],
        12: [0,13],
        14: [0,15]
    }
   
    
    tag = 'SS'
    tag_ = 'PS'
    fit_out = test_NPoint(tag,corr_gv,prior=priors)
    fit_ = test_NPoint_snk(tag_,corr_gv,prior=priors)
    print(fit_)
    print(fit_out)
    plot.plot_correlators(corr_gv,t_plot_max=20)
    # plot.plot_effective_mass(corr_gv,fit=fit_out,show_fit=False)
    fit_3pt = test_NPoint_3pt('gA_PS',ydict,fit_out,fit_)
    print(fit_3pt)
    # test_BaseTimes()
   
def test_main():
    """Test the main function."""
    corr.main()

def test_BaseTimes():
    """Test TimeContainer.BaseTimes class."""
    tdata = range(100)
    tmin = 1
    tmax = 50
    nt = 200
    tp = -1

    times = cf.TimeContainer(tdata)
    print(times)
    assert times.tmin == 5, "Expected default tmin=5"
    assert times.tmax == len(tdata) - 1, "Expected default tmax=len(tdata)-1"
    assert times.nt == len(tdata), "Expected default nt=len(tdata)"
    assert times.tp == times.nt, "Expected default tp=nt"

    with pytest.raises(ValueError):
        times = cf.TimeContainer(tdata, tmin=-1)

    with pytest.raises(ValueError):
        times = cf.TimeContainer(tdata, tmax=len(tdata)+1)

def test_NPoint(tag,data,prior): #prior
    """Test cf.C_2pt and cf.C_3pt class."""
    # print(data[tag].shape)
    nt = data[tag].shape
    # print(corr)
    data_ = data.pop(tag)
    c2_src = cf.C_2pt(tag, data_)
    
    # print(len(c2_src.meff(avg=True)))
    model =get_two_point_model(c2_src)
    # t_start = c2_src.times.tmin 
    # t_end = c2_src.times.tmax
    Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
    nstates = Nstates(n=1, no=0)
    # prior = priors.MesonPriorPDG(nstates, 'pi',a_fm = .09)
    fitter = C_2pt_Analysis(c2_src)
    # fit = fitter.run_fit()

    fit = fitter.run_fit()

    # c2.__setitem__
    c2_src[0] = 1.0
    print(fit)
    # # Figures
    _ = plot.plot_correlators(data,t_plot_max=20)
    _ = plot.plot_effective_mass(data, 1, 16)
    return c2_src

def test_NPoint_snk(tag,data,prior):
    # tag = 'PS'
    nt = data[tag].shape
    data_ = data.pop(tag)
    c2_snk = cf.C_2pt(tag, data_)
    # print(c2_snk)
    assert len(c2_snk) == nt[0],\
        "Unexpected len(c2_snk)"
    assert len(c2_snk[:]) == nt[0],\
        "Unexpected len(c2_snk[:])"
    Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
    nstates = Nstates(n=1, no=0)
    # n=nstates.n, no=nstates.no
    # prior = priors.MesonPriorPDG(nstates, 'pi',a_fm = .09)
    
    # Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
    # nstates = Nstates(n=1, no=0)
    fitter = C_2pt_Analysis(c2_snk)
    fit = fitter.run_fit()
    print(fit)
    _ = plot.plot_correlators(data,t_plot_max=20)
    _ = plot.plot_effective_mass(data, 1, 16)
    return c2_snk

    
    

def test_NPoint_3pt(tag,data,c2_src,c2_snk):
    # nt = data[tag].shape
    # print(nt)
    # ds = {key: val for key, val in data.items()}
    c3 = cf.C_3pt(tag, data)
    # print(c3.ydict)
    # avg = c3.avg(m_src=c2_src.mass, m_snk=c2_snk.mass)
    # prior = priors.vmatrix(nstates)
    Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
    nstates = Nstates(n=1, no=0)
    avg = c3.avg(m_src=c2_src.mass, m_snk=c2_snk.mass)
    fitter = C_3pt_Analysis(c3)
    fit = fitter.run_sequential_fits(nstates)
    print(fit)
    return c3

def get_two_point_model(two_point, osc=True):
    tag = two_point.tag
    a_pnames = (f'{tag}:a', f'{tag}:ao')
    b_pnames = (f'{tag}:a', f'{tag}:ao')
    dE_pnames = (f'{tag}:dE', f'{tag}:dEo')

    if not osc:
        a_pnames = a_pnames[0]
        b_pnames = b_pnames[0]
        dE_pnames = dE_pnames[0]

    model = corrfitter.Corr2(
        datatag=tag,
        tp=two_point.times.tp,
        tmin=two_point.times.tmin,
        tmax=two_point.times.tmax,
        tdata=two_point.times.tdata,
        a=a_pnames,
        b=b_pnames,
        dE=dE_pnames,
        s=(1.0, -1.0)
    )
    return model

def get_three_point_model(t_snk, tfit, tdata, nstates, tags=None, constrain=False):
    """Gets a model for a 3pt function."""
    if tags is None:
        tags = Tags(src='SS', snk='PS') #this should be state+_sink
    src = tags.src
    snk = tags.snk
    # if max value of t_range of fit exceeds that of data, abort 
    if max(tfit) > max(tdata):
        LOGGER.error('Caution: max(tfit) exceeds max(tdata)')
        LOGGER.error('Restrict max(tfit) to max(tdata)')
        raise ValueError('Error: invalid tfit.')
    if tfit.size and np.all(np.isin(tfit, tdata)):
        a_pnames = (f'{src}:a', f'{src}:ao')
        b_pnames = (f'{snk}:a', f'{snk}:ao')
        dEa_pnames = (f'{src}:dE', f'{src}:dEo')
        dEb_pnames = (f'{snk}:dE', f'{snk}:dEo')
        vnn = 'Vnn'
        von = 'Von'
        vno = 'Vno'
        voo = 'Voo'

        if nstates.no == 0:
            a_pnames = a_pnames[0]
            dEa_pnames = dEa_pnames[0]
            von = None
            voo = None
        if nstates.mo == 0:
            b_pnames = b_pnames[0]
            dEb_pnames = dEb_pnames[0]
            vno = None
            voo = None
        # if constrain:
        #     _Model = models.ConstrainedCorr3 #TODO hold all fit params fixed except for matrix elements V__ 
       
        # Create 3-pt model to fit data from Lepage's library 
        model = corrfitter.Corr3(
            datatag=t_snk, T=t_snk, tdata=tdata, tfit=tfit,
            # Amplitudes in src 2-pt function
            a=a_pnames,
            # Amplitudes in snk 2-pt function
            b=b_pnames,
            # Energies in src 2-pt function
            dEa=dEa_pnames,
            # Energies in src 2-pt function
            dEb=dEb_pnames,
            # sign factors in src 2-pt function
            sa=(1.0, -1.0),
            # sign factors in snk 2-pt function
            sb=(1.0, -1.0),
            # connect src decay --> snk decay
            Vnn=vnn,
            # connect src decay --> snk oscillating
            Vno=vno,
            # connect src oscillating --> snk decay
            Von=von,
            # connect src oscillating --> snk oscillating
            Voo=voo
        )
    else:
        # Empty tfit -- no model
        model = None

    return model

def get_model(ds, tag, nstates, constrain=False):
    """Gets a corrfitter model"""
    # TODO implement constrain 
    if isinstance(ds[tag], cf.C_2pt):
        osc = bool(nstates.no) if tag == ds.tags.src else bool(nstates.mo)
        return get_two_point_model(ds[tag], osc)
    if isinstance(tag, int):
        t_snk = tag
        tdata = ds.c3.times.tdata
        return get_three_point_model(t_snk, ds.tfit[t_snk], tdata, nstates,
                                     constrain=constrain)

def count_nstates(params, key_map=None, tags=''):
    """
    Count the number of states used in fit.
    Default behavior assumes names 'SS' and 'PS' for src and
    snk, respectively.
    """
    if tags is None:
        tags = Tags(src='SS', snk='PS')
    src = tags.src
    snk = tags.snk
    if key_map is None:
        key_map = {
            'n': f'{src}:dE', 'no': f'{src}:dEo',
            'm': f'{snk}:dE', 'mo': f'{snk}:dEo',
        }
    kwargs = {key1: len(params.get(key2, [])) for key1, key2 in key_map.items()}
    return Nstates(**kwargs)

def compute_yfit(ds, params):
    """
    Computes the model values "yfit" from fit params to be compared with
    data stored in "ds". Should share same structure.
    Args:
        ds:  correlated data
        params: dict of fit parameters
    Returns:
        yfit: dict
    """
    nstates = count_nstates(params)
    yfit = {}
    for tag in ds:
        model = get_model(ds, tag, nstates)
        tdata = np.array(model.tdata, dtype=float)
        yfit[tag] = model.fitfcn(t=tdata, p=params)
    return yfit

class C_2pt_Analysis(object):
    """
    A basic fitter class for two-point correlation functions.
    Args:
        c2: corr_functions.C_2pt object
    Returns:
        lsqfit fit object
    """
    def __init__(self, c2):
        self.tag = c2.tag
        self.c2 = c2
        self.prior = None 
        self.fitter = None
        self._nstates = None
        self._fit = None

    def run_fit(self, nstates=Nstates(1, 0), prior=None, **fitter_kwargs):
        """
        Run the fit.
        Args:
            nstates: tuple (n_decay, n_osc) specifying the number of decaying
                and oscillating states, respectively. Defaults to (1,0).
            prior: dict from label file. Default is None, for which the fitter
                tries to constuct a prior itself. TODO 
        """
        self._nstates = nstates
        if prior is None:
            prior = priors.MesonPrior(
                nstates.n, nstates.no, amps=['a', 'ao'],
                tag=self.tag, ffit=self.c2.fastfit,
                extend=True)
        #     ) #TODO this should try a meson/baryon prior class first then resort to
        #       dict of priors in label file 
        self.prior = prior
        # Model construction infers the fit times from c2
        model = get_two_point_model(self.c2, bool(nstates.no))
        self.fitter = corrfitter.CorrFitter(models=model)
        data = {self.tag: self.c2}
        fit = self.fitter.lsqfit(data=data, prior=prior, p0=None, **fitter_kwargs)
        self._fit = fit
        # fit.show_plots()

        # if fit.failed:
        #     fit = None
        return fit

class C_3pt_Analysis(object):
    '''
    Perform simult. fits to 2pt, 3pt fcns
    Output:
        form factors: gA, gV  
    '''
    def __init__(self, ds, positive_ff=True):

        self.ds = ds #gvar dataset
        self.positive_ff = positive_ff #forces form factor to be positive
        self.prior = None
        self.fits = {} #fit dict to fill 
        self.fitter = None
        self.r = None

    def run_sequential_fits(
            self, nstates, tmin_override=None,
            prior=None, chain=False, constrain=False,
            **fitter_kwargs):
        """
        First runs two-point functions, Then runs the simult fit.
        """
        self.prior = prior
        if tmin_override is not None:
            if tmin_override.src is not None:
                self.ds.c2_src.times.tmin = tmin_override.src
            if tmin_override.snk is not None:
                self.ds.c2_snk.times.tmin = tmin_override.snk
        self.fit_two_point(
            nstates=nstates,
            **fitter_kwargs)
        self.fit_three_point(
            nstates=nstates,
            chain=chain,
            constrain=constrain,
            **fitter_kwargs)

    def mass(self, tag):
        """Gets the mass/energy of the ground state from full fit."""
        params = self.fits['full'].p
        return params[f'{tag}:dE'][0]

    @property
    def m_src(self):
        """Gets the mass/energy of the "source" ground state."""
        src_tag = self.ds.tags.src
        return self.mass(src_tag)

    @property
    def m_snk(self):
        """Gets the mass/energy of the "snk" ground state."""
        snk_tag = self.ds.tags.snk
        return self.mass(snk_tag)

    @property
    def matrix_element(self):
        """Gets the matrix element Vnn[0, 0] needed for the form factor."""
        if self.fits['full'] is not None:
            return self.fits['full'].p['Vnn'][0, 0]

    @property
    def r_prior(self):
        src_tag = self.ds.tags.src
        m_src = gv.mean(self.prior[f'{src_tag}:dE'][0])
        matrix_element = self.prior['Vnn'][0, 0]
        return convert_vnn_to_ratio(m_src, matrix_element)

    def fit_two_point(self, nstates, **fitter_kwargs):
        """Run the fits of two-point functions."""
        for tag in self.ds.c2:
            _nstates = nstates
            if tag == self.ds.tags.snk:
              _nstates = Nstates(n=nstates.m, no=nstates.mo)
            
            fit = C_2pt_Analysis(self.ds.c2[tag]).\
                run_fit(_nstates, **fitter_kwargs)
            if fit is None:
                LOGGER.warning('Fit failed for two-point function %s.', tag)
            else:
                pass
            # TODO: handle possible re-running if fit fails initially
                #TODO generate prior from fit to be reused 
                # self.prior.update(
                #     update_with=fit.p, width=width, fractional_width=fractional_width)
            self.fits[tag] = fit #tags are of form state+_sink

    def fit_simult(self, nstates, chain=False, constrain=False, **fitter_kwargs):
        """Run the simultaneous fit of 2- and 3-point functions for form factor."""

        # Handle prior
        prior = fitter_kwargs.get('prior')
        if prior is not None:
            self.prior = prior
        else:
            prior = self.prior
        fitter_kwargs['prior'] = prior

        # Handle models
        models_list = []
        for tag in self.ds:
            model = get_model(self.ds, tag, nstates,constrain)
            if model is not None:
                models_list.append(model)

        # Abort if too few models found; There should be a model corresponding to each key in dataset
        if len(models_list) != len(set(self.ds.keys())):
            self.fitter = None
            fit = None
            LOGGER.warning('Insufficient models found. Skipping joint fit.')
            return

        # Run fit
        self.fitter = corrfitter.CorrFitter(models=models_list) #call Lepage's library
        if chain:
            _lsqfit = self.fitter.chained_lsqfit #chained lsqfit 
        else:
            _lsqfit = self.fitter.lsqfit
        fit = _lsqfit(data=self.ds, **fitter_kwargs)
        self.fits['full'] = fit
        if fit.failed:
            LOGGER.warning('Full joint fit failed.')
        else:
            self.ds.set_masses(fit.p['proton_SS:dE'][0],
                               fit.p['proton_PS:dE'][0])

       
        vnn = fit.p['Vnn'][0, 0]
        self.r = convert_vnn_to_ratio(self.m_src, vnn)
    #TODO routine to convert fit into plaintext, possible use for bootstrapping

        # TODO call plotting commans here 

        # def plot_results():

        # def plot_E_summary():

        # def plot_z_summary():

        # def plot_gA():

        # def plot_gV():

class Ratio:
    '''
    Analysis of :math:`\< \frac{C3pt}{C2pt} \>`
    Returns:
    matrix elements eg. form factors(gA,gV) at each T
    '''
    def __init__(self,x_params,fit_params):
        self.x_params   = x_params 
        self.fit_params = p
        self._fit = None
        fcn = cf.CorrFunction()
    def get_ratio_model(self):
        result = {}
        for t_snk, t in self.fit_params.items():

            result[t_snk] = fcn.exp_open(self.x_params, self.fit_params)

        return result 

    def make_fit_data(self,tmin_src,tmin_snk):
        # keys for dicts: T values, values: Ratio(T) 
        y_fit = {}
        x_fit = {}
        pt3_data_start = 2 #tmin_src
        pt3_data_end = 15 # T+1- tmin_snk
        linspace_num = 100
        #tsep values are keys, range is value
        gA_tsep = []
        gA_tau = []
        for t in range(pt3_data_start, pt3_data_end): # tsep and tau to calc fit values of no transition 
            if t % 2 == 0: # only when tsep is even, tau = tsep/2
                gA_tsep.append(t)
                gA_tau.append(t/2) # tau = tsep/2



    def make_prior(self):
        prior = {}

    def __call__(self):
        x_fit, y_fit = self.make_fit_data(tmin_src, tmin_snk)
        prior = self.make_prior()
        model = self.get_ratio_model()
        fit = lsqfit.nonlinear_fit(data=(x_fit, y_fit), fcn=model, prior=prior, **fitter_kwargs)
        self._fit = fit #this just fills empty fit with generated fit, lsqfit technicality? 
        return fit 




# class SequentialFit_out:
#     def __init__(self):
#         self.src = None
#         self.snk = None
#         self.ratio = None
#         self.direct = None

#     def __iter__(self):
#         for fit in [self.src, self.snk]:#, self.ratio, self.direct]:
#             yield fit

#     def asdict(self):
#         return self.__dict__

# class SequentialFitter:
#     """
#     Run a sequential set of fits in order to determine a matrix element / form factor.
#     Args:
#         data: gvar dataset from load_data.gv_data
#         a_fm: lattice spacing in fm
#     Notes:
#     ------
#     The sequence of fits is
#         1) Fit the "source" 2pt function
#         2) Fit the "sink" 2pt function
#         3) Fit the ratio Rbar of 3pt and 2pt function
#         4) Fit the spectral decompostion directly
#         - simult fit with same covariance matricx (block diagonl mat)
#     """
#     def __init__(self, data, a_fm,prior):
#         self.data = data
#         self.a_fm = a_fm
#         self.prior = prior
#         self.fits = SequentialFit_out()
#         self.r_ratio = None
#         self.r_direct = None


#     def get_source_fit(self, n, no, **fitter_kwargs):
#         """ Fits the source 2pt function. """
#         tag = 'SS'
#         c2 = self.data.c2_src

#         c2.tag = tag
#         nstates = Nstates(n=n, no=no)
#         prior = self.prior
#         fitter = TwoPointAnalysis(c2)
#         fit = fitter.run_fit(nstates, prior=prior, **fitter_kwargs)
#         self.fits.src = fit

#     def get_sink_fit(self, m, mo, **fitter_kwargs):
#         """ Fits the sink 2pt function. """
#         tag = 'PS'
#         c2 = self.data.c2_snk

#         c2.tag = tag
#         nstates = Nstates(n=m, no=mo)
#         prior = self.prior
#         _Times = collections.namedtuple('Times', ['tmin_src', 'tmin_snk', 't_step'])
#         times = _Times(tmin_src, tmin_snk, t_step)

#         fitter = TwoPointAnalysis(c2)
#         fit = fitter.run_fit(nstates, prior=prior, **fitter_kwargs)
#         self.fits.snk = fit

#     def get_ratio_fit(self, n, m, tmin_src, tmin_snk, t_step, **fitter_kwargs):




#         # if args.fit_2pt:
#         # if args.fit_3pt:
#         # if args.fit_ratio:
#         # if ars.fit_sequential:
#     # def __call__(self, nstates, times, p2_boost=None, **fitter_kwargs):
#     #     """
#     #     Runs the sequential fits.
#     #     Args:
#     #         TODO
#     #     Returns:
#     #         TODO
#     #     """
#     #     # Set times once and for all
#     #     self.data.c2_src.times.tmin = times.tmin_src
#     #     self.data.c2_src.times.tmax = times.tmax_src
#     #     self.data.c2_snk.times.tmin = times.tmin_snk
#     #     self.data.c2_snk.times.tmax = times.tmax_snk

#     #     self.run_source(n=nstates.n, no=nstates.no, p2_boost=p2_boost, **fitter_kwargs)

#     #     self.run_sink(m=nstates.m, mo=nstates.mo, **fitter_kwargs)

#     # def summarize(self):
#     #     """ Print a summary of the results"""

#     #     print(" Source Fit ".center(80, "#"))
#     #     print(self.fits.src.format(maxline=False))

#     #     print(" Sink Fit ".center(80, "#"))
#     #     print(self.fits.snk.format(maxline=False))

if __name__ == '__main__':
    main()
