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

import fitter.corr_functions as cf 
import fitter.load_data as ld
import fitter.plotting as plot
import fitter.fastfit as prelim #LePage's preliminary corrfitter to generate p0
import fitter.priors as priors 

SVDCUT = 0.002
Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
def main():
    parser = argparse.ArgumentParser(
        description='Perform analysis of two-point correlation function')
    parser.add_argument('fit_params',    help='input file to specify fit')
    parser.add_argument('--fit',         default=True, action='store_true',
                        help=            'do fit? [%(default)s]')
    parser.add_argument('-b', '--block', default=1, type=int,
                        help=            'specify bin/blocking size in terms of saved configs')
    parser.add_argument('--uncorr_corrs', default=False, action='store_true',
                        help=            'uncorrelate different correlation functions? [%(default)s]')
    parser.add_argument('--uncorr_all',  default=False, action='store_true',
                        help=            'uncorrelate all snk,src for each correlation function? [%(default)s]')
    parser.add_argument('--states',      nargs='+',
                        help=            'specify states to fit?')                        

    args = parser.parse_args()
    # if args.save_figs and not os.path.exists('figures'):
    #     os.makedirs('figures')
    # print(args)
    # add path to the input file and load it
    sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
    fp = importlib.import_module(
        args.fit_params.split('/')[-1].split('.py')[0])

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

    x,y,n_states,priors  = ld.make_fit_params(fp=fp,states=states,gv_data=gv_data)
    fit_funcs = cf.FitCorr()

    fit_lst, p0, x_fit, y_fit = fit_funcs.get_fit(priors=priors, states=states,x=x,y=y)
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
   
    mass = prelim.FastFit(gv_data['proton_PS'])
    print(mass)
    # y = {}

    
    
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
    tag = 'proton_SS'
    tag_ = 'proton_PS'
    tag_3pt = {'gA_SS','gA_PS'}
    # print(data[tag].shape)
    
    fit_out = test_NPoint(tag,gv_data,prior=priors)
    print(fit_out)
    fit_ = test_NPoint_snk(tag_,gv_data,prior=priors)
    print(fit_)
    # plot.plot_effective_mass(corr_gv,fit=fit_out,show_fit=True)
    test_NPoint_3pt('gA_SS',gv_data)
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

def test_NPoint(tag,data,prior):
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
    fitter = C_2pt_Analysis(c2_src,prior)
    # fit = fitter.run_fit()

    fit = fitter.run_fit()


    assert len(c2_src) == nt[0],\
        "Unexpected len(c2)"
    assert len(c2_src[:]) == nt[0],\
        "Unexpected len(c2[:])"
    assert c2_src.times.tmax == (len(c2_src) - 1),\
        "Unexpected c2.times.tmax"
    assert len(c2_src.meff(avg=False)) == (len(c2_src) - 2),\
        "Unexpected len(c2_src.meff())"

    # c2.__setitem__
    c2_src[0] = 1.0

    # # Figures
    _ = plot.plot_correlators(data,t_plot_max=20)
    _ = plot.plot_effective_mass(data, 1, 16)
    # _ = c2_src.plot_corr(avg=False)
    # _ = c2_src.plot_meff(avg=False)
    # _ = c2_src.plot_meff(avg=True)
    return fit

def test_NPoint_snk(tag,data,prior):
    # tag = 'PS'
    nt = data[tag].shape
    data_ = data.pop(tag)
    c2_snk = cf.C_2pt(tag, data_,skip_fastfit=False)
    print(c2_snk)
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
    fitter = C_2pt_Analysis(c2_snk,prior)
    fit = fitter.run_fit()
    print(fit)

    # _ = plot.plot_correlators(data,t_plot_max=20)
    # _ = plot.plot_effective_mass(data, 1, 16)

def test_NPoint_3pt(tag,data):
    # nt = data[tag].shape
    # print(nt)
    ds = {key: val for key, val in data.items()}
    c3 = cf.C_3pt(tag, ds)
    print(c3)
    # avg = c3.avg(m_src=c2_src.mass, m_snk=c2_snk.mass)
    # prior = priors.vmatrix(nstates)
    Nstates = collections.namedtuple('NStates', ['n', 'no', 'm', 'mo'], defaults=(1, 0, 0, 0))
    nstates = Nstates(n=1, no=0)
   
    fitter = ThreePointAnalysis(c3)
    fit = fitter.run_sequential_fits(nstates)
    print(fit)

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
        tags = dataset.Tags(src='SS', snk='PS')
    src = tags.src
    snk = tags.snk
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
        if constrain:
            _Model = models.ConstrainedCorr3
        else:
            _Model = models.Corr3
        model = _Model(
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
        tags = dataset.Tags(src='SS', snk='PS')
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
        ds: FormFactorDataset with the data
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
        c2: C_2pt object
    """
    def __init__(self, c2,prior):
        self.tag = c2.tag
        self.c2 = c2
        self.prior = prior 
        self.fitter = None
        self._nstates = None
        self._fit = None

    def run_fit(self, nstates=Nstates(1, 0), prior=None, **fitter_kwargs):
        """
        Run the fit.
        Args:
            nstates: tuple (n_decay, n_osc) specifying the number of decaying
                and oscillating states, respectively. Defaults to (1,0).
            prior: BasicPrior object. Default is None, for which the fitter
                tries to constuct a prior itself. TODO 
        """
        self._nstates = nstates
        if prior is None:
            prior = priors.MesonPrior(
                nstates.n, nstates.no, amps=['a', 'ao'],
                tag=self.tag, ffit=self.c2.fastfit,
                extend=True
            )
        self.prior = prior
        # Model construction infers the fit times from c2
        model = get_two_point_model(self.c2, bool(nstates.no))
        self.fitter = corrfitter.CorrFitter(models=model)
        data = {self.tag: self.c2}
        fit = self.fitter.lsqfit(data=data, prior=prior, p0=prior.p0, **fitter_kwargs)
        # fit = serialize.SerializableNonlinearFit(fit)
        self._fit = fit
        fit.show_plots()

        # if fit.failed:
        #     fit = None
        return fit

class C_3pt_Analysis(object):
    '''
    Extract form factors from simult. fits to 2pt, 3pt fcns  
    '''
    def __init__(self, ds, positive_ff=True):

        self.ds = ds #dataset
        self.positive_ff = positive_ff
        self.prior = None
        self.fits = {}
        self.fitter = None
        self.r = None

    def run_sequential_fits(
            self, nstates, tmin_override=None,
            width=0.1, fractional_width=False,
            prior=None, chain=False, constrain=False,
            **fitter_kwargs):
        """
        Runs sequential fits.
        First runs two-point functions, Then runs the simult fit.
        """
        if prior is None:
            self.prior = priors.FormFactorPrior(
                nstates,
                self.ds,
                positive_ff=self.positive_ff)
        else:
            self.prior = prior
        if tmin_override is not None:
            if tmin_override.src is not None:
                self.ds.c2_src.times.tmin = tmin_override.src
            if tmin_override.snk is not None:
                self.ds.c2_snk.times.tmin = tmin_override.snk
        self.fit_two_point(
            nstates=nstates,
            width=width,
            fractional_width=fractional_width,
            **fitter_kwargs)
        self.fit_form_factor(
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
        """Fetches the matrix element Vnn[0, 0] needed for the form factor."""
        if self.fits['full'] is not None:
            return self.fits['full'].p['Vnn'][0, 0]

    @property
    def r_prior(self):
        src_tag = self.ds.tags.src
        m_src = gv.mean(self.prior[f'{src_tag}:dE'][0])
        matrix_element = self.prior['Vnn'][0, 0]
        return convert_vnn_to_ratio(m_src, matrix_element)

    def fit_two_point(self, nstates, width=0.1, fractional_width=False, **fitter_kwargs):
        """Run the fits of two-point functions."""
        for tag in self.ds.c2:
            _nstates = nstates
            if tag == self.ds.tags.snk:
              _nstates = Nstates(n=nstates.m, no=nstates.mo)
            # TODO: handle possible re-running if fit fails initially
            # In the other code, I reset the priors on dEo to match dE
            fit = TwoPointAnalysis(self.ds.c2[tag]).\
                run_fit(_nstates, **fitter_kwargs)
            if fit is None:
                LOGGER.warning('Fit failed for two-point function %s.', tag)
            else:
                pass
                # self.prior.update(
                #     update_with=fit.p, width=width, fractional_width=fractional_width)
            self.fits[tag] = fit

    def fit_form_factor(self, nstates, chain=False, constrain=False, **fitter_kwargs):
        """Run the joint fit of 2- and 3-point functions for form factor."""

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

        # Abort if too few models found
        if len(models_list) != len(set(self.ds.keys())):
            self.fitter = None
            fit = None
            LOGGER.warning('Insufficient models found. Skipping joint fit.')
            return

        # Run fit
        self.fitter = cf.CorrFitter(models=models_list)
        if chain:
            _lsqfit = self.fitter.chained_lsqfit
        else:
            _lsqfit = self.fitter.lsqfit
        fit = _lsqfit(data=self.ds, **fitter_kwargs)
        self.fits['full'] = fit
        if fit.failed:
            LOGGER.warning('Full joint fit failed.')
        else:
            self.ds.set_masses(fit.p['SS:dE'][0],
                               fit.p['PS:dE'][0])

       
        vnn = fit.p['Vnn'][0, 0]
        self.r = convert_vnn_to_ratio(self.m_src, vnn)

class SequentialFitResult:
    def __init__(self):
        self.src = None
        self.snk = None
        self.ratio = None
        self.direct = None

    def __iter__(self):
        for fit in [self.src, self.snk]:#, self.ratio, self.direct]:
            yield fit

    def asdict(self):
        return self.__dict__

class SequentialFitter:
    """
    Run a sequential set of fits in order to determine a matrix element / form factor.
    Args:
        data: dataset.FormFactorDataset
        a_fm: approximate lattice spacing in fm
    Notes:
    ------
    The sequence of fits is
        1) Fit the "source" 2pt function
        2) Fit the "sink" 2pt function
        3) Fit the ratio Rbar of 3pt and 2pt function
        4) Fit the spectral decompostion directly
        - simult fit with same covariance matricx (block diagonl mat)
    """
    def __init__(self, data, a_fm):
        self.data = data
        self.a_fm = a_fm
        self.fits = SequentialFitResult()
        self.r_ratio = None
        self.r_direct = None


    def run_source(self, n, no, p2_boost=None, **fitter_kwargs):
        """ Fits the source 2pt function. """
        tag = 'SS'
        c2 = self.data.c2_src

        c2.tag = tag
        nstates = Nstates(n=n, no=no)
        prior = priors.MesonPrior(
                nstates.n, nstates.no, amps=['a', 'ao'],
                tag=c2.tag, ffit=self.c2.fastfit,
                extend=True)

        # Boost the energies as necessary
        # if p2_boost:
        #     prior[f"{tag}:dE"] = bayes_prior.boost(prior[f"{tag}:dE"], p2_boost)

        fitter = TwoPointAnalysis(c2)
        fit = fitter.run_fit(nstates, prior=prior, **fitter_kwargs)
        self.fits.src = fit

    def run_sink(self, m, mo, **fitter_kwargs):
        """ Fits the sink 2pt function. """
        tag = 'PS'
        c2 = self.data.c2_snk

        c2.tag = tag
        nstates = Nstates(n=m, no=mo)
        prior = priors.MesonPrior(
                nstates.n, nstates.no, amps=['a', 'ao'],
                tag=c2.tag, ffit=self.c2.fastfit,
                extend=True)
        _Times = collections.namedtuple('Times', ['tmin_src', 'tmin_snk', 't_step'])
        times = _Times(tmin_src, tmin_snk, t_step)

        fitter = TwoPointAnalysis(c2)
        fit = fitter.run_fit(nstates, prior=prior, **fitter_kwargs)
        self.fits.snk = fit

    def __call__(self, nstates, times, p2_boost=None, **fitter_kwargs):
        """
        Runs the sequential fits.
        Args:
            TODO
        Returns:
            TODO
        """
        # Set times once and for all
        self.data.c2_src.times.tmin = times.tmin_src
        self.data.c2_src.times.tmax = times.tmax_src
        self.data.c2_snk.times.tmin = times.tmin_snk
        self.data.c2_snk.times.tmax = times.tmax_snk

        self.run_source(n=nstates.n, no=nstates.no, p2_boost=p2_boost, **fitter_kwargs)

        self.run_sink(m=nstates.m, mo=nstates.mo, **fitter_kwargs)

    def summarize(self):
        """ Print a summary of the results"""

        print(" Source Fit ".center(80, "#"))
        print(self.fits.src.format(maxline=False))

        print(" Sink Fit ".center(80, "#"))
        print(self.fits.snk.format(maxline=False))

if __name__ == '__main__':
    main()
