import fitter.lanl_plotting as plot
import fitter.corr_functions as cf
import fitter.lanl_load_data as ld
import fitter.bootstrap as bs 
import sys 
import lsqfit
import gvar as gv
import importlib
import h5py as h5 
import copy
import os
import sys
import pathlib
import random
import matplotlib.pyplot as plt
import argparse
import numpy as np
np.seterr(invalid='ignore')
#import tables as h5


# user libs
#sys.path.append('util')


def main():
    # do sys.argv routine here 
    parser = argparse.ArgumentParser(
        description='Perform analysis of two-point correlation function')
    parser.add_argument('fit_params',    help='input file to specify fit')
    parser.add_argument('--fit',         default=True, action='store_true',
                        help=            'do fit? [%(default)s]')
    parser.add_argument('--svdcut',      type=float, help='add svdcut to fit')
    parser.add_argument('--svd_test',    default=True, action='store_false',
                        help=            'perform gvar svd_diagnosis? [%(default)s]')
    parser.add_argument('--svd_nbs',     type=int, default=50, help='number of BS samples for estimating SVD cut [%(default)s]')
    parser.add_argument('--fold',        default=True, action='store_false',
                        help=            'fold data about T/2? [%(default)s]')
    parser.add_argument('-b', '--block', default=1, type=int,
                        help=            'specify bin/blocking size in terms of saved configs')
    parser.add_argument('--eff',         default=True, action='store_true',
                        help=            'plot effective mass and z_eff data? [%(default)s]')
    parser.add_argument('--scale',       default=None, nargs=2,
                        help=            'add right axis with physical scale specified by input value [scale, units]')
    parser.add_argument('--stability',   nargs='+',
                        help=            'specify states to perform t_min and n_state sweep')
    parser.add_argument('--es_stability', default=False, action='store_true',
                        help=            'plot excited state stability? [%(default)s]')
    parser.add_argument('--states',      nargs='+',
                        help=            'specify states to fit?')
    parser.add_argument('-v','--verbose', default=False,action='store_true',
                        help=            'add verbosity [%(default)s]')
    parser.add_argument('--verbose_fit', default=False, action='store_true',
                        help=            'print y vs f(x,p) also? [%(default)s]')
    parser.add_argument('--save_figs',   default=False, action='store_true',
                        help=            'save figs? [%(default)s]')
    parser.add_argument('--bs',          default=False, action='store_true',
                        help=            'run bootstrap fit? [%(default)s]')
    parser.add_argument('--Nbs',         type=int, default=20,
                        help=            'specify the number of BS samples to compute [%(default)s]')
    parser.add_argument('--bs_seed',     default=None,
                        help=            'set a string to seed the bootstrap - None will be random [%(default)s]')
    parser.add_argument('--bs_write',    default=True, action='store_false',
                        help=            'write bs results to file? [%(default)s]')
    parser.add_argument('--bs_results',  default='bs_results/spectrum_bs.h5',
                        help=            'set file to write bootstrap results [%(default)s]')
    parser.add_argument('--overwrite',   default=False, action='store_true',
                        help=            'overwrite existing bootstrap results? [%(default)s]')
    parser.add_argument('--bs_path',     default='spec',
                        help=            'specify path in h5 file for bs results')
    parser.add_argument('--uncorr_corrs', default=False, action='store_true',
                        help=            'uncorrelate different correlation functions? [%(default)s]')
    parser.add_argument('--uncorr_all',  default=False, action='store_true',
                        help=            'uncorrelate all snk,src for each correlation function? [%(default)s]')
    parser.add_argument('--interact',    default=False, action='store_true',
                        help=            'open IPython instance after to interact with results? [%(default)s]')
    parser.add_argument('--gui',    default=False, action='store_true',
                        help=            'open dashboard for analyzing fit. Must be used together with fit flag. [%(default)s]')

    args = parser.parse_args()
    if args.save_figs and not os.path.exists('figures'):
        os.makedirs('figures')
    print(args)
    # add path to the input file and load it
    sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
    fp = importlib.import_module(
        args.fit_params.split('/')[-1].split('.py')[0])

    

    # can only uncorrelate all or sets of corrs
    if args.uncorr_all and args.uncorr_corrs:
        sys.exit('you can only select uncorr_corrs or uncorr_all')
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
    print(gv_data,data_cfg)

    x,y,n_states,priors= plot.eff_plots.make_fit_params(fp=fp,states=states,gv_data=gv_data)
    print(x,y,n_states,priors)
    
    # # if args.eff:
    # #     plt.ion()
    # #     effective = plot.eff_plots()
    # #     effective.make_eff_plots(states=states,fp=fp,x_fit=None,priors=priors,gv_data=gv_data,fit=None, scale=args.scale,show_fit=False,save_figs=args.save_figs)
        
    # # # set up svdcut if added
    # # if args.svdcut is not None:
    # #     svdcut = args.svdcut
    # #     has_svd = True
    # # else:
    # #     try:
    # #         svdcut = fp.svdcut
    # #         has_svd = True
    # #     except:
    # #         has_svd = False

    # # if args.stability:
    # #     plot.eff_plots.make_stability_plot(
    # #     states=states,x=x,fp=fp, priors=priors, gv_data=gv_data, stability=args.stability, 
    # #     scale = args.scale, svd_test=args.svd_test, data_cfg = data_cfg,n_states=n_states, 
    # #     svd_nbs=args.svd_nbs, es_stability=args.es_stability,save_figs=args.save_figs)


    # # if args.fit:
    # #     fit_funcs = cf.FitCorr()
    # #     p0, x_fit, y_fit = fit_funcs.get_fit(priors=priors, states=states,x=x,y=y)

    # #     if args.svd_test:
    # #         data_chop = dict()
    # #         for d in y:
    # #             if d in x_fit:
    # #                 if d in x_fit and 'mres' not in d:
    # #                     data_chop[d] = data_cfg[d][:,x_fit[d]['t_range']]
    # #                 if 'mres' in d and len(d.split('_')) > 1:
    # #                     data_chop[d] = data_cfg[d][:,x_fit[d.split('_')[0]]['t_range']]
    # #         svd_test = gv.dataset.svd_diagnosis(data_chop, nbstrap=args.svd_nbs)
    # #         svdcut = svd_test.svdcut
    # #         has_svd = True
    # #         if args.svdcut is not None:
    # #             print('    s.svdcut = %.2e' %svd_test.svdcut)
    # #             print(' args.svdcut = %.2e' %args.svdcut)
    # #             use_svd = input('   use specified svdcut instead of that from svd_diagnosis? [y/n]\n')
    # #             if use_svd in ['y','Y','yes']:
    # #                 svdcut = args.svdcut
    # #     if has_svd:
    # #          fit = lsqfit.nonlinear_fit(data=(x_fit, y_fit), prior=priors, p0=p0, fcn=fit_funcs.fit_function,
    # #                                    svdcut=svdcut)
    # #     else:
    # #         fit = lsqfit.nonlinear_fit(
    # #             data=(x_fit, y_fit), prior=priors, p0=p0, fcn=fit_funcs.fit_function)
    # #     if args.verbose_fit:
    # #         print(fit.format(maxline=True))
    # #     else:
    # #         print(fit)

    # #     # if args.mres_avg:
    # #     #     plot.eff_plots.plot_mres(ax, dsets, key)

            


    # #     if args.gui:
    # #         from lsqfitgui import run_server
    # #         run_server(fit, name="c51 Two-Point Fitter")

    # #     x_plot = copy.deepcopy(x_fit)
    # #     #generate eff mass plot with fit overlay
    # #     if args.eff:
    # #         effective.make_eff_plots(states, fp, x_fit=x_fit, fit=fit,gv_data=gv_data, priors=priors, 
    # #                             scale=args.scale,show_fit=True,save_figs=args.save_figs)

    # #     # run bootstrapping utility 
    # #     # NOTE: we should abstract this into boostrap.py 
    # #     if args.bs:
    # #         # make sure results dir exists
    # #         if args.bs_write:
    # #             if not os.path.exists('bs_results'):
    # #                 os.makedirs('bs_results')
    # #         if len(args.bs_results.split('/')) == 1:
    # #             bs_file = 'bs_results/'+args.bs_results
    # #         else:
    # #             bs_file = args.bs_results
    # #         # check if we already wrote this dataset
    # #         if args.bs_write:
    # #             have_bs = False
    # #             if os.path.exists(bs_file):
    # #                 #with h5.open_file(bs_file,'r') as f5:
    # #                 with h5py.File(bs_file, 'r') as f5:
    # #                     if args.bs_path in f5:
    # #                         if len(f5[args.bs_path]) > 0 and not args.overwrite:
    # #                             have_bs = True
    # #                             print(
    # #                                 'you asked to write bs results to an existing dset and overwrite =', args.overwrite)
    # #         else:
    # #             have_bs = False
    # #         if not have_bs:
    # #             print('beginning Nbs=%d bootstrap fits' % args.Nbs)
    # #             import bootstrap as bs

    # #             # let us use the fit posterior to set the initial guess for bs loop
    # #             p0_bs = dict()
    # #             for k in fit.p:
    # #                 p0_bs[k] = fit.p[k].mean

    # #             if not args.bs_seed and 'bs_seed' not in dir(fp):
    # #                 tmp = input('you have not passed a BS seed nor is it defined in the input file\nenter a seed or hit return for none')
    # #                 if not tmp:
    # #                     bs_seed = None
    # #                 else:
    # #                     bs_seed = tmp
    # #             elif 'bs_seed' in dir(fp):
    # #                 bs_seed = fp.bs_seed
    # #             if args.bs_seed:
    # #                 if args.verbose:
    # #                     print('WARNING: you are overwriting the bs_seed from the input file')
    # #                 bs_seed = args.bs_seed

    # #             # make BS data
    # #             corr_bs = {}
    # #             for k in data_cfg:
    # #                 corr_bs[k] = bs.bs_corrs(data_cfg[k], Nbs=args.Nbs,
    # #                                             seed=bs_seed, return_mbs=True)
    # #             # make BS list for priors
    # #             p_bs_mean = dict()
    # #             for k in priors:
    # #                 p_bs_mean[k] = bs.bs_prior(args.Nbs, mean=priors[k].mean,
    # #                                         sdev=priors[k].sdev, seed=bs_seed+'_'+k)

    # #             # set up posterior lists of bs results
    # #             post_bs = dict()
    # #             for k in fit.p:
    # #                 post_bs[k] = []

    # #             for bs in range(args.Nbs):
    # #                 sys.stdout.write('%4d / %d\r' % (bs, args.Nbs))
    # #                 sys.stdout.flush()

    # #                 ''' all gvar's created in this switch are destroyed at restore_gvar [they are out of scope] '''
    # #                 gv.switch_gvar()

    # #                 bs_data = dict()
    # #                 for k in corr_bs:
    # #                     bs_data[k] = corr_bs[k][bs]
    # #                 bs_gv = gv.dataset.avg_data(bs_data)
    # #                 #import IPython; IPython.embed()
    # #                 if any(['mres' in k for k in bs_gv]):
    # #                     bs_tmp = {k:v for (k,v) in bs_gv.items() if 'mres' not in k}
    # #                     for k in [key for key in bs_gv if 'mres' in key]:
    # #                         mres = k.split('_')[0]
    # #                         if mres not in bs_tmp:
    # #                             bs_tmp[mres] = bs_gv[mres+'_MP'] / bs_gv[mres+'_PP']
    # #                     bs_gv = bs_tmp
    # #                 y_bs = {k: v[x_fit[k]['t_range']]
    # #                         for (k, v) in bs_gv.items() if k in states}
    # #                 p_bs = dict()
    # #                 for k in p_bs_mean:
    # #                     p_bs[k] = gv.gvar(p_bs_mean[k][bs], priors[k].sdev)
    # #                 # do the fit
    # #                 if has_svd:
    # #                     fit_bs = lsqfit.nonlinear_fit(data=(x_fit, y_bs), prior=p_bs, p0=p0_bs,
    # #                                                   fcn=fit_funcs.fit_function, svdcut=svdcut)
    # #                 else:
    # #                     fit_bs = lsqfit.nonlinear_fit(data=(x_fit, y_bs), prior=p_bs, p0=p0_bs,
    # #                                                   fcn=fit_funcs.fit_function)

    # #                 for r in post_bs:
    # #                     post_bs[r].append(fit_bs.p[r].mean)

    # #                 ''' end of gvar scope used for bootstrap '''
    # #                 gv.restore_gvar()

    # #             for r in post_bs:
    # #                 post_bs[r] = np.array(post_bs[r])
    # #             if args.bs_write:
    # #                 # write the results
    # #                 with h5py.File(bs_file, 'a') as f5:
    # #                     try:
    # #                         f5.create_group(args.bs_path)
    # #                     except Exception as e:
    # #                         print(e)
    # #                     for r in post_bs:
    # #                         if len(post_bs[r]) > 0:
    # #                             if r in f5[args.bs_path]:
    # #                                 del f5[args.bs_path+'/'+r]
    # #                             f5.create_dataset(
    # #                                 args.bs_path+'/'+r, data=post_bs[r])

    # #             print('DONE')

    # #     if args.svd_test:
    # #         fig = plt.figure('svd_diagnosis', figsize=(7, 4))
    # #         svd_test.plot_ratio(show=True)
    # # if args.interact:
    # #     import IPython; IPython.embed()

    # # plt.ioff()
    # # plt.show()

    # # if args.fit and args.gui:
    # #     from lsqfitgui import run_server
    # #     run_server(fit, name="c51 Two-Point Fitter")



if __name__ == "__main__":
    main()
