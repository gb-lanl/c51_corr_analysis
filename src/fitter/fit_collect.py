import lsqfit
import numpy as np
import gvar as gv

class fitter(object):

    def __init__(self, n_states, t_range,
                pt2_nstates, pt3_nstates, sum_nstates,
                nucleon_corr_gv=None, axial_num_gv=None, vector_num_gv=None):

        self.n_states = n_states
        self.t_range = t_range
        # self.prior = prior
        self.nucleon_corr_gv = nucleon_corr_gv
        self.axial_num_gv = axial_num_gv
        self.vector_num_gv = vector_num_gv
        self.fit = None
        self.prior = None
        self.pt2_nstates = pt2_nstates
        self.pt3_nstates = pt3_nstates
        self.sum_nstates = sum_nstates


    def get_fit(self):
        if self.fit is not None:
            return self.fit
        else:
            return self._make_fit()
    def get_prior(self):
        if self.prior is not None:
            return self.prior
        else:
            return self._make_prior()

    def get_data(self):
        return self._make_data()
    


    def get_energies(self):

        # Don't rerun the fit if it's already been made
        if self.fit is not None:
            temp_fit = self.fit
        else:
            temp_fit = self.get_fit()

        max_n_states = np.max([self.n_states[key] for key in list(self.n_states.keys())])
        output = gv.gvar(np.zeros(max_n_states))
        output[0] = temp_fit.p['E0']
        for k in range(1, max_n_states):
            output[k] = output[0] + np.sum([(temp_fit.p['dE'][j]) for j in range(k)], axis=0)
        return output

    def _make_fit(self):
        # LOGIC FOR SIMULTANEOUS FITS
        # This is the convoluted way we use MultiFitter
        # Essentially: first we create a model (which is a subclass of MultiFitter)
        # Then we make a fitter using the models
        # Finally, we make the fit with our two sets of correlators

        models = self._make_models_simult_fit()

        fitter = lsqfit.MultiFitter(models=models)
        fit = fitter.lsqfit(data=self._make_data(), prior=self._make_prior())
        self.fit = fit
        return fit

    def _make_models_simult_fit(self):
        models = np.array([])

        # if self.nucleon_corr_gv is not None:
        #     for sink in list(self.nucleon_corr_gv.keys()):
        #         param_keys = {
        #             'E0'      : 'E0',
        #             # 'E1'      : 'E1',
        #             # 'E2'      : 'E2',
        #             # 'E3'      : 'E3',
        #             'log(dE)' : 'log(dE)'+sink,
        #             'z0'      : 'z_'+sink,
        #         }
        #         models = np.append(models,
        #                 baryon_model(datatag="nucleon_"+sink,
        #                 t=list(range(self.t_range['corr'][0], self.t_range['corr'][1])),
        #                 param_keys=param_keys, n_states=self.n_states['corr']))

        if self.axial_num_gv is not None:
            for keys in list(self.axial_num_gv.keys()):
                param_keys = {
                    'E0'      : 'E0',
                    'log(dE)' : 'log(dE)',
                    'd'       : 'd_A_'+keys,
                    'A3_nm'    : 'A3_nm',
                    'z'      : 'z_'+keys,
                }
                models = np.append(models,
                        fh_num_model(datatag="axial_num_"+keys,
                        t=list(range(self.t_range['gA'][0], self.t_range['gA'][1])),
                        param_keys=param_keys, n_states=self.n_states['gA']))

        if self.vector_num_gv is not None:
            for keys in list(self.vector_num_gv.keys()):
                
                param_keys = {
                    'E0'      : 'E0',
                    'log(dE)' : 'log(dE)',
                    'd'       : 'd_V_'+keys,
                    'V4_nm'    : 'V4_nm',
                    'z'      : 'z_'+keys,
                }
                models = np.append(models,
                        fh_num_model(datatag="vector_num_"+keys,
                        t=list(range(self.t_range['gV'][0], self.t_range['gV'][1])),
                        param_keys=param_keys, n_states=self.n_states['gV']))

        return models

    # data array needs to match size of t array
    def _make_data(self):
        data = {}
        if self.nucleon_corr_gv is not None:
            for sink in list(self.nucleon_corr_gv.keys()):
                data["nucleon_"+sink] = self.nucleon_corr_gv[sink][self.t_range['corr'][0]:self.t_range['corr'][1]]

        if self.axial_num_gv is not None:
            for ens in list(self.axial_num_gv.keys()):
                data["axial_num_"+ens] = self.axial_num_gv[ens][self.t_range['gA'][0]:self.t_range['gA'][1]]

        if self.vector_num_gv is not None:
            for ens in list(self.vector_num_gv.keys()):
                data["vector_num_"+ens]  = self.vector_num_gv[ens][self.t_range['gV'][0]:self.t_range['gV'][1]]

        return data

    def _make_prior(self):

        prior =  {} 

        prior['E0'] = gv.gvar(0.50, 0.05)
        prior['z0'] = gv.gvar(0.00034, 0.00034)

        prior['log(dE1)'] = gv.gvar(-1.25, 0.5)
        prior['z1'] = gv.gvar(0, 0.00025)

        prior['log(dE2)'] = gv.gvar(-1.25, 0.5)
        prior['z2'] = gv.gvar(0, 0.00015)

        max_nstate = max(np.array([self.pt2_nstates, self.pt3_nstates, self.sum_nstates]))

        if max_nstate > 3:
            for i in range(3, max_nstate):
                prior['log(dE'+str(i)+')'] = gv.gvar(-1.25, 0.5)
                prior['z'+str(i)] = gv.gvar(0, 0.00015)

        prior['log(dE'+str(max_nstate-1)+')'] = gv.gvar(-1.25, 0.5*5) # garbage can

        prior['log(E_sum)'] = gv.gvar(-1.25, 0.5*5) 
        prior['z_sum'] = gv.gvar(0, 0.00015) 

        prior['A3_00'] = gv.gvar(1.2, 0.2)
        prior['V4_00'] = gv.gvar(1.0, 0.2)

        ff_nstates = max(np.array([self.pt3_nstates, self.sum_nstates]))

        for i in range(ff_nstates):
            for j in range(ff_nstates):
                if i+j >= 1:
                    if j < i:  
                        prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                        prior['V4_'+str(j)+str(i)] = gv.gvar(0, 1)

                    elif j == i:
                        prior['A3_'+str(j)+str(i)] = gv.gvar(0, 1)
                        prior['V4_'+str(j)+str(i)] = gv.gvar(1, 0.2)

        for i in range(self.sum_nstates-1):
            prior['sum_A3_'+str(i)] = gv.gvar(0, 1)
            prior['sum_V4_'+str(i)] = gv.gvar(0, 1)

        prior['sum_A3_'+str(self.sum_nstates-1)] = gv.gvar(0, 1)
        prior['sum_V4_'+str(self.sum_nstates-1)] = gv.gvar(1, 0.2)

        return prior
        # resized_prior = {}

        # max_n_states = np.max([self.n_states[key] for key in list(self.n_states.keys())])
        # for key in list(prior.keys()):
        #     if key == 'g_A_nm':
        #         resized_prior[key] = prior[key][:self.n_states['gA'], :self.n_states['gA']]
        #     elif key == 'g_V_nm':
        #         resized_prior[key] = prior[key][:self.n_states['gV'], :self.n_states['gV']]
        #     elif key in ['d_A_dir', 'd_A_smr']:
        #         resized_prior[key] = prior[key][:self.n_states['gA']]
        #     elif key in ['d_V_dir', 'd_V_smr']:
        #         resized_prior[key] = prior[key][:self.n_states['gV']]
        #     else:
        #         resized_prior[key] = prior[key][:max_n_states]

        # new_prior = resized_prior.copy()
        # new_prior['E0'] = resized_prior['E'][0]
        # # Don't need this entry
        # new_prior.pop('E', None)

        # # We force the energy to be positive by using the log-normal dist of dE
        # # let log(dE) ~ eta; then dE ~ e^eta
        # new_prior['log(dE)'] = gv.gvar(np.zeros(len(resized_prior['E']) - 1))
        # for j in range(len(new_prior['log(dE)'])):
        #     #excited_state_energy = p[self.mass] + np.sum([np.exp(p[self.log_dE][k]) for k in range(j-1)], axis=0)

        #     # Notice that I've coded this s.t.
        #     # the std is determined entirely by the excited state
        #     # dE_mean = gv.mean(resized_prior['E'][j+1] - resized_prior['E'][j])
        #     # dE_std = gv.sdev(resized_prior['E'][j+1])
        #     temp = gv.gvar(resized_prior['E'][j+1]) - gv.gvar(resized_prior['E'][j])
        #     temp2 = gv.gvar(resized_prior['E'][j+1])
        #     temp_gvar = gv.gvar(temp.mean,temp2.sdev)
        #     new_prior['log(dE)'][j] = np.log(temp_gvar)

        # return new_prior



# This class is needed to instantiate an object for lsqfit.MultiFitter
# There's probably a better way to do it, but the documentation is lacking
# Used for particles that obey fermi-dirac statistics
class baryon_model(lsqfit.MultiFitterModel):
    def __init__(self, datatag, t, param_keys, n_states):
        super(baryon_model, self).__init__(datatag)
        # variables for fit
        self.t = np.array(t)
        self.n_states = n_states

        # keys (strings) used to find the wf_overlap and energy in a parameter dictionary
        self.param_keys = param_keys

    def fitfcn(self, p, t=None):
        if t is None:
            t = self.t

        z = p[self.param_keys['z']]
        E0 = p[self.param_keys['E0']]
        log_dE = p[self.param_keys['log(dE)']]

        output = z[0] * np.exp(-E0 * t)
        for j in range(1, self.n_states):
            excited_state_energy = E0 + np.sum([np.exp(log_dE[k]) for k in range(j)], axis=0)
            output = output + z[j] * np.exp(-excited_state_energy * t)
        return output

    # The prior determines the variables that will be fit by multifitter --
    # each entry in the prior returned by this function will be fitted
    def buildprior(self, prior, mopt=None, extend=False):
        # Extract the model's parameters from prior.
        return prior

    def builddata(self, data):
        # Extract the model's fit data from data.
        # Key of data must match model's datatag!
        return data[self.datatag]

    def fcn_effective_mass(self, p, t=None):
        if t is None:
            t=self.t

        return np.log(self.fitfcn(p, t) / self.fitfcn(p, t+1))

# This class is needed to instantiate an object for lsqfit.MultiFitter
# There's probably a better way to do it, but the documentation is lacking
# Used for particles that obey fermi-dirac statistics
class fh_num_model(lsqfit.MultiFitterModel):
    def __init__(self, datatag, t, param_keys, n_states):
        super(fh_num_model, self).__init__(datatag)
        # variables for fit
        self.t = np.array(t)
        self.n_states = n_states

        # keys (strings) used to find the wf_overlap and energy in a parameter dictionary
        self.param_keys = param_keys

    def fitfcn(self, p, t=None):
        if t is None:
            t = self.t


        z = p[self.param_keys['z']]
        E0 = p[self.param_keys['E0']]
        log_dE = np.append(-np.inf, p[self.param_keys['log(dE)']])
        d = p[self.param_keys['d']]
        g_nm = p[self.param_keys['g_nm']]

        E = np.array([np.sum([np.exp(log_dE[j]) for j in range(n+1)]) for n in range(self.n_states)]) + E0

        output = 0
        for n in range(self.n_states):
            for m in range(self.n_states):
                if n == m:
                    #if m > n: g_nm[n, m] = g_nm[m, n]
                    output += ((t-1)*z[n]*g_nm[n, m] + d[n]) * np.exp(-E[n] * t)
                else:
                    pass
                    E_n = E[n]
                    E_m = E[m]
                    dE_nm = E_n - E_m
                    dE_mn = -dE_nm

                    output += (z[n]*g_nm[n, m]) * ((np.exp(dE_nm/2 - E_n*t) - np.exp(dE_mn/2 - E_m*t)) /
                                                    (np.exp(dE_mn/2) - np.exp(dE_nm/2)))
        return output

    # The prior determines the variables that will be fit by multifitter --
    # each entry in the prior returned by this function will be fitted
    def buildprior(self, prior, mopt=None, extend=False):
        # Extract the model's parameters from prior.

        return prior

    def builddata(self, data):
        # Extract the model's fit data from data.
        # Key of data must match model's datatag!
        return data[self.datatag]
