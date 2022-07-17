import lsqfit
import numpy as np
import gvar as gv

class fitter(object):

    def __init__(self, n_states, prior, t_range,
                 nucleon_corr_gv=None, axial_fh_num_gv=None, vector_fh_num_gv=None):

        self.n_states = n_states
        self.t_range = t_range
        self.prior = prior
        self.nucleon_corr_gv = nucleon_corr_gv
        self.axial_fh_num_gv = axial_fh_num_gv
        self.vector_fh_num_gv = vector_fh_num_gv
        self.fit = None
        self.prior = self._make_prior(prior)


    def get_fit(self):
        if self.fit is not None:
            return self.fit
        else:
            return self._make_fit()

    def get_energies(self):

        # Don't rerun the fit if it's already been made
        if self.fit is not None:
            temp_fit = self.fit
        else:
            temp_fit = self.get_fit()

        max_n_states = np.max([self.n_states[key] for key in self.n_states.keys()])
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
        fit = fitter.lsqfit(data=self._make_data(), prior=self.prior)
        self.fit = fit
        return fit

    def _make_models_simult_fit(self):
        models = np.array([])

        if self.nucleon_corr_gv is not None:
            for sink in self.nucleon_corr_gv.keys():
                param_keys = {
                    'E0'      : 'E0',
                    'log(dE)' : 'log(dE)',
                    'z'      : 'z_'+sink,
                }
                models = np.append(models,
                           baryon_model(datatag="nucleon_"+sink,
                           t=range(self.t_range['corr'][0], self.t_range['corr'][1]),
                           param_keys=param_keys, n_states=self.n_states['corr']))

        if self.axial_fh_num_gv is not None:
            for sink in self.axial_fh_num_gv.keys():
                param_keys = {
                    'E0'      : 'E0',
                    'log(dE)' : 'log(dE)',
                    'd'       : 'd_A_'+sink,
                    'g_nm'    : 'g_A_nm',
                    'z'      : 'z_'+sink,
                }
                models = np.append(models,
                           fh_num_model(datatag="axial_fh_num_"+sink,
                           t=range(self.t_range['gA'][0], self.t_range['gA'][1]),
                           param_keys=param_keys, n_states=self.n_states['gA']))

        if self.vector_fh_num_gv is not None:
            for sink in self.vector_fh_num_gv.keys():
                param_keys = {
                    'E0'      : 'E0',
                    'log(dE)' : 'log(dE)',
                    'd'       : 'd_V_'+sink,
                    'g_nm'    : 'g_V_nm',
                    'wf'      : 'wf_'+sink,
                }
                models = np.append(models,
                           fh_num_model(datatag="vector_fh_num_"+sink,
                           t=range(self.t_range['gV'][0], self.t_range['gV'][1]),
                           param_keys=param_keys, n_states=self.n_states['gV']))

        return models

    # data array needs to match size of t array
    def _make_data(self):
        data = {}
        if self.nucleon_corr_gv is not None:
            for sink in self.nucleon_corr_gv.keys():
                data["nucleon_"+sink] = self.nucleon_corr_gv[sink][self.t_range['corr'][0]:self.t_range['corr'][1]]

        if self.axial_fh_num_gv is not None:
            for sink in self.axial_fh_num_gv.keys():
                data["axial_fh_num_"+sink] = self.axial_fh_num_gv[sink][self.t_range['gA'][0]:self.t_range['gA'][1]]

        if self.vector_fh_num_gv is not None:
            for sink in self.vector_fh_num_gv.keys():
                data["vector_fh_num_"+sink] = self.vector_fh_num_gv[sink][self.t_range['gV'][0]:self.t_range['gV'][1]]

        return data

    def _make_prior(self, prior):
        resized_prior = {}

        max_n_states = np.max([self.n_states[key] for key in self.n_states.keys()])
        for key in prior.keys():
            if key == 'g_A_nm':
                resized_prior[key] = prior[key][:self.n_states['gA'], :self.n_states['gA']]
            elif key == 'g_V_nm':
                resized_prior[key] = prior[key][:self.n_states['gV'], :self.n_states['gV']]
            elif key in ['d_A_PS', 'd_A_SS']:
                resized_prior[key] = prior[key][:self.n_states['gA']]
            elif key in ['d_V_PS', 'd_V_SS']:
                resized_prior[key] = prior[key][:self.n_states['gV']]
            else:
                resized_prior[key] = prior[key][:max_n_states]

        new_prior = resized_prior.copy()
        new_prior['E0'] = resized_prior['E'][0]

        # Don't need this entry
        new_prior.pop('E', None)

        # We force the energy to be positive by using the log-normal dist of dE
        # let log(dE) ~ eta; then dE ~ e^eta
        new_prior['log(dE)'] = gv.gvar(np.zeros(len(resized_prior['E']) - 1))
        for j in range(len(new_prior['log(dE)'])):
            #excited_state_energy = p[self.mass] + np.sum([np.exp(p[self.log_dE][k]) for k in range(j-1)], axis=0)

            # Notice that I've coded this s.t.
            # the std is determined entirely by the excited state
            dE_mean = gv.mean(resized_prior['E'][j+1] - resized_prior['E'][j])
            dE_std = gv.sdev(resized_prior['E'][j+1])

            new_prior['log(dE)'][j] = np.log(gv.gvar(dE_mean, dE_std))

        return new_prior

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

        wf = p[self.param_keys['wf']]
        E0 = p[self.param_keys['E0']]
        log_dE = p[self.param_keys['log(dE)']]

        output = wf[0] * np.exp(-E0 * t)
        for j in range(1, self.n_states):
            excited_state_energy = E0 + np.sum([np.exp(log_dE[k]) for k in range(j)], axis=0)
            output = output + wf[j] * np.exp(-excited_state_energy * t)
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


        wf = p[self.param_keys['wf']]
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
                    output += ((t-1)*wf[n]*g_nm[n, m] + d[n]) * np.exp(-E[n] * t)
                else:
                    pass
                    E_n = E[n]
                    E_m = E[m]
                    dE_nm = E_n - E_m
                    dE_mn = -dE_nm

                    output += (wf[n]*g_nm[n, m]) * ((np.exp(dE_nm/2 - E_n*t) - np.exp(dE_mn/2 - E_m*t)) /
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
