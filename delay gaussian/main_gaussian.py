import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import math
import emcee
import dynesty
import h5py
#import matplotlib
import matplotlib.pyplot as plt
import pylab as pl
import gzip
import time
import csv
import corner
import random
from schwimmbad import MPIPool

import bhmergerrate
import bhmergernumber

#load samples of uncertainty parameters in stellar mass function
with open('../samples_of_parameters/stellar_mass_function/Mstar.csv',
          newline='') as f:
    read = csv.reader(f)
    Mstar = list(read)
    Mstar = [[float(y) for y in x] for x in Mstar]
    Mstar = np.array(Mstar)

with open('../samples_of_parameters/stellar_mass_function/Phi1star.csv',
          newline='') as f:
    read = csv.reader(f)
    phi1star = list(read)
    phi1star = [[float(y) for y in x] for x in phi1star]
    phi1star = np.array(phi1star)

with open('../samples_of_parameters/stellar_mass_function/Phi2star.csv',
          newline='') as f:
    read = csv.reader(f)
    phi2star = list(read)
    phi2star = [[float(y) for y in x] for x in phi2star]
    phi2star = np.array(phi2star)

with open('../samples_of_parameters/stellar_mass_function/alpha1.csv',
          newline='') as f:
    read = csv.reader(f)
    alpha1_s = list(read)
    alpha1_s = [[float(y) for y in x] for x in alpha1_s]
    alpha1_s = np.array(alpha1_s)

with open('../samples_of_parameters/stellar_mass_function/alpha2.csv',
          newline='') as f:
    read = csv.reader(f)
    alpha2_s = list(read)
    alpha2_s = [[float(y) for y in x] for x in alpha2_s]
    alpha2_s = np.array(alpha2_s)


#load samples of uncertainty parameters in galaxy merger rate per galaxy
with open('../samples_of_parameters/galaxy_merger_per_galaxy/M0.csv',
          newline='') as f:
    read = csv.reader(f)
    M0 = list(read)
    M0 = [[float(y) for y in x] for x in M0]
    M0 = np.array(M0).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/A0.csv',
          newline='') as f:
    read = csv.reader(f)
    A0 = list(read)
    A0 = [[float(y) for y in x] for x in A0]
    A0 = np.array(A0).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/eta.csv',
          newline='') as f:
    read = csv.reader(f)
    eta = list(read)
    eta = [[float(y) for y in x] for x in eta]
    eta = np.array(eta).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/alpha0.csv',
          newline='') as f:
    read = csv.reader(f)
    alpha0_g = list(read)
    alpha0_g = [[float(y) for y in x] for x in alpha0_g]
    alpha0_g = np.array(alpha0_g).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/alpha1.csv',
          newline='') as f:
    read = csv.reader(f)
    alpha1_g = list(read)
    alpha1_g = [[float(y) for y in x] for x in alpha1_g]
    alpha1_g = np.array(alpha1_g).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/beta0.csv',
          newline='') as f:
    read = csv.reader(f)
    beta0_g = list(read)
    beta0_g = [[float(y) for y in x] for x in beta0_g]
    beta0_g = np.array(beta0_g).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/beta1.csv',
          newline='') as f:
    read = csv.reader(f)
    beta1_g = list(read)
    beta1_g = [[float(y) for y in x] for x in beta1_g]
    beta1_g = np.array(beta1_g).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/gamma.csv',
          newline='') as f:
    read = csv.reader(f)
    gamma_g = list(read)
    gamma_g = [[float(y) for y in x] for x in gamma_g]
    gamma_g = np.array(gamma_g).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/delta0.csv',
          newline='') as f:
    read = csv.reader(f)
    delta0_g = list(read)
    delta0_g = [[float(y) for y in x] for x in delta0_g]
    delta0_g = np.array(delta0_g).ravel()

with open('../samples_of_parameters/galaxy_merger_per_galaxy/delta1.csv',
          newline='') as f:
    read = csv.reader(f)
    delta1_g = list(read)
    delta1_g = [[float(y) for y in x] for x in delta1_g]
    delta1_g = np.array(delta1_g).ravel()


##########################################

# set truth value
alpha_true = -1/16
Beta_true = 0.8
sigma_true = 0.8

# load posterior data of GW events
data_list = h5py.File('../GW_data_Gaussian.h5', 'r')

data1 = data_list.get('z')
data_z = np.array(data1)

data2 = data_list.get('Log10Mbh')
data_mass = np.array(data2)

data={'mass':data_mass,'redshift':data_z}
for key in data:
    data[key] = np.array(data[key])

event_num = int(len(data['mass']))
event_data = []
for ii in range (event_num):
    eventii_data = np.transpose(  np.array([ data['mass'][ii], data['redshift'][ii] ]) )
    event_data.append(list(eventii_data))


def log_likelihood(theta, params_galaxy, params_stellarmassfun):
    alpha, Beta, sigma = theta
    samples_posterior = 10 ###
    prob_ii = []

    MergerNumber = bhmergernumber.mergernumber(theta, params_galaxy, params_stellarmassfun, 100, 5000)
    #print(MergerNumber*4)
    
    if (abs(MergerNumber*4 - event_num) < 40):

        for data_event_ii in event_data: 
            
            samp_data_event_ii = random.sample(data_event_ii, samples_posterior)

            prob_sum = 0
            for jj in range(0, samples_posterior):
                R_value = bhmergerrate.mergerrate(samp_data_event_ii[jj][0], samp_data_event_ii[jj][1], theta, params_galaxy, params_stellarmassfun, 100)
                prob_sum = prob_sum + R_value
            prob_ii.append( prob_sum )
        
        ln_bayes_factors_per_event =  -np.log( samples_posterior ) + np.log( prob_ii )
        ln_l = np.sum( ln_bayes_factors_per_event  ) - event_num * np.log( MergerNumber ) + np.log( (MergerNumber*4)** event_num * np.exp(- MergerNumber*4) )
        
    elif (abs(MergerNumber*4 - event_num) >= 40):
        ln_l = - np.inf
    return ln_l


# averaged Likelihood L({d}|Lambda) over params_galaxy and params_stellarmassfun
def log_averaged_likelihood(theta):
    log_l_ii = [] 
    for ii in range(10):  ###
        #np.random.seed(123)
        n = np.random.randint(0,10000,size=15)
        params_phi = np.array([ Mstar[:,n[0]], phi1star[:,n[1]], phi2star[:,n[2]], alpha1_s[:,n[3]], alpha2_s[:,n[4]] ])
        params_Rg  = np.array([ M0[n[5]], A0[n[6]], eta[n[7]], alpha0_g[n[8]], alpha1_g[n[9]], beta0_g[n[10]], beta1_g[n[11]], gamma_g[n[12]], delta0_g[n[13]], delta1_g[n[14]] ])
##        params_phi = bestfit_stellarmass
##        params_Rg  = bestfit_galaxymerger
        
        log_l = log_likelihood(theta, params_Rg, params_phi)
        if (np.isnan(log_l)!=True):
            log_l_ii.append( log_l )
    likelihood_value = np.exp(log_l_ii)
    average = np.sum(likelihood_value)/len(likelihood_value)
    return np.log(average)

##start = time.time()
##result =log_averaged_likelihood( np.array([-0.25,0.8,0.8]) )
##end = time.time()
##serial_time = end - start
##print("interpolation took {0:.4f} seconds".format(serial_time))
##print("the result is")
##print(result)


############################################
### prior of hyperparameter Lambda
def log_prior(theta):
    alpha, Beta, sigma = theta
    if -3 < alpha < 3 and  0.001 < Beta < 3 and 0.001 < sigma < 4:
        return 0.0
    return -np.inf

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_averaged_likelihood(theta)



############### parallel jobs and sample with emcee
with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    # Initialize the walkers
    np.random.seed(42)
    initial = np.array([alpha_true, Beta_true, sigma_true]) + 0.05 * np.random.randn(40, 3)*np.array([alpha_true, Beta_true, sigma_true])
    pos = initial
    nwalkers, ndim = pos.shape
    nsteps = 10000
    # save data
    filename = "data_hyperparam.h5"
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
    # Initialize the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool, backend=backend)
    sampler.run_mcmc(pos, nsteps, progress=True)


fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["\N{GREEK SMALL LETTER ALPHA}","\N{GREEK SMALL LETTER BETA}","\N{GREEK SMALL LETTER SIGMA}"]

for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
plt.savefig('para_steps.eps')


# plot posteriors
flat_samples = sampler.get_chain(discard=2000, flat=True, thin=5)

fig = corner.corner(
    flat_samples, labels=labels, truths=[alpha_true, Beta_true, sigma_true], smooth=True, quiet=True, quantiles=(0.16, 0.84), levels=(0.68, 0.954, ), plot_datapoints=False
)
plt.tight_layout()
plt.savefig('delay_parameters.eps')
