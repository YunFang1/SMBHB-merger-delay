import numpy as np
import math
import emcee
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import pylab as pl
import h5py
import corner

# set truth value
alpha_true = -1/16
Beta_true = 0.8
sigma_true = 0.8

# read postreiors from MCMC run with averaged likelihood
reader = emcee.backends.HDFBackend("data_hyperparam.h5")
flat_samples_A = reader.get_chain(discard=1000, flat=True, thin=5)

# read postreiors from MCMC run using likelihood assuming best fits
Flat_samples1 = np.loadtxt('array_alpha.txt')
Flat_samples2 = np.loadtxt('array_beta.txt')
Flat_samples3 =np.loadtxt('array_sigma.txt')
flat_samples_B = np.transpose(np.array([Flat_samples1, Flat_samples2, Flat_samples3]))

### plot posteriors
labels = ["\N{GREEK SMALL LETTER ALPHA}","\N{GREEK SMALL LETTER BETA}","\N{GREEK SMALL LETTER SIGMA}"]
limits = [(-2.0, 1.2), (0.0,1.6), (0.0,2.0)]


fig = corner.corner(
    flat_samples_A, labels=labels, truths=[alpha_true, Beta_true, sigma_true], truth_color = 'black',
    color='deeppink', fill_contours=False, range=limits, quiet=True, plot_datapoints=False,
    weights=np.ones(len(flat_samples_A))/len(flat_samples_A), smooth=True,
    quantiles=(0.16, 0.84), levels=(0.68, 0.954, ), label_kwargs={"fontsize": 16})

corner.corner(flat_samples_B, fig=fig, labels=labels, color='blue', fill_contours=False,
              range=limits, quiet=True, plot_datapoints=False,
              weights=np.ones(len(flat_samples_B))/len(flat_samples_B), smooth=True,
              quantiles=(0.16, 0.84), levels=(0.68,0.954, ), label_kwargs={"fontsize": 16})

for ax in fig.get_axes():
    ax.tick_params(axis='both', labelsize=12)
plt.tight_layout()
plt.savefig('delay_parameters.eps')



