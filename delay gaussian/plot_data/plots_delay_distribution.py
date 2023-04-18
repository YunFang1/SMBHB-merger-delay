import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import math
import emcee
import h5py
import matplotlib
import matplotlib.pyplot as plt
from scipy import special


#set initial data
alpha_true = -1/16
beta_true = 0.8
sigma_true = 0.8

# read postreiors from MCMC run with averaged likelihood
reader = emcee.backends.HDFBackend("data_hyperparam.h5")
flat_samples_A = reader.get_chain(discard=1000, flat=True, thin=3)

# read postreiors from MCMC run using likelihood assuming best fits
Flat_samples1 = np.loadtxt('array_alpha.txt')
Flat_samples2 = np.loadtxt('array_beta.txt')
Flat_samples3 =np.loadtxt('array_sigma.txt')

flat_samples_B = np.transpose(np.array([Flat_samples1[40000:(len(Flat_samples1)-1)],
                                        Flat_samples2[40000:(len(Flat_samples1)-1)],
                                        Flat_samples3[40000:(len(Flat_samples1)-1)]] ))




################# plot P(t) #######################

alpha_post_A = flat_samples_A[:,0]
beta_post_A = flat_samples_A[:,1]
sigma_post_A = flat_samples_A[:,2]

alpha_post_B = flat_samples_B[:,0]
beta_post_B = flat_samples_B[:,1]
sigma_post_B = flat_samples_B[:,2]

def P_delay(Log10Mbh, t, alpha, Beta, sigma):
    mu = Beta * 10**( alpha * (Log10Mbh - 6) )
    Pdelay = np.exp(- 0.5 * (t - mu)**2 / sigma**2 ) / ( (2 *np.pi)**0.5 * sigma ) * 2/( 1 + special.erf( mu / ( np.sqrt(2) * sigma ) ) )
    return Pdelay

t=np.linspace(0,6,100)
logmass=6.3

prob_distri_A=[]
for ii in range(len(t)):
    prob_distri_A.append( P_delay(logmass, t[ii], alpha_post_A, beta_post_A, sigma_post_A) )

prob_distri_B=[]
for ii in range(len(t)):
    prob_distri_B.append( P_delay(logmass, t[ii], alpha_post_B, beta_post_B, sigma_post_B) )



colours = {}
colours['distribute'] = '#1f77b4'
colours['true'] = '#ff7f0e'


# plot it!
fig, ax = plt.subplots(1)

prob_distri = prob_distri_A
p16_A, p50_A, p84_A = {}, {}, {}
p16_A = np.percentile(prob_distri,16,axis=1)
p50_A = np.percentile(prob_distri,50,axis=1)
p84_A = np.percentile(prob_distri,84,axis=1)

prob_distri = prob_distri_B
p16_B, p50_B, p84_B = {}, {}, {}
p16_B = np.percentile(prob_distri,16,axis=1)
p50_B = np.percentile(prob_distri,50,axis=1)
p84_B = np.percentile(prob_distri,84,axis=1)


my_plot_A = ax.fill_between(t, p16_A, p84_A, facecolor='deeppink', alpha=0.8, interpolate=True)
my_plot_B = ax.fill_between(t, p16_B, p84_B, facecolor='blue', alpha=0.65, interpolate=True)
ax.plot(t, P_delay(logmass, t, alpha_true, beta_true, sigma_true), color = 'orange', alpha = 1, linewidth = 2.0)


ax.set_xlabel('t [$\mathrm{Gyr}$] ', fontsize=15)
ax.set_ylabel('$P_{\mathrm{delay}}$',fontsize=15)
plt.xlim(0, 3)
plt.ylim(1e-1, 1)
plt.tick_params(axis='both', labelsize=13)
plt.title('$M= 2 * 10^6 \mathrm{M}_\odot$', fontsize=15)
#plt.tight_layout()
plt.show()
#plt.savefig('confidence_delay_distribution.eps')


