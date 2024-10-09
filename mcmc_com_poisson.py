import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from mpmath import mpf, log, exp, loggamma, fsum, mp
from src.bounding_pairs_mp import bounding_pairs_mp as bpm
import math

def log_Z_term(theta: tuple, k: int):
    if k == 1:
        return mpf(0)
    elif (k == 2):
        return theta[0]
    else:
        return theta[1] * ((mpf(k)-1) * theta[0] - loggamma(mpf(k)))

# Log-likelihood function: Poisson distribution for count data with frequencies
def log_likelihood_freq(values, freqs, theta1, theta2):
    """
    Log-likelihood function for COM-Poisson distribution with frequencies.
    values: unique observed values
    freqs: corresponding frequencies for each unique value
    """
    theta1 = mpf(theta1)
    theta2 = mpf(theta2)
    
    mp.dps = 400
    terms = [0] * len(values)

    log_Z = bpm(log_Z_term, (log(theta1), theta2), 10**6, mpf(0), mpf(2)**mpf(-52), initial_k=1)[1]
    for i in range(len(values)):
        terms[i] = freqs[i] * (theta2 * (values[i] * theta1 - loggamma(values[i] + 1)) - log_Z)

    return fsum(terms)

# Prior for theta1: Gamma distribution (log)
def log_gamma_prior(theta1, alpha=0.01, beta=0.01):
    return stats.gamma.logpdf(theta1, alpha, scale=1/beta)

# Prior for theta2: Normal distribution with mean 0 (log)
def log_normal_prior(theta2, mean=0, sigma=1.0):
    return stats.norm.logpdf(theta2, mean, sigma)

# Log-posterior distribution
def log_posterior_freq(values, freqs, theta1, theta2, alpha=0.01, beta=0.01, sigma=1.0):
    """
    Log-posterior is the sum of log-likelihood + log-prior(theta1) + log-prior(theta2).
    """
    log_likelihood_value = log_likelihood_freq(values, freqs, theta1, theta2)
    log_prior_theta1 = log_gamma_prior(theta1, alpha, beta)
    log_prior_theta2 = log_normal_prior(theta2, mean=0, sigma=sigma)
    
    return float(log_likelihood_value) + log_prior_theta1 + log_prior_theta2

# MCMC with Metropolis-Hastings algorithm (log version)
def mcmc_log_freq(values, freqs, iterations=1000, theta1_init=0.1, theta2_init=0.1, alpha=0.1, beta=0.1, sigma=1.0):
    theta1_samples = [theta1_init]
    theta2_samples = [theta2_init]
    
    for _ in range(iterations):
        # Propose new theta1 and theta2
        theta1_proposed = np.random.gamma(alpha, beta)
        theta2_proposed = np.abs(np.random.normal(0, sigma))

        # Compute log-posterior for current and proposed
        log_posterior_current = log_posterior_freq(values, freqs, theta1_samples[-1], theta2_samples[-1], alpha, beta, sigma)
        log_posterior_proposed = log_posterior_freq(values, freqs, theta1_proposed, theta2_proposed, alpha, beta, sigma)

        # Acceptance probability (log version)
        log_acceptance_prob = log_posterior_proposed - log_posterior_current
        acceptance_prob = exp(log_acceptance_prob)  # Convert back to probability space

        # Accept or reject the proposal
        if np.random.rand() < acceptance_prob:
            theta1_samples.append(theta1_proposed)
            theta2_samples.append(theta2_proposed)
        else:
            theta1_samples.append(theta1_samples[-1])
            theta2_samples.append(theta2_samples[-1])
    
    return np.array(theta1_samples), np.array(theta2_samples)

# Simulating some value-frequency data (unique values with frequencies)
np.random.seed(42)
values = [i for i in range(31)]  # Unique values (example counts)
freqs = [514, 503, 457, 423, 326, 233, 195, 139, 101, 77, 56, 40, 37, 22, 9, 7, 10, 9, 3, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1]  # Frequencies of those values

# Running the MCMC with log-likelihood
theta1_samples, theta2_samples = mcmc_log_freq(values, freqs)

# Posterior mean estimates for theta1 and theta2
theta1_mean = np.mean(theta1_samples)
theta2_mean = np.mean(theta2_samples)
log_Z = bpm(log_Z_term, (log(mpf(theta1_mean)), theta2_mean), 10**6, mpf(0), mpf(2)**mpf(-52), initial_k=1)[1]

# Generate Poisson distribution based on fitted lambda
x_values = np.arange(0, max(values)+1)
y_values = np.array([exp((theta2_mean * (x * theta1_mean - loggamma(x + 1)) - log_Z)) for x in x_values])

# Plot the data and the fitted distribution
plt.figure(figsize=(12, 6))

# Bar plot of observed data (value-frequency pairs)
plt.bar(values, freqs/np.sum(freqs), width=0.6, color='g', alpha=0.6, label='Observed Data (Freqs)')

# Fitted Poisson distribution
plt.plot(x_values, y_values, 'r-', lw=2, label=f'Fitted Poisson\n mu = {theta1_mean:.6f},\n nu = {theta2_mean:.6f}')

# Titles and labels
plt.title('Observed Data (Frequencies) and Fitted COM-Poisson Distribution (Log-Likelihood)', fontsize=14)
plt.xlabel('Value', fontsize=12)
plt.ylabel('Probability Density', fontsize=12)
plt.legend()

plt.show()
