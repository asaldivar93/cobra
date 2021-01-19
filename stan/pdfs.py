import matplotlib.pyplot as plot
import numpy as np

from scipy.special import gamma, beta
from numpy import exp


def gamma_pdf(x, a, b):
    return b**a * x**(a - 1) * exp(-b * x) / gamma(a)


def beta_pdf(x, a, b):
    return x**(a - 1) * (1 - x)**(b - 1) / beta(a, b)


x = np.linspace(0, 1, 100)
mu_0 = np.median(np.linspace(1, 17, 17))
b = float(1)
a = mu_0 * b
plot.plot(x, beta_pdf(x, 200, 200))
