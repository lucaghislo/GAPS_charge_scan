import numpy as np
from scipy import special


def erf_function(x, mu, sigma):
    return (1 / 2) * (1 + special.erf((x - mu) / (sigma * np.sqrt(2))))
