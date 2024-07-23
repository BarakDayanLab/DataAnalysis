import numpy as np
import scipy.constants as const
from scipy.stats import poisson

## Calculation of acceleration due to scattering ##

Number_of_photon_scattered_per_second = 3000
Lambda = 780e-9
k = 2 * const.pi / Lambda
delta_p = const.hbar * k * 3000
m_rb87 = 1.5e-25
acceleration = delta_p / m_rb87

## Calculation for PNSA, QBER ##

mu_COW = 2
x = 1 - poisson.cdf(k=0, mu=mu_COW)


## Calculation of Magnetic Field / magnetic dipole moment ##

mu_0 = 4 * const.pi * 1e11  # vacuum permeability [G * m / A] - gauss metre per ampere
B_0 =
