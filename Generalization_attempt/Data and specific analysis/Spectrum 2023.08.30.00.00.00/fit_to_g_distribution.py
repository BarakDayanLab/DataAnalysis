import numpy as np
from scipy.stats import poisson
from scipy.optimize import curve_fit
from Archive.Helpers.Helpers import cartesian_product, fit_zero
from matplotlib import pyplot as plt
import pandas as pd

r = np.linspace(0, 1500, 100)
y = np.linspace(-600, 600, 100)
ry = cartesian_product(r, y)
p_in_rect = np.array(0.5, )
spatial_distribution = 1 / (r.shape[0] * y.shape[0]) * np.ones((r.shape[0] * y.shape[0]))
spatial_distribution = spatial_distribution / np.sum(spatial_distribution)
k_extrinsic = np.array(8, )
k_intrinsic = np.array(7, )
gamma = np.array(6.065 / 2, )
h = np.array(2, )
alpha_squared = np.array(1, )
num_photons_for_detection = 1  # inclusive
g0 = np.array(40.0, )
sigma_y = np.array(3000.0, )
c_decay_r = np.array(0.1 * (2 * np.pi) / 780.0, )


def transmission_differentiable(f, g1, g2, gamma, k_ex, k_i, h):
    k = k_i + k_ex
    C = (g1 ** 2 + g2 ** 2) / (2 * (k + 1j * f) * (gamma + 1j * f))
    g1g1 = np.power(np.abs(g1), 2) / (np.power(np.abs(g1), 2) + np.power(np.abs(g2), 2))
    g1g2 = g1 * g2 / (np.power(np.abs(g1), 2) + np.power(np.abs(g2), 2))

    t_0 = 1 + 2 * 1j * k_ex * (f - 1j * (k_i + k_ex)) / \
          (np.power((f - 1j * (k_i + k_ex)), 2) - np.power(h, 2))
    r_0 = 2 * k_ex * h / (np.power((1j * f + (k_i + k_ex)), 2) + np.power(h, 2))

    z = np.power(np.abs((2 * (k + 1j * f) * k_ex / (np.power((k + 1j * f), 2) + np.power(h, 2))) * g1g1 *
                        (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * f), 2) + 2 * C)) + t_0), 2) + \
        np.power(np.abs(r_0 * (2 * k_ex / (k + 1j * f)) * g1g2 *
                        (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * f), 2) + 2 * C))), 2)
    return z


def transmission_g0(f, k_ex, k_i, h):
    k = k_i + k_ex
    C = 0
    g1g1 = 0
    g1g2 = 0

    t_0 = 1 + 2 * 1j * k_ex * (f - 1j * (k_i + k_ex)) / \
          (np.power((f - 1j * (k_i + k_ex)), 2) - np.power(h, 2))
    r_0 = 2 * k_ex * h / (np.power((1j * f + (k_i + k_ex)), 2) + np.power(h, 2))

    z = np.power(np.abs((2 * (k + 1j * f) * k_ex / (np.power((k + 1j * f), 2) + np.power(h, 2))) * g1g1 *
                        (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * f), 2) + 2 * C)) + t_0), 2) + \
        np.power(np.abs(r_0 * (2 * k_ex / (k + 1j * f)) * g1g2 *
                        (2 * C / (1 + np.power(h, 2) / np.power((k + 1j * f), 2) + 2 * C))), 2)
    return z


def differentiable_g_profile(ry, g0, sigma_y, c_decay_r):
    r, y = ry[:, 0], ry[:, 1]
    ret = g0 * np.exp(-c_decay_r * r) * np.exp(- y ** 2 / (2 * sigma_y ** 2))
    return ret


def fit_func(f, p_in_rect, g0, sigma_y, c_decay_r, alpha_squared, k_extrinsic, k_intrinsic, gamma, h):
    spatial_distribution_ = p_in_rect * spatial_distribution
    spatial_distribution_ = np.concatenate((spatial_distribution_, [1 - p_in_rect]), axis=0)

    g = differentiable_g_profile(ry, g0, sigma_y=sigma_y, c_decay_r=c_decay_r)
    # g = np.concatenate((g, np.array([0])), dim=0)

    # RESHAPES
    f = np.array(f).reshape((1, -1))
    g = g.reshape((-1, 1))
    spatial_distribution_ = spatial_distribution_.reshape((-1, 1))

    T = transmission_differentiable(f=f, g1=g, g2=0, gamma=gamma, k_ex=k_extrinsic, k_i=k_intrinsic, h=h)
    T_f0 = transmission_differentiable(f=np.array(0), g1=g, g2=0, gamma=gamma, k_ex=k_extrinsic, k_i=k_intrinsic, h=h)

    g0 = 0
    T_0 = transmission_g0(f=f, k_ex=k_extrinsic, k_i=k_intrinsic, h=h)
    T_0_f0 = transmission_g0(f=0, k_ex=k_extrinsic, k_i=k_intrinsic, h=h)
    g = np.concatenate((g, np.array([[g0]])), axis=0)
    T = np.concatenate((T, T_0), axis=0)
    T_f0 = np.concatenate((T_f0, [T_0_f0]), axis=0)

    poisson_lambda = T * alpha_squared
    poisson_lambda_f0 = T_f0 * alpha_squared

    probability_of_detection_f0 = poisson.sf(num_photons_for_detection - 1, poisson_lambda_f0)
    prob_squared = probability_of_detection_f0 ** 2

    averaged_selected_poisson_lambda = \
        np.sum(poisson_lambda * prob_squared * spatial_distribution_, axis=0, keepdims=True) / \
        np.sum(prob_squared * spatial_distribution_, axis=0, keepdims=True)

    return averaged_selected_poisson_lambda.flatten()


transmission_onres = 0.01
k_extrinsic = fit_zero(
    lambda k_ex: transmission_onres - transmission_g0(f=0, k_ex=k_ex, k_i=k_intrinsic, h=h).flatten(), k_intrinsic + 1).flatten()


data = pd.read_csv('Archive/Data and specific analysis/Spectrum 2023.08.30/output.csv')
x = data['x'].to_numpy()
y = data['y'].to_numpy()
popt, pcov = curve_fit(lambda *x: fit_func(*x, k_extrinsic, k_intrinsic, gamma, h), x, y,
                       p0=[0.7, 40.0, 300.0, (2 * np.pi / 780), 1.0], bounds=([0, 0, 0, 0, 0], [1, 1000, 1000, 1000, 10]))

dense_x = np.linspace(x[0], x[-1], 500)
plt.scatter(x, y, label='data')
plt.plot(dense_x, fit_func(dense_x, *popt, k_extrinsic, k_intrinsic, gamma, h), 'r-', label='fit')
plt.legend()
plt.xlabel('Detuning Frequency (MHz)')
plt.ylabel('Counts after passing test, normalized by counts for no cavity')
plt.title('Counts after passing test vs detuning frequency')
plt.ylim([0, 1.1 * np.max(y)])
plt.show()

print(popt)
print(pcov)
