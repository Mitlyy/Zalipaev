import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots
from numpy import loadtxt
from scipy.fftpack import *

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


THz = 1e+12
f0 = 1 * THz
tau_p = 0.5/ THz
z = 1e-3
d = 0.5e-3
data = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)
omega = np.linspace(5e+12, 8e+12, data[1].shape[0])
E = data[1]
t = data[0]
c = 3e+8
F = fft(E)

def T_(w_):
    dt = t[1]-t[0]
    a = []
    for i in range(t.shape[0]):
        a.append(1/F[i] * sum(np.exp(1j * w_[i] * (t - z / c)) * E) * dt)
    a = np.asarray(a)
    # plt.plot(omega, a.real)
    # plt.plot(omega, a.imag)
    ans = np.asarray(a)
    return ans

# def eps(w_):
# 	w_pl = 1.5 * THz # THz
# 	G = 0.5 * THz # THz
# 	w1 = 2 * 3.14 * THz  # THz
# 	eps_ = 1 + w_pl ** 2 * ((w1 ** 2 - w_ ** 2) / ((w1 ** 2 - w_ ** 2) ** 2 + w_ ** 2 * G ** 2) +
# 	                        + 1j * (G * w_) / (((w1 ** 2 - w_ ** 2) ** 2) + w_ ** 2 * G ** 2))
#
# 	k0 = w_ / 3e+8
# 	n = np.sqrt(eps_)
# 	T = 4 * n / ((n+1)**2 * np.exp(-1j*k0*n*d)-(n-1)**2 * np.exp(1j*k0*n*d))
# 	# plt.plot(omega, T.real)
# 	# plt.plot(omega, T.imag)
# 	return n, T, k0
#
#
# plt.show()
# n, T, k0 = eps(omega)
def newton():
	def f(n_):
		g = 4 * n_ / ((n_+1)**2 * np.exp(-1j*k0*n_*d)-(n_-1)**2 * np.exp(1j*k0*n_*d)) - T
		return g
	def df(n_, h = 0.001):
		dg = (f(n_ + h) - f(n_)) / h
		return dg

	def iter(func, i = 0, xn = 1):
		xn_ = xn - f(xn) / df(xn)
		if i == 5:
			return xn_
		else:
			i += 1
			ans = iter(func, i, xn = xn_)
			return ans
		return xn_
	ans = iter(func = f(n), xn = 1)
	return ans
n = newton()

# plt.plot(omega, (n**2).imag, "*--")
# plt.plot(omega, (n**2).real, "*--")
plt.show()