import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots
from numpy import loadtxt
from scipy.fftpack import *

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])

omega = np.linspace(0, 9e+12, 1024)
THz = 1e+12
f0 = 1 * THz
tau_p = 0.5/ THz
omega_0 = 2 * np.pi * f0
z = 1e-3
d = 0.5e-3
data = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)
E = data[1]
t = data[0]

c = 3e+8
def F(w_):
    return 2 * np.sqrt(np.pi) * tau_p * np.exp(-(w_-omega_0)**2 * tau_p**2)
# plt.plot(omega, F(omega))
# plt.show()
F = F(omega)

def T(w_):
    dt = t[1]-t[0]
    a = []
    for i in range(t.shape[0]):
        a.append(1/F[i] * sum(np.exp(1j * w_[i] * (t - z / c)) * E) * dt)
    ans = np.asarray(a)/max(a)
    return ans


def eps(w_):
	w_pl = 1.5 * THz # THz
	G = 0.5 * THz # THz
	w1 = 2 * 3.14 * THz  # THz
	eps_ = 1 + w_pl ** 2 * ((w1 ** 2 - w_ ** 2) / ((w1 ** 2 - w_ ** 2) ** 2 + w_ ** 2 * G ** 2) +
	                        + 1j * (G * w_) / (((w1 ** 2 - w_ ** 2) ** 2) + w_ ** 2 * G ** 2))
	# plt.plot(omega, eps_.real)
	# plt.plot(omega, eps_.imag)
	k0 = w_ / 3e+8
	n = np.sqrt(eps_)
	T = 4 * n / ((n+1)**2 * np.exp(-1j*k0*n*d)-(n-1)**2 * np.exp(1j*k0*n*d))
	return n, T, k0


n, T_, k0 = eps(omega)
# T = T(omega)
# plt.plot(omega, F)
plt.plot(omega, T_)
plt.show()
def newton():
	k0 = omega / 3e+8
	def f(n_):
		g = 4 * n_ / ((n_+1)**2 * np.exp(-1j*k0*n_*d)-(n_-1)**2 * np.exp(1j*k0*n_*d)) - T_
		return g
	def df(n_, h = 0.001):
		dg = (f(n_ + h) - f(n_)) / h
		return dg
	def iter(i = 0, xn = 1):
		xn_ = xn - f(xn) / df(xn)
		# print(xn_[0])
		if i == 4:
			return xn_
		else:
			i += 1
			ans = iter(i, xn = xn_)
			return ans
		return xn_
	ans = iter(0, xn = 1)
	return ans

n = newton()
print(n.shape)
plt.plot(omega, abs(n.imag), "*--")
plt.plot(omega, abs(n.real), "*--")
plt.show()