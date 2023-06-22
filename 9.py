import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots
from numpy import loadtxt
from scipy.fftpack import *

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])



omega = np.linspace(0, 5e+12, 1024)
wpl = 0.5e12
G = 0.1e12
w1 = 2*np.pi*1e12
F1 = 1
f0 = 1e12
tp = 0.5e-12
z = 1e-3
d = 0.5e-3
w0 = 2*np.pi*f0
c = 3e8
k0 = w0/c


eps = 1e-7
eps1 = 1e-9
dc = 10 / 3
f0 = 1.0
tp = 0.5
delf = 0.6
fmin = f0 - delf
fmax = f0 + delf
z = 1.0
omj = 2 * np.pi * 1.3
ompl = 1.0
gam = 0.5


data = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)
E = data[1]
t = data[0] * 1e-12
t_ = np.linspace(-10, 10, 1024)

c = 3e+8


def func(x):
	return np.exp(-x**2 / (4 * tp**2)) * np.cos(w0 * x)


# F = fft(func(t_))
# plt.plot(t, E)
# plt.show()

# plt.plot(omega, F)
# plt.show()


def T(w_):
    dt = t[1]-t[0]
    a = []
    for i in range(t.shape[0]):
        a.append(sum(np.exp(1j * w_[i] * (t - z / c)) * E) * dt)
    ans = np.asarray(a)
    return ans
T = T(omega)


plt.plot(omega, T.imag)
plt.plot(omega, T.real)

plt.show()

def eps(w_):
	w_pl = 1.5 * THz # THz
	G = 0.5 * THz # THz
	w1 = 1 * THz  # THz
	eps_ = 1 + w_pl ** 2 * ((w1 ** 2 - w_ ** 2) / ((w1 ** 2 - w_ ** 2) ** 2 + w_ ** 2 * G ** 2) +
	                        + 1j * (G * w_) / (((w1 ** 2 - w_ ** 2) ** 2) + w_ ** 2 * G ** 2))
	# plt.plot(omega, eps_.real)
	# plt.plot(omega, eps_.imag)
	k0 = w_ / 3e+8
	n = np.sqrt(eps_)
	T = 4 * n / ((n+1)**2 * np.exp(-1j*k0*n*d)-(n-1)**2 * np.exp(1j*k0*n*d))
	return n, T, k0


# n, T_, k0 = eps(omega)
# T = T(omega)
# plt.plot(omega, T)
# plt.plot(omega, T)
# plt.show()
def newton():
	k0 = omega / 3e+8
	def f(n_):
		g = 4 * n_ / ((n_+1)**2 * np.exp(-1j*k0*n_*d)-(n_-1)**2 * np.exp(1j*k0*n_*d)) - T
		return g
	def df(n_, h = 0.001):
		dg = (f(n_ + h) - f(n_)) / h
		return dg
	def iter(i = 0, xn = 1):
		xn_ = xn - f(xn) / df(xn)
		# print(xn_[0])
		if i == 3:
			return xn_
		else:
			i += 1
			ans = iter(i, xn = xn_)
			return ans
		return xn_
	ans = iter(0, xn = 1)
	return ans

# n_ = newton()
#
# plt.plot(omega, abs(n_.imag), "*")
# plt.plot(omega, abs(n_.real), "*")
# plt.plot(omega, abs(n.imag))
# plt.plot(omega, abs(n.real))
# plt.show()