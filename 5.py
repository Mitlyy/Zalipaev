import numpy as np
from matplotlib import pyplot as plt
import math
import scienceplots
from scipy.fftpack import *
from numpy import loadtxt
plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


data = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)
E = data[1]

def T(w_):
    F = fft(E)
	ans = 1/F * np.exp(1j * w(t - z/c)) * E
    return ans


# plt.plot(omega, eps(omega).real)
# plt.plot(omega, eps(omega).imag)
plt.show()

n, T, k0 = eps(omega)

def newton():

	def f(n_):
		g = 4 * n_ / ((n_+1)**2 * np.exp(-1j*k0*n_*d)-(n_-1)**2 * np.exp(1j*k0*n_*d)) - T
		return g

	def df(n_, h = 0.001):
		dg = (f(n_ + h) - f(n_)) / h
		return dg

	def iter(func, i = 0, xn = 1):
		xn_ = xn - f(xn) / df(xn)
		print(xn_[0])
		if i == 4:
			return xn_
		else:
			i += 1
			ans = iter(func, i, xn = xn_)
			return ans
		return xn_

	ans = 	iter(func = f(n), xn = 1)
	return ans
n = newton()


plt.plot(omega, (n**2).imag, "*--")
plt.plot(omega, (n**2).real, "*--")

plt.show()
# plt.plot(omega, eps(omega).imag)
