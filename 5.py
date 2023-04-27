import numpy as np
from matplotlib import pyplot as plt
import math
import scienceplots
from scipy.fftpack import *
from numpy import loadtxt
plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


data = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)
E = data[1]
t = data[0]
z = 3e-3
c = 3e+8
omega = np.linspace(6.27e+12, 6.3e+12, t.shape[0])
n = omega*0 + 10
k0 = omega/c
d = 0.5e-3

F = fft(E)
def T(w_):
    dt = t[1]-t[0]
    a = []
    for i in range(t.shape[0]):
        a.append(1/F[i] * sum(np.exp(1j * w_[i] * (t - z / c)) * E) * dt)
    ans = np.asarray(a)
    return ans

def newton():

    def f(n_):
        g = 4 * n_ / ((n_+1)**2 * np.exp(-1j*k0*n_*d)-(n_-1)**2 * np.exp(1j*k0*n_*d)) - T(omega)
        plt.plot(omega, g.real)
        return g

    def df(n_, h = 0.001):
        dg = (f(n_ + h) - f(n_)) / h
        return dg

    def iter(func, i = 0, xn = 1):
        xn_ = xn - f(xn) / df(xn)
        print(xn_)
        if i == 2:
            return xn_
        else:
            i += 1
            ans = iter(func, i, xn = xn_)
            return ans
        return xn_

    ans = iter(func = f(1), xn = n)
    return ans
# n = newton()


def eps():
    omega1 = 2 * np.pi * 1e+12
    G_1 = 0.5e+12
    omega_pl = 1.5e+12
    b = omega1**2 - omega**2
    return 1 + omega_pl**2 * (b/(b**2 + omega**2 * G_1) + 1j * (G_1* omega)/(b**2 + omega**2 * G_1))
n = np.sqrt(eps())

f = 4 * n / ((n+1)**2 * np.exp(-1j*k0*n*d)-(n-1)**2 * np.exp(1j*k0*n*d))

plt.plot(omega, f.real)
# plt.plot(omega, n.imag)
plt.show()
# plt.plot(omega, (n**2).imag, "*--")
# plt.plot(omega, (n**2).real, "*--")

plt.show()
# plt.plot(omega, eps(omega).imag)
