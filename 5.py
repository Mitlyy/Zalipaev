import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from numpy import loadtxt
from math import *
from scipy.fftpack import *
import scipy.integrate as integrate
import scipy.special as special

plt.style.use(['science', 'notebook', 'grid'])

fig, ax = plt.subplots(2, 2, figsize=(8,8))
data = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)
E_tr = data[1]
t = data[0]


eps = 1e-7
eps1 = 1e-9
dc = 10 / 3
f0 = 1.0
taup = 0.5
delf = 0.6
fmin = f0 - delf
fmax = f0 + delf
z = 1.0
ompl = 1.0
gam = 0.5
tz = t - (z/3)*10
omega = np.linspace(fmin, fmax, 1024)

E_in = 2 * np.cos(2 * np.pi * f0 * tz-4) * np.exp(-(tz-4)**2 / 4 / taup**2)

ax[0, 0].plot(t, E_in)

F = fftshift(fft(E_in))
# ax[1, 0].plot(omega, abs(F))
ax[0, 1].plot(t, E_tr)



def T(w_):
    dt = t[1]-t[0]
    a = []
    for i in range(t.shape[0]):
        a.append(sum(np.exp(2*np.pi*1j*w_[i]*(t-z/3*10)) * E_tr * dt))
    ans = np.asarray(a)
    return ans

T_ = T(omega)
ax[1, 0].plot(omega, T_.real)
ax[1, 0].plot(omega, T_.imag)




eps = 1e-7
eps1 = 1e-9
dc = 10 / 3
f0 = 1.0
taup = 0.5
delf = 0.6
fmin = f0 - delf
fmax = f0 + delf
z = 1.0
ompl = 1.0
gam = 0.5
maxt = 300
th = 1 / 20
tin = -1.
maxf=512
New=20
fh=2*delf/maxf
n=1
fa1=[]
maxt=1024
maxtd2=maxt//2
Tfm = np.zeros(1024) * 1j
epsm = np.zeros(1024) * 1j
for mf in range(1, maxf+2):
    f=fmin+(mf-1)*fh
    fa1.append(f)
    sum4=0
    sum2=0
    t=tin+th
    for m in range(1, maxtd2):
        m4=m*2-1
        sum4=sum4+E_tr[m4-1]*np.exp(2*np.pi*1j*f*(t-z/3*10))
        t=t+th
        m2=m*2
        sum2=sum2+E_tr[m2-1]*np.exp(2*np.pi*1j*f*(t-z/3*10))
        t=t+th
    Tf= (4*sum4+2*sum2) * th/3 * exp((2*pi*taup)**2*(f-f0)**2)/(2*sqrt(pi)*taup)
    Tfm[mf-1] = Tf

# ax[1, 1].plot(omega, Tfm.imag)
# ax[1, 1].plot(omega, Tfm.real)
def newton():
	def f(n_):
		g = 4 * n_ / ((n_+1)**2 * np.exp(-1j*n_* 2 * pi * omega * dc)-(n_-1)**2 * np.exp(1j*n_* 2 * pi * omega * dc)) - T_
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

n = newton()**2
print(n)
ax[1,1].plot(omega[300:], abs(n.imag)[300:], "--")
ax[1,1].plot(omega[300:], abs(n.real)[300:], "--")


plt.show()

