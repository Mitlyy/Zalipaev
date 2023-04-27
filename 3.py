import numpy as np
from matplotlib import pyplot as plt
import math
import scienceplots
from scipy.fftpack import *

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


A1 = 4
A2 = 2
w1 = 3
w2 = 5
Ndots = 100

fig, ax = plt.subplots(2,1,figsize=(5,5))
x = np.linspace(-1, 1, Ndots)

def f(x_):
	f = A1 * np.cos(w1 * x_) + A2 * np.sin(w2 * x_)
	return f


def fourier(i):
	def a_n(n):
		ans = f(x) * np.cos(x * math.pi * n) * 1/50
		return ans.sum()

	def b_n(n):
		ans = f(x) * np.sin(x * math.pi * n) * 1/50
		return ans.sum()

	S_N = np.array([a_n(n) * np.cos(math.pi * n * x) + b_n(n) * np.sin(math.pi * n * x) for n in range(1,i)])

	return sum(S_N) + a_n(0)/2

ff = fourier(64)

ax[0].plot(x,ff, 'o', color = 'orange')
ax[0].plot(x,f(x), color = 'brown')


t = np.linspace(0, math.pi, 200)
x = np.cos(t)

def chebishev(N):
	a_n0 = (f(np.cos(t))).sum() / 200

	def a_nn(n):
		ans = f(np.cos(t)) * np.cos(n * t)
		return ans.sum()/200 * 2

	S_N = np.array([a_nn(n) * np.cos(t*n) for n in range(1, N)])

	return sum(S_N) + a_n0

fft = chebishev(64)



x = x[2:fft.shape[0]-2]
fft = fft[2:fft.shape[0]-2]
# ax[1].plot(x,fft,"o", color = 'orange')
# ax[1].plot(x,f(x), color = 'brown')

ax[1].plot(x,fft,"o", color = 'orange')
ax[1].plot(x,f(x), color = 'orange')
plt.show()