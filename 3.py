import numpy as np
from matplotlib import pyplot as plt
import math
import scienceplots
from scipy.fftpack import *

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


A1 = 10
A2 = 2.6
w1 = 1.3
w2 = 23.6
N = 500

fig, ax = plt.subplots(2,1,figsize=(5,5))
x = np.linspace(-1, 1, N)

def f(x_):
	f = A1 * np.cos(w1 * x_) + A2 * np.sin(w2 * x_)
	return f




def fourier(i):
	def a_n(n):
		ans = f(x) * np.cos(x * math.pi * n)/N * 2
		return ans.sum()

	def b_n(n):
		ans = f(x) * np.sin(x * math.pi * n)/N * 2
		return ans.sum()

	S_N = np.array([a_n(n) * np.cos(math.pi * n * x) + b_n(n) * np.sin(math.pi * n * x) for n in range(1,i)])

	return sum(S_N) + a_n(0)/2

ff = fourier(500)

ax[0].plot(x,ff, '*', color = 'orange')
ax[0].plot(x,f(x), color = 'brown')


t = np.linspace(0, math.pi, 500)
x = np.cos(t)

def chebishev(N):
	a_n0 = (f(x) / N).sum()

	def a_nn(n):
		ans = f(x) * np.cos(n * t) / N * 2
		return ans.sum()

	S_N = np.array([a_nn(n) * np.cos(n * t) for n in range(1, N)])

	return sum(S_N) + a_n0

fft = chebishev(200)




ax[1].plot(x,ff-f(x), color = 'orange')
# ax[1].plot(x,f(x), color = 'brown')

plt.show()