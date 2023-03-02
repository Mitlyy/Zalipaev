import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


_ = np.linspace(-2, 2, 100)
x, y = np.meshgrid(_,_)
z = 5 * x**2 + 1/4 * x*y + 2 * y**2 - 4 * x - 5 * y + 5

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(x,y,z, cmap='coolwarm',
                       linewidth=0, antialiased=False)
ax.view_init(elev=10, azim=50)
plt.show()


def newtone():

	def f(x = 0, y = 0):
		g = 5 * x**2 + 1/4 * x*y + 2 * y**2 + 4 * x + 5 * y + 5
		return g

	def df(x, y = 0, h = 1e-5):
		dg_x = (f(x + h) - f(x)) / h
		dg_y = (f(x, y + h) - f(x, y)) / h
		return dg_x, dg_y

	def df2(x, y = 0,  h = 1e-5):
		dg2x = (df(x + h)[0] - df(x)[0]) / h
		dg2y = (df(x,y + h)[1] - df(y)[1]) / h
		return dg2x, dg2y

	def iterx(n = 0, xn = 0):
		xn_ = xn - df(xn)[0] / df2(xn)[0]
		# if abs(xn_ - xn) < 1e-5:
		if n > 400:
			return xn_
		else:
			n += 1
			ans = iterx(n, xn = xn_)
			return ans
		return xn_

	def itery(n = 0, yn = 0, x = 0):
		yn_ = yn - df(x, yn)[1] / df2(x, yn)[1]
		if abs(yn_ - yn) < 1e-3:
			return yn_
		else:
			n += 1
			ans = itery(n, yn = yn_, x = x)
			return ans
		return yn_
	x_ = iterx()
	y_ = itery()
	print([x_,y_, f(x_, y_)])
newtone()

