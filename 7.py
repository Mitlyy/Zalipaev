import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])

dots = 4000
_ = np.linspace(-1, 1, dots)

x, y = np.meshgrid(_,_)

func = 2 * x**2 + y**2 + 4
func = np.select([(x ** 2 + y ** 2 <= 1), True], [func, 0])
rn = np.random.random((dots,dots)) * 6
rn = np.select([rn < func, True], [1, 0])

def one(f, d = dots):
	ans = 4/d**2 * sum(sum(f))
	ans2 = 4*6 * sum(sum(rn))/dots**2
	print(np.around(ans, 3))
	print(np.around(ans2, 3))
one(func)

# plt.contourf(x, y, func)
# plt.colorbar()
# plt.show()