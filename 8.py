import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots
from numpy import loadtxt
plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


lines = loadtxt(r"C:\Users\User\Desktop\epsdata.txt", comments="#", delimiter=" ", unpack=True)



def F(w):
	f = 2 * np.sqrt(np.pi) * tau_p * np.exp(-((w-w0)*tau_p)**2)
	return f
print(lines[0], lines[1].real)
plt.plot(lines[0], lines[1].real)
plt.show()