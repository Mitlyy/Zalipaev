import numpy as np
from matplotlib import pyplot as plt
import math
from matplotlib.animation import ArtistAnimation
import scienceplots

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])

mat = np.matrix([[10, 3, 0], [3, 15, 1], [0, 1, 7]])
vect = np.array([2, 12, 5])
test = np.array([-29 / 977, 748 / 977, 591 / 977])



def jacobi(A, f=0, psi = 1e-5, x = np.array([0.1, 0.1, 0.1])):
	D = np.matrix(np.zeros(A.shape))
	for n in range(A.shape[0]):
		D[n, n] = A[n, n]

	def iter1(x_n, n = 0):
		n += 1
		x_n_ = np.asarray(D.I.dot(f) - (D.I.dot(A - D)).dot(x_n))[0]
		if (x_n_ - test).all() > psi:
			x_n_ = iter1(x_n_, n)
			return x_n_
		else:
			return x_n_, n

	ans = iter1(x)
	return ans

def seidel(A, f, psi = 1e-5, x = np.array([0.1, 0.1, 0.1])):
	D = np.matrix(np.zeros(A.shape))
	L = np.matrix(np.zeros(A.shape))
	U = np.matrix(np.zeros(A.shape))
	for n in range(A.shape[0]):
		for m in range(A.shape[1]):
			if n < m:
				U[n, m] = A[n, m]
			elif n > m:
				L[n, m] = A[n, m]
			else:
				D[n, m] = A[n, m]

	def iter2(x_n, n = 0):
		n += 1
		x_n_ = np.asarray((L+D).I.dot(f) - ((L+D).I.dot(U)).dot(x_n))[0]
		if (x_n_ - test).all() > psi:
			x_n_ = iter2(x_n_, n)
			return x_n_
		else:
			return x_n_, n
	ans2 = iter2(x)
	return ans2

def sor(A, f, psi = 1e-2, x = np.array([0.1, 0.1, 0.1]), w = 0.99):
	D = np.matrix(np.zeros(A.shape))
	L = np.matrix(np.zeros(A.shape))
	U = np.matrix(np.zeros(A.shape))
	for n in range(A.shape[0]):
		for m in range(A.shape[1]):
			if n < m:
				U[n, m] = A[n, m]
			elif n > m:
				L[n, m] = A[n, m]
			else:
				D[n, m] = A[n, m]
	def iter3(x_n, n=0):
		n += 1
		x_n_ = np.asarray((L+D/w).I.dot(f) - ((L+D/w).I.dot(U+(w-1)/w*D)).dot(x_n))[0]
		if n > 300:
			return x_n_

		elif (x_n_ - test).all() > psi:
			x_n_ = iter3(x_n_, n)
			return x_n_
		else:
			return x_n_, n
	ans3 = iter3(x)
	return ans3

#
# Y = []
# x = np.linspace(1, 1.1, 100)
# for i in x:
# 	Y.append(sor(mat, f = vect, w = i)[1])
#
# plt.plot(x,Y)
# plt.show()

print(jacobi(mat, f = vect))
print(seidel(mat, f = vect))
print(sor(mat, f = vect, w = 1.018))
