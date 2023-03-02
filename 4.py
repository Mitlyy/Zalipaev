import numpy as np
from matplotlib import pyplot as plt
import math
import scienceplots
from scipy.fftpack import *

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])
graph = input("График нужен?(y , n):")
dot = int(input("Примерный корень(1):"))

space = int(input("Рандом число(1):"))
A = np.matrix([[16, 3, 2], [3, 5, 1], [2, 1, 10]])

x = np.linspace(0, 20, 100)


def f(x):
	g = -x ** 3 + 31 * x ** 2 - 276 * x + 686
	return g


def df(x, h = 0.001):
	dg = (f(x + h) - f(x)) / h
	return dg


def iter(n = 0, xn = 0, x0 = 0):
	xn_ = xn - f(xn) / df(x0)
	if abs(xn_ - xn) < 1e-3:
		return xn_
	else:
		n += 1
		ans = iter(n, xn = xn_)
		return ans
	return xn_


print("Lambda:", '{0:.3}'.format(iter(xn = dot, x0 = space)))
# plt.plot(x, iter(func = f(x)), "*")
if graph == "y":
	plt.plot(x, f(x))
	plt.show()


def vectors(a):
	mat = np.matrix([[16, 3, 2], [3, 5, 1], [2, 1, 10]]) - a * np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
	vect = np.array([0, 0, 0])
	test = np.array([-0.3, 0.016, 1])

	def jacobi(A, f = 0, psi = 1e-5, x = np.array([0.1, 0.1, 0.1])):
		D = np.matrix(np.zeros(A.shape))
		for n in range(A.shape[0]):
			D[n, n] = A[n, n]

		def iter(x_n, n = 0):
			n += 1
			x_n_ = np.asarray(D.I.dot(f) - (D.I.dot(A - D)).dot(x_n))[0]

			if n > 400:
				print("Метод Якоби не сошелся в пределах psi")
				return x_n_, n

			elif all(i > psi for i in abs(x_n_ - x_n)):
				x_n_ = iter(x_n_, n)
				return x_n_
			else:
				return x_n_, n

		ans = iter(x)
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

		def iter(x_n, n = 0):
			n += 1
			x_n_ = np.asarray((L + D).I.dot(f) - ((L + D).I.dot(U)).dot(x_n))[0]
			if n > 400:
				print("Метод Зейделя не сошелся в пределах psi")
				return x_n_, n

			elif all(i > psi for i in abs(x_n_ - x_n)):
				x_n_ = iter(x_n_, n)
				return x_n_
			else:
				return x_n_, n

		ans2 = iter(x)
		return ans2

	def sor(A, f, psi = 1e-4, x = np.array([0.1, 0.1, 0.1]), w = 0.99):
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

		def iter(x_n, n = 0):
			n += 1
			x_n_ = np.asarray((L + D / w).I.dot(f) - ((L + D / w).I.dot(U + (w - 1) / w * D)).dot(x_n))[0]

			if n > 400:
				print("Метод SOR не сошелся в пределах psi")
				return x_n_, n

			elif all(i > psi for i in abs(x_n_ - x_n)):
				x_n_ = iter(x_n_, n)
				return x_n_
			else:
				return x_n_, n

		ans3 = iter(x)
		return ans3

	x, n = jacobi(mat, f = vect)
	print(np.around(x / x[2], 2), "Iter = ", n)
	x, n = seidel(mat, f = vect)
	print(np.around(x / x[2], 2), "Iter = ", n)
	x, n = sor(mat, f = vect)
	print(np.around(x / x[2], 2), "Iter = ", n)


vectors(iter(xn = dot, x0 = space))
