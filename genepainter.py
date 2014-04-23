import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import time

class Painting(object):
	pass

class Stroke(object):
	def __init__(self, Nx, Ny, points=2, radius=1):
		pass

def stroke(x, y, n):
	if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
		raise Exception()

	t = np.linspace(0, 1, x.size)

	sx = scipy.interpolate.interp1d(t, x, kind='cubic')
	sy = scipy.interpolate.interp1d(t, y, kind='cubic')

	st = np.linspace(0, 1, n)

	return sx(st), sy(st)

def min_distance(sx, sy, px, py):
	dx = sx - px
	dy = sy - py
	d = dx ** 2 + dy ** 2
	return math.sqrt(np.amin(d))

def render(image):
	shape = image.shape
	Nx, Ny = shape[0], shape[1]

	x = np.array([5, 10, 15, 20, 5, 5])
	y = np.array([5, 5, 20, 15, 10, 30])

	sx, sy = stroke(x, y, 100)
	print min_distance(sx, sy, 10, 10)

	plt.plot(sx, sy, label='spline')

	for xi in xrange(Nx):
		for yi in xrange(Ny):
			d = min_distance(sx, sy, xi, yi)
			if d < 3.0:
				image[yi, xi, :] = np.array([3, 3, 3])

	plt.imshow(image, interpolation='none')
	plt.show()

if __name__ == "__main__":
	print 'lol'

	image_size = (128, 128, 3)
	image = np.zeros(image_size, dtype=np.float32)

	render(image)
