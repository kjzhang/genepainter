import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate


class Drawable(object):
	pass

class Line(Drawable):
	pass

def stroke(x, y, n):
	if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
		raise Exception()

	t = np.linspace(0, 1, x.size)

	sx = scipy.interpolate.interp1d(t, x, kind='cubic')
	sy = scipy.interpolate.interp1d(t, y, kind='cubic')

	st = np.linspace(0, 1, n)

	return sx(st), sy(st)

def render(image):
	x = np.array([5, 10, 15, 20, 5, 5])
	y = np.array([5, 5, 20, 15, 10, 30])

	sx, sy = stroke(x, y, 100)
	plt.plot(sx, sy, label='spline')

	plt.imshow(image)
	plt.show()

if __name__ == "__main__":
	print 'lol'

	image_size = (64, 64, 3)
	image = np.zeros(image_size, dtype=np.float64)

	render(image)
