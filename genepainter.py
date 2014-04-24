import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
import time

from PIL import Image
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.lines import Line2D

class GenePainter(object):
	def __init__(self, source):
		self.source = source

	def render(self):
		output = np.zeros(self.source.shape, dtype=np.float32)
		Ny, Nx = output.shape[0], output.shape[1]

		#x = np.array([5, 10, 15, 20, 5, 5])
		#y = np.array([5, 5, 20, 15, 10, 30])

		x = np.array(np.random.random(4) * 128, dtype=np.float32)
		y = np.array(np.random.random(4) * 128, dtype=np.float32)

		sx, sy = spline(x, y, 1000)

		t = time.time()
		for yi in xrange(Ny):
			for xi in xrange(Nx):
				d = min_distance(sx, sy, xi, yi)
				if d < 10.: # radius
					output[yi, xi, :] = np.array([1, 1, 0, 0.5])
		print time.time() - t


		w, h = 256, 256
		r = RendererAgg(256, 256, 72)
		arr = np.frombuffer(r.buffer_rgba(), np.uint8)
		arr.shape = r.height, r.width, -1
		t = np.linspace(0, 2*np.pi, 100)
		x = np.sin(2*t) * w*0.45 + w*0.5
		y = np.cos(3*t) * h*0.45 + h*0.5

		t = time.time()
		line = Line2D(sx, sy, linewidth=50, color=(1.0, 0.0, 0.0), alpha=0.3)
		line.draw(r)
		print time.time() - t
		# plt.imsave("test.png", arr)

		# t = time.time()
		# for _ in xrange(100):
		# 	plt.plot(sx, sy, label='spline', linewidth=10, aa=False, solid_capstyle="round")
		# print time.time() - t

		plt.imshow(arr, interpolation='none')
		plt.show()

	def score(self, image):
		return np.linalg.norm(self.source - image, 2)

def spline(x, y, n):
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

def read_image(file):
	image_raw = Image.open(file)
	image_raw.load()
	# return np.array(image_raw, dtype=np.float32)
	image_rgb = Image.new('RGB', image_raw.size)
	image_rgb.paste(image_raw, None)
	return np.array(image_rgb, dtype=np.float32)

class GeneImage(object):
	def __init__(self):
		self.strokes = []

class GeneStroke(object):
	def __init__(self, x, y, radius):
		if x.ndim != 1 or y.ndim != 1 or x.size != y.size:
			raise Exception()

		self.x = x
		self.y = y
		self.radius = radius

if __name__ == "__main__":

	# source = read_image('ML129.png')
	source = np.zeros((256, 256, 4), dtype=np.float32)

	p = GenePainter(source)
	p.render()
