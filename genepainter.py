import math
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.interpolate
import time

from PIL import Image
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.lines import Line2D

epsilon = 10
objects = 10
population_size = 10
mutation_rate = 0.5
min_iter = 10
num_spline_points = 4

class GeneStroke(object):
	def __init__(self, sx, sy, color, alpha, width):
		if sx.ndim != 1 or sy.ndim != 1 or sx.size != sy.size:
			raise Exception()

		self.sx = sx
		self.sy = sy
		self.color = color
		self.alpha = alpha
		self.width = width

	def spline(self, resolution=1000, interpolation='cubic'):
		t = np.linspace(0, 1, self.sx.size)
		sx = scipy.interpolate.interp1d(t, self.sx, kind=interpolation)
		sy = scipy.interpolate.interp1d(t, self.sy, kind=interpolation)
		st = np.linspace(0, 1, resolution)
		return sx(st), sy(st)

	def render(self, r):
		sx, sy = self.spline()
		stroke = Line2D(sx, sy, color=self.color, alpha=self.alpha, lw=self.width, solid_capstyle="round")
		stroke.draw(r)

	@classmethod
	def random(cls, shape):
		ymax, xmax = shape
		pointsx = np.random.random(num_spline_points) * xmax
		pointsy = np.random.random(num_spline_points) * ymax

		color = tuple(np.random.random(3))
		alpha = np.random.uniform(0.25, 1)
		width = random.randint(0, min(ymax, xmax) / 4)

		return cls(pointsx, pointsy, color, alpha, width)

class GeneImage(object):
	def __init__(self, shape):
		self.shape = shape
		self.strokes = []

		for _ in xrange(10):
			self.strokes.append(GeneStroke.random(self.shape))

	def render(self):
		w, h = self.shape
		r = RendererAgg(w, h, 72)
		arr = np.frombuffer(r.buffer_rgba(), np.uint8)
		arr.shape = r.height, r.width, -1

		t = time.time()
		for stroke in self.strokes:
			stroke.render(r)
		print time.time() - t

		return np.array(arr, dtype=np.float32) / 255

class GenePainter(object):
	def __init__(self, source):
		self.source = source

	def fitness(self, image):
		return np.linalg.norm(self.source - image, 2)

	def paint(self):
		pass

def min_distance(sx, sy, px, py):
	dx = sx - px
	dy = sy - py
	d = dx ** 2 + dy ** 2
	return math.sqrt(np.amin(d))

def read_image(file):
	image_raw = Image.open(file)
	image_raw.load()
	return np.array(image_raw, dtype=np.float32)

if __name__ == "__main__":

	# source = read_image('ML129.png')
	source = np.zeros((256, 256, 4), dtype=np.float32)

	p = GeneImage((256, 256))
	image = p.render()

	plt.imshow(image, interpolation='none')
	plt.show()
