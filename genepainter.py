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
population_size = 100
mutation_rate = 0.5
min_iter = 10
num_spline_points = 4

def read_image(file, use_alpha=False):
	image_raw = Image.open(file)
	image_raw.load()
	image = np.array(image_raw, dtype=np.float32) / 255
	return image if use_alpha else rgba_to_rgb(image)

def rgba_to_rgb(image_raw):
	shape = (image_raw.shape[0], image_raw.shape[1], 3)

	R = image_raw[:, :, 0]
	G = image_raw[:, :, 1]
	B = image_raw[:, :, 2]
	A = image_raw[:, :, 3]

	image = np.zeros(shape, dtype=np.float32)
	image[:, :, 0] = np.multiply(A, R) + ((1.0 - A) * 1.0)
	image[:, :, 1] = np.multiply(A, G) + ((1.0 - A) * 1.0)
	image[:, :, 2] = np.multiply(A, B) + ((1.0 - A) * 1.0)

	return image

class GeneStroke(object):
	def __init__(self, sx, sy, color, alpha, width):
		if sx.ndim != 1 or sy.ndim != 1 or sx.size != sy.size:
			raise Exception()

		self.sx = sx
		self.sy = sy
		self.color = color
		self.alpha = alpha
		self.width = width

	def spline(self, resolution=100, interpolation='cubic'):
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
		self.fitness = float('inf')

		for _ in xrange(10):
			self.strokes.append(GeneStroke.random(self.shape))

	def render(self):
		w, h = self.shape
		r = RendererAgg(w, h, 72)
		arr = np.frombuffer(r.buffer_rgba(), np.uint8)
		arr.shape = r.height, r.width, -1

		for stroke in self.strokes:
			stroke.render(r)

		image_raw = np.float32(arr, dtype=np.float32) / 255
		return rgba_to_rgb(image_raw)

	def update_fitness(self, source):
		image = self.render()
		self.fitness = np.sum(np.square(source - image))

class GenePainter(object):
	def __init__(self, source):
		self.source = source
		self.shape = source.shape[0], source.shape[1]
		self.population = []

	def paint(self):

		for _ in xrange(population_size):
			image = GeneImage(self.shape)
			image.update_fitness(self.source)
	    	self.population.append(image)

		for i in xrange(10):
			print i

			self.population = sorted(self.population, key=lambda image: image.fitness)
			print [image.fitness for image in self.population]

			min_error = self.population[0].fitness

			for _ in xrange(population_size):
				image = GeneImage(self.shape)
				image.update_fitness(self.source)
				self.population.append(image)

			self.population = sorted(self.population, key=lambda image: image.fitness)
			self.population = self.population[:population_size]

		plt.imshow(self.population[0].render(), interpolation='none')
		plt.show()

if __name__ == "__main__":

	source = read_image('ML257.png')

	p = GenePainter(source)

	i = GeneImage((256, 256))
	image = i.render()

	p.paint()
