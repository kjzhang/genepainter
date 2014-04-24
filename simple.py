import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
import math
import random

filename = 'sample_small.jpg'
epsilon = 10
objects = 10
population_size = 10


class RectangleEncoding:
	def __init__(self, p1, p2, color, width):
		self.p1 = p1
		self.p2 = p2
		self.color = color
		self.width = width


def objective_function(current_rendition, actual):
	return np.sum(np.square(current_rendition - actual))


def point_in_square(encoding, x, y):
	(x1, y1) = encoding.p1
	(x2, y2) = encoding.p2
	m = 1.0*(y2-y1)/(x2-x1)
	const = y1 - m*x1

	a = m
	b = -1
	c = const

	dist = abs(a*x+b*y+c)/math.sqrt(a*a+b*b)
	xint = (b*(b*x - a*y)-a*c)/(a*a+b*b)
	yint = (a*(-b*x+a*y)-b*c)/(a*a+b*b)

	return dist < encoding.width/2 and xint >= min(x1,x2) and xint <= max(x2,x1)

def generate_random_encoding(shape):
	ymax, xmax, _ = shape
	red = random.randint(0,255)
	blue = random.randint(0,255)
	green = random.randint(0,255)
	x1 = random.randint(0, xmax)
	y1 = random.randint(0, ymax)
	x2 = random.randint(0, xmax)
	y2 = random.randint(0, ymax)
	width = random.randint(0, min(ymax,xmax)/2)
	return RectangleEncoding((x1,y1), (x2,y2), [red,green,blue],  width)


"""
object has 2 points, color, and width
-----------------------------------------

"""
def render_encoding(encodings, shape, dtype):
	aggregate = np.zeros(shape=shape, dtype=np.int32)
	weights = np.zeros(shape=shape, dtype=np.int32)
	row, col, channel = shape

	for i in xrange(row):
		for j in xrange(col):
				y = row - i
				x = j
				for encoding in encodings:
					if point_in_square(encoding, x, y):
						for k in xrange(channel):
							aggregate[i,j,k] = encoding.color[k]
							weights[i,j] += 1

	np.seterr(divide='ignore')
	return np.divide(aggregate,weights).astype(np.uint8)

img = misc.imread(filename)
population = []
for i in xrange(population_size):
	dna_strand = [generate_random_encoding(img.shape)]
	population.append(dna_strand)

min_error = 999999
for strand in population:
	error = objective_function(render_encoding(strand, img.shape, img.dtype), img)
	if error < min_error:
		min_error = error


for i in xrange(1,objects):

    """
    repeat process until satisfied
    """
    max_num_iter = 1000
    epsilon = 10
    num_iter = 0
    error_change = 9999
    while num_iter < max_num_iter and error_change < epsilon:

		#generate new candidates via mutation?
		new_candidates = gen_candidates(population)

		fitness = []
		for strand in new_candidates:
			error = objective_function(render_encoding(strand, img.shape, img.dtype), img)
			fitness.append(error)



	error_diff = 100

	curr_error = objective_function(render_encoding(curr_list), img)

		lowest_error = 9999999999
		for _ in xrange(max_iters):
			mutation = curr_list[:-1]
			mutation.append(generate_random_encoding(img.shape))
			new_error = objective_function(render_encoding(mutation), img)
			if new_error < lowest_error:
				lowest_error = new_error


"""
print img.shape
stuff = generate_random_encoding(img.shape)
print stuff.p1, stuff.p2, stuff.width, stuff.color
"""

plt.imshow(render_encoding([stuff],  img.shape, img.dtype))
plt.show()


