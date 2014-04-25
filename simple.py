import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
import math
import random
import copy
import genepainter

filename = 'sample_small.jpg'
img = misc.imread(filename)
epsilon = 10
objects = 10
population_size = 10
mutation_rate = 0.5
min_iter = 10
num_spline_points = 3

class StrokeEncoding:
    def __init__(self, pointsx, pointsy, color, alpha, width):
        self.pointsx = pointsx
        self.pointsy = pointsy
        self.color = color
        self.alpha = alpha
        self.width = width

class DNAStrand:
    def __init__(self, strand, fitness):
        self.strand = strand
        self.fitness = fitness


def objective_function(current_rendition, actual):
    fitness = np.sum(np.square(current_rendition - actual))
    return fitness

def generate_random_encoding(shape):
    ymax, xmax, _ = shape
    pointsx = np.random.random(num_spline_points) * xmax
    pointsy = np.random.random(num_spline_points) * ymax

    color = tuple(np.random.random(3))
    alpha = random.random()
    width = random.randint(0, min(ymax,xmax)/2)
    return StrokeEncoding(pointsx, pointsy, color, alpha, width)

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

def gen_candidates(population):
    new_candidates = []
    for dna_strand in population:
        new_strand = copy.deepcopy(dna_strand.strand)
        for i in xrange(len(new_strand)):
            if random.random() <= mutation_rate:
                new_strand[i] = generate_random_encoding(img.shape)
        new_error = objective_function(render_encoding(new_strand, img.shape, img.dtype), img)
        new_dna_strand = DNAStrand(new_strand, new_error)
        new_candidates.append(new_dna_strand)
    return new_candidates


population = []
for i in xrange(population_size):
    strand = []
    error = 0
    dna_strand = DNAStrand(strand, error)
    population.append(dna_strand)


for i in xrange(1,objects):

    print "add object ", i

    # add a gene
    for dna_strand in population:
        dna_strand.strand.append(generate_random_encoding(img.shape))
        dna_strand.fitness = objective_function(render_encoding(dna_strand.strand, img.shape, img.dtype), img)

    population = sorted(population, key=lambda dna_strand: dna_strand.fitness)
    min_error = population[0].fitness

    """
    repeat process until satisfied
    """
    max_num_iter = 1000
    epsilon = 10
    num_iter = 0
    error_change = 9999
    while (num_iter < max_num_iter and error_change >= epsilon) or num_iter <= min_iter:

        print "iter"

        #generate new candidates via mutation?
        new_candidates = gen_candidates(population)

        compete_population = population + new_candidates

        # pick the best
        surviving_children = sorted(compete_population, key=lambda dna_strand: dna_strand.fitness)[:population_size]

        #do some crossovers
        #new_candidates = gen_crossovers(population)

        #compete_population = population + new_candidates

        # pick the best
        #surviving_children = sorted(compete_population, key=lambda dna_strand: dna_strand.fitness)[:population_size]

        new_best_error = surviving_children[0].fitness

        error_change = abs(min_error - new_best_error)
        print error_change
        min_error =new_best_error
        population = surviving_children
        num_iter += 1

    


best_encoding = population[0].strand


"""
print img.shape
stuff = generate_random_encoding(img.shape)
print stuff.p1, stuff.p2, stuff.width, stuff.color
"""

plt.imshow(render_encoding(best_encoding,  img.shape, img.dtype))
plt.show()


