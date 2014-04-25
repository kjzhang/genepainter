import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
import math
import random
import copy
import genepainter

filename = 'sample_small.jpg'
img = misc.imread(filename)

#population size and strokes per strand
num_strokes = 50
population_size = 10
min_spline_points = 2
max_spline_points = 3

#evolving parameters
add_mutation_rate = 0.2
change_or_delete_mutation_rate = 0.07
# delete from 0 to 0.5, 0.5 to 1 is change
conditional_delete_mutation_rate = 0.5
num_children = 25

#crossover rates
crossover_rate = 0.2
num_crossovers = 25

#iterations
min_iter = 100
max_num_iter = 1000
epsilon = 10
num_iter = 0
error_change = 9999


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
    num_spline_points = random.randint(min_spline_points, max_spline_points)
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

    return None

def gen_candidates(population):
    new_candidates = []

    for _ in xrange(num_children):
        random_child_index = random.randint(0, population_size)
        new_strand = copy.deepcopy(population[random_child_index].strand)

        delete_indexes = []

        for i in xrange(len(new_strand)):
            rand = random.random()
            if rand <= change_or_delete_mutation_rate:
                #rescale rand
                rand = rand/change_or_delete_mutation_rate
                if rand <= conditional_delete_mutation_rate:
                    delete_indexes.append(i)
                else:
                    new_strand[i] = generate_random_encoding(img.shape)

        #delete
        for i in xrange(len(delete_indexes)-1, -1, -1):
            del new_strand[i]

        if random.random() <= add_mutation_rate:
            new_strand.append(generate_random_encoding(img.shape))

        new_error = objective_function(render_encoding(new_strand, img.shape, img.dtype), img)
        new_dna_strand = DNAStrand(new_strand, new_error)
        new_candidates.append(new_dna_strand)
    return new_candidates

def crossover(population):
    new_candidates = []

    for _ in xrange(num_crossovers):
        rand1 = random.randint(0, population_size-1)
        rand2 = random.randint(0, population_size-1)
        if rand1 != rand2:
            strand_a = population[rand1].strand
            strand_b = population[rand2].strand
            new_strand = copy.deepcopy(strand_a)
            for i in xrange(len(strand_a)):
                if i > len(strand_b):
                    break
                if random.random() <= crossover_rate:
                    new_strand[i] = strand_b[i]
            new_error = objective_function(render_encoding(new_strand, img.shape, img.dtype), img)
            new_dna_strand = DNAStrand(new_strand, new_error)
            new_candidates.append(new_dna_strand)
    return new_candidates


population = []
for i in xrange(population_size):
    strand = []
    for _ in xrange(num_strokes):
        strand.append(generate_random_encoding(img.shape))
    error = objective_function(render_encoding(strand, img.shape, img.dtype), img)
    dna_strand = DNAStrand(strand, error)
    population.append(dna_strand)



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


