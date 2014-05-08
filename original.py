import math
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.interpolate
import time
import copy
import sys

from PIL import Image
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.lines import Line2D

from genepainter import dump_dna

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
min_iter = 50000
max_num_iter = 10000
epsilon = 10
num_iter = 0

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
    def __init__(self, shape, sx, sy, color, alpha, width):
        if sx.ndim != 1 or sy.ndim != 1 or sx.size != sy.size:
            raise Exception()

        self.shape = shape
        self.sx = sx
        self.sy = sy
        self.color = color
        self.alpha = alpha
        self.width = width

    def spline(self, resolution=1000, interpolation='cubic'):
        if self.sx.size <= 3:
            return self.sx, self.sy

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

        num_spline_points = random.randint(min_spline_points, max_spline_points)

        pointsx = np.random.random(num_spline_points) * xmax
        pointsy = np.random.random(num_spline_points) * ymax
        color = tuple(np.random.random(3))
        alpha = np.random.uniform(0.25, 1)
        width = random.randint(0, min(ymax, xmax) / 4)

        return cls(shape, pointsx, pointsy, color, alpha, width)

class DNAImage(object):
    def __init__(self, shape, strokes=None):
        self.shape = shape
        self.fitness = float('inf')

        if strokes is None:
            self.strokes = []
            for _ in xrange(num_strokes):
                self.strokes.append(GeneStroke.random(self.shape))
        else:
            self.strokes = strokes

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

    def mutate(self):
        new_strokes = copy.deepcopy(self.strokes)

        delete_indexes = []

        for i in xrange(len(new_strokes)):
            rand = random.random()
            if rand <= change_or_delete_mutation_rate:
                #rescale rand
                rand = rand/change_or_delete_mutation_rate
                if rand <= conditional_delete_mutation_rate:
                    delete_indexes.append(i)
                else:
                    new_strokes[i] = GeneStroke.random(self.shape)

        #delete
        for i in xrange(len(delete_indexes)-1, -1, -1):
            del new_strokes[i]

        if random.random() <= add_mutation_rate:
            new_strokes.append(GeneStroke.random(self.shape))
        new_dna = DNAImage(self.shape, new_strokes)
        return new_dna


class GenePainter(object):
    def __init__(self, source):
        self.source = source
        self.shape = source.shape[0], source.shape[1]
        self.population = []

    def paint(self):

        for _ in xrange(population_size):
            dna = DNAImage(self.shape)
            dna.update_fitness(self.source)
            self.population.append(dna)

        num_iter = 0
        error_change = 9999
        min_error = 999999

        while (num_iter < max_num_iter and error_change >= epsilon) or num_iter <= min_iter:

            print num_iter

            #generate new candidates via mutation?
            new_candidates = []
            for _ in xrange(num_children):
                rand_int = random.randint(0, population_size-1)
                mutation = self.population[rand_int].mutate()
                mutation.update_fitness(self.source)
                new_candidates.append(mutation)

            compete_population = self.population + new_candidates

            # pick the best
            surviving_children = sorted(compete_population, key=lambda dna: dna.fitness)[:population_size]

            #do some crossovers
            #new_candidates = gen_crossovers(population)

            #compete_population = population + new_candidates

            # pick the best
            #surviving_children = sorted(compete_population, key=lambda dna_strand: dna_strand.fitness)[:population_size]

            new_best_error = surviving_children[0].fitness

            error_change = min_error - new_best_error
            print error_change
            min_error = new_best_error
            self.population = surviving_children
            num_iter += 1

            if error_change > 0.0:
                for z in xrange(10):
                    dump_dna(self.population[0].strokes, 'original-' + sys.argv[1] + '-' + str(num_iter) + '-' + str(z) + str(np.random.random()) + '.txt')

        plt.imshow(self.population[0].render(), interpolation='none')
        plt.show()

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: %s image" % sys.argv[0])
        sys.exit(-1)

    source = read_image(sys.argv[1])

    p = GenePainter(source)

    p.paint()
