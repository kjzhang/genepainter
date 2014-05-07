import math
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.interpolate
import time
import copy
import datetime

from PIL import Image
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.lines import Line2D

# images per population
population_size = 100
population_interval = 1

# bias of how long a newly created gene lasts
dna_life_bias = 1

# migration
migrate_count = 20


### STROKE PARAMETERS

# strokes per image
strokes_start = 10
strokes_min = 10
strokes_max = 50

# points per stroke
min_spline_points = 2
max_spline_points = 3

# alpha start
alpha_min = 0.25
alpha_max = 1.0

# mutation, stroke level
p_spline_point_add = 0.1
p_spline_point_edit = 0.5
p_spline_point_delete = 0.5

### IMAGE PARAMETERS

# mutation, image level
p_stroke_add = 0.
p_stroke_edit = 0
p_stroke_delte = 0.5

### POPULATION PARAMETERS

# mutation, population level
p_crossover = 0

#evolving parameters
add_mutation_rate = 0.3
change_or_delete_mutation_rate = 0.2
# delete from 0 to 0.5, 0.5 to 1 is change
conditional_delete_mutation_rate = 0.5
num_children = 50

#crossover rates
crossover_rate = 0.2
num_crossovers = 25

#iterations
min_num_iter = 50000
max_num_iter = 100000
num_iter = 0
epsilon = 10


def load_dna(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    strokes = []
    for line in lines:
        info = line.split(';')
        shape_tmp = info[0].split(':')[1][1:-1].split(',')
        shape = (int(shape_tmp[0][:-1]), int(shape_tmp[1][:-1]))

        sx_tmp = info[1].split(':')[1][1:-1].split(' ')
        sx = []
        for str_rep in sx_tmp:
            if len(str_rep) > 0:
                sx.append(float(str_rep))
        sx = np.array(sx)

        sy_tmp = info[2].split(':')[1][1:-1].split(' ')
        sy = []
        for str_rep in sy_tmp:
            if len(str_rep) > 0:
                sy.append(float(str_rep))
        sy = np.array(sy)

        color_tmp = info[3].split(':')[1][1:-1].split(',')
        color = tuple([ float(val) for val in color_tmp])

        alpha = float(info[4].split(':')[1])
        width = int(info[5].split(':')[1])
        strokes.append(GeneStroke(shape, sx, sy, color, alpha, width))

    return strokes

def dump_dna(dna, filename=None):
    if filename is None:
        filename = 'dna' + str((datetime.datetime.now() - datetime.datetime(1970,1,1)).total_seconds()) + '.txt'
    f = open(filename, 'w')
    for stroke in dna:
        f.write('shape:' + str(stroke.shape) + ';sx:' + str(stroke.sx) + ';sy:' + str(stroke.sy) + ';color:' + str(stroke.color) + ';alpha:' + str(stroke.alpha) + ';width:' + str(stroke.width) + "\n")
    f.close()

def output_name(iteration, id):
    pass    

def read_image(file, use_alpha=False):
    image_raw = Image.open(file)
    image_raw.load()
    image = np.array(image_raw, dtype=np.float32) / 255
    return image if use_alpha else rgba_to_rgb(image)

def rgba_to_rgb(image_raw, background=(1.0, 1.0, 1.0)):
    """ Convert an RGBA image to an RGB image.
    """

    shape = (image_raw.shape[0], image_raw.shape[1], 3)

    R = image_raw[:, :, 0]
    G = image_raw[:, :, 1]
    B = image_raw[:, :, 2]
    A = image_raw[:, :, 3]

    BG_R, BG_G, BG_B = background

    image = np.zeros(shape, dtype=np.float32)
    image[:, :, 0] = np.multiply(A, R) + ((1.0 - A) * BG_R)
    image[:, :, 1] = np.multiply(A, G) + ((1.0 - A) * BG_G)
    image[:, :, 2] = np.multiply(A, B) + ((1.0 - A) * BG_B)

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

    def mutate(self):
        # spline point add
        if np.random.random() < p_spline_point_add and self.sx.size < max_spline_points:
            px, py = self.random_point()
            self.sx.append(px)
            self.sy.append(py)

        # spline point delete
        if np.random.random() < p_spline_point_delete and self.sx.size > min_spline_points:
            pass

    def mutate_value(self, center, radius, min_value, max_value):
        valid = False
        while not valid:
            new_value = np.random.normal(center, radius)
            if new_value >= min_value and new_value <= max_value:
                valid=True
        return new_value

    def random_color(self, center=None, radius=None):
        if center is not None and radius is not None:
            old_R, old_G, old_B = center
            new_R = self.mutate_value(old_R, radius, 0.0, 1.0)
            new_G = self.mutate_value(old_G, radius, 0.0, 1.0)
            new_B = self.mutate_value(old_B, radius, 0.0, 1.0)
        else:
            new_R = np.random.random()
            new_G = np.random.random()
            new_B = np.random.random()

        return new_R, new_G, new_B

    def random_point(self, center=None, radius=None):
        ymax, xmax = self.shape

        if center is not None and radius is not None:
            cx, cy = center
            px = self.mutate_value(cx, radius, 0.0, xmax)
            py = self.mutate_value(cy, radius, 0.0, ymax)
        else:
            px = np.random.random() * ymax
            py = np.random.random() * xmax

        return px, py

    @classmethod
    def random(cls, shape):
        ymax, xmax = shape

        num_spline_points = np.random.randint(min_spline_points, max_spline_points + 1)
        pointsx = np.random.random(num_spline_points) * xmax
        pointsy = np.random.random(num_spline_points) * ymax

        color = tuple(np.random.random(3))
        alpha = np.random.uniform(alpha_min, alpha_max)

        # TODO: set width limits
        width = random.randint(0, min(ymax, xmax) / 4)

        return cls(shape, pointsx, pointsy, color, alpha, width)

class DNAImage(object):
    def __init__(self, shape, strokes=None, age=0):
        self.shape = shape

        if strokes is None:
            self.strokes = []
            for _ in xrange(strokes_start):
                self.strokes.append(GeneStroke.random(self.shape))
        else:
            self.strokes = strokes

        self.age = age

        self.fitness = float('inf')

    def render(self):
    	""" Renders the image from the set of strokes and returns an RGB image
    		with dimensions self.shape.
    	"""
        w, h = self.shape
        r = RendererAgg(w, h, 72)
        arr = np.frombuffer(r.buffer_rgba(), np.uint8)
        arr.shape = r.height, r.width, -1

        for stroke in self.strokes:
            stroke.render(r)

        image_raw = np.float32(arr, dtype=np.float32) / 255
        return rgba_to_rgb(image_raw)

    def update_fitness(self, source):
    	""" Update the fitness score of the image, the MSE between the current
    		set of strokes and the source image.
    	"""
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

    def paint_hill_climbing(self):
        dna = DNAImage(self.shape)
        dna.update_fitness(self.source)
        self.population.append(dna)

    def paint(self):

        for _ in xrange(population_size):
            dna = DNAImage(self.shape)
            dna.update_fitness(self.source)
            self.population.append( (0, dna) )

        num_iter = 0
        error_change = 9999
        min_error = 999999

        while num_iter <= min_num_iter or (num_iter < max_num_iter and error_change >= epsilon):

            print num_iter

            #generate new candidates via mutation?
            new_candidates = []
            for _ in xrange(num_children):
                rand_int = random.randint(0, population_size-1)
                mutation = self.population[rand_int][1].mutate()
                mutation.update_fitness(self.source)
                new_candidates.append((0,mutation))

            # increase the age of the old population
            for i in xrange(population_size):
                self.population[i] = (self.population[i][0]+1, self.population[i][1])

            compete_population = self.population + new_candidates

            younglings = []
            aged = []
            for (age, dna) in compete_population:
                if age < dna_life_bias:
                    younglings.append((age,dna))
                else:
                    aged.append((age,dna))

            for _ in xrange(migrate_count):
                dna = DNAImage(self.shape)
                dna.update_fitness(self.source)
                younglings.append( (0, dna) )


            # pick the best out of the old children
            strong_parents = sorted(aged, key=lambda (age, dna): dna.fitness)[:population_size]

            alive = strong_parents + younglings
            alive = sorted(alive, key=lambda (age, dna): dna.fitness)


            #do some crossovers
            #new_candidates = gen_crossovers(population)

            #compete_population = population + new_candidates

            # pick the best
            #surviving_children = sorted(compete_population, key=lambda dna_strand: dna_strand.fitness)[:population_size]


            new_best_error = alive[0][1].fitness

            error_change = min_error - new_best_error
            print error_change
            min_error = new_best_error
            self.population = alive

            #print [dna.fitness for dna in self.population]
            num_iter += 1

            if num_iter > 1000 and error_change > 0:
                for z in xrange(10):
                    dump_dna(self.population[0][1].strokes, str(num_iter) + '-' + str(z) + '.txt')

        for z in xrange(10):
            dump_dna(self.population[0][1].strokes, str(z) + '.txt')

        plt.imshow(self.population[0][1].render(), interpolation='none')
        plt.show()

if __name__ == "__main__":

    source = read_image('ML257.png')
    #(a,b,c) = source.shape
    #dump_dna(DNAImage((a,b)).strokes)

    p = GenePainter(source)
    p.paint()

    """
    image = DNAImage((256, 256), strokes=[])
    stroke = GeneStroke.random((256, 256))
    image.strokes.append(stroke)

    plt.imshow(image.render(), interpolation='none')
    plt.show()

    stroke.color = stroke.random_color(stroke.color, 0.25)

    plt.imshow(image.render(), interpolation='none')
    plt.show()
    """
