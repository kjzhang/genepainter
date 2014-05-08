import sys
import matplotlib.pyplot as plt
from genepainter import DNAImage, load_dna

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: %s dump.txt" % sys.argv[0])
        sys.exit(-1)

    strokes = load_dna(sys.argv[1])
    dna = DNAImage(strokes[0].shape, strokes)
    plt.imshow(dna.render(), interpolation='none')
    plt.show()
