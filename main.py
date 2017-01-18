# We will need some things from several places
from __future__ import division, absolute_import, print_function
import sys
if sys.version_info < (3,):
    range = xrange
import os
from pylab import *  # for plotting
from numpy.random import *  # for random sampling
seed(42)
from parse import *


A, D, G = cPickle.load(open('data/disgenet.pickle'))
disgenet = load_graph('data/disgenet_graph.xml.gz')


def degree_plot():
    # Let's plot its in-degree distribution
    in_hist = vertex_hist(g, "in")

    y = in_hist[0]
    err = sqrt(in_hist[0])
    err[err >= y] = y[err >= y] - 1e-2

    figure(figsize=(6,4))
    errorbar(in_hist[1][:-1], in_hist[0], fmt="o", yerr=err,
            label="in")
    gca().set_yscale("log")
    gca().set_xscale("log")
    gca().set_ylim(1e-1, 1e5)
    gca().set_xlim(0.8, 1e3)
    subplots_adjust(left=0.2, bottom=0.2)
    xlabel("$k_{in}$")
    ylabel("$NP(k_{in})$")
    tight_layout()
    savefig("price-deg-dist.pdf")
    savefig("price-deg-dist.png")
    pass

if __name__ == "__main__":
    graph_draw(disgenet, output='disgenet_graph.pdf')



