"""
basic-coloring.py
Abhi Sinha
October 08, 2019

Implementations of naïve algorithms for general graph coloring.
"""
import networkx as nx


# Greedy coloring
def firstColor(colors):
    ncolors = len(colors) + 1
    count = [0] * ncolors
    for c in colors:
        if c < ncolors:
            count[c] += 1
    for c in range(ncolors):
        if count[c] == 0:
            return c


def greedyColor(G: nx.Graph(), ordering, modify=False) -> dict:
    color = dict()
    for v in ordering:
        color[v] = firstColor([color[w] for w in G.neighbors(v) if w in color])
    if modify:
        for v in G.nodes():
            G.node[v]["color"] = color[v]
    return color
