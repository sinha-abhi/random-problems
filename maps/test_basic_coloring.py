from basic_coloring import greedyColor
import networkx as nx

# TODO: make unit tests
G = nx.Graph()
for i in range(6):
    G.add_node(i + 1)
G.add_edges_from([(1, 6), (1, 2), (1, 5), (2, 4), (2, 3), (3, 4)])
print(greedyColor(G, G.nodes(), modify=True))
print(G.nodes(data=True))
