import networkx as nx
from random import seed, choices
seed(1)

net = nx.read_graphml('metaNetwork.graphml')
print(len(net.edges()))

def shuffleNet (net):
    net2 = nx.double_edge_swap(net.to_undirected(), nswap=len(net.edges())/2, max_tries=len(net.edges())*100)
    net2 = net2.to_directed()
    signal = dict(zip(net2.edges(), choices(["a", "i"],k=len(net2.edges()))))
    nx.set_edge_attributes(net2, values = signal, name = "signal")
    return net2

for i in range(0,101):
    fileName = "shuffled_metaNetwork_" + str(i) + ".graphml"
    nx.write_graphml(shuffleNet(net), fileName)
    print(fileName)

