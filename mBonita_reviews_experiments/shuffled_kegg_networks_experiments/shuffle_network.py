import networkx as nx
from random import seed, choices
from glob import glob
from pathway_analysis_setup import *

seed(1)

# shuffle metanetwork
#net = nx.read_graphml('metaNetwork.graphml')
#print(len(net.edges()))

def shuffleNet (net):
    net2 = nx.double_edge_swap(net.to_undirected(), nswap=len(net.edges())/2, max_tries=len(net.edges())*100)
    net2 = net2.to_directed()
    signal = dict(zip(net2.edges(), choices(["a", "i"],k=len(net2.edges()))))
    nx.set_edge_attributes(net2, values = signal, name = "signal")
    return net2

#for i in range(0,101):
#    fileName = "shuffled_metaNetwork_" + str(i) + ".graphml"
#    nx.write_graphml(shuffleNet(net), fileName)
#    print(fileName)

# shuffle 5 kegg networks
# selected KEGG networks:

nets = tuple(glob("graphmls/hsa*.graphml"))
print(nets)
for net in nets:
    for i in range(6, 51):
        print(net)
        print(str(net[:-8]))
        fileName = "shuffled_" + str(net[:-8].replace("graphmls/", "")) + "_" + str(i)+ ".graphml"
        net2 = nx.read_graphml(net)
        nx.write_graphml(shuffleNet(net2), fileName)
        print(fileName)


ss, geneDict, cvDict = readFpkmData("concatenated_datasets.csv", ",")  # read in data
nets = tuple(glob("*.graphml"))
print(nets)
for net in nets:
    retrieveGraph_customGraph(geneDict, net) # generate gpickles needed for pathway analysis  


