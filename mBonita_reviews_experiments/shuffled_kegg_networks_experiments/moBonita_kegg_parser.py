from bioservices.kegg import KEGG
import networkx as nx
import requests
import re
import string
from bs4 import BeautifulSoup
import itertools
import argparse as argparse
import pandas as pd
import csv
import pickle

def parseKEGGdict():
    """makes a dictionary to convert ko numbers from KEGG into real gene names
    #this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two."""
    converterDict = {}
    pathwayFile = requests.get("http://rest.kegg.jp/get/br:ko00001", stream=True)
    for line in pathwayFile.iter_lines():
        line = line.decode("utf-8")
        if (
            len(line) > 1 and line[0] == "D"
        ):  # lines which begin with D translate kegg codes to gene names
            # to split into kegg code, gene names
            converter = re.split(r"\s+", re.split(r";", line)[0])
            converterDict[converter[1].upper()] = converter[2].upper()
    return converterDict


def expand_groups(node_id, groups):
    """
    node_id: a node ID that may be a group
    groups: store group IDs and list of sub-ids
    return value: a list that contains all group IDs deconvoluted
    """
    node_list = []
    if node_id in groups.keys():
        for component_id in groups[node_id]:
            node_list.extend(expand_groups(component_id, groups))
    else:
        node_list.extend([node_id])
    return node_list


def readKEGG(lines, graph, KEGGdict, hsaDict):
    # read all lines into a bs4 object using libXML parser
    soup = BeautifulSoup("".join(lines), "xml")
    groups = {}  # store group IDs and list of sub-ids
    id_to_name = {}  # map id numbers to names

    for entry in soup.find_all("entry"):
        entry_split = entry["name"].split(":")
        if len(entry_split) > 2:
            if entry_split[0] == "hsa" or entry_split[0] == "ko":
                if entry_split[0] == "hsa":
                    useDict = hsaDict
                else:
                    useDict = KEGGdict
                nameList = []
                entry_name = ""
                namer = entry_split.pop(0)
                namer = entry_split.pop(0)
                namer = namer.split()[0]
                entry_name = (
                    entry_name + useDict[namer]
                    if namer in useDict.keys()
                    else entry_name + namer
                )
                for i in range(len(entry_split)):
                    nameList.append(entry_split[i].split()[0])
                for namer in nameList:
                    entry_name = (
                        entry_name + "-" + useDict[namer]
                        if namer in useDict.keys()
                        else entry_name + "-" + namer
                    )
                entry_type = entry["type"]
            else:
                entry_name = entry["name"]
                entry_type = entry["type"]
        else:
            if entry_split[0] == "hsa":
                entry_name = entry_split[1]
                entry_type = entry["type"]
                entry_name = (
                    hsaDict[entry_name] if entry_name in hsaDict.keys() else entry_name
                )
            elif entry_split[0] == "ko":
                entry_name = entry_split[1]
                entry_type = entry["type"]
                entry_name = (
                    KEGGdict[entry_name]
                    if entry_name in KEGGdict.keys()
                    else entry_name
                )
            elif entry_split[0] == "path":
                entry_name = entry["name"]
                entry_type = "path"
            else:
                entry_name = entry["name"]
                entry_type = entry["type"]

        entry_id = entry["id"]
        entry_name = re.sub(",", "", entry_name)
        id_to_name[entry_id] = entry_name

        if entry_type == "group":
            group_ids = []
            for component in entry.find_all("component"):
                group_ids.append(component["id"])
            groups[entry_id] = group_ids
        else:
            graph.add_node(entry_name, name=entry_name, type=entry_type)

    for relation in soup.find_all("relation"):
        (color, signal) = ("black", "a")

        relation_entry1 = relation["entry1"]
        relation_entry2 = relation["entry2"]
        relation_type = relation["type"]

        subtypes = []

        for subtype in relation.find_all("subtype"):
            subtypes.append(subtype["name"])

        if (
            ("activation" in subtypes)
            or ("expression" in subtypes)
            or ("glycosylation" in subtypes)
            or ("binding/association" in subtypes)
            or ("compound" in subtypes)
            or ("phosphorylation" in subtypes)
            or ("indirect effect" in subtypes)
        ):
            color = "red"
            signal = "a"
            interaction = "a"
        elif (
            ("inhibition" in subtypes)
            or ("repression" in subtypes)
            or ("dephosphorylation" in subtypes)
            or ("dissociation" in subtypes)
            or ("ubiquitination" in subtypes)
        ):
            color = "blue"
            signal = "i"
            interaction = "i"
        else:
            print("color not detected. Signal assigned to activation arbitrarily")
            print(subtypes)
            signal = "a"
            interaction = "a"

        entry1_list = expand_groups(relation_entry1, groups)
        entry2_list = expand_groups(relation_entry2, groups)

        for (entry1, entry2) in itertools.product(entry1_list, entry2_list):
            node1 = id_to_name[entry1]
            node2 = id_to_name[entry2]
            graph.add_edge(
                node1,
                node2,
                color=color,
                subtype="/".join(subtypes),
                type=relation_type,
                signal=signal,
                interaction=interaction,
            )
    """
    for node in graph.nodes():
        if graph.degree(node)==0:
            graph.remove_node(node)

    if len(geneList) > 0:
        # 3. remove nodes with no input data
        removeNodeList = [x for x in list(graph.nodes()) if not x in geneList]
        for rm in removeNodeList:
            for start in graph.predecessors(rm):
                for finish in graph.successors(rm):
                    edge1 = graph.get_edge_data(start, rm)["signal"]
                    edge2 = graph.get_edge_data(rm, finish)["signal"]
                    inhCount = 0
                    if edge1 == "i":
                        inhCount = inhCount + 1
                    if edge2 == "i":
                        inhCount = inhCount + 1
                    if inhCount == 1:
                        graph.add_edge(start, finish, signal="i")
                    else:
                        graph.add_edge(start, finish, signal="a")
            graph.remove_node(rm)
    """
    return graph

def readFpkmData(dataName, delmited):
    """read rpm or fpkm data into format necessary for simulations and pathway analysis"""
    with open(dataName) as csvfile:
        data = []
        reader = csv.reader(csvfile, delimiter=delmited)
        for row in reader:
            data.append(row)
    sampleList = []
    geneDict = {}
    for j in range(1, len(data[1])):
        sampleList.append({})
    for i in range(1, len(data)):
        tempDatalist = []
        for j in range(1, len(data[i])):
            tempDatalist.append(float(data[i][j]))
        maxdata = max(tempDatalist)
        if maxdata == 0:
            maxdata = 1.0
        geneDict[str.upper(data[i][0])] = [
            temperDataPoint / maxdata for temperDataPoint in tempDatalist
        ]
        for j in range(0, len(data[i]) - 1):
            sampleList[j][str.upper(data[i][0])] = float(data[i][j + 1]) / maxdata
    return sampleList, geneDict

def find_pathways_kegg(
    geneList=[], preDefList=[], writeGraphml=True, organism="hsa", minimumOverlap=1
):

    """
    geneList = the list of genes included in dataset
    writeGraphml = whether or not to write out a graphml (usually true)
    organism = organism code from kegg. Eg human = 'hsa', mouse = 'mus'
    """
    koDict = parseKEGGdict()  # parse the dictionary of ko codes
    try:  # try to retrieve and parse the dictionary containing organism gene names to codes conversion
        url = requests.get("http://rest.kegg.jp/list/" + organism, stream=True)
        # reads KEGG dictionary of identifiers between numbers and actual protein names and saves it to a python dictionary
        aliasDict = {}
        orgDict = {}
        for line in url.iter_lines():
            line = line.decode("utf-8")
            line_split = line.split("\t")
            k = line_split[0].split(":")[1]
            nameline = line_split[1].split(";")
            name = nameline[0]
            if "," in name:
                nameline = name.split(",")
                name = nameline[0]
                for entry in range(1, len(nameline)):
                    aliasDict[nameline[entry].strip()] = name.upper()
            orgDict[k] = name
    except:
        print("Could not get library: " + organism)
    k = KEGG()  # read KEGG from bioservices
    k.organism = organism
    if len(preDefList) == 0:
        pathwayList = list(k.pathwayIds)
    else:
        pathwayList = list(preDefList)
    print(pathwayList)
    pathwayDict = {}
    for x in pathwayList:
        x = x.replace("path:", "")
        code = str(x)
        code = re.sub(
            "[a-zA-Z]+", "", code
        )  # eliminate org letters - retain only numbers from KEGG pathway codes
        origCode = code
        coder = str("ko" + code)  # add ko
        graph = nx.DiGraph()  # open a graph object
        # get ko pathway
        for code in [coder]:
            print(code)
            try:
                url = requests.get(
                    "http://rest.kegg.jp/get/" + code + "/kgml", stream=True
                )
                text = [line.decode("utf-8") for line in url.iter_lines()]
            except:
                print("could not read code: " + code)
                continue
            # read kegg
            graph = readKEGG(text, graph, koDict, orgDict)
        coder = str(organism + origCode)  # set up with org letters
        # get org pathway
        text = []
        for code in [coder]:
            try:
                url = requests.get(
                    "http://rest.kegg.jp/get/" + code + "/kgml", stream=True
                )
                text = []
                for line in url.iter_lines():
                    line = line.decode("utf-8")
                    text.append(line)
            except:
                print("could not read code: " + code)
                continue
            # read kegg
            graph = readKEGG(text, graph, koDict, orgDict)

        # remove complexes and rewire components
        removeNodeList = [x for x in graph.nodes() if "-" in x]
        for rm in removeNodeList:
            for start in graph.predecessors(rm):
                edge1 = graph.get_edge_data(start, rm)["signal"]
                if edge1 == "i":
                    for element in rm.split("-"):
                        graph.add_edge(start, element, signal="i")
                else:
                    for element in rm.split("-"):
                        graph.add_edge(start, element, signal="a")
            for finish in graph.successors(rm):
                edge2 = graph.get_edge_data(rm, finish)["signal"]
                if edge2 == "i":
                    for element in rm.split("-"):
                        graph.add_edge(element, finish, signal="i")
                else:
                    for element in rm.split("-"):
                        graph.add_edge(element, finish, signal="a")
            graph.remove_node(rm)
        # remove dependence of nodes on complexes that include that node
        for node in list(graph.nodes()):
            predlist = graph.predecessors(node)
            for pred in predlist:
                if "-" in pred:
                    genes = pred.split("-")
                    flag = True
                    for gene in genes:
                        if not gene in predlist:
                            flag = False
                    if flag:
                        graph.remove_edge(pred, node)

        # remove self edges
        for edge in list(graph.edges()):
            if edge[0] == edge[1]:
                graph.remove_edge(edge[0], edge[1])
        # check to see if there is a connected component, simplify graph and print if so
        allNodes = set(graph.nodes())
        if len(geneList) > 0:
            test = len(allNodes.intersection(geneList))
            removeNodeList = [x for x in list(graph.nodes()) if not x in geneList]
            for rm in removeNodeList:
                graph.remove_node(rm)
            print("Pathway: ", x, " Overlap: ", test, " Edges: ", len(graph.edges()))
        else:
            print("Gene list not provided, returning graph with all original nodes")
            print(
                "Pathway: ",
                x,
                " Nodes: ",
                len(graph.nodes()),
                " Edges: ",
                len(graph.edges()),
            )
            minimumOverlap = len(graph.nodes())
            test = minimumOverlap + 1
        nx.write_graphml(graph, coder + "_before.graphml")
        if (
            test > minimumOverlap and len(graph.edges()) > 0
        ):  # if there are at least minimumOverlap genes shared between the network and the genes in the dataset
            # nx.write_graphml(graph,coder+'_processed.graphml') # write graph out

            nx.write_graphml(graph, coder + ".graphml")
            nx.write_gpickle(graph, coder + ".gpickle")
            print(
                "nodes: ",
                str(len(graph.nodes())),
                ",   edges:",
                str(len(graph.edges())),
            )
            pathwayDict[code] = graph
            # save the removed nodes and omics data values for just those nodes in the particular pathway
            pathwaySampleList = [{} for q in range(len(geneDict[list(graph.nodes())[0]]))]
            for noder in list(graph.nodes()):
                for jn in range(len(pathwaySampleList)):
                    pathwaySampleList[jn][noder] = geneDict[noder][jn]
                pickle.dump(pathwaySampleList, open(coder + "_sss.pickle", "wb"))
    return pathwayDict

def processMetaNetwork(metaNetwork, geneList):
    
    removeNodeList = [x for x in list(metaNetwork.nodes()) if not x in geneList]
    for rm in removeNodeList:
        for start in metaNetwork.predecessors(rm):
            for finish in metaNetwork.successors(rm):
                edge1 = metaNetwork.get_edge_data(start, rm)["signal"]
                edge2 = metaNetwork.get_edge_data(rm, finish)["signal"]
                inhCount = 0
                if edge1 == "i":
                    inhCount = inhCount + 1
                if edge2 == "i":
                    inhCount = inhCount + 1
                if inhCount == 1:
                    metaNetwork.add_edge(start, finish, signal="i")
                else:
                    metaNetwork.add_edge(start, finish, signal="a")
        metaNetwork.remove_node(rm)
    metaNetwork.remove_nodes_from(
        [
            n
            for n in metaNetwork
            if n
            not in set(max(nx.connected_components(metaNetwork.to_undirected()), key=len))
        ]
    )
    return metaNetwork

def providedNetworks(coder, geneDict):
    geneList=geneDict.keys()
    graph = nx.read_graphml(coder)
    print(coder)
    graph = processMetaNetwork(graph, geneList)
    print(
    "Provided Network: ",
    coder,
    " Nodes: ",
    len(graph.nodes()),
    " Edges: ",
    len(graph.edges()),
    )
    minimumOverlap = len(graph.nodes())
    test = minimumOverlap + 1
    #nx.write_graphml(graph, coder + "_before.graphml")
    if (
        test > minimumOverlap and len(graph.edges()) > 0
    ):  # if there are at least minimumOverlap genes shared between the network and the genes in the dataset
        nx.write_graphml(graph,coder+'_processed.graphml') # write graph out
        #nx.write_graphml(graph, coder + ".graphml")
        nx.write_gpickle(graph, coder + ".gpickle")
        print(
            "nodes: ",
            str(len(graph.nodes())),
            ",   edges:",
            str(len(graph.edges())),
        )
    # save the removed nodes and omics data values for just those nodes in the particular pathway
    pathwaySampleList = [{} for q in range(len(geneDict[list(graph.nodes())[0]]))]
    for noder in list(graph.nodes()):
        for jn in range(len(pathwaySampleList)):
            pathwaySampleList[jn][noder] = geneDict[noder][jn]
        pickle.dump(pathwaySampleList, open(coder + "_sss.pickle", "wb"))

if __name__ == "__main__":
    # read in options
    parser = argparse.ArgumentParser()
    parser.set_defaults(
        sep=",", org="hsa", pathways="None", usePredefined = "False"
    )
    parser.add_argument(
        "-sep",
        "--sep",
        metavar="separator",
        help="How are columns in datafile delimited?",
    )
    parser.add_argument(
        "-t", action="store_const", const="\t", dest="sep", help="Specify a tab-delimiter"
    )
    parser.add_argument(
        "-org", "--org", metavar="org", help="Three-letter organism code to search KEGG"
    )
    parser.add_argument(
        "-paths",
        "--paths",
        dest="pathways",
        help="File with list of pathways to be considered each on one line",
    )
    parser.add_argument("-data", "--data", help="Delimited data file with columns = samples and rows = genes")
    parser.add_argument('-usePredefined', '--usePredefined', help = "use graphml files in the current working directory")
    results = parser.parse_args()
    dataName = results.data
    org = results.org
    paths = results.pathways
    usePredefined = results.usePredefined
    sss, geneDict = readFpkmData(dataName, results.sep)  # read in data
    pickle.dump( sss, open( 'sss.pickle', "wb" ) ) # save data in correct format for runs
    if paths == "None" and usePredefined == "False":
        find_pathways_kegg(
            geneList=geneDict.keys(),
            preDefList=[],
            organism=org
        )
    else:
        if paths != "None" and usePredefined == "False":
            with open(paths, "r") as inputfile:
                lines = inputfile.readlines()
            pathList = []
            for line in lines:
                for element in line.split(","):
                    pathList.append(element.strip())
            find_pathways_kegg(
                geneList=geneDict.keys(),
                preDefList=pathList,
                organism=org
            )
        else:
            if usePredefined == "True":
                from glob import glob
                graphs = tuple(glob("*.graphml"))
                for coder in graphs:
                    #print(nx.get_edge_attributes(nx.read_graphml(coder), "signal"))
                    providedNetworks(coder, geneDict)


"""
##########TESTS##########

import pandas as pd
import networkx as nx

data = pd.read_csv("phospho_LSP1.csv")
data.head()
geneList = data.Gene.tolist()
print("Genes in Dataset: ", len(geneList))
graphs = find_pathways_kegg(
    geneList=geneList,
    preDefList=[
        "hsa04010",
        "hsa04062",
        "hsa04064",
        "hsa04066",
        "hsa04150",
        "hsa04151",
        "hsa04370",
        "hsa04625",
        "hsa04514",
        "hsa04630",
        "hsa04668",
        "hsa04670",
        "hsa04810",
    ],
    writeGraphml=True,
    organism="hsa",
)
metaNetwork = nx.compose_all(list(graphs.values()))
# metaNetwork.remove_nodes_from([n for n in metaNetwork if n not in set(geneList)])
removeNodeList = [x for x in list(metaNetwork.nodes()) if not x in geneList]
for rm in removeNodeList:
    for start in metaNetwork.predecessors(rm):
        for finish in metaNetwork.successors(rm):
            edge1 = metaNetwork.get_edge_data(start, rm)["signal"]
            edge2 = metaNetwork.get_edge_data(rm, finish)["signal"]
            inhCount = 0
            if edge1 == "i":
                inhCount = inhCount + 1
            if edge2 == "i":
                inhCount = inhCount + 1
            if inhCount == 1:
                metaNetwork.add_edge(start, finish, signal="i")
            else:
                metaNetwork.add_edge(start, finish, signal="a")
    metaNetwork.remove_node(rm)
metaNetwork.remove_nodes_from(
    [
        n
        for n in metaNetwork
        if n
        not in set(max(nx.connected_components(metaNetwork.to_undirected()), key=len))
    ]
)
origGraph = nx.read_graphml("gephi_edited_hif1Agraph_rules.graphml")
origGraph.remove_nodes_from([n for n in origGraph if n not in set(geneList)])
# origGraph.remove_nodes_from([n for n in origGraph if n not in set(max(nx.connected_components(origGraph.to_undirected()), key=len))])
print("Genes in metaNetwork: ", len(metaNetwork))
print("Genes in origGraph: ", len(origGraph))
overlap = set(metaNetwork.nodes()).intersection(set(origGraph.nodes()))
print("Overlap: ", len(overlap))
# print(["GRK5", "GRK5" in overlap, "CXCR4", "CXCR4" in overlap, "SYK", "SYK" in overlap, "STAT3", "STAT3" in overlap, "PLCG2", "PLCG2" in overlap, "CD19", "CD19" in overlap, "PDPK1", "PDPK1" in overlap, "BLNK", "BLNK" in overlap, "MTOR", "MTOR" in overlap, "CSNK2B", "CSNK2B" in overlap, "PRKCG", "PRKCG" in overlap, "PHLPP1", "PHLPP1" in overlap, "MAPKAPK2", "MAPKAPK2" in overlap, "LSP1", "LSP1" in overlap ])
# print(sum([origNode in metaNetwork.nodes() for origNode in origGraph.nodes()]))
for origNode in origGraph.nodes():
    if not origNode in metaNetwork.nodes():
        print(origNode)

######

find_pathways_kegg(geneList=["IL6", "IL6R", "STAT3","TLR4","RELA","IFNG","NFKB1","ERBB2"], preDefList = ["hsa04010"], writeGraphml=True,  organism="hsa")

origGraph = nx.read_graphml("hsa04010_orig.graphml")
print(len(origGraph.nodes()))
print(len(origGraph.edges()))
origNodes = list(origGraph.nodes())
origNodes.sort()
print(origNodes)
print("*********")
newGraph = nx.read_graphml("hsa04010_processed.graphml")
print(len(newGraph.nodes()))
print(len(newGraph.edges()))
newNodes = list(newGraph.nodes())
newNodes.sort()
print(newNodes)
print("*********")

#print(set(list(origGraph.edges())).intersection(set(list(newGraph.edges()))))
#print(set(origNodes).intersection(newNodes))
print(len(set(origNodes).intersection(newNodes)))

print(set(origNodes).difference(newNodes))
print(len(set(origNodes).difference(newNodes)))
"""
