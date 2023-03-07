# import necessary modules
import argparse as argparse
import operator
import networkx as nx
import pickle
from ctypes import *
import glob
from simulation import paramClass, modelClass, NPsync
from utils import genInitValueList, setupEmptyKOKI, writeModel
from GA import GAsearchModel, localSearch

# calculate importance scores
def calcImportance(individual, params, model, sss, knockoutLists, knockinLists, boolC):
    importanceScores = []  # holder for impact scores
    # loop over all nodes
    for node in range(len(model.nodeList)):
        SSEs = []  # add SSE across samples
        nodeValues = [sss[j][model.nodeList[node]] for j in range(0, len(sss))]
        for j in range(
            0, len(sss)
        ):  # knock each node out and in to observe differences
            ss = sss[j]
            initValues = list(model.initValueList[j])
            knockerOuter = list(knockoutLists[j])
            knockerOuter.append(node)  # knock out node
            boolValues1 = NPsync(
                individual,
                model,
                params.cells,
                initValues,
                params,
                knockerOuter,
                knockinLists[j],
                boolC,
            )
            knockerInner = list(knockinLists[j])
            knockerInner.append(node)  # knock in node
            boolValues2 = NPsync(
                individual,
                model,
                params.cells,
                initValues,
                params,
                knockoutLists[j],
                knockerInner,
                boolC,
            )
            # find difference between knockout and knockin
            SSE = 0
            for i in range(0, len(model.nodeList)):
                SSE += (boolValues1[i] - boolValues2[i]) ** 2
            SSEs.append(SSE)
        importanceScores.append(sum(SSEs))
    return importanceScores


if __name__ == "__main__":
    import time

    start_time = time.time()

    # read in arguments from shell scripts
    parser = argparse.ArgumentParser()
    parser.add_argument("graph")
    parser.add_argument("iterNum", nargs='?', const="1",type=str)
    results = parser.parse_args()
    graphName = results.graph
    if results.iterNum is None:
        iterNum = int(1)
    else:
        iterNum = results.iterNum
        iterNum = int(iterNum)
    name = graphName[:-8] + "_" + str(iterNum)
    with open(graphName, 'rb') as f:
        graph = pickle.load(f)
    graph= graph.to_directed()
    ruleFiles = glob.glob(name+"_rules.txt")
    if len(ruleFiles) < 1:
        # read in C function to run simulations
        updateBooler = cdll.LoadLibrary("./simulator.so")
        boolC = updateBooler.syncBool

        # load data
        sampleList = pickle.Unpickler(open(graphName[:-8] + "_sss.pickle", "rb")).load()

        # set up parameters of run, model
        params = paramClass()
        model = modelClass(graph, sampleList, False)
        model.updateCpointers()

        storeModel = [
            (model.size),
            list(model.nodeList),
            list(model.individualParse),
            list(model.andNodeList),
            list(model.andNodeInvertList),
            list(model.andLenList),
            list(model.nodeList),
            dict(model.nodeDict),
            list(model.initValueList),
        ]

        # put lack of KOs, initial values into correct format
        knockoutLists, knockinLists = setupEmptyKOKI(len(sampleList))
        newInitValueList = genInitValueList(sampleList, model)
        model.initValueList = newInitValueList

        # find rules by doing GA then local search
        model1, dev, bruteOut = GAsearchModel(
            model, sampleList, params, knockoutLists, knockinLists, name, boolC
        )  # run GA
        bruteOut1, equivalents, dev2 = localSearch(
            model1, bruteOut, sampleList, params, knockoutLists, knockinLists, boolC
        )  # run local search

        # output results
        storeModel3 = [
            (model.size),
            list(model.nodeList),
            list(model.individualParse),
            list(model.andNodeList),
            list(model.andNodeInvertList),
            list(model.andLenList),
            list(model.nodeList),
            dict(model.nodeDict),
            list(model.initValueList),
        ]
        outputList = [bruteOut1, dev, storeModel, storeModel3, equivalents, dev2]
        pickle.dump(outputList, open(name + "_local1.pickle", "wb"))  # output rules

        # calculate importance scores and output
        scores1 = calcImportance(
            bruteOut1, params, model1, sampleList, knockoutLists, knockinLists, boolC
        )
        pickle.dump(scores1, open(name + "_scores1.pickle", "wb"))

        # write rules
        with open(name + "_rules.txt", "w") as text_file:
            text_file.write(writeModel(bruteOut1, model1))
        print(("--- %s seconds ---" % (time.time() - start_time)))
