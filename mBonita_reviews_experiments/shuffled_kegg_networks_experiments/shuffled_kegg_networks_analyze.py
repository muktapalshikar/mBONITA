import networkx as nx
import pandas as pd
import numpy as np
import os
import pickle
from pathway_analysis_score_pathways import *
import re
import scipy.stats

# finds pathways that should be compared
def findPathwayList_noReps():
    pathways = []
    codes = []
    # looks for pathways that have gpickles generated originally
    for file in os.listdir("gpickles"):
        if file.endswith(".gpickle"):
            if os.path.isfile("pickles/" + file[:-8] + "_1_local1.pickle"):
                codes.append(file[:-8])
            else:
                print((file[:-8] + " has no output"))
    #print(codes)
    # for each of these pathways, we find the output of the rule determination and scoring procedures and put them together.
    for code in codes:
        pathVals = []
        rules = []
        [
            bruteOut1,
            dev,
            storeModel,
            storeModel3,
            equivalents,
            dev2,
        ] = pickle.Unpickler(
            open("pickles/" + code + "_" + "1_local1.pickle", "rb")
        ).load()
        model = modelHolder(storeModel3)
        pathVals.append(
            pickle.Unpickler(
                open("pickles/" + code + "_" + "1_scores1.pickle", "rb")
            ).load()
        )
        rules.append(writeModel(bruteOut1, model))
        #print(pathVals)
        #graph = nx.read_gpickle("gpickles/" + code + ".gpickle")
        graph = nx.read_graphml("graphmls/" + code + ".graphml")
        ImportanceVals = {}  # average importance vals over trials
        i=0
        for node in range(len(storeModel[1])):
            ImportanceVals[storeModel[1][node]] = pathVals[i][node]
        # add nodes removed during network simplification back in
        #print(ImportanceVals)
        pathways.append([code, ImportanceVals, rules, graph])
    return pathways
pathList = findPathwayList_noReps() # identify pathways under consideration


pathImportances_complete = []
for i in range(0, len(pathList)): # iterate over pathways
    pathway = pathList[i]
    pathImportances = pd.DataFrame.from_dict(pathway[1], orient = 'index')
    pathImportances['Pathway'] = pathway[0]
    pathImportances[0] = pathImportances[0]/max(pathImportances[0].values)
    if len(pathImportances_complete) == 0:
        pathImportances_complete = pathImportances
    else:
        pathImportances_complete = pd.concat([pathImportances_complete , pathImportances])

pathImportances_complete.columns = ['ImportanceScore', 'Pathway']
pathImportances_complete = pathImportances_complete.reset_index()
pathImportances_complete.columns = ['Gene', 'ImportanceScore', 'Pathway']
pathImportances_complete.to_csv("ImportanceScores.csv")

pathImportances_complete['OriginalNetwork'] = [path.replace("shuffled_","") for path in pathImportances_complete.Pathway]
pathImportances_complete.OriginalNetwork = [re.sub("_\d+","", path) for path in pathImportances_complete.OriginalNetwork]
pathImportances_complete['IsShuffled'] = [str("shuffled" in path) for path in pathImportances_complete.Pathway]

importanceScorePvals = {}
importanceScorePvalsDF = pd.DataFrame()

for path in set(pathImportances_complete.OriginalNetwork):
    importanceScorePvals[path] = {}
    for gene in set(pathImportances_complete.Gene[pathImportances_complete.OriginalNetwork==path]):
        importanceScorePvals[path][gene] = {}
        importanceScores = pathImportances_complete[pathImportances_complete.OriginalNetwork==path]
        importanceScores = importanceScores[importanceScores.Gene==gene]
        score = importanceScores.ImportanceScore[importanceScores.IsShuffled=="False"].tolist()[0]
        importanceScores = importanceScores['ImportanceScore'].tolist()
        meanScore = np.mean(importanceScores)
        stdDev = np.std(importanceScores)
        zscore = (score - meanScore) / stdDev
        pval = scipy.stats.norm.sf(zscore)
        importanceScorePvals[path][gene]['zscore'] = zscore
        importanceScorePvals[path][gene]['pvalue'] = pval
        importanceScorePvalsDF = importanceScorePvalsDF.append(pd.Series([path, gene, score, zscore, pval]), ignore_index=True)
#print(pd.DataFrame.from_dict(importanceScorePvals, orient='columns'))
importanceScorePvalsDF.columns = ['OriginalNetwork', 'Gene', 'ImportanceScore', 'Zscore', 'Pvalue']
print(importanceScorePvalsDF)
print(scipy.stats.pearsonr(importanceScorePvalsDF.ImportanceScore, importanceScorePvalsDF.Pvalue))
print(scipy.stats.pearsonr(importanceScorePvalsDF.ImportanceScore, abs(importanceScorePvalsDF.Zscore)))
print(np.mean(importanceScorePvalsDF.Pvalue))



