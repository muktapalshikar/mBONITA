import networkx as nx
import pandas as pd
import numpy as np
import os
import pickle
from pathway_analysis_score_pathways import *
import re
import scipy.stats
import matplotlib.pyplot as plt
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
importanceScorePvalsDF.columns = ['Original Network', 'Gene', 'ImportanceScore', 'Zscore', 'Pvalue']
importanceScorePvalsDF["NegativeLog10Pvalue"] = [(-1)*np.log10(i) for i in  importanceScorePvalsDF.Pvalue]
print(importanceScorePvalsDF)
importanceScorePvalsDF.to_csv("importanceScorePvalsDF.csv")
print(scipy.stats.spearmanr(importanceScorePvalsDF.ImportanceScore, importanceScorePvalsDF.Pvalue))
print(scipy.stats.spearmanr(importanceScorePvalsDF.ImportanceScore, abs(importanceScorePvalsDF.Zscore)))
print(np.mean(importanceScorePvalsDF.Pvalue))
print(importanceScorePvalsDF[importanceScorePvalsDF.ImportanceScore==1])
print(importanceScorePvalsDF[importanceScorePvalsDF.ImportanceScore>=0.75])
print(importanceScorePvalsDF[importanceScorePvalsDF.ImportanceScore>=0.5])
threshold = (-1)*np.log10(0.1)
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context("paper")
fig, ax = plt.subplots(figsize=(6,4.5))
graph = sns.scatterplot(data=importanceScorePvalsDF, x="ImportanceScore", y="NegativeLog10Pvalue", style="Original Network", hue="Original Network", palette="Set2", s=40)
graph.axhline(threshold, linestyle="--", color='black')
graph.axvline(0.75, color='black', linestyle="--" )
graph.set_xlabel("mBONITA Importance Score" , size = 14)
graph.set_ylabel("-log10 (p-value)" , size = 14)
#print(importanceScorePvalsDF.ImportanceScore[importanceScorePvalsDF.ImportanceScore==1].tolist()[0])
plt.text(x=1-0.1,y=importanceScorePvalsDF.NegativeLog10Pvalue[importanceScorePvalsDF.ImportanceScore==1].tolist()[0], s=importanceScorePvalsDF.Gene[importanceScorePvalsDF.ImportanceScore==1].tolist()[0], size='large')
plt.text(x=1-0.1,y=importanceScorePvalsDF.NegativeLog10Pvalue[importanceScorePvalsDF.ImportanceScore==1].tolist()[1], s=importanceScorePvalsDF.Gene[importanceScorePvalsDF.ImportanceScore==1].tolist()[1], size='large')
plt.text(x=1-0.1,y=importanceScorePvalsDF.NegativeLog10Pvalue[importanceScorePvalsDF.ImportanceScore==1].tolist()[2]-0.5, s=importanceScorePvalsDF.Gene[importanceScorePvalsDF.ImportanceScore==1].tolist()[2], size='large')
plt.text(x=1-0.1,y=importanceScorePvalsDF.NegativeLog10Pvalue[importanceScorePvalsDF.ImportanceScore==1].tolist()[3], s=importanceScorePvalsDF.Gene[importanceScorePvalsDF.ImportanceScore==1].tolist()[3], size='large')
plt.text(x=1-0.1,y=importanceScorePvalsDF.NegativeLog10Pvalue[importanceScorePvalsDF.ImportanceScore==1].tolist()[4], s=importanceScorePvalsDF.Gene[importanceScorePvalsDF.ImportanceScore==1].tolist()[4], size='large')
fig = graph.get_figure()
fig.savefig("supplementary_figure_2.pdf")