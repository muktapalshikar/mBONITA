# %% [markdown]
# # Figure 4: moBonita identifies mechanisms of hypoxia-mediated chemotaxis in RAMOS B cells 
# 
# * (a) Reduction of ERS – show ERS sizes of all three datasets on one plot 
# * (b) Highly modulated nodes from the dataset – modulation score vs expression scatterplot 
# * (c ) Pathway analysis from all three datasets on one plot 
# * (d ) Subnetworks of influential genes from PKNs 

# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import statsmodels.stats.multitest
import os
import re
import requests
from ast import literal_eval
from glob import glob
from numpy import arange
import deap 

def getPathwayName(hsaURL):
    """Use KEGG API to get the readable name using the KEGG code (eg. hsa00010 corresponds to the glycolysis pathway)"""
    fileReg = re.compile("NAME\s+(\w+.*)")
    pathwayFile = requests.get("http://rest.kegg.jp/get/" + hsaURL, stream=True)
    for line in pathwayFile.iter_lines():
        line = line.decode("utf-8")
        result = fileReg.match(line)
        if result:
            return result.group(1)
    return hsaURL

# %% [markdown]
# # Collect the *local1.pickle* files in the specified directory (change directory variable below)

# %% [markdown]
# directory = "/gpfs/fs2/scratch/mpalshik/multiomics_networks_2022/BONITA_experiments/"
# outputfiles = []
# for root, dirs, files in os.walk(directory):
#     for f in files:
#         if f.endswith("_local1.pickle"):
#             outputfiles.append(os.path.join(root, f))
# print(len(outputfiles), outputfiles[0:5])

# %% [markdown]
# # Open local1.pickle files and process the information into a single dataframe
# **One row in the dataframe contains information for one node. The dataframe has the following columns:**
#  - Network name - readable, descriptive KEGG network name
#  - Method name - subfolder of the main directory in which the pickle was found
#  - andNodeList - indices of parent nodes
#  - andNodeInvertList - a bitstring encoding the activation and inhibition edges. True implies that the edge from the corresponding parent node in the andNodeList is an inhibitory edge
#  - ruleLengths - length (ie, size) of the ERS for the node
#  - equivs - bitstring representation of the equivalent rule set
#  - plainRules - plain text representation of the rules in the ERS
#  - randomERSIndividual - random individual from the ERS
#  - minLocalSearchError - lowest error for the rules tried for each node

# %%
import pickle
from utils import *
from networkConstructor import *
from simulation import *
from pathway_analysis_score_pathways import *
from pathway_analysis_score_nodes import *
from GA import *
# from pathway_analysis_setup import *
import glob

def findEnds2(model, node, indiv):
    """ find the end of a node in the bitstring """
    node = model.nodeList.index(node)
    if node == len(model.nodeList) - 1:
        end1 = len(indiv)
    else:
        end1 = model.individualParse[node + 1]
    start1 = model.individualParse[node]
    return start1, end1

""" directory = "/gpfs/fs2/scratch/mpalshik/multiomics_networks_2022/BONITA_experiments/"
fileReg = re.compile('.*hsa(\d+)\_.+\.pickle')
seriesIndices = ["networkName", "methodName", "nodeList", "andNodeList", "andNodeInvertList", "ruleLengths", "equivs", "plainRules", "randomERSIndividual", "minLocalSearchError"]
df = pd.DataFrame(columns = seriesIndices)
i = 1
counter = 1
fileReg2 = re.compile(re.escape(directory) + r"(\w+.*)")
outputfiles = glob.glob("/gpfs/fs2/scratch/mpalshik/multiomics_networks_2022/BONITA_experiments/"+"*/pickles/*local1.pickle")

print(len(outputfiles))

for f in outputfiles:
    if (counter % 20 == 0):
        print(counter)
    counter = counter + 1
    getMethod=fileReg2.match(f)
    if getMethod:
        methodName=getMethod.group(1)
    else:
        methodName="N.A."
    result = fileReg.match(f)
    networkName = getPathwayName('hsa'+result.group(1))
    print(f)
    outputList = pickle.load(open(f, mode = "rb"))
    print(len(outputList))
    bruteOut1,dev,storeModel, storeModel3, equivalents, dev2 = [outputList[k] for k in range(len(outputList))]
    randomERSIndividual = bruteOut1 #random individual from the ERS
    minLocalSearchError = dev2 #lowest error for the rules tried for each node
    #equivalents = ERS for the network
    #storeModel = model from the GA
    minGAErrors = dev # minimum errors returned by the GA
    model1=modelHolder(storeModel3) #model from the local search
    for node in range(0,len(model1.nodeList)):
        plainRules=[]
        start1,end1=findEnds2(model1, model1.nodeList[node], equivalents[node])
        ers=equivalents[node] # find the bitstring for just this node
        #inEdges=findInEdges(model1, model1.nodeList.index(model1.nodeList[node]))
        for rule in ers:
            plainRules.append(writeNode(model1.nodeList.index(model1.nodeList[node]), rule, model1))
        ruleLengths=len(ers)
        ersAllNodes=plainRules
        s = pd.Series([networkName, methodName, model1.nodeList[node], model1.andNodeList[node], model1.andNodeInvertList[node], ruleLengths, str(ers), plainRules, randomERSIndividual, minLocalSearchError[node]], index = seriesIndices)
        df.loc[i] = s
        i = i + 1
df['methodName'] = df['methodName'].str.extract(r'(\w+)\/', expand=False)
#df['indegree'] = [len(set([item for sublist in literal_eval(i) for item in sublist])) for i in df.andNodeList]
CSVfile = "local1Data.csv"
dfCSVFile = open(CSVfile, mode = "w")
df.to_csv(dfCSVFile)
dfCSVFile.close()
 """

# %%
df = pd.read_csv("local1Data.csv", index_col=0)
df['indegree'] = [
    len(set([item for sublist in literal_eval(i) for item in sublist]))
    for i in df.andNodeList
]
df = df.replace('concatenated', 'Concatenated')
df.head()

# %%
sns.set_theme(context='talk', style='ticks', rc={'figure.figsize': (1, 1)})
g = sns.displot(
    data=df[df.indegree >= 3],
    x="ruleLengths",
    #multiple="stack",
    row="methodName",
    row_order = ["Proteomics", "Transcriptomics", "Phosphoproteomics", "Concatenated"], 
    hue="methodName",
    palette="colorblind",
    legend=False,
    stat="count",
    facet_kws={
        'sharey': False,
        'sharex': True
    })
g.set_axis_labels("Size of ERS for\nnodes with in-degree >= 3", "Count")
g.set_titles("{row_name}")
g.figure.figsize = (1, 1)
g.savefig("figure3a.png", dpi=300)
g.savefig("figure3a.pdf")
plt.clf()
# %%
sns.set_theme(context='talk', style='ticks', rc={'figure.figsize': (1, 1)})
g = sns.catplot(data=df[df.indegree >= 3],
                x="ruleLengths",
                y="methodName",
                hue="methodName",
                palette="colorblind",
                legend=False)
g.set_axis_labels("Size of ERS for\nnodes with in-degree >= 3", "Dataset")
#g.set_titles("{col_name}")

# %%
transcriptomics = df[df['methodName'].str.contains('Transcriptomics')]
print(transcriptomics.shape)
phosphoproteomics = df[df['methodName'].str.contains('Phosphoproteomics')]
print(phosphoproteomics.shape)
proteomics = df[df['methodName'].str.match(
    'Proteomics'
)]  #note difference between contains and match - DO NOT CHANGE THIS
print(proteomics.shape)

# %%
# Get optimized networks
def maxERS(inDegree):
    if inDegree == 0:
        maxERS = 1
    else:
        inDegree = min(inDegree, 3)
        if inDegree == 2:
            maxERS = 15
        else:
            maxERS = 127
    return maxERS


df.loc[:, 'maxERS'] = [maxERS(i) for i in df.indegree]
df.loc[:, 'hasReducedERS'] = df.loc[:, 'ruleLengths'] < df.loc[:, 'maxERS']

# %% [markdown]
# # Compare importance scores
# 

# %% [markdown]

# nodeTable = pd.DataFrame(index=["temp"],
#                          columns=[
#                              "andNode", "Display Name", "IS", "name", "RA",
#                              "selected", "shared name", "Dataset", "Contrast",
#                              "Network"
#                          ])
# for file in glob.glob(
#         "/gpfs/fs2/scratch/mpalshik/multiomics_networks_2022/BONITA_experiments/*/*/*rules.graphml"
# ):
#     temp_df1 = pd.DataFrame(index=["temp"],
#                             columns=[
#                                 "andNode", "Display Name", "IS", "name", "RA","Dataset",
#                                 "Contrast", "Network"
#                             ])
#     temp_G = nx.read_graphml(file)
#     raw_table = list(temp_G.nodes.data())
#     for node in raw_table:
#         data = dict(node[1])
#         temp_df2 = pd.DataFrame(data, index=[0])
#         temp_df1 = pd.concat((temp_df1, temp_df2), axis=0)
#     temp_df1 = temp_df1.iloc[1:]
#     temp_df1 = temp_df1[temp_df1["andNode"] == 0]
#     temp_df1.drop(temp_df1.columns[[0, 3, 5, 6]], axis=1, inplace=True)
#     temp_df1.index = arange(1, len(temp_df1) + 1)
#     filename = filepath = ''
#     parsed1 = file.split('/')
#     temp_df1["Dataset"] = parsed1[7]
#     temp_df1["Contrast"] = parsed1[8]
#     temp_df1["Network"] = parsed1[9].replace("_rules.graphml", "")
#     #parsed2 = parsed1[3:]
#     #parsed3 = parsed2[0].split(
#     #    '_')[0][:2] + '_' + parsed2[1] + '_' + parsed2[2][:-14]
#     #parsed4 = parsed3.replace('percent', '')
#     #filename = parsed4.replace('-', '_')
#     #filepath = '/'.join(parsed1[:3]) + "/Node_Tables/"
#     #destination = filepath + filename
#     #print(filename)
#     #temp_df1.to_csv(destination)
#     nodeTable = pd.concat([nodeTable, temp_df1])
# nodeTable = nodeTable.dropna(axis=0, how='all')
# nodeTable.to_csv("node_table.csv")
# nodeTable

# %%
nodeTable = pd.read_csv("node_table.csv", index_col = 0)
nodeTable = nodeTable.replace('concatenated', 'Concatenated')
# %%
nt2 = nodeTable[['Display Name', 'IS', 'Dataset', 'Network']].drop_duplicates()
nt2 = nt2.pivot_table(values="IS", index = ["Display Name", "Network"], columns = "Dataset")
nt2 = nt2.dropna(axis = 0, how = 'any')
nt2

# %%
sns.set_context("paper")
g = sns.pairplot(nt2,
            x_vars = ["Proteomics", "Transcriptomics", "Phosphoproteomics", "Concatenated"],
                 y_vars = ["Proteomics", "Transcriptomics", "Phosphoproteomics", "Concatenated"],
             diag_kind='hist',
             kind='reg',
             plot_kws={
                 'scatter_kws': {
                     'alpha': 0.25,
                     'color': "#005AB5"
                 },
                 'line_kws': {
                     'color': '#DC3220'
                 }
             }, diag_kws = {'color': '#005AB5'},
             corner=False)
g.savefig("relation_between_importance_scores.png", dpi = 300)
g.savefig("relation_between_importance_scores.pdf")

# %%
nt2.corr(method="pearson").to_csv("pearson_importance_score.csv")
nt2.corr(method="spearman").to_csv("spearman_importance_score.csv") # better to choose spearman

# %%
commonNodesRules = df[df.nodeList.isin(
    nt2.index.get_level_values('Display Name').tolist())]
netid_to_name = {
    key: value
    for (key, value) in
    zip(set(nt2.index.get_level_values('Network')),
        [getPathwayName(f) for f in set(nt2.index.get_level_values('Network'))])
}
nt2['networkName'] = [netid_to_name[i] for i in nt2.index.get_level_values('Network')]
commonNodesRules = commonNodesRules[commonNodesRules.networkName.isin(
    nt2.networkName)].reset_index(drop=True)
commonNodesRules.head()

# %%
cnr_transcriptomics = commonNodesRules[
    commonNodesRules['methodName'].str.contains('Transcriptomics')]
print(cnr_transcriptomics.shape)
cnr_phosphoproteomics = commonNodesRules[
    commonNodesRules['methodName'].str.contains('Phosphoproteomics')]
print(cnr_phosphoproteomics.shape)
cnr_proteomics = commonNodesRules[commonNodesRules['methodName'].str.match(
    'Proteomics'
)]  #note difference between contains and match - DO NOT CHANGE THIS
print(cnr_proteomics.shape)
cnr_concatenated = commonNodesRules[
    commonNodesRules['methodName'].str.contains('Concatenated')]
print(cnr_concatenated.shape)

tallyDF_final = commonNodesRules
tallyDF_final["percentOverlap"] = [np.nan for i in range(0, tallyDF_final.shape[0])]
print(tallyDF_final.shape)

def percentOverlap2(equivs, exactCompare=False):
    """Compare n lists to get a percentage overlap compared to the "universe" (i.e. the union of all n lists)"""
    overlap = set()
    universe = set()
    for l in equivs:
        rule = set([tuple(i) for i in l])
        universe = universe | rule
        if len(overlap) == 0:
            overlap = rule
        else:
            overlap = rule & overlap
    if len(universe) == 0:
        return(np.nan)
    else:
        if exactCompare:
            percentOverlap = float(len(overlap)) / len(universe) *100
        else:
            if len(overlap) > 0:
                return(float(1))
            else:
                return(float(0))
    return(percentOverlap)

def getNodeMinimalRule(ruleSet, andNodeList):
    numOrRules = [sum(ruleSet[j]) for j in range(len(ruleSet))]
    deciderVar = max(numOrRules)  # maximal or rules
    maxOrRules = ruleSet[
        numOrRules.index(deciderVar)
    ]  # rules with maximum or terms
    maxUpstreamNodes = 0
    minimalRule = []
    for orRule in [maxOrRules]:
        if sum(orRule) > 0:
            numUpstreamNodes = [
                andNodeList[orTerm]
                for orTerm in range(len(orRule))
                if orRule[orTerm] == 1
            ]
        else:
            minimalRule = orRule
            continue
        numUpstreamNodes = [len(element) for element in numUpstreamNodes]
        numUpstreamNodes = sum(numUpstreamNodes)
        if numUpstreamNodes > maxUpstreamNodes:
            maxUpstreamNodes = numUpstreamNodes
            minimalRule = orRule
        else:
            maxUpstreamNodes = maxUpstreamNodes
            minimalRule = minimalRule
    return minimalRule

for network in list(set(commonNodesRules.networkName))[0:15]:
    print("############")
    print(network.upper())
    print("######\n")
    net_trans = cnr_transcriptomics[
        cnr_transcriptomics['networkName'].str.contains(network, regex=False)]
    net_phosph = cnr_phosphoproteomics[
        cnr_phosphoproteomics['networkName'].str.contains(network,
                                                          regex=False)]
    net_prot = cnr_proteomics[cnr_proteomics['networkName'].str.contains(
        network, regex=False)]
    allNodes = pd.concat(
        [net_trans.nodeList, net_phosph.nodeList, net_prot.nodeList])
    for node in set(allNodes):
        tempTrans = pd.DataFrame(net_trans[net_trans["nodeList"] == node])
        tempProt = pd.DataFrame(net_prot[net_prot["nodeList"] == node])
        tempPhosph = pd.DataFrame(net_phosph[net_phosph["nodeList"] == node])
        #print(tempTrans.shape, tempPhosph.shape, tempProt.shape)
        """
        percentOverlap = percentOverlap2([
            tempTrans.plainRules.tolist(),
            tempProt.plainRules.tolist(),
            tempPhosph.plainRules.tolist()
        ])
        if percentOverlap != 1:
            print(node, network)
            print([
                set(tempTrans.plainRules.tolist()),
                set(tempProt.plainRules.tolist()),
                set(tempPhosph.plainRules.tolist())
            ], percentOverlap)
        """
        minimalTrans = getNodeMinimalRule(
            literal_eval(tempTrans.equivs.iloc[0]),
            literal_eval(tempTrans.andNodeList.iloc[0]))
        minimalTrans = literal_eval(tempTrans.plainRules.iloc[0])[literal_eval(
            tempTrans.equivs.iloc[0]).index(minimalTrans)]
        minimalProt = getNodeMinimalRule(
            literal_eval(tempProt.equivs.iloc[0]),
            literal_eval(tempProt.andNodeList.iloc[0]))
        minimalProt = literal_eval(tempProt.plainRules.iloc[0])[literal_eval(
            tempProt.equivs.iloc[0]).index(minimalProt)]
        minimalPhosph = getNodeMinimalRule(
            literal_eval(tempPhosph.equivs.iloc[0]),
            literal_eval(tempPhosph.andNodeList.iloc[0]))
        minimalPhosph = literal_eval(tempPhosph.plainRules.iloc[0])[literal_eval(
            tempPhosph.equivs.iloc[0]).index(minimalPhosph)]
        print(node)
        print("Trans: ", minimalTrans, "\nProt:", minimalProt, "\nPhosph:",
              minimalPhosph)
        print("######\n")

# # %%
minimalTrans = "RB1*=(CDKN1A) or (CDK2 and CDKN1A) or (CDK2)"
minimalProt = "RB1*=(GRB2) or (TP53 and CDK2 and GRB2) or (TP53 and CDK2) or (TP53 and GRB2) or (TP53) or (CDK2 and GRB2) or (CDK2)"
minimalPhosph = "RB1*=(not PTK2B) or (not TP53 and  not MAP3K7 and  not PTK2B) or (not TP53 and  not MAP3K7) or (not TP53) or (not MAP3K7 and  not PTK2B) or (not MAP3K7)"

# simple intersection
def makeRuleList(rule):
    line=rule
    line=line.strip()
    matcher=re.search("(.*)\*=(.*)", line)
    node=matcher.group(1).strip()
    rule=matcher.group(2).strip()
    preds=re.sub("and|not|or|\(|\)", "", rule)
    preds=preds.strip()
    preds=re.sub("\s+", ",", preds)
    preds=tuple(preds.split(","))
    rule=rule.split(" or ")
    return rule, preds
    
transRule, transPreds = makeRuleList(minimalTrans)
protRule, protPreds = makeRuleList(minimalProt)
phosphRule, phosphPreds = makeRuleList(minimalPhosph)

print(set(transRule),set(protRule),set(phosphRule))
print(set(transRule) & set(protRule) & set(phosphRule))

allPreds = [transPreds, protPreds, phosphPreds]
allRules = [transRule, protRule, phosphRule]
lenAllPreds = [len(i) for i in allPreds]
maxPreds = allPreds[lenAllPreds.index(max(lenAllPreds))]
maxRule = allRules[lenAllPreds.index(max(lenAllPreds))]

# a few cases arise:
# case 1: the largest set of predictors is a superset of the others. In this case the minimal rule set with the highest number of predictors is the accepted rule (should this be lowest??)
if len([set(i).issubset(set(maxPreds)) for i in allPreds]) == len(maxPreds):
    finalRule = maxRule
# case 2: the sets of predictors are all disjoint. In this case -
# a. create a list of all predictors
# b. find a set of rules that includes all predictors - iteratively add an or term to the final rule
pooledPreds = set([x for xs in [transPreds,phosphPreds,protPreds] for x in xs])
# intersections between preds from each experiment and pooled predictors








# # %%
# from utils import Get_expanded_network
# Get_expanded_network([literal_eval(t)[0] for t in net_trans.plainRules.tolist()], equal_sign="*=").nodes


# # %%
# for network in list(set(commonNodesRules.networkName)):
#     net_trans = cnr_transcriptomics[cnr_transcriptomics['networkName'].str.contains(network, regex=False)]
#     net_phosph = cnr_phosphoproteomics[cnr_phosphoproteomics['networkName'].str.contains(network, regex=False)]
#     net_prot = cnr_proteomics[cnr_proteomics['networkName'].str.contains(network, regex=False)]
#     print(net_trans.shape)
#     for node in set(df.nodeList):
#         tempTrans = pd.DataFrame(net_trans[net_trans["nodeList"] == node])
#         tempProt = pd.DataFrame(net_prot[net_prot["nodeList"] == node])
#         tempPhosph = pd.DataFrame(net_phosph[net_phosph["nodeList"] == node])
#         if (net_trans.shape[0] >= 0 and net_prot.shape[0] >= 0 and net_phosph.shape[0]>= 0):
#             percentOverlap = percentOverlap2([tempTrans.plainRules.tolist(), tempProt.plainRules.tolist(), tempPhosph.plainRules.tolist()])
#             temp2 = [t1 & t2 for t1, t2 in zip(list(tallyDF_final["networkName"]== network), list(tallyDF_final["nodeList"] == node))]
#             tallyDF_final.loc[temp2, "percentOverlap"] = percentOverlap

# # %%
# tallyDF_final.head()
# tallyDF_final.to_csv("tallyDF.csv")
# """Aggregate by network and show average percent overlap"""
# temp1 = pd.DataFrame(tallyDF_final.groupby(["networkName"]).mean())
# temp1.head()


