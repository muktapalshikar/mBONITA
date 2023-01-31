from pathway_analysis_score_pathways import findPathwayList, retrievePathKey
import pandas as pd
from numpy import count_nonzero
#import argparse
import numpy as np
from scipy.stats import norm as norm

"""
# do pathway analysis! store in one folder for each comparison
def analyze_pathways_mBonita(contrastName, conditionName, dataName, delimited):
    data = pd.read_csv(dataName, sep = delimited, index_col = 0).T
    contrasts = pd.read_csv(contrastName, sep = delimited, header = None)
    conditionMatrix = pd.read_csv(conditionName, sep = delimited, index_col=0)
    # check if 'Dataset' column exists in conditionMatrix
    if not 'Dataset' in conditionMatrix.columns:
        print("Condition matrix must have a column 'Dataset' indicating the original omics dataset")
        return None
    elif not 'Condition' in conditionMatrix.columns:
        print("Condition matrix must have a column 'Condition' indicating the condition of the sample.")
        return None
    else:
        groupedData = data.groupby([conditionMatrix['Dataset'], conditionMatrix['Condition']], axis = 0)
        meanAbundance = groupedData.mean()
        CV = data.groupby(conditionMatrix['Dataset'], axis = 0).std()
        score = data.groupby(conditionMatrix['Dataset']).agg(lambda x: count_nonzero(x)).agg(lambda x: count_nonzero(x))
        print(CV)
        print(meanAbundance)
        print(score)
        pathList = findPathwayList() # identify pathways under consideration
        #for i in range(0, len(contrasts)): # iterate over contrasts
        i = 0
        print(i)
        condition1 = contrasts.iloc[i,0]
        condition2 = contrasts.iloc[i,1]
        #print(condition1)
        #print(condition2)
        RA = pd.DataFrame.from_dict(meanAbundance[meanAbundance.index.get_level_values('Condition') == condition1].values - meanAbundance[meanAbundance.index.get_level_values('Condition') == condition2].values)
        RA.index = meanAbundance[meanAbundance.index.get_level_values('Condition') == condition1].index.get_level_values('Dataset')
        RA.columns = meanAbundance.columns
        print(RA)
        nodeModulation = RA * CV
        print(nodeModulation)
        nodeModulation = nodeModulation.sum()* score
        #for pathway in pathList: # iterate over pathways
        pathway = pathList[0]
        pathImportances = pd.DataFrame.from_dict(pathway[1], orient = 'index')
        print(pathImportances)
        print(nodeModulation.loc[pathImportances.index])
        nodeModulation_pathway = pathImportances * nodeModulation.loc[pathImportances.index]
        print(nodeModulation_pathway)
        print(nodeModulation_pathway.values.sum())
        # To-do:
        # calculate z-score
        # calculate p-values
        # print out pvalues table
        # print out graphs like Bonita

if __name__ == "__main__":
    
    # load arguments from user
    parser = argparse.ArgumentParser(prog="BONITA")
    parser.set_defaults(sep=",")
    parser.add_argument(
        "-sep",
        "--sep",
        metavar="separator",
        help="How are columns in datafile specified",
    )
    parser.add_argument(
        "-t", action="store_const", const="\t", dest="sep", help="Tab delimited?"
    )
    parser.add_argument("data")

    parser.add_argument("matrix")
    parser.add_argument("diffName")

    results = parser.parse_args()
    
    contrastName = "contrasts.csv"
    conditionName="concatenated_conditions.csv"
    dataName = "concatenated_datasets.csv"
    delimited = ","

    analyze_pathways_mBonita(contrastName = "contrasts.csv", conditionName="concatenated_conditions.csv", dataName = "concatenated_datasets.csv", delimited = ",")

"""
###TEST###
#module load anaconda3/2020.11
#activate BONITA
from pathway_analysis_score_pathways import findPathwayList
import pandas as pd
import numpy as np
contrastName = "contrasts.csv"
conditionName="concatenated_conditions.csv"
dataName = "concatenated_datasets.csv"
delimited = ","
data = pd.read_csv(dataName, sep = delimited, index_col = 0).T
contrasts = pd.read_csv(contrastName, sep = delimited, header = None)
conditionMatrix = pd.read_csv(conditionName, sep = delimited, index_col=0)
groupedData = data.groupby([conditionMatrix['Dataset'], conditionMatrix['Condition']], axis = 0)
meanAbundance = groupedData.mean()
CV = data.groupby(conditionMatrix['Dataset'], axis = 0).std()
CV2 = pd.DataFrame(CV).T
pd.DataFrame(CV2).to_csv("StandardDeviation.csv")

score = data.groupby(conditionMatrix['Dataset']).agg(lambda x: np.count_nonzero(x)).agg(lambda x: np.count_nonzero(x))
score2 = pd.DataFrame(score)
score2.columns = ["EvidenceScore"]
pd.DataFrame(score2).to_csv("EvidenceScore.csv")
pvaluesDF = pd.DataFrame()
pathList = findPathwayList() # identify pathways under consideration
pathDict = retrievePathKey()
nodeMod_complete = pd.DataFrame.from_dict({'nodeModulation':[np.nan], 'Contrast' : np.nan, 'Pathway': np.nan}, orient = 'columns')
RA_complete = pd.DataFrame()

pathImportances_complete = pd.DataFrame()
for i in range(0, len(pathList)): # iterate over pathways
    pathway = pathList[i]
    pathImportances = pd.DataFrame.from_dict(pathway[1], orient = 'index')
    pathImportances['Pathway'] = pathway[0]
    if len(RA_complete) == 0:
        pathImportances_complete = pathImportances
    else:
        pathImportances_complete = pd.concat([pathImportances_complete , pathImportances])
pathImportances_complete.columns = ['ImportanceScore', 'Pathway']
pathImportances_complete = pathImportances_complete.reset_index()
pathImportances_complete.to_csv("ImportanceScores.csv")

for j in range(0, len(contrasts)): # iterate over contrasts
    condition1 = contrasts.iloc[j,0]
    condition2 = contrasts.iloc[j,1]
    #print(condition1)
    #print(condition2)
    RA = pd.DataFrame.from_dict(meanAbundance[meanAbundance.index.get_level_values('Condition') == condition1].values - meanAbundance[meanAbundance.index.get_level_values('Condition') == condition2].values).abs()
    RA.index = meanAbundance[meanAbundance.index.get_level_values('Condition') == condition1].index.get_level_values('Dataset')
    RA.columns = meanAbundance.columns
    RA2 = RA.T
    RA2['Contrast'] = str(condition1)+"_vs_"+str(condition2)
    if len(RA_complete) == 0:
        RA_complete = RA2
    else:
        RA_complete = pd.concat([RA2, RA_complete])
    nodeModulation = RA * CV
    nodeModulation = nodeModulation.sum() * score
    z_scores = {}
    for i in range(0, len(pathList)): # iterate over pathways
        pathway = pathList[i]
        pathImportances = pd.DataFrame.from_dict(pathway[1], orient = 'index')
        nodeModulation_pathway = pathImportances * pd.DataFrame(nodeModulation.loc[set(nodeModulation.index).intersection(set(pathImportances.index))])
        nodeModPathway = nodeModulation_pathway.values.sum()
        nodeModulation_pathway['Contrast'] = str(condition1)+"_vs_"+str(condition2)
        nodeModulation_pathway['Pathway'] = str(pathway[0])
        nodeModulation_pathway.columns = ['nodeModulation', 'Contrast', 'Pathway']
        #pd.DataFrame(nodeModulation, columns = ["nodeModulation"]).to_csv(str(condition1)+"_vs_"+str(condition2)+"_"+str(pathway[0])+"_nodeModulation.csv")
        nodeMod_complete = pd.concat([nodeMod_complete, nodeModulation_pathway])

nodeMod_complete = nodeMod_complete.reset_index().dropna()
#nodeMod_complete.columns = ['Gene','nodeModulation', 'Contrast, Pathway']
nodeMod_complete.to_csv("nodeModulation.csv", index=True)
RA_complete.to_csv("RelativeAbundance.csv")
"""
        #Calculate z-score
        randomScores = []
        for i in range(1000):
            nm_temp = pd.DataFrame(RA.T.sample(n = len(pathImportances.index))).reset_index(drop=True) * pd.DataFrame(CV.T.sample(n = len(pathImportances.index))).reset_index(drop=True)
            nm_temp = nm_temp.sum() * pd.DataFrame(score.T.sample(n = len(pathImportances.index))).reset_index(drop=True)
            nm_temp = (pd.DataFrame(pathImportances.T.reset_index(drop=True)) * nm_temp).values.sum()
            randomScores.append(nm_temp)
        meaner = np.mean(randomScores)
        stdev = np.std(randomScores)
        if stdev == 0:
            zscore = 0
        else:
            zscore = (score - meaner) / stdev
        z_scores[pathway[0]] = zscore
    pvals = norm.sf(list(z_scores.values()))  # calculate p value
    pvalDict = {}
    for i in range(0, len(pathList)):
        pathway = pathList[i]
        pvalDict[pathway[0]] = pvals[i]
    print(pvalDict)
    temp = pd.DataFrame.from_dict(pvalDict, orient = 'index')
    print(temp)
    temp['Contrast'] = ' vs '.join([condition1, condition2])
    temp = temp.reset_index(drop = False)
    print(temp)
    temp['-log10Pvalue'] = np.nan
    print(temp.head())
    temp.columns = ['Pathway', 'P-value', 'Contrast', '-log10Pvalue']
    temp['-log10Pvalue'] = (-1)*np.log10(temp['P-value'])
    if len(pvaluesDF) == 0:
        pvaluesDF = temp
    else:
        pvaluesDF = pd.concat([pvaluesDF, temp])
    #print(temp)
pvaluesDF['Pathway Name'] = [pathDict[i[3:]] for i in pvaluesDF.Pathway]
pvaluesDF.to_csv("pvals_temp.csv", index = False)
"""