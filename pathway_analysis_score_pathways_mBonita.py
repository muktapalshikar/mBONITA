from pathway_analysis_score_pathways import findPathwayList
import pandas as pd
from numpy import count_nonzero
#import argparse


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
        print(condition1)
        print(condition2)
        RA = pd.DataFrame(meanAbundance[meanAbundance.index.get_level_values('Condition') == condition1].values - meanAbundance[meanAbundance.index.get_level_values('Condition') == condition2].values)
        RA.index = meanAbundance[meanAbundance.index.get_level_values('Condition') == condition1].index.get_level_values('Dataset')
        RA.columns = meanAbundance.columns
        print(RA)
        nodeModulation = RA * CV
        nodeModulation = nodeModulation.sum()* score
        #for pathway in pathList: # iterate over pathways
        pathway = pathList[0]
        print(pathway[0])
        pathImportances = pd.DataFrame.from_dict(pathway[1], orient = 'index')
        nodeModulation_pathway = pathImportances * nodeModulation.loc[pathImportances.index]
        print(nodeModulation_pathway)
        print(nodeModulation_pathway.values.sum())

        # To-do:
        # calculate z-score
        # calculate p-values
        # print out pvalues table
        # print out graphs like Bonita

if __name__ == "__main__":
    """
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
    """
    contrastName = "contrasts.csv"
    conditionName="concatenated_conditions.csv"
    dataName = "concatenated_datasets.csv"
    delimited = ","

    analyze_pathways_mBonita(contrastName = "contrasts.csv", conditionName="concatenated_conditions.csv", dataName = "concatenated_datasets.csv", delimited = ",")