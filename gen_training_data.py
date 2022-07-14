import pandas as pd
import numpy as np
from random import random

def readData(datafiles = {}):
    """
    Read csv training data files.
    datafiles: dictionary where keys = labels (protein, proteomics, mRNA, etc), values = paths to training data files
    returns:
        a concatenated dataset with labels
    """
    finaldf = pd.DataFrame()
    for key, value in zip(datafiles.keys(), datafiles.values()):
        df = pd.read_csv(value, index_col="Gene")
        df.loc['Type'] = key
        if len(finaldf) == 0:
            finaldf = df
        else:
            finaldf = pd.concat([finaldf, df], axis = 1).fillna(0)
    return finaldf

def splitData(concatDF, trainingLabel, targetLabel):
    """
        splits the output of readData into target and input datasets
        Input:
            concatDF: output of readData
            trainingLabel: values of concatDF['Type'] representing the input data
            targetLabel: values of concatDF['Type'] representing the target data
        Output:
            a dictionary of two dataframes - target and input
    """
    splitDF = {'input': pd.DataFrame(), 'target': pd.DataFrame()}
    splitDF['input'] = concatDF.loc[:, concatDF.isin(trainingLabel).any()]
    splitDF['target'] = concatDF.loc[:, concatDF.isin(targetLabel).any()]
    return splitDF

def fakeSingleCells(column, numberOfCells = 10):
    """
    Takes a bulk sample, rescales, and generates a set of single cells
    """
    column = pd.Series(column)
    column_scaled = column/max(column)
    fsc = {}
    cell = 0
    while cell < numberOfCells:
        fsc[str(cell)] = [random() < node for node in column_scaled]
        cell = cell + 1
    fsc = pd.DataFrame.from_dict(fsc)
    fsc.index = column.index
    return fsc

def makeFSCdataset(df, numberOfCells = 10):
    """
    Apply fakeSingleCells over all samples in the given dataset
    """
    alldata = pd.DataFrame()
    for col in df.columns:
        temp = [not i for i in df.index.isin(['Type'])]
        fsc = fakeSingleCells(df.loc[temp, col], numberOfCells=5)
        type = df.loc['Type', col]
        fsc.index = pd.MultiIndex.from_tuples([(type, str(i), col) for i in fsc.index], names=["Type", "Entity", "Sample"])
        if len(alldata) == 0:
            alldata = fsc
        else:
            alldata = pd.concat([alldata, fsc], axis = 0).fillna(0)
        alldata.columns = [str(i) for i in range(0, len(alldata.columns))]
    return alldata

def experimentPartOneWrapper():
    concatDF = readData(datafiles = {"proteins": "bonita_proteomics.csv", "phosphoproteins": "bonita_phosphoproteomics.csv","mRNA": "bonita_transcriptomics.csv"})
    print("concatdf:", concatDF.shape)
    splitDF = splitData(concatDF, ["mRNA", "phosphoproteins"], ["proteins"])
    print("input", splitDF['input'].shape)
    print("target", splitDF['target'].shape)
    # generate input data
    inputdata = makeFSCdataset(splitDF['input'], numberOfCells = 10)
    print(inputdata)
    # generate target data
    targetdata = makeFSCdataset(splitDF['target'], numberOfCells = 10)
    print(targetdata)
    return concatDF, splitDF, inputdata, targetdata

if __name__ == "__main__":
    concatDF, splitDF, inputdata, targetdata = experimentPartOneWrapper()
