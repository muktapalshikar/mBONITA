import pandas as pd
import numpy as np
from random import random

def genInitValueList(newSampleList, model):
    """use dictionaries of values at each node for each sample to construct initial value list for the model"""
    newInitValueList = []
    for j in range(0, len(newSampleList)):
        newInitValueList.append([])
    for j in range(0, len(model.nodeList)):
        for k in range(0, len(newSampleList)):
            ss = newSampleList[k]
            if model.nodeList[j] in newSampleList[0]:
                newInitValueList[k].append(ss[model.nodeList[j]])
            else:
                newInitValueList[k].append(0.5)
    return newInitValueList

def genEBNInitValues(individual, model, sampleProbs):
    """init value generator for EBNs"""
    return [True if (random()<sampleProbs[node]) else False for node in range(0,len(sampleProbs))]
    #initValues = np.zeros(500, dtype=np.intc, order="C")
    #for node in range(0, len(sampleProbs)):
    #    if random() < sampleProbs[node]:
    #        initValues[node] = 1
    #return initValues
#initValues = genEBNInitValues(individual, model, sampleProbs)  # get initial values for all nodes

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
        print(alldata.shape)
        print(col)
        temp = [not i for i in df.index.isin(['Type'])]
        fsc = fakeSingleCells(df.loc[temp, col], numberOfCells=5)
        fsc.loc['Type'] = df.loc['Type']
        if len(alldata) == 0:
            alldata = fsc
        else:
            alldata = pd.concat([alldata, df], axis = 1).fillna(0)
    return alldata

if __name__ == "__main__":
    concatDF = readData(datafiles = {"proteins": "bonita_proteomics.csv", "phosphoproteins": "bonita_phosphoproteomics.csv","mRNA": "bonita_transcriptomics.csv"})
    print("concatdf:", concatDF.shape)
    splitDF = splitData(concatDF, ["proteins", "phosphoproteins"], ["mRNA"])
    print("input", splitDF['input'].shape)
    print("target", splitDF['target'].shape)
    alldata = makeFSCdataset(splitDF['input'], numberOfCells = 10)
    print(alldata.columns)
