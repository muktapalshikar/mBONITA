import pandas as pd
from random import random
import networkx as nx
from tensorflow.keras import models, layers, utils, backend as K
import tensorflow.keras
import matplotlib.pyplot as plt


def readData(datafiles={}):
    """
    Read csv training data files.
    datafiles: dictionary where keys = labels (protein, proteomics, mRNA, etc), values = paths to training data files
    returns:
        a concatenated dataset with labels
    """
    finaldf = pd.DataFrame()
    for key, value in zip(datafiles.keys(), datafiles.values()):
        df = pd.read_csv(value, index_col="Gene")
        df.loc["Type"] = key
        if len(finaldf) == 0:
            finaldf = df
        else:
            finaldf = pd.concat([finaldf, df], axis=1).fillna(0)
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
    splitDF = {"input": pd.DataFrame(), "target": pd.DataFrame()}
    splitDF["input"] = concatDF.loc[:, concatDF.isin(trainingLabel).any()]
    splitDF["target"] = concatDF.loc[:, concatDF.isin(targetLabel).any()]
    return splitDF


def fakeSingleCells(column, numberOfCells=10):
    """
    Takes a bulk sample, rescales, and generates a set of single cells
    """
    column = pd.Series(column)
    column_scaled = column / max(column)
    fsc = {}
    cell = 0
    while cell < numberOfCells:
        fsc[str(cell)] = [random() < node for node in column_scaled]
        cell = cell + 1
    fsc = pd.DataFrame.from_dict(fsc)
    fsc.index = column.index
    return fsc


def makeFSCdataset(df, numberOfCells):
    """
    Apply fakeSingleCells over all samples in the given dataset
    """
    alldata = pd.DataFrame()
    for col in df.columns:
        temp = [not i for i in df.index.isin(["Type"])]
        fsc = fakeSingleCells(df.loc[temp, col], numberOfCells=numberOfCells)
        type = df.loc["Type", col]
        #fsc.index = pd.MultiIndex.from_tuples(
        #    [(type, str(i), col) for i in fsc.index], names=["Type", "Entity", "Sample"]
        #)
        fsc.index = pd.MultiIndex.from_tuples(
            [(type, str(i)) for i in fsc.index], names=["Type", "Entity"]
        )
        if len(alldata) == 0:
            alldata = fsc
        else:
            alldata = pd.concat([alldata, fsc], axis=1).fillna(0)
        alldata.columns = [str(i) for i in range(0, len(alldata.columns))]
    return alldata


def experimentPartOneWrapper():
    concatDF = readData(
        datafiles={
            "proteins": "bonita_proteomics.csv",
            "phosphoproteins": "bonita_phosphoproteomics.csv",
            "mRNA": "bonita_transcriptomics.csv",
        }
    )
    prot_common_conditions = [
        "norm_19-A",
        "norm_19-B",
        "norm_19-C",
        "norm_1-A",
        "norm_1-B",
        "norm_1-C",
        "norm_1-1ug-A",
        "norm_1-1ug-B",
        "norm_1-1ug-C",
    ]
    mrna_common_conditions = [
        "Ramos_19O2_NoCyclo_1",
        "Ramos_19O2_NoCyclo_2",
        "Ramos_19O2_NoCyclo_3",
        "Ramos_1O2_NoCyclo_1",
        "Ramos_1O2_NoCyclo_2",
        "Ramos_1O2_NoCyclo_3",
        "Ramos_1O2_PlusCyclo_1",
        "Ramos_1O2_PlusCyclo_2",
        "Ramos_1O2_PlusCyclo_3",
    ]
    phosph_common_conditions = [
        "Gr1_F1: 126, Sample, 1",
        "Gr1_F2: 126, Sample, 1",
        "Gr1_F3: 126, Sample, 1",
        "Gr1_F1: 128C, Sample, 4",
        "Gr1_F2: 128C, Sample, 4",
        "Gr1_F3: 128C, Sample, 4",
        "Gr2_F1: 128N, Control, 781",
        "Gr2_F2: 128N, Control, 781",
        "Gr2_F3: 128N, Control, 781",
    ]
    concatDF = concatDF.loc[:,  prot_common_conditions+mrna_common_conditions+phosph_common_conditions]
    #print("concatdf:", concatDF.shape)
    splitDF = splitData(concatDF, ["mRNA", "phosphoproteins"], ["proteins"])
    #print("input", splitDF["input"].shape)
    #print("target", splitDF["target"].shape)
    # generate input data
    inputdata = makeFSCdataset(splitDF["input"], numberOfCells=5)
    #print(inputdata)
    # generate target data
    targetdata = makeFSCdataset(splitDF["target"], numberOfCells=10)
    #print(targetdata)
    return concatDF, splitDF, inputdata, targetdata


if __name__ == "__main__":
    # prepare data
    concatDF, splitDF, inputdata, targetdata = experimentPartOneWrapper()
    # prepare network
    testNet = nx.read_graphml("consensus_net.graphml")
    # get node
    testNode = dict(
        filter(
            lambda elem: elem[1] == max(dict(testNet.in_degree).values()),
            dict(testNet.in_degree).items(),
        )
    )
    testNode = list(testNode.keys())[0]
    # get upstream nodes
    upstream = testNet.in_edges(testNode)
    upstream = [i[0] for i in upstream]
    # subset input data
    testInput = inputdata.iloc[
        inputdata.index.get_level_values("Entity").isin(upstream)
    ].T.astype('int')
    print(testInput.shape)
    # subset target data
    testTarget = targetdata.iloc[
        targetdata.index.get_level_values("Entity").isin([testNode])
    ].T.astype('int')
    print(testTarget.shape)
    print(testTarget)
    print(testInput)
    
    model = models.Sequential(
        name="Perceptron",
        layers=[
            layers.Dense(  # a fully connected layer
                name="dense",
                input_dim=len(
                    testInput.columns
                ),  # number of features = number of upstream nodes from the PKN
                units=1,  # and 1 node because we want 1 output
                activation="linear",  # f(x)=x
            )
        ],
    )
    model.summary()

    # define metrics
    def R2(y, y_hat):
        ss_res = K.sum(K.square(y - y_hat))
        ss_tot = K.sum(K.square(y - K.mean(y)))
        return 1 - ss_res / (ss_tot + K.epsilon())

    # compile the neural network
    model.compile(
        optimizer="adam",
        loss="mean_absolute_error",
        metrics=[
            tensorflow.keras.metrics.AUC(),
            tensorflow.keras.metrics.FalsePositives(),
            tensorflow.keras.metrics.FalseNegatives(),
            tensorflow.keras.metrics.Accuracy(),
        ],
    )

    X = testInput
    y = testTarget
    training = model.fit(
        x=X, y=y, batch_size=2, epochs=10, shuffle=True, verbose=0, validation_split=0.3
    )

    # plot
    metrics = [
        k for k in training.history.keys() if ("loss" not in k) and ("val" not in k)
    ]
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(15, 3))

    ## training
    ax[0].set(title="Training")
    ax11 = ax[0].twinx()
    ax[0].plot(training.history["loss"], color="black")
    ax[0].set_xlabel("Epochs")
    ax[0].set_ylabel("Loss", color="black")
    for metric in metrics:
        ax11.plot(training.history[metric], label=metric)
    ax11.set_ylabel("Score", color="steelblue")
    ax11.legend()

    ## validation
    ax[1].set(title="Validation")
    ax22 = ax[1].twinx()
    ax[1].plot(training.history["val_loss"], color="black")
    ax[1].set_xlabel("Epochs")
    ax[1].set_ylabel("Loss", color="black")
    for metric in metrics:
        ax22.plot(training.history["val_" + metric], label=metric)
    ax22.set_ylabel("Score", color="steelblue")
    plt.savefig("temp.png")
    