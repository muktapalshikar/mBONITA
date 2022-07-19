from math import floor
from re import L
import pandas as pd
from random import random
import networkx as nx
from sklearn.metrics import explained_variance_score
from tensorflow.keras import models, layers, utils, backend as K
import tensorflow.keras
import tensorflow as tf
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import shap
#import random
from random import random, sample, seed
import numpy as np

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
    #targetdata = pd.concat([targetdata, makeFSCdataset(splitDF["target"], numberOfCells=5)])
    #print(targetdata)
    return concatDF, splitDF, inputdata, targetdata

def explainer_shap(model, X_names, X_instance, X_train=None, task="classification", top=10):
    '''
    Use shap to build an explainer.
    :parameter
        :param model: model instance (after fitting)
        :param X_names: list
        :param X_instance: array of size n x 1 (n,)
        :param X_train: array - if None the model is simple machine learning, if not None then it's a deep learning model
        :param task: string - "classification", "regression"
        :param top: num - top features to display
    :return
        dtf with explanations
    '''
    ## create explainer
    ### machine learning
    if X_train is None:
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_instance)
    ### deep learning
    else:
        explainer = shap.DeepExplainer(model, data=X_train[:100])
        shap_values = explainer.shap_values(X_instance.reshape(1,-1))[0].reshape(-1)

    ## plot
    ### classification
    if task == "classification":
        shap.decision_plot(explainer.expected_value, shap_values, link='logit', feature_order='importance',
                           features=X_instance, feature_names=X_names, feature_display_range=slice(-1,-top-1,-1))
    ### regression
    else:
        shap.waterfall_plot(explainer.expected_value[0], shap_values, 
                            features=X_instance, feature_names=X_names, max_display=top)

def R2(y, y_hat):
    ss_res = K.sum(K.square(y - y_hat))
    ss_tot = K.sum(K.square(y - K.mean(y)))
    return 1 - ss_res / (ss_tot + K.epsilon())

def experimentPartTwoWrapper():
    """Simple perceptron"""
    # prepare data
    concatDF, splitDF, inputdata, targetdata = experimentPartOneWrapper()
    # prepare network
    testNet = nx.read_graphml("large_hif1Agraph.graphml")
    # get node
    testNode = dict(
        filter(
            #lambda elem: elem[1] == max(dict(testNet.in_degree).values()),
            lambda elem: elem[1] >= 3,
            dict(testNet.in_degree).items(),
        )
    )
    import random
    answers = {}
    testNodes = list(testNode.keys())
    #testNode = testNode[random.sample(range(0, len(testNode)), 1)[0]]
    for testNode in testNodes:
        print(testNode)
        # get upstream nodes
        upstream = testNet.in_edges(testNode, data = True)
        print(upstream)
        signal = [i[2]['signal'] for i in upstream]
        signal = [-1 if x == "i" else 1 for x in signal]
        signal.extend(signal)
        upstream = [i[0] for i in upstream]
        # subset input data
        testInput = inputdata.iloc[
            inputdata.index.get_level_values("Entity").isin(upstream)
        ].T.astype('int')#.mul(signal)
        print(testInput.shape)
        # subset target data
        testTarget = targetdata.iloc[
            targetdata.index.get_level_values("Entity").isin([testNode])
        ].T.astype('int')
        print(testTarget.shape)
        print(testTarget)
        print(testInput)
        """
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
        """
        model = models.Sequential(name="DeepNN", layers=[
        ### hidden layer 1
        layers.Dense(name="h1", input_dim=len(
                        testInput.columns
                    ),
                    units=int(round((len(
                        testInput.columns
                    )+1)/2)), 
                    activation='linear'),
        layers.Dropout(name="drop1", rate=0.2),
        
        ### hidden layer 2
        layers.Dense(name="h2", units=int(round((len(
                        testInput.columns
                    )+1)/4)), 
                    activation='relu'),
        layers.Dropout(name="drop2", rate=0.2),
        
        ### layer output
        layers.Dense(name="output", units=1, activation='sigmoid')
        ])
        model.summary()
        


        # compile the neural network
        model.compile(
            optimizer="adam",
            loss="mean_absolute_error",
            metrics=[
                #tensorflow.keras.metrics.AUC(),
                #tensorflow.keras.metrics.FalsePositives(),
                #tensorflow.keras.metrics.FalseNegatives(),
                tensorflow.keras.metrics.BinaryAccuracy(),
                R2,
                #tensorflow.keras.metrics.Recall(),
            ],
        )

        n_samples = len(testInput.index)
        trainingSamples = random.sample(range(0, n_samples), floor(3*n_samples/4))
        testSamples = list(set(range(0,n_samples)).difference(set(trainingSamples)))
        X = testInput.iloc[trainingSamples]
        y = testTarget.iloc[trainingSamples]
        training = model.fit(
            x=X, y=y, batch_size=2, epochs=16, shuffle=True, verbose=0, validation_split=0.3
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
        colors = ["blue", "red", "orange"]
        i = 0
        for metric in metrics:
            color = colors[i]
            i = i + 1
            ax11.plot(training.history[metric], label=metric, color = color)
        ax11.set_ylabel("Score", color="red")
        ax11.legend()

        ## validation
        i = 0
        ax[1].set(title="Validation")
        ax22 = ax[1].twinx()
        ax[1].plot(training.history["val_loss"], color="black")
        ax[1].set_xlabel("Epochs")
        ax[1].set_ylabel("Loss", color="black")
        for metric in metrics:
            ax22.plot(training.history["val_" + metric], label=metric, color=colors[i])
            i = i + 1
        ax22.set_ylabel("Score", color="red")
        ax11.legend()
        plt.savefig("temp.png")
        plt.close()

        ## explainer_shap(model, upstream, X, X_train=X, task="regression", top=10)
        from numpy import argmax
        test_predictions = argmax(model.predict(testInput.iloc[testSamples]),axis=-1).flatten()

        a = plt.axes(aspect='equal')
        plt.scatter(testTarget.iloc[testSamples], test_predictions)
        plt.xlabel('True Values')
        plt.ylabel('Predictions')
        lims = [0, 1]
        plt.xlim(lims)
        plt.ylim(lims)
        plt.savefig('temp2.png')

        print(test_predictions)
        #print(testTarget.iloc[testSamples])
        #print(len(testTarget.iloc[testSamples]))
        answer = testTarget.iloc[testSamples]
        answer = answer.iloc[:,0]
        answer = answer.to_list()
        print(answer)
        print(len(answer))
        print(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)]))

        answers[testNode] = [sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer), test_predictions, answer]

    #print(answers)

    for tn in answers.keys():
        if sum(answers[tn][2]) > 0:
            print(answers[tn])

def xorTutorial():
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
    # input X vector
    X = [[0, 0], [0, 1], [1, 0], [1, 1]]
    # output Y vector
    Y = [[0], [1], [1], [0]]
    
    # Placeholders for input and output
    x = tf.placeholder(tf.float32, shape=[4,2])
    y = tf.placeholder(tf.float32, shape=[4,1])
    
    # W matrix
    W1 = tf.Variable([[1.0, 0.0], [1.0, 0.0]], shape=[2,2])
    W2 = tf.Variable([[0.0], [1.0]], shape=[2,1])
    
    # Biases
    B1 = tf.Variable([0.0, 0.0], shape=[2])
    B2 = tf.Variable([0.0], shape=1)
    
    # Hidden layer and outout layer
    output =tf.sigmoid(tf.matmul(tf.sigmoid(tf.matmul(x, W1) + B1), W2) + B2)
    
    # error estimation
    e = tf.reduce_mean(tf.squared_difference(y, output))
    train = tf.train.GradientDescentOptimizer(0.1).minimize(e)
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    
    for i in range (100001):
        error = sess.run(train, feed_dict={x: X, y: Y})
        if i % 10000 == 0:
            print('\nEpoch: ' + str(i))
            print('\nError: ' + str(sess.run(e, feed_dict={x: X, y: Y})))
            for el in sess.run(output, feed_dict={x: X, y: Y}):
                print('    ',el)
    sess.close()
    
    print ("Complete")

def three_vars_gateDeepNN(i):
    seed(i)
    import numpy as np
    # input X vector
    #X1 = pd.DataFrame(np.repeat(np.array([[0, 0], [0, 1], [1, 0], [1, 1]]), 100, axis = 0))
    a = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(a)
    b = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(b)
    d = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(d)
    b = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(b)
    c = [ int(int(j) & int(i) | int(k)) for i,j,k in zip(a, b, d)] # B & A | D
    X1 = pd.DataFrame({'A': a, 'B': b, 'D': d})
    print(X1)
    # output Y vector
    #Y1 = pd.DataFrame(np.repeat(np.array([[0], [1], [1], [0]]), 100, axis =  0))
    Y1 = pd.DataFrame({'C': c})
    print(Y1)
    model = models.Sequential(name="DeepNN", layers=[
    ### hidden layer 1
    layers.Dense(name="h1", input_dim=3,
                units=int(2), 
                activation='tanh'),
    layers.Dropout(name="drop1", rate=0.2),
    ### hidden layer 2
    layers.Dense(name="h2", units=int(2), 
                activation='sigmoid'),
    layers.Dropout(name="drop2", rate=0.2),
    ### layer output
    layers.Dense(name="output", units=1, activation='sigmoid')
    ])
    model.summary()

    # compile the neural network
    model.compile(
        optimizer="adam",
        loss="mean_absolute_error",
        metrics=[
            #tensorflow.keras.metrics.AUC(),
            #tensorflow.keras.metrics.FalsePositives(),
            #tensorflow.keras.metrics.FalseNegatives(),
            tensorflow.keras.metrics.BinaryAccuracy(),
            #R2,
            #tensorflow.keras.metrics.Recall(),
        ],
    )
    n_samples = len(X1)
    trainingSamples = sample(range(0, n_samples), floor(3*n_samples/4))
    testSamples = list(set(range(0,n_samples)).difference(set(trainingSamples)))
    X = X1#.iloc[trainingSamples]
    y = Y1#.iloc[trainingSamples]
    training = model.fit(
        x=X, y=y, batch_size=10, epochs=1000, shuffle=False, verbose=0, validation_split=0.3
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
    colors = ["blue", "red", "orange"]
    i = 0
    for metric in metrics:
        color = colors[i]
        i = i + 1
        ax11.plot(training.history[metric], label=metric, color = color)
    ax11.set_ylabel("Score", color="red")
    ax11.legend()

    ## validation
    i = 0
    ax[1].set(title="Validation")
    ax22 = ax[1].twinx()
    ax[1].plot(training.history["val_loss"], color="black")
    ax[1].set_xlabel("Epochs")
    ax[1].set_ylabel("Loss", color="black")
    for metric in metrics:
        ax22.plot(training.history["val_" + metric], label=metric, color=colors[i])
        i = i + 1
    ax22.set_ylabel("Score", color="red")
    ax11.legend()
    plt.savefig("3_variable_gate_"+str(i)+".png")
    plt.close()

    ## explainer_shap(model, upstream, X, X_train=X, task="regression", top=10)
    answer = Y1.iloc[testSamples]
    answer = answer.iloc[:,0]
    answer = answer.to_list()
    from numpy import argmax
    test_predictions = model.predict(X1.iloc[testSamples]) #argmax(model.predict(X1.iloc[testSamples]),axis=-1).flatten()
    test_predictions = np.round(test_predictions)
    plt.axes(aspect='equal')
    plt.title("Proportion true predictions: " + str(round(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer), 2)))
    #plt.scatter(Y1.iloc[testSamples], test_predictions)
    plt.hist([int(int(i) == int(j)) for i,j in zip(answer, test_predictions)])
    plt.xlabel('True Values')
    plt.ylabel('Predictions')
    lims = [0, 1]
    plt.xlim(lims)
    #plt.ylim(lims)
    plt.savefig("3_variable_gate_predictions_"+str(i)+".png")

    print(test_predictions)
    #print(testTarget.iloc[testSamples])
    #print(len(testTarget.iloc[testSamples]))
    print(answer)
    print(len(answer))
    print(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer))

def five_vars_gateDeepNN(i):
    seed(i)
    import numpy as np
    # input X vector
    #X1 = pd.DataFrame(np.repeat(np.array([[0, 0], [0, 1], [1, 0], [1, 1]]), 100, axis = 0))
    a = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(a)
    b = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(b)
    d = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(d)
    e = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(e)
    f = np.hstack((np.ones(150), np.zeros(150)))
    np.random.shuffle(f)
    c = [ int(int(j) & int(i) | int(k) | int(l) & int(m)) for i,j,k,l,m in zip(a, b, d, e, f)] # B & A | D
    X1 = pd.DataFrame({'A': a, 'B': b, 'D': d, 'E' : e, 'F': f})
    print(X1)
    # output Y vector
    #Y1 = pd.DataFrame(np.repeat(np.array([[0], [1], [1], [0]]), 100, axis =  0))
    Y1 = pd.DataFrame({'C': c})
    print(Y1)
    model = models.Sequential(name="DeepNN", layers=[
    ### hidden layer 1
    layers.Dense(name="h1", input_dim=5,
                units=int(3), 
                activation='tanh'),
    layers.Dropout(name="drop1", rate=0.2),
    ### hidden layer 2
    layers.Dense(name="h2", units=int(2), 
                activation='sigmoid'),
    layers.Dropout(name="drop2", rate=0.2),
    ### layer output
    layers.Dense(name="output", units=1, activation='sigmoid')
    ])
    model.summary()

    # compile the neural network
    model.compile(
        optimizer="adam",
        loss="mean_absolute_error",
        metrics=[
            #tensorflow.keras.metrics.AUC(),
            #tensorflow.keras.metrics.FalsePositives(),
            #tensorflow.keras.metrics.FalseNegatives(),
            tensorflow.keras.metrics.BinaryAccuracy(),
            #R2,
            #tensorflow.keras.metrics.Recall(),
        ],
    )
    n_samples = len(X1)
    trainingSamples = sample(range(0, n_samples), floor(3*n_samples/4))
    testSamples = list(set(range(0,n_samples)).difference(set(trainingSamples)))
    X = X1#.iloc[trainingSamples]
    y = Y1#.iloc[trainingSamples]
    training = model.fit(
        x=X, y=y, batch_size=10, epochs=1000, shuffle=False, verbose=0, validation_split=0.3
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
    colors = ["blue", "red", "orange"]
    i = 0
    for metric in metrics:
        color = colors[i]
        i = i + 1
        ax11.plot(training.history[metric], label=metric, color = color)
    ax11.set_ylabel("Score", color="red")
    ax11.legend()

    ## validation
    i = 0
    ax[1].set(title="Validation")
    ax22 = ax[1].twinx()
    ax[1].plot(training.history["val_loss"], color="black")
    ax[1].set_xlabel("Epochs")
    ax[1].set_ylabel("Loss", color="black")
    for metric in metrics:
        ax22.plot(training.history["val_" + metric], label=metric, color=colors[i])
        i = i + 1
    ax22.set_ylabel("Score", color="red")
    ax11.legend()
    plt.savefig("5_variable_gate_"+str(i)+".png")

    plt.close()

    ## explainer_shap(model, upstream, X, X_train=X, task="regression", top=10)
    answer = Y1.iloc[testSamples]
    answer = answer.iloc[:,0]
    answer = answer.to_list()
    test_predictions = model.predict(X1.iloc[testSamples]) #argmax(model.predict(X1.iloc[testSamples]),axis=-1).flatten()
    test_predictions = np.round(test_predictions)
    #plt.axes(aspect='equal')
    plt.title("Proportion true predictions: " + str(round(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer),2)))
    #plt.scatter(Y1.iloc[testSamples], test_predictions)
    plt.hist([int(int(i) == int(j)) for i,j in zip(answer, test_predictions)])
    plt.xlabel('True Values')
    plt.ylabel('Predictions')
    lims = [0, 1]
    plt.xlim(lims)
    #plt.ylim(lims)
    plt.savefig("5_variable_gate_predictions_"+str(i)+".png")

    print(test_predictions)
    #print(testTarget.iloc[testSamples])
    #print(len(testTarget.iloc[testSamples]))
    print(answer)
    print(len(answer))
    print(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer))

def ten_vars_gateDeepNN(i):
    seed(i)
    import numpy as np
    # input X vector
    cols = [chr(x) for x in range(97, 123)]
    cols = cols[0:10]
    data = {}
    for letter in cols:
        print(letter)
        randomData = np.hstack((np.ones(150), np.zeros(150)))
        np.random.shuffle(randomData)
        data[letter] = randomData
        X1 = pd.DataFrame(data)
    print(X1.columns)
    z = [ int(int(a1) & int(b1) | int(c1) | int(d1) & int(e1) | (not int(f1)) | int (g1) | int(h1) | int(i1) & int(j1) ) for a1, b1, c1, d1, e1, f1, g1, h1, i1, j1 in zip(X1['a'], X1['b'], X1['c'], X1['d'], X1['e'], X1['f'], X1['g'], X1['h'], X1['i'], X1['j'])] # B & A | D
            
    print(X1)
    # output Y vector
    #Y1 = pd.DataFrame(np.repeat(np.array([[0], [1], [1], [0]]), 100, axis =  0))
    Y1 = pd.DataFrame({'Z': z})
    print(Y1)
    model = models.Sequential(name="DeepNN", layers=[
    ### hidden layer 1
    layers.Dense(name="h1", input_dim=10,
                units=int(3), 
                activation='tanh'),
    layers.Dropout(name="drop1", rate=0.2),
    ### hidden layer 2
    layers.Dense(name="h2", units=int(2), 
                activation='sigmoid'),
    layers.Dropout(name="drop2", rate=0.2),
    ### layer output
    layers.Dense(name="output", units=1, activation='sigmoid')
    ])
    model.summary()

    # compile the neural network
    model.compile(
        optimizer="adam",
        loss="mean_absolute_error",
        metrics=[
            #tensorflow.keras.metrics.AUC(),
            #tensorflow.keras.metrics.FalsePositives(),
            #tensorflow.keras.metrics.FalseNegatives(),
            tensorflow.keras.metrics.BinaryAccuracy(),
            #R2,
            #tensorflow.keras.metrics.Recall(),
        ],
    )
    n_samples = len(X1)
    trainingSamples = sample(range(0, n_samples), floor(3*n_samples/4))
    testSamples = list(set(range(0,n_samples)).difference(set(trainingSamples)))
    X = X1#.iloc[trainingSamples]
    y = Y1#.iloc[trainingSamples]
    training = model.fit(
        x=X, y=y, batch_size=10, epochs=1000, shuffle=False, verbose=0, validation_split=0.3
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
    colors = ["blue", "red", "orange"]
    i = 0
    for metric in metrics:
        color = colors[i]
        i = i + 1
        ax11.plot(training.history[metric], label=metric, color = color)
    ax11.set_ylabel("Score", color="red")
    ax11.legend()

    ## validation
    i = 0
    ax[1].set(title="Validation")
    ax22 = ax[1].twinx()
    ax[1].plot(training.history["val_loss"], color="black")
    ax[1].set_xlabel("Epochs")
    ax[1].set_ylabel("Loss", color="black")
    for metric in metrics:
        ax22.plot(training.history["val_" + metric], label=metric, color=colors[i])
        i = i + 1
    ax22.set_ylabel("Score", color="red")
    ax11.legend()
    plt.savefig("10_variable_gate_"+str(i)+".png")

    plt.close()

    ## explainer_shap(model, upstream, X, X_train=X, task="regression", top=10)
    answer = Y1.iloc[testSamples]
    answer = answer.iloc[:,0]
    answer = answer.to_list()
    test_predictions = model.predict(X1.iloc[testSamples]) #argmax(model.predict(X1.iloc[testSamples]),axis=-1).flatten()
    test_predictions = np.round(test_predictions)
    #plt.axes(aspect='equal')
    plt.title("Proportion true predictions: " + str(round(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer),2)))
    #plt.scatter(Y1.iloc[testSamples], test_predictions)
    plt.hist([int(int(i) == int(j)) for i,j in zip(answer, test_predictions)])
    plt.xlabel('True Values')
    plt.ylabel('Predictions')
    lims = [0, 1]
    plt.xlim(lims)
    #plt.ylim(lims)
    plt.savefig("10_variable_gate_predictions_"+str(i)+".png")

    print(test_predictions)
    #print(testTarget.iloc[testSamples])
    #print(len(testTarget.iloc[testSamples]))
    print(answer)
    print(len(answer))
    print(sum([1 if i == j else 0 for i, j in zip(test_predictions,answer)])/len(answer))

if __name__ == "__main__":
    # prepare data
    #concatDF, splitDF, inputdata, targetdata = experimentPartOneWrapper()
    # prepare network
    #testNet = nx.read_graphml("large_hif1Agraph.graphml")
    #experimentPartTwoWrapper()
    for i in range(0,3):
        print(i)
        ten_vars_gateDeepNN(i)
        five_vars_gateDeepNN(i)
        three_vars_gateDeepNN(i)