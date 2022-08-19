from gen_training_data import *
model = experimentPartFiveWrapper()

finalInput = pd.read_csv("finalInput.csv", index_col = 0, header = [0,1])
finalTarget = pd.read_csv("finalTarget.csv", index_col = 0)

n_samples = len(finalInput.index)
trainingSamples = sample(range(0, n_samples), floor(3 * n_samples / 4))
testSamples = list(set(range(0, n_samples)).difference(set(trainingSamples)))

answer = finalTarget.iloc[[0]]
#answer = answer.iloc[:, 0]
answer = answer.to_list()
    
test_predictions = np.round(model.predict(np.asarray(finalInput.iloc[0]).astype('float32')))
print(
    "Proportion true predictions: "
    + str(
        round(
            sum([1 if i == j else 0 for i, j in zip(test_predictions, answer)])
            / len(answer),
            2,
        )
    )
)