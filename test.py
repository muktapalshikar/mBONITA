# test drive some tools

from tensorflow.keras import models, layers, utils, backend as K
import matplotlib.pyplot as plt

# import shap
import numpy as np

model = models.Sequential(
    name="Perceptron",
    layers=[
        layers.Dense(  # a fully connected layer
            name="dense",
            input_dim=3,  # with 3 features as the input
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
model.compile(optimizer="adam", loss="mean_absolute_error", metrics=[R2])


X = np.random.rand(100, 3)
y = np.random.choice([1, 0], size=100)  # train/validation
training = model.fit(
    x=X, y=y, batch_size=2, epochs=10, shuffle=True, verbose=0, validation_split=0.3
)

# plot
metrics = [k for k in training.history.keys() if ("loss" not in k) and ("val" not in k)]
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

"""
Use shap to build an a explainer.
:parameter
    :param model: model instance (after fitting)
    :param X_names: list
    :param X_instance: array of size n x 1 (n,)
    :param X_train: array - if None the model is simple machine learning, if not None then it's a deep learning model
    :param task: string - "classification", "regression"
    :param top: num - top features to display
:return
    dtf with explanations

def explainer_shap(model, X_names, X_instance, X_train=None, task="classification", top=10):
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

i = 1
explainer_shap(model, 
               X_names=range(0,3), 
               X_instance=X[i], 
               X_train=X, 
               task="regression",
               top=10)
"""
