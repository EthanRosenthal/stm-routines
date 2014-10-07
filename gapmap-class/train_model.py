import scipy.io as scio
from sklearn import linear_model
import sklearn.cross_validation as cross_val
import numpy as np

# Things to try:
# sklearn.preprocessing.PolynomialFeatures for generating higher order features.
# sklearn.learning_curve
# normalize, scale, etc...

def split_data(x_train, y_train):
    """
    Given training data cropped from the original dataset by create_training_set.py, split this data up into training, cross-validation, and test data.

    INPUTS:
    x_train = Features cropped from original dataset
    y_train = Labels manually inputed from x_train

    OUTPUTS:
    new_x_train = New training data randomly selected from x_train
    new_x_crossval = Cross-validation samples from x_train
    new_x_test = Test samples from x_train
    new_y_train = Training labels
    new_y_crossval = Cross-validation labels
    new_y_test = Testing labels
    """
    new_x_train, new_x_test, new_y_train, new_y_test \
     = cross_val.train_test_split(x_train,
                                  y_train,
                                  test_size=0.3,
                                  random_state=53)
    new_x_crossval, new_x_test, new_y_crossval, new_y_test \
     = cross_val.train_test_split(new_x_test,
                                  new_y_test,
                                  test_size=0.5,
                                  random_state=41)
    return new_x_train, new_x_crossval, new_x_test, new_y_train, \
            new_y_crossval, new_y_test

def sweep_lambda(x_train, x_crossval, y_train, y_crossval, lamb_array):
    """
    Train model for an array of regularization values (lamb_array) and return prediction score from cross-validation data.

    INPUTS:
    x_train, x_crossval = Training and cross-validation features
    y_train, y_crossval = Training and cross-validation labels
    lamb_array = Array of lambda regularization values

    OUTPUTS:
    scores = array of prediction probabilities corresponding to each regularization value in lamb_array
    """
    scores = np.zeros(len(lamb_array))
    i = 0
    for lamb in lamb_array:
        model = train_regression(x_train, y_train, lamb)
        scores[i] = model.score(x_crossval, y_crossval)
        i += 1

    return scores

def sweep_sample_size(x_train,
                      x_crossval,
                      y_train,
                      y_crossval,
                      lamb,
                      sample_sizes):
    """
    Train model for an array of regularization values (lamb_array) and return prediction score from cross-validation data.

    INPUTS:
    x_train, x_crossval = Training and cross-validation features
    y_train, y_crossval = Training and cross-validation labels
    lamb = Lambda regularization value
    sample_sizes = Integer array with different values for number of samples to use in training the model.

    OUTPUTS:
    scores = array of prediction probabilities corresponding to training the model on each sample size in sample_sizes
    """
    scores = np.zeros(len(sample_sizes))
    i = 0
    for sz in sample_sizes:
        sz = int(sz)
        model = train_regression(x_train[:sz], y_train[:sz], lamb)
        scores[i] = model.score(x_crossval, y_crossval)
        i+=1
    return scores

def train_regression(x, y, lamb):
    """
    Training logistic regression model on samples, x, and labels, y, using regularization parameter "lamb".
    """
    regress_model = linear_model.LogisticRegression(penalty='l2',
                                                    C=1./lamb)
    regress_model.fit(x, y)
    return regress_model

def create_gapmap_class(raw, regress_model):
    """
    Classify regions of the original dataset as SDW or SC based upon trained regression model.
    INPUTS:
    raw = original dI/dV dataset
    regress_model = Model trained by train_regression

    OUTPUTS:
    gapmap_class = Map of gaps classified by gap-type (0 for SDW, 1 for SC)
    """
    px1, px2, nenergy = raw.shape
    raw_test = np.reshape(raw, (px1*px2, nenergy)) # Flatten map
    gapmap_class = regress_model.predict(raw_test)
    gapmap_class = np.reshape(gapmap_class, (px1, px2)) # Unflatten map
    return gapmap_class


def main(x_train, y_train):
    x_train, x_crossval, x_test, y_train, y_crossval, y_test = \
    split_data(x_train, y_train)

    model = train_regression(x_train, y_train, 0.3)

    cross_val_pred = model.predict(x_crossval)
    score = model.score(x_crossval, y_crossval)
    return model