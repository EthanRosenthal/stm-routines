import scipy.io as scio
from sklearn import linear_model
import sklearn.cross_validation as cross_val
import numpy as np

# Things to try:
# sklearn.preprocessing.PolynomialFeatures for generating higher order features.
# sklearn.learning_curve
# normalize, scale, etc...

def split_data(x_train, y_train):
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


def run_regression(x, y, lamb):
    regress_model = linear_model.LogisticRegression(penalty='l2',
                                                    C=lamb)
    regress_model.fit(x, y)
    return regress_model

def create_gapmap_class(raw, regress_model):
    px1, px2, nenergy = raw.shape
    raw_test = np.reshape(raw, (px1*px2, nenergy)) # Flatten map
    gapmap_class = regress_model.predict(raw_test)
    gapmap_class = np.reshape(gapmap_class, (px1, px2))
    return gapmap_class


def main(x_train, y_train):
    x_train, x_crossval, x_test, y_train, y_crossval, y_test = \
    split_data(x_train, y_train)

    model = run_regression(x_train, y_train, 2)

    cross_val_pred = model.predict(x_crossval)
    score = model.score(x_crossval, y_crossval)
    return model