import scipy.io as scio
import sklearn.cross_validation as cross_val
import numpy as np
import matplotlib.pyplot as plt

def shape_input(raw):
    """
    Takes an input dI/dV map "raw" and splits off a small portion to be used to build a training set. For a 256 x 256 pixel map, a test_size of 0.997 leaves 196 pixels to be used in the training set.
    Contrary to conventional training and test data, this training data will be split into training, cross validation, and test data. The alorithm trained with this data will then be used, unsupervised, on the rest of the map.

    INPUTS:
    raw = 3-D dI/dV map

    OUTPUTS:
    x_train = The portion of the map that has been randomly chosen as training data. This variable is two dimensional. It can be thought of as an array of row vectors where each row corresponds to a single dI/dV curve.
    x_test = The portion of the map that has been randomly left to be used as test data.
    """
    px1, px2, nenergy = raw.shape
    raw_test = np.reshape(raw, (px1*px2, nenergy)) # Flatten map
    x_train, x_test = cross_val.train_test_split(raw_test,
                                                test_size=0.997,
                                                random_state=43)
    return (x_train, x_test)

def build_training_set(x_train, V):
    """
    Takes a set of training data and allows the user to interactively classify each dI/dV curve.

    INPUTS:
    x_train = training data outputted from shape_input()

    OUTPUTS:
    y_train = classifications corresponding to x_train
    """
    px, nenergy = x_train.shape
    y_train = np.zeros(px)
    print 'Classify as 1 for superconductivity or \n 2 for spin density wave: \n'
    # For each curve in x_train
    for ctr in range(px):
        plt.clf()
        plt.plot(V, x_train[ctr, :])
        plt.show()
        y_train[ctr] = get_training_classification()
        print str(ctr) + '/' + str(px)
    return y_train

def get_training_classification():
    """
    Have user input a value for y_train. This routine makes sure that the user only inputs a 1 or a 2.

    INPUTS:
    none

    OUTPUTS:
    y_train = value of 1 or 2 which corresponds to whether or not the curve is indiciative of superconductivity or spin density wave, respectively.
    """
    y_train = raw_input('INPUT:\t')
    try:
        y_train = int(y_train)
    except ValueError:
        pass # Do nothing

    if y_train in [1, 2]:
        return y_train
    else:
        print '\nOnly enter a 1 or a 2'
        get_training_classification()


def main():
    raw = scio.loadmat('data/raw.mat')
    raw = raw['raw']
    V = np.linspace(-25, 25, 51)
    x_train, rest_of_map = shape_input(raw)
    y_train = build_training_set(x_train, V)

    print 'Change SDW 2\'s into 0\'s\n'

    y_train[y_train == 2] == 0

    return x_train, y_train, rest_of_map