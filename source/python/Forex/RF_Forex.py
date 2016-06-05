from zigzag import peak_valley_pivots, max_drawdown, compute_segment_returns, pivots_to_modes
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from sklearn import preprocessing
from datetime import datetime
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
import operator
import re
from dateutil import parser
import time
import sys
from IPython.html.widgets import FloatProgress
from IPython.display import display
from time import sleep

def plot_quotes(quotes, dataframe): 
    """
    Plot all columns in the daframe, shared axis x
    """   
    f, axr = pyplot.subplots(len(quotes), sharex=True, figsize=(15,10))
    for i, ax in enumerate(axr):
        dataframe.iloc[:,i].plot(ax=axr[i])
        axr[i].set_ylabel(quotes[i])



def create_binary_zig_zag_class(X, target_quote, pts_variation, plot=True, lasttrend=False):
    
    pivots = peak_valley_pivots(X[target_quote].values, pts_variation, -pts_variation)
    spivots = pd.Series(pivots, index=X.index)
    spivots = spivots[pivots != 0] # just the pivots point, when it changes to up or down    

    # remove the last segment of zigzag from training, can not be used it
    # is a misleading trend?? is it?
    if lasttrend == True:
        spivots.drop(spivots.tail(1).index, inplace=True) # remove the last point it is not real??
    
    if plot:
        f, axr = plt.subplots(2, sharex=True, figsize=(15,4))
        f.subplots_adjust(hspace=0)
        X[target_quote].plot(ax=axr[0])
        X.loc[spivots.index, target_quote].plot(ax=axr[0], style='.r-')

    # calculate average interval between up's and downś 
    deltas = pd.Series([spivots.index[i+1]-spivots.index[i] for i in range(spivots.index.size-1)])
    
    mean = np.sum(deltas)/spivots.index.size
    deltamin = np.min(deltas)
    deltamax = np.max(deltas)
        
    #X['du'] = pd.Series(np.zeros(X.index.size)*np.nan, index=X.index)
    # make class of up (1) or down (0) trends from the zigzap indicator
    for i in range(spivots.index.size-1):
        X.loc[ ((X.index >= spivots.index[i]) & 
              (X.index < spivots.index[i+1])), 'du' ] = spivots[i]

    # du is column 2
    #X.iat[-1, 2] = X.iat[-2, 2] # replicate in the end ups or downs

    # turn in [0-1]
    # 1 is up trend, 0 is down trend
    X['du'] = -X['du']
    X.loc[ X['du'] < 0, 'du'] = 0 # where is smaller than zero put 0
    
    
    if plot:
        X['du'].plot(style='-b', ax=axr[1])
        plt.ylim(-0.15, 1.15)
        plt.ylabel('up=1 down=0')
    
    # variations in time between ups and downs of zig/zag
    print("delta mean: ", mean, " max: ", deltamax, " min: ", deltamin)
    return (mean, deltamax, deltamin), spivots


def make_shift_binary_class(X, shift, plot=True):
    """
    Also mark as NaN the last n samples shifted
    """

    # last not NaN Sample
    # for controlling plot, plot 50 samples before for seeing the displacement
    # caused by the shift
    last_sample = X['du'].dropna(inplace=False).index[-50-shift] 

    if plot:
        f, axr = plt.subplots(1, sharex=True, figsize=(15,1))
        X['du'][last_sample:].plot(style='.-b', ax=axr)
        axr.set_ylim(-0.15, 1.15)
    
    X['du_shifted'] = X['du'].shift(-shift) #shift = 1 # shift -1 sample to the left

    # SET NAN last n shifted samples That cannot 
    # be used since they come from FUTURE VALUES we dont know yet.
    X.loc[ X.tail(shift).index, 'du_shifted'] = np.nan


    if plot:
        X['du_shifted'][last_sample:].plot(style='.-r', ax=axr)
        axr.set_ylim(-0.15, 1.15)


    du_shifted = X['du_shifted'] # classification shifted class
    # wont be inserted on this table
    X.drop('du_shifted', axis='columns', inplace=True)
    # wont be needed
    X.drop('du', axis='columns', inplace=True)

    return du_shifted

def MACD(y, a=26, b=12):
    return pd.Series.ewm(y, span=b).mean() - pd.Series.ewm(y, span=a).mean() #12-26

def RSI(y, windown=14):
    dy = y.diff()
    dy.iat[0] = dy.iat[1]
    u = dy.apply(lambda x: x if (x > 0) else 0) # uptrend 0 with where it goes down
    d = dy.apply(lambda x: -x if (x < 0) else 0) # downtred 0 with where it goes up
    # simple exponential moving average
    # pd.Series.ewm(df[quote], span=5, adjust=True, min_periods=0).mean()
    rs = pd.Series.ewm(u, span=windown).mean()/pd.Series.ewm(d, span=windown).mean()
    return 100. - (100./(1.+rs))

# def RSx(y, windown=14):
#     dy = y.diff()
#     dy.iat[0] = dy.iat[1]
#     u = dy.apply(lambda x: x if (x > 0) else 0) # uptrend 0 with where it goes down
#     d = dy.apply(lambda x: -x if (x < 0) else 0) # downtred 0 with where it goes up
#     # simple exponential moving average
#     rs = u.cumsum()
#     return 100. - (100./(1.+rs))

def create_indicators(df, target_variable):
    # wont be needed
    df.drop(target_variable, axis='columns', inplace=True)
    quotes = df.columns
    for quote in quotes:
        df['ema_5 '+quote] = pd.Series.ewm(df[quote], span=5, adjust=True, min_periods=0).mean()
        df['ema_10 '+quote] = pd.Series.ewm(df[quote], span=10, adjust=True, min_periods=0).mean()
        df['ema_15 '+quote] = pd.Series.ewm(df[quote], span=15, adjust=True, min_periods=0).mean()
        df['dema_5 '+quote] = df[quote] - df['ema_5 '+quote]
        df['dema_10 '+quote] = df[quote] - df['ema_10 '+quote]
        df['dema_15 '+quote] = df[quote] - df['ema_15 '+quote]
        df['macd '+quote] = MACD(df[quote])
        df['rsi_5 '+quote] = RSI(df[quote], 5)
        df['rsi_20 '+quote] = RSI(df[quote], 20)
        df['rsi_30 '+quote] = RSI(df[quote], 30)
    for i in range(len(quotes)-1):
        df[quotes[i]+'+'+quotes[i+1]] = df[quotes[i]]+df[quotes[i+1]]

def calculate_corr_and_remove(df, serie_binary_class, verbose=True):
    """
    calculate correlation with binary classification serie
    collumns less then 5% are removed and wont be used on
    random forest
    """
    df['target_variable_dup'] = pd.Series(serie_binary_class, index=df.index)
    # drop NANs cannot be used for correlation
    dfcpy = df.dropna(inplace=False)
    corr = dfcpy.corr().ix['target_variable_dup', :-1]
    corr = corr.apply(lambda x: np.abs(x)) # make all correlations positive
    corr.sort_values(ascending=False, inplace=True)
    if verbose:
        print(corr.head(5))
        print(corr.tail(5))
    df.drop('target_variable_dup', axis='columns', inplace=True)
    df.drop(corr.index[corr < 0.05], axis='columns', inplace=True)

def performRFClass(X_train, y_train, X_test, y_test, algorithm):
    """
    Random Forest Binary Classification
    """

    if algorithm == 'RF':
        clf = RandomForestClassifier(n_estimators=700, n_jobs=-1)
    else:
        clf = ExtraTreesClassifier(n_estimators=700, n_jobs=-1)

    #print(len(X_train))
    clf.fit(X_train, y_train)
    
    accuracy = clf.score(X_test, y_test)
    
    return accuracy


def prepareDataForClassification(dataset):
    """
    generates categorical to be predicted column, 
    attach to dataframe 
    and label the categories
    """
    le = preprocessing.LabelEncoder()
    
    dataset['UpDown'] = dataset['target_variable_dup']
    dataset.UpDown = le.fit(dataset.UpDown).transform(dataset.UpDown) #create classes for randomForest
    
    features = dataset.columns[0:-2] # last two are the target variable. IGNORE THEM
    X = dataset[features]    
    y = dataset.UpDown    
    
    return X, y


def performTimeSeriesCV(X_train, y_train, number_folds):
    """
    Given X_train and y_train (the test set is excluded from the Cross Validation),
    number of folds, the ML algorithm to implement and the parameters to test,
    the function acts based on the following logic: it splits X_train and y_train in a
    number of folds equal to number_folds. Then train on one fold and tests accuracy
    on the consecutive as follows:
    - Train on fold 1, test on 2
    - Train on fold 1-2, test on 3
    - Train on fold 1-2-3, test on 4
    ....
    Returns mean of test accuracies.
    """
 
    print('Size train set: ', X_train.shape)
    
    # k is the size of each fold. It is computed dividing the number of 
    # rows in X_train by number_folds. This number is floored and coerced to int
    k = int(np.floor(float(X_train.shape[0]) / number_folds))
    print('Size of each fold: ', k)
    
    # initialize to zero the accuracies array. It is important to stress that
    # in the CV of Time Series if I have n folds I test n-1 folds as the first
    # one is always needed to train
    accuracies = np.zeros(number_folds-1)
 
    # loop from the first 2 folds to the total number of folds    
    for i in range(2, number_folds + 1):
        print('')
        
        # the split is the percentage at which to split the folds into train
        # and test. For example when i = 2 we are taking the first 2 folds out 
        # of the total available. In this specific case we have to split the
        # two of them in half (train on the first, test on the second), 
        # so split = 1/2 = 0.5 = 50%. When i = 3 we are taking the first 3 folds 
        # out of the total available, meaning that we have to split the three of them
        # in two at split = 2/3 = 0.66 = 66% (train on the first 2 and test on the
        # following)
        split = float(i-1)/i
        
        # example with i = 4 (first 4 folds):
        #      Splitting the first       4        chunks at          3      /        4
        print('Splitting the first ' + str(i) + ' chunks at ' + str(i-1) + '/' + str(i))
        
        # as we loop over the folds X and y are updated and increase in size.
        # This is the data that is going to be split and it increases in size 
        # in the loop as we account for more folds. If k = 300, with i starting from 2
        # the result is the following in the loop
        # i = 2
        # X = X_train[:(600)]
        # y = y_train[:(600)]
        #
        # i = 3
        # X = X_train[:(900)]
        # y = y_train[:(900)]
        # .... 
        X = X_train[:(k*i)]
        y = y_train[:(k*i)]
        print('Size of train + test: ', X.shape) # the size of the dataframe is going to be k*i
 
        # X and y contain both the folds to train and the fold to test.
        # index is the integer telling us where to split, according to the
        # split percentage we have set above
        index = int(np.floor(X.shape[0] * split))
        
        # folds used to train the model        
        X_trainFolds = X[:index]        
        y_trainFolds = y[:index]
        
        # fold used to test the model
        X_testFold = X[(index + 1):(index + 4)]
        y_testFold = y[(index + 1):(index + 4)]
        print('size of samples to test ', len(y_testFold))
        print('first index of tested sample ', index+1)        
        
        # i starts from 2 so the zeroth element in accuracies array is i-2. performClassification() is a function which takes care of a classification problem. This is only an example and you can replace this function with whatever ML approach you need.
        accuracies[i-2] = performRFClass(X_trainFolds, y_trainFolds, X_testFold, y_testFold, ' ')


def performCV_analysis(X_train, y_train, windown, nanalysis=-1, nvalidate=2, algorithm='ET'):    
    """
    Cross-Validate the model    

    train the model using a sliding window of size windown

    the training is progressive cumulative

    and will produce nanalysis (array) of size number of window movements
    """

    # nvalidate  number of samples to use for validation
    pshifts = round(X_train.index.size-windown+1-nvalidate) # possible shifts N-windown+1    

    print('Size train set: ', X_train.shape)
    print('samples in each window: ', windown)

    if nanalysis < 0:
        # number of samples default, equal number of 
        # windows inside data
        nanalysis = round((X_train.index.size-nvalidate)/windown)
        print('number of analysis: ', nanalysis)
    elif nanalysis >= pshifts:
        print("Error")
        return

    if algorithm == 'RF': # classification model
        clf = RandomForestClassifier(n_estimators=700, n_jobs=-1)
    else:
        clf = ExtraTreesClassifier(n_estimators=700, n_jobs=-1)

    step = int(round(pshifts/nanalysis)) # window shift step in samples
    #diff = pshifts-
    shifts = range(0, pshifts, step)
    accuracies = np.zeros(len(shifts)) # store the result of cross-validation
    iaccuracies = np.zeros(len(shifts)) # index of the last sample in the window

    f = FloatProgress(min=0, max=nanalysis)
    display(f)

    # sliding window with step sample of shift     
    for j, i in enumerate(shifts): # shift, classify and cross-validate
        f.value = j # counter of analysis

        # train using the window samples
        X_trainFolds = X_train[i:i+windown]
        y_trainFolds = y_train[i:i+windown]
        # test using the samples just after the window       
        X_testFold = X_train[i+windown+1:i+windown+nvalidate]
        y_testFold = y_train[i+windown+1:i+windown+nvalidate]    

        clf.fit(X_trainFolds, y_trainFolds)        
        accuracies[j] = clf.score(X_testFold, y_testFold)
        iaccuracies[j] = i+windown

    return (iaccuracies, accuracies), clf


def performCV_shift(X_train, y_train, windown, algorithm, ntest):    
    print('Size train set: ', X_train.shape)
    print('samples in each window: ', windown)
    ntest = ntest  # number of samples to estimate
    nshifts = len(y_train)-windown+1-ntest # possible shifts N-windown+1    
    print('number of shifts', nshifts)

    accuracies = np.zeros(nshifts)  
    
    f = FloatProgress(min=0, max=nshifts)
    display(f)

    if algorithm == 'RF':
        clf = RandomForestClassifier(n_estimators=700, n_jobs=-1)
    else:
        clf = ExtraTreesClassifier(n_estimators=700, n_jobs=-1)

    # sliding window with one sample of shift
    for i in range(nshifts): # shift, classify and estimate/test
        f.value = i
        # train using the window samples
        X_trainFolds = X_train[i:i+windown]
        y_trainFolds = y_train[i:i+windown]
        # test using the samples just after the window       
        X_testFold = X_train[i+windown+1:i+windown+ntest]
        y_testFold = y_train[i+windown+1:i+windown+ntest]    

        #print(len(X_trainFolds), len(y_trainFolds))        
        #print(len(X_testFold), len(y_testFold))       
        #print(len(X_train))
        clf.fit(X_trainFolds, y_trainFolds)        
        accuracies[i] = clf.score(X_testFold, y_testFold)
        #accuracies[i] = performRFClass(X_trainFolds, y_trainFolds, X_testFold, y_testFold, algorithm)

    return accuracies, clf
