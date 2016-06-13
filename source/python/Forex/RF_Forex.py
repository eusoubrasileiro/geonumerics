from zigzag import peak_valley_pivots, max_drawdown, compute_segment_returns, pivots_to_modes
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from IPython.html.widgets import FloatProgress
from IPython.display import display


def plot_quotes(quotes, dataframe): 
    """
    Plot all columns in the daframe, shared axis x
    """   
    f, axr = plt.subplots(len(quotes), sharex=True, figsize=(15,10))
    for i, ax in enumerate(axr):
        dataframe.iloc[:,i].plot(ax=axr[i])
        axr[i].set_ylabel(quotes[i])

def create_binary_zig_zag_class_sm(X, target_quote, span=10, smooth=3, plot=True):
    sd = X[target_quote]
    # using ewm or simple moving average
    #sdsm = sd.ewm(10).mean() # exponential moving average
    sdsm = sd.rolling(span).mean() # simple moving average
    sdsm.fillna(method='bfill', inplace=True)
    sddiff = sd-sdsm # signal change points crossings
    def updown(x):
        if x > 0:
            return 1
        elif x < 0:
            return 0
        else:
            return np.nan
    supdown = pd.Series(np.append(np.array(list(map(updown, (np.diff(np.sign(sddiff)))))), np.nan), 
                        index=sd.index)
    spivots = supdown.dropna()
    supdown.fillna(method='ffill', inplace=True) # fill nans with last valid value
    # remove spikes of one/n sample works only for odd sample windows, mediam works only for odd samples window
    supddown_smooth = supdown.rolling(smooth, center=True).median()

    if plot:
        fig, arx = plt.subplots(2, sharex=True, figsize=(18,6))
        sd[supdown.index].plot(ax=arx[0], style='--b', lw=0.5) # plot just the intersection points
        fig.subplots_adjust(hspace=0)
        sdsm.plot(ax=arx[0], style='g', lw=0.5)
        sd.plot(ax=arx[0], style='-', lw=1.0)
#        supdown.plot(ax=arx[1], style='-b', lw=0.5)          
        supddown_smooth.plot(ax=arx[1], style='-g', lw=1) 
#    supddown_smooth.index.size == supdown.index.size

    X['du'] = supddown_smooth

    return spivots

def create_binary_zig_zag_class(X, target_quote, pts_variation, plot=True, lasttrend=False):
    
    pivots = peak_valley_pivots(X[target_quote].values, pts_variation, -pts_variation)
    spivots = pd.Series(pivots, index=X.index)
    spivots = spivots[pivots != 0] # just the pivots point, when it changes to up or down    

    # remove the last segment of zigzag from training, can not be used it
    # is a misleading trend?? is it?
    if not lasttrend:
        spivots.drop(spivots.tail(1).index, inplace=True) # remove the last point it is not real??

    X['du'] = spivots
    
    if plot:
        f, axr = plt.subplots(2, sharex=True, figsize=(15,4))
        f.subplots_adjust(hspace=0)
        X[target_quote].plot(ax=axr[0])
        X.loc[spivots.index, target_quote].plot(ax=axr[0], style='.r-')


    X['du'].fillna(method='ffill', inplace=True) # fill nans with last valid value

    # turn in [0-1]
    # 1 is up trend, 0 is down trend
    X['du'] = -X['du']
    X.loc[ X['du'] < 0, 'du'] = 0 # where is smaller than zero put 0
    
    
    if plot:
        X['du'].plot(style='-b', ax=axr[1])
        plt.ylim(-0.15, 1.15)
        plt.ylabel('up=1 down=0')
    
    return spivots


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
    # wont be needed, and will create problemns?
    df.drop(target_variable, axis='columns', inplace=True)
    quotes = df.columns
    for quote in quotes:
        df['ema_2 '+quote] = pd.Series.ewm(df[quote], span=2, adjust=True, min_periods=0).mean()
        df['ema_3 '+quote] = pd.Series.ewm(df[quote], span=3, adjust=True, min_periods=0).mean()        
        df['ema_5 '+quote] = pd.Series.ewm(df[quote], span=5, adjust=True, min_periods=0).mean()
        df['ema_10 '+quote] = pd.Series.ewm(df[quote], span=10, adjust=True, min_periods=0).mean()
        df['ema_15 '+quote] = pd.Series.ewm(df[quote], span=15, adjust=True, min_periods=0).mean()
        df['dema_2 '+quote] = df[quote] - df['ema_2 '+quote]
        df['dema_3 '+quote] = df[quote] - df['ema_3 '+quote]
        df['dema_5 '+quote] = df[quote] - df['ema_5 '+quote]
        df['dema_10 '+quote] = df[quote] - df['ema_10 '+quote]
        df['dema_15 '+quote] = df[quote] - df['ema_15 '+quote]
        df['macd '+quote] = MACD(df[quote])
        df['macd_tiny '+quote] = MACD(df[quote], 7, 3)
        df['macd_mean '+quote] = MACD(df[quote], 13, 6)
        df['macd_double '+quote] = MACD(df[quote], 52, 24)
        df['rsi_2 '+quote] = RSI(df[quote], 2)
        df['rsi_3 '+quote] = RSI(df[quote], 3)
        df['rsi_5 '+quote] = RSI(df[quote], 5)
        df['rsi_10 '+quote] = RSI(df[quote], 30)
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
    # random forest binary classifier
        clf = RandomForestClassifier(n_estimators=700, n_jobs=-1)
    else:
    # extra tree binary classifier
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


def performQuickTrainigAnalysis(X_train, y_train, nvalidate, clfmodel=None, windown=60, nanalysis=-1):    
    """
    Cross-Validate the model    

    train the model using a sliding window of size windown

    the training is not progressive cumulative, a new tree every time

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
    elif nanalysis >= pshifts:
        print("Error")
        return

    print('number of analysis: ', nanalysis)

    if clfmodel is None : # in the absence create one
        clfmodel = ExtraTreesClassifier(n_estimators=700, n_jobs=-1)

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

        clfmodel.fit(X_trainFolds, y_trainFolds)        
        accuracies[j] = clfmodel.score(X_testFold, y_testFold)
        iaccuracies[j] = i+windown

    print("percent of data used ", 100*float(iaccuracies[-1]/X_train.index.size))

    return (iaccuracies, accuracies), clfmodel


def performRandomQuickTrainigAnalysis(X_train, y_train, nvalidate, windown=60, nanalysis=-1):    
    """
    Cross-Validate the model    

    train the model using a sliding window of size windown

    the training is not progressive cumulative, a new tree every time

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
    elif nanalysis >= pshifts:
        print("Error")
        return

    print('number of analysis: ', nanalysis)

    clfmodel = ExtraTreesClassifier(n_estimators=700, n_jobs=-1)

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
        
        # create a random begin for the window
        i = np.random.randint(0, X_train.index.size-windown-nvalidate)

        # train using the window samples
        X_trainFolds = X_train[i:i+windown]
        y_trainFolds = y_train[i:i+windown]
        # test using the samples just after the window       
        X_testFold = X_train[i+windown+1:i+windown+nvalidate]
        y_testFold = y_train[i+windown+1:i+windown+nvalidate]    

        clfmodel.fit(X_trainFolds, y_trainFolds)        
        accuracies[j] = clfmodel.score(X_testFold, y_testFold)
        iaccuracies[j] = i+windown

    return (iaccuracies, accuracies), clfmodel