import urllib.request
import re
import datetime
import time
import pickle
import sys

def record_quotes_real_time(quotes=['EUR/USD','AUD/USD'], sample_interval=15, wbproxy=None):
    #columns = ['time', 'AUD/USD' etc..]
    # one colum per currency forex pair ONLY AVERAGE (ASK + BID PRICE)/2 stored
    data_askquotes = [] # list of dict

    # if Proxy enabled
    if wbproxy is not None:
        proxy_support = urllib.request.ProxyHandler({'http': wbproxy})
        opener = urllib.request.build_opener(proxy_support)
        urllib.request.install_opener(opener)

    # requesting page TrueFx
    login ='u=aflopes7&p=gig1684&q=ozrates&c='+','.join(quotes)+'&f=csv&s=y'
    pgtrfx = 'http://webrates.truefx.com/rates/connect.html?'

    # first connection create session
    response = urllib.request.urlopen(pgtrfx+login)
    session = response.read().decode().strip()

    # to read the txt response (only float values)
    rx = re.compile("[-+]?\d*\.\d+|\d+", re.VERBOSE)

    i = 0
    while True:
        try:
            # after connection
            # request the cote, it is page response in bytes that needs to be converted to 'str'
            pgreturn = urllib.request.urlopen(pgtrfx+'id='+session).read().decode()
            pgreturn = pgreturn.split('\n') # each line is a quote (currency pair)
            quotes_ask = [None] * len(quotes) # empty list of ask price for each cote (currency pair)
            for p in range(len(quotes)): # for each cote read/store its respective line
                # ignore first decimal number is timestamp and high open etc...
                bid, bidpt, ask, askpt = list(map(float, rx.findall(pgreturn[p])))[1:5]
                bid, ask = bid+bidpt*0.00001, ask+askpt*0.00001
                # price is divided in two parts points + ticks
                quotes_ask[p] = (bid+ask)/2 # like IQ option average between two
            dit = dict(zip(['time'] + quotes, [datetime.datetime.now()] + quotes_ask ))
            data_askquotes.append(dit)
            with open('real_time_quotes.pickle', 'wb') as handle:
                pickle.dump(data_askquotes, handle)
        except: # reconnect dont increment
            if wbproxy is not None:
                proxy_support = urllib.request.ProxyHandler({'http': wbproxy})
                opener = urllib.request.build_opener(proxy_support)
                urllib.request.install_opener(opener)
            response = urllib.request.urlopen(pgtrfx+login)
            session = response.read().decode().strip()
        else: # increment and wait time delay
            time.sleep(sample_interval)
            i = i+1

import pandas as pd

def read_real_time_quotes(name='real_time_quotes.pickle'):

    with open(name, 'rb') as handle:
        data = pickle.load(handle)
    dfdata = pd.DataFrame(data)
    dfdata = dfdata.set_index(pd.DatetimeIndex(dfdata['time']))
    dfdata.drop('time', axis='columns', inplace=True)  # not needed anymore
    return dfdata

quotes=['EUR/USD', 'AUD/USD', 'GBP/USD', 'EUR/GBP',
        'USD/CHF', 'USD/CAD', 'EUR/CHF', 'EUR/AUD',
        'AUD/JPY', 'USD/NOK', 'GBP/JPY']


if __name__ == "__main__":
    """
    python 3 background process to record real time quotes

    Optional: pass proxy string as argument if needed

    python3 Forex_Real_Time.py http://usr:password@proxy:port

    """

    wbproxy = None
    if len(sys.argv) > 1:
        wbproxy = sys.argv[1]

    # every 30 s, 120 samples per hour x N hours
    dt=30;
    record_quotes_real_time(quotes=quotes, sample_interval=dt, wbproxy=wbproxy)