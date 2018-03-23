###########################################################################
######################## EISAGWGI VIVLIOTHIKWN ############################
###########################################################################

import os.path

import re

import numpy as np # vivliothiki statistikis analysis
from numpy import ndarray
import numpy.testing


import numpy.core
import numpy.core.multiarray
from numpy import pi

import setuptools

import matplotlib
import matplotlib.pyplot as pyp
import matplotlib.pylab as plt

import pandas as pd

from scipy.linalg import toeplitz

import sklearn
from sklearn import metrics
from sklearn.preprocessing import normalize
from sklearn.metrics.cluster import mutual_info_score

from scipy import ndimage
from scipy.stats import chisquare
from scipy import stats
from scipy.linalg import det


import statsmodels.tsa


###########################################################################
#################### VOITHITIKES SYNARTISEIS ##############################
###########################################################################

def grangercausalitytests(x, maxlag, addconst=True, verbose=True):
    """four tests for granger non causality of 2 timeseries

    all four tests give similar results
    `params_ftest` and `ssr_ftest` are equivalent based on F test which is
    identical to lmtest:grangertest in R

    Parameters
    ----------
    x : array, 2d
        data for test whether the time series in the second column Granger
        causes the time series in the first column
    maxlag : integer
        the Granger causality test results are calculated for all lags up to
        maxlag
    verbose : bool
        print results if true

    Returns
    -------
    results : dictionary
        all test results, dictionary keys are the number of lags. For each
        lag the values are a tuple, with the first element a dictionary with
        teststatistic, pvalues, degrees of freedom, the second element are
        the OLS estimation results for the restricted model, the unrestricted
        model and the restriction (contrast) matrix for the parameter f_test.

    Notes
    -----
    TODO: convert to class and attach results properly

    The Null hypothesis for grangercausalitytests is that the time series in
    the second column, x2, does NOT Granger cause the time series in the first
    column, x1. Grange causality means that past values of x2 have a
    statistically significant effect on the current value of x1, taking past
    values of x1 into account as regressors. We reject the null hypothesis
    that x2 does not Granger cause x1 if the pvalues are below a desired size
    of the test.

    The null hypothesis for all four test is that the coefficients
    corresponding to past values of the second time series are zero.

    'params_ftest', 'ssr_ftest' are based on F distribution

    'ssr_chi2test', 'lrtest' are based on chi-square distribution

    References
    ----------
    http://en.wikipedia.org/wiki/Granger_causality
    Greene: Econometric Analysis

    """
    from scipy import stats
    from statsmodels.tsa.tsatools import lagmat2ds
    from statsmodels.tools.tools import add_constant
    from statsmodels.regression.linear_model import OLS

    x = np.asarray(x)

    if x.shape[0] <= 3 * maxlag + int(addconst):
        raise ValueError("Insufficient observations. Maximum allowable "
                         "lag is {0}".format(int((x.shape[0] - int(addconst)) /
                                                 3) - 1))

    resli = {}

    for mlg in range(1, maxlag + 1):
        result = {}
#        if verbose:
#            print('\nGranger Causality')
#            print('number of lags (no zero)', mlg)
        mxlg = mlg

        # create lagmat of both time series
        dta = lagmat2ds(x, mxlg, trim='both', dropex=1)

        #add constant
        if addconst:
            dtaown = add_constant(dta[:, 1:(mxlg + 1)], prepend=False)
            dtajoint = add_constant(dta[:, 1:], prepend=False)
        else:
            raise NotImplementedError('Not Implemented')
            #dtaown = dta[:, 1:mxlg]
            #dtajoint = dta[:, 1:]

        # Run ols on both models without and with lags of second variable
        res2down = OLS(dta[:, 0], dtaown).fit()
        res2djoint = OLS(dta[:, 0], dtajoint).fit()

        #print results
        #for ssr based tests see:
        #http://support.sas.com/rnd/app/examples/ets/granger/index.htm
        #the other tests are made-up

        # Granger Causality test using ssr (F statistic)
        fgc1 = ((res2down.ssr - res2djoint.ssr) /
                res2djoint.ssr / mxlg * res2djoint.df_resid)
#        if verbose:
#            print('ssr based F test:         F=%-8.4f, p=%-8.4f, df_denom=%d,'
#                   ' df_num=%d' % (fgc1,
#                                    stats.f.sf(fgc1, mxlg,
#                                               res2djoint.df_resid),
#                                    res2djoint.df_resid, mxlg))
        result['ssr_ftest'] = (fgc1,
                               stats.f.sf(fgc1, mxlg, res2djoint.df_resid),
                               res2djoint.df_resid, mxlg)

        # Granger Causality test using ssr (ch2 statistic)
        fgc2 = res2down.nobs * (res2down.ssr - res2djoint.ssr) / res2djoint.ssr
#        if verbose:
#            print('ssr based chi2 test:   chi2=%-8.4f, p=%-8.4f, '
#                   'df=%d' % (fgc2, stats.chi2.sf(fgc2, mxlg), mxlg))
        result['ssr_chi2test'] = (fgc2, stats.chi2.sf(fgc2, mxlg), mxlg)

        #likelihood ratio test pvalue:
        lr = -2 * (res2down.llf - res2djoint.llf)
#        if verbose:
#            print('likelihood ratio test: chi2=%-8.4f, p=%-8.4f, df=%d' %
#                   (lr, stats.chi2.sf(lr, mxlg), mxlg))
        result['lrtest'] = (lr, stats.chi2.sf(lr, mxlg), mxlg)

        # F test that all lag coefficients of exog are zero
        rconstr = np.column_stack((np.zeros((mxlg, mxlg)),
                                   np.eye(mxlg, mxlg),
                                   np.zeros((mxlg, 1))))
        ftres = res2djoint.f_test(rconstr)
#        if verbose:
#            print('parameter F test:         F=%-8.4f, p=%-8.4f, df_denom=%d,'
#                   ' df_num=%d' % (ftres.fvalue, ftres.pvalue, ftres.df_denom,
#                                    ftres.df_num))
        result['params_ftest'] = (np.squeeze(ftres.fvalue)[()],
                                  np.squeeze(ftres.pvalue)[()],
                                  ftres.df_denom, ftres.df_num)

        resli[mxlg] = (result, [res2down, res2djoint, rconstr])

    return resli

###########   AMOIVAIA PLIROFORIA & ELEGXOS  ##############

def entropy_gaussian(C):
    '''
    Entropy of a gaussian variable with covariance matrix C
    '''
    if np.isscalar(C): # C is the variance :  FALSE ARA PROSPERNAEI
        return .5*(1 + np.log(2*pi)) + .5*np.log(C)
    else:
        n = C.shape[0] # dimension : EDW MPAINEI KAI YPOLOGIZEI
        return .5*n*(1 + np.log(2*pi)) + .5*np.log(abs(det(C)))  # ayto pairnei kai vazei gia n to size tou cov matrix = 2

def mutual_information_2d(x, y, sigma=1, normalized=False):
    """
    Computes (normalized) mutual information between two 1D variate from a
    joint histogram.
    Parameters
    ----------
    x : 1D array
        first variable
    y : 1D array
        second variable
    sigma: float
        sigma for Gaussian smoothing of the joint histogram
    Returns
    -------
    nmi: float
        the computed similariy measure
    """
    EPS = np.finfo(float).eps

    bins = (256, 256)

    jh = np.histogram2d(x, y, bins=bins)[0]

    # smooth the jh with a gaussian filter of given sigma
    ndimage.gaussian_filter(jh, sigma=sigma, mode='constant',
                                 output=jh)

    # compute marginal histograms
    jh = jh + EPS
    sh = np.sum(jh)
    jh = jh / sh
    s1 = np.sum(jh, axis=0).reshape((-1, jh.shape[0]))
    s2 = np.sum(jh, axis=1).reshape((jh.shape[1], -1))

    # Normalised Mutual Information of:
    # Studholme,  jhill & jhawkes (1998).
    # "A normalized entropy measure of 3-D medical image alignment".
    # in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
    if normalized:
        mi = ((np.sum(s1 * np.log(s1)) + np.sum(s2 * np.log(s2)))
                / np.sum(jh * np.log(jh))) - 1
    else:
        mi = ( np.sum(jh * np.log(jh)) - np.sum(s1 * np.log(s1))
               - np.sum(s2 * np.log(s2)))

    return mi


def test_mutual_information_2d(timeseries1,timeseries2):
    # Mutual information between two correlated gaussian variables
    # Entropy of a 2-dimensional gaussian variable
    n = 5000
    rng = np.random.RandomState(0)
    #P = np.random.randn(2, 2)
    P = np.array([timeseries1,timeseries2])   # 2X72
    C = np.dot(P, P.T) # Pollaplasiazei (.T anastrofos) ---> C=P*PT 2X2  COV matrix
    U = rng.randn(2, n) # dimiourgei pinaka 2X5000 times apo normal distribution times
    Z = np.dot(P.T, U)  #pollaplasiazei toys pinakes PT,U & dinei to apotelesma 71X5000
    X = Z[:, 0] # pairnei tin prwti grammi toy Z (72X5000)
    X = X.reshape(len(X), 1) # ton kanei pinaka stili  ----> 72X1
    Y = Z[:, 1] #pairnei ti deyteri grammi tou Z 
    Y = Y.reshape(len(Y), 1) #antistoixa  ----> 72X1
    # in bits
    MI_est = mutual_information_2d(X.ravel(), Y.ravel()) # Ftiaxnei ton 72X72
    MI_th = (entropy_gaussian(C[0, 0])
             + entropy_gaussian(C[1, 1])
             - entropy_gaussian(C)
            )   # Entropy of a gaussian variable with covariance matrix C

    #Symfwna me tis symeiwseis toy i MI_th einai i mutual information kai einai symmetriki
    #Otan ayti einai poly konta sto 0 kai katw apo to p=0.01 tote einai Aneksartita kai vazei 0
    #Oso megalyteri einai i timi toys toso megalyteri einai kai i sxesi toys kai vaze 1
    #Ixy=H(X)+H(Y)-H(X,Y)   :   H(X)=1/2log[(2*pi*e)^n*|C|] opoy C o COV matrix, kai n oi dyastaseis toy
    #                           H(X)=1/2[log(2*pi*e)^n+log|C|]
    #                           H(X)=1/2*n*[log(2*pi*e)]+1/2*log(det(C))]
    #                           H(X)=0.5*n*[loge+log(2*pi)]+0.5log(det(C))]
    #                           H(X)=0.5*n*[1+log(2*pi)]+0.5*log(det(C))]

    #Parousiasi 4 - Selida 9  Kougioymtzis
    #Typos entropias : tema1bwp.pdf -- selida : 10/24
    #Synthiki syndesis : coverch8.pdf -- selida 25/35
    
    return  (MI_est, MI_th)  #Epistrefei tin ektimitea mut_info & tin entropia
    # Our estimator should undershoot once again: it will undershoot more
    # for the 2D estimation that for the 1D estimation
    #np.testing.assert_array_less(MI_est, MI_th)
    #np.testing.assert_array_less(MI_th, MI_est + .2)

#####################################################################
################### ELEGXOS APALOIFIS TASIS #########################
                 

def test_stationarity(timeseries):

    from statsmodels.tsa.stattools import adfuller
    
    #Determing rolling statistics
    rolmean = pd.rolling_mean(timeseries, window=12)
    rolstd = pd.rolling_std(timeseries, window=12)

    #Plot rolling statistics:
    orig = plt.plot(timeseries, color='blue',label='Original')
    mean = plt.plot(rolmean, color='red', label='Rolling Mean')
    std = plt.plot(rolstd, color='black', label = 'Rolling Std')
    plt.legend(loc='best')
    plt.title('Rolling Mean & Standard Deviation')
    plt.show(block=False)
    
    #Perform Dickey-Fuller test:
    print('Results of Dickey-Fuller Test :')
    dftest = adfuller(timeseries, autolag='AIC')
    dfoutput = pd.Series(dftest[0:4], index=['Test Statistic','p-value','#Lags Used','Number of Observations Used'])
    for key,value in dftest[4].items():
        dfoutput['Critical Value (%s)'%key] = value
    print(dfoutput) 

    # The test statistic is smaller than the critical values of 1% so with 99% evidence we can say it is a stationary timeserie
    # P-value is also 0.002 ~ 0 so we can surely reject the null ypothesis of the Test
    # antistoixa kai gia alla parathyra pou elegxthikan
    # To diastima epilogis tasis einai 12mines = etisia tasi
    

#######################################################################
#################### CROSS CORRELATION WITH LAG  ######################
                 
def lagcorr(x,y,lag=None,verbose=True):
    '''
    Compute lead-lag correlations between 2 time series.

    <x>,<y>: 1-D time series.
    <lag>: lag option, could take different forms of <lag>:
          if 0 or None, compute ordinary correlation and p-value;
          if positive integer, compute lagged correlation with lag
          upto <lag>;
          if negative integer, compute lead correlation with lead
          upto <-lag>;
          if pass in an list or tuple or array of integers, compute 
          lead/lag correlations at different leads/lags.

    Note: when talking about lead/lag, uses <y> as a reference.
    Therefore positive lag means <x> lags <y> by <lag>, computation is
    done by shifting <x> to the left hand side by <lag> with respect to
    <y>.
    Similarly negative lag means <x> leads <y> by <lag>, computation is
    done by shifting <x> to the right hand side by <lag> with respect to
    <y>.

    Return <result>: a (n*2) array, with 1st column the correlation 
    coefficients, 2nd column correpsonding p values.

    Currently only works for 1-D arrays.
    '''
    import numpy
    from scipy.stats import pearsonr

    if len(x)!=len(y):
        raise('Input variables of different lengths.')

    #--------Unify types of <lag>-------------
    if numpy.isscalar(lag):
        if abs(lag)>=len(x):
            raise('Maximum lag equal or larger than array.')
        if lag<0:
            lag=-numpy.arange(abs(lag)+1)
        elif lag==0:
            lag=[0,]
        else:
            lag=numpy.arange(lag+1)    
    elif lag is None:
        lag=[0,]
    else:
        lag=numpy.asarray(lag)

    #-------Loop over lags---------------------
    result=[]
#    if verbose:
#        print ('\n#<lagcorr>: Computing lagged-correlations at lags:',lag)

    for ii in lag:
        if ii<0:
            result.append(pearsonr(x[:ii],y[-ii:]))
        elif ii==0:
            result.append(pearsonr(x,y))
        elif ii>0:
            result.append(pearsonr(x[ii:],y[:-ii]))

    result=numpy.asarray(result)

    return result   #result=(Pearson correlation coefficient,p-value)
                    # cor=[-1,1] --> -1=antistrofi metavoli (thetiki sysx)
                    #            --> +1=analogi metavoli    (arnitiki sysx)
                    #            --> 0 =den paratireitai sysxetisi


###########################################################################
######################## KYRIOS PROGRAMMA #################################
###########################################################################

#### EISAGWGI TWN DEDOMENWN APO TA ARXEIA & APOTHIKEYSI TOYS SE LISTES ####

os.chdir('/home/filippos/Python/kugiu')
with open('WBCommodityUSDnames.dat','r',encoding='utf-8') as n:
    names=n.read().split('\n')
    del names[-1]
prices=[]
with open('WBCommodityUSDprices.dat','r',encoding='utf-8') as p:
    prices=p.read().split('\n')
    prices.pop(-1)

############# DIORTHWSI LATHWN ###############

for i,temp in enumerate(prices):
    prices[i]=str(temp).split('   ')
for i in range(len(prices)):
    prices[i].remove('')

############# Lista 50 onomatwn : names[metoxi]
ln=len(names)
############## Lista 374imerwn X 50 : prices[imera][metoxi]
li=len(prices)

############### Epilogi 20 metoxwn kai ayksousa taksinomisi
ep=[5,6,7,8,9,10,11,30,31,32,12,14,16,17,18,23,24,25,26,27]
#ep.sort()  proairetika
d={}
stocks=[]
print('Oi 20 deiktes agathwn pou epilexthikan einai :')
for i in ep:
    d[i]=names[i]
    stocks.append(names[i])
    print(d[i])


######## DIMIOYRGIA TWN PARATHYRWN TVN XRONOSEIRWN ANA METOXI #########

####### Epilogi Diastimatos : 6 xronia
years=6
################# Ypologismos parathyrwn
w6=li//(years*12)

start=[]
stop=[]
for i in range(1,w6+1):
    start.append((i-1)*(years*12)+1)
    stop.append(i*(years*12))

###### DIMIOYRGIA LISTAS TIMWN P KAI LEXIKOY PDICT= {'ONOMA':PARATHYRO TIMWN }
p=[]  #p[parathyro][metoxi][imera]
pdict={}
for k in range(0,w6): 
    y=[]
    for e in ep:
        x=[]
        for i in range(start[k],stop[k]+1): 
            x.append(float(prices[i-1][e]))
        y.append(x)
    p.append(y) #H lista me tis 15 xronoseires ana metoxi
    pdict[e]=p  #Lexiko me onoma kai ti lista gia na kseroume poia einai
'''
########## PLOT TWN METOXWN TOY 1ou PARATHYROY KAI OLWN TWN XRONOSEIRWN
for i in range(20) : pyp.plot(p[0][i])
pyp.show()
'''

###########################################################################
######################### PRE WHITENING ###################################
###########################################################################

prewhite=[]

for q in range(w6):

    ########## Ypologismos kai aferesi tis tasis

    prewhite_temp=[]

    for w in range(20):

        ########## Diadikasia Pre - whitening  #############

        #transformation :
        #       Penalising higher values more than smaller values
        #       log, square root, cube root

        ts_log = np.log(p[q][w])

        '''
        pyp.plot(ts_log) # log transform 

        #technics to model the trend :
        #       Aggregation (OK), Smoothing , Polynomial Fitting
        '''

        moving_avg = pd.rolling_mean(ts_log,2)

        #petontas tis teleytaies 2 times afou ypologizei mesous orous ana 2 times

        '''
        pyp.plot(ts_log)
        pyp.plot(moving_avg, color='red')
        pyp.show()
        '''
        
        #Observation of first 2 means
        ts_log_moving_avg_diff = ts_log - moving_avg
        ts_log_moving_avg_diff[0:2]

        #droping first 1 means pou einai Nan
        ts_log_moving_avg_diff=pd.Series(ts_log_moving_avg_diff).dropna()

        ts_log_moving_avg_diff.dropna(inplace=True)

        ts_log_moving_avg_diff=np.asarray(ts_log_moving_avg_diff)

        prewhite_temp.append(ts_log_moving_avg_diff)

    #   test_stationarity(prewhite[][])
    #       den trexei apeytheias alla mono toy sto compiler
    #       Dokimazontas gia meses times ana 2,3,4,5,6... kataligw oti logw megethous xronoseiras
    #       o pio katalilos mesos einai ana 2..vasi tou p-value pou  kai toy critival value pou
    #       se sxesi me to critical test vgainei synexws katw apo 5%

    prewhite.append(prewhite_temp) #ana parathyro

#############################################################################################
####                     cross correlation gia kathe parathyro :                        #####
####                     ------  XWRIS & ME prolefkansi  ------                             #####
####                             Pearson2_cor = coef                                    #####
####                             & Pearson2_p = p-test                                  #####
####                           & Ysterisi lag = 0   kai   1                             #####
#############################################################################################

Parathyro_undirected={}

################  LAG 0

Pearson60_cor=[]
Pearson60_cor_white=[]
Pearson60_p=[]

ccvx=[]
ccvy=[]
ccv=[]


for s in range(0,w6): # ana parathyro xronoseiras
    (y1,y2)=([],[])
    for j in range(20):
        (x1,x2)=([],[])
        for i in range(20): 
            (m,l)=([],[])
            m=lagcorr(prewhite[s][j],prewhite[s][i])   #######  me prewhite  #######
            l=lagcorr(p[s][j],p[s][i])                 ####### xwris prewhite ######
            ccvx=m[0][0]
            if m[0][1]<=0.01 : x1.append(int(1))       # 1 gia statistika simantiki >(a=0.01)
            elif m[0][1]>0.01: x1.append(int(0))       # 0 gia statistika asimanti <(a=0.01)
            if l[0][1]<=0.01 : x2.append(int(1))
            elif l[0][1]>0.01 : x2.append(int(0))
    #       x[j] =  oi syndeseis tou j me ola ta alla            
        y1.append(x1)         #  Apothieyei oles tis syndeseis tou parathyrou se pinaka
        y2.append(x2)
        ccvy.append(ccvx)   #   Apothikeyei tis cross cov gia elegxo
                            
    Pearson60_cor_white.append(y1) # Lista pinakwn syndesewm ana parathyro
    Pearson60_cor.append(y2)

    ccv.append(ccvy)
    Parathyro_undirected[str(s)]=y

################ LAG 1

Parathyro_directed={}

Pearson61_cor=[]
Pearson61_cor_white=[]
Pearson61_p=[]

ccvx=[]
ccvy=[]
ccv=[]


for s in range(0,w6): # ana parathyro xronoseiras
    (y1,y2)=([],[])
    for j in range(20):
        (x1,x2)=([],[])
        for i in range(20): 
            (m,l)=([],[])
            m=lagcorr(prewhite[s][j],prewhite[s][i],lag=1)   #######  me prewhite  #######
            l=lagcorr(p[s][j],p[s][i],lag=1)                 ####### xwris prewhite ######
            ccvx=m[0][0]
            if m[1][1]<=0.01 : x1.append(int(1))       # 1 gia statistika simantiki >(a=0.01)
            elif m[1][1]>0.01: x1.append(int(0))       # 0 gia statistika asimanti <(a=0.01)
            if l[1][1]<=0.01 : x2.append(int(1))
            elif l[1][1]>0.01 : x2.append(int(0))
    #       x[j] =  oi syndeseis tou j me ola ta alla            
        y1.append(x1)         #  Apothieyei oles tis syndeseis tou parathyrou se pinaka
        y2.append(x2)
        ccvy.append(ccvx)   #   Apothikeyei tis cross cov gia elegxo
                            
    Pearson61_cor_white.append(y1) # Lista pinakwn syndesewm ana parathyro
    Pearson61_cor.append(y2)

    ccv.append(ccvy)
    Parathyro_directed[str(s)]=y


##########################################################################
##### PLOT MEMONOMENWN & SYGKENTRWTIKWN PARATHYRWN GIA PREWHITENED #######
#####              1.   LAG = 0   KAI   LAG = 1                    #######
#####        2. XWRIS PROLEFKANSI   KAI    ME PROLEFKANSI          #######
#####         3.  SYGKENTROTIKWN GIA TIN KATHE PERIPTWSI           #######
##########################################################################


####################       LAG = 0     #######################

############## PLOT OLWN TWN PARATHYRWN XWRIS PREWHITE

for i in range(0,w6):
    matrix = Pearson60_cor[i]  # symmetrikos pinakas 0-1

    fig=plt.figure()

    ax=fig.add_subplot(1,1,1)
    ax.set_aspect('equal')

    plt.imshow(matrix,interpolation='nearest',cmap=plt.cm.Reds)
    plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
    plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
    fig.canvas.set_window_title(str(i+1)+'o Parathyro CROSS CORRELATION xwris PREWHITE me LAG=0')
    plt.title('Starts:'+str(start[i])+'-  Stops:'+str(stop[i]))
    plt.colorbar()
    plt.show()

############## SYGKENTRWTIKO XWRIS PREWHITEN ME LAG=0
x=[[0]*20]*20
for i in range (0,w6):
    x+=np.asarray(Pearson60_cor[i])

meg=np.amax(x)
norm_matrix=x.tolist()/meg

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(norm_matrix,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro diasysxetisis me LAG = 0')
plt.title('Sygkentrwtikes normalized syndeseis i--j XWRIS PREWHITE xwris YSTERISI sta {} parathyra'.format(w6))
plt.colorbar()
plt.show()

# Krataw mono osous parousiasan panw apo ta misa parathyra syndeseis
y=(norm_matrix>0.5).astype(int)

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(y,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro diasysxetisis me LAG = 0')
plt.title('Osoi emfanisan syndesi i--j se perissotera apo ta misa parathyra')
plt.colorbar()
plt.show()




############## PLOT OLWN TWN PARATHYRWN ME PREWHITE ME LAG = 0

for i in range(0,w6):
    matrix = Pearson60_cor_white[i]  # symmetrikos pinakas 0-1

    fig=plt.figure()

    ax=fig.add_subplot(1,1,1)
    ax.set_aspect('equal')

    plt.imshow(matrix,interpolation='nearest',cmap=plt.cm.Reds)
    plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
    plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
    fig.canvas.set_window_title(str(i+1)+'o Parathyro CROSS CORRELATION me PREWHITE me LAG=0')
    plt.title('Starts:'+str(start[i])+'-  Stops:'+str(stop[i]))
    plt.colorbar()
    plt.show()


############## SYGKENTRWTIKO ME PREWHITE ME LAG=0

x=[[0]*20]*20
for i in range (0,w6):
    x+=np.asarray(Pearson60_cor_white[i])

meg=np.amax(x)
norm_matrix=x.tolist()/meg

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(norm_matrix,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro CROSS CORRELATION ME PREWHITE me LAG = 0')
plt.title('Sygkentrwtikes normalized syndeseis i--j XWRIS YSTERISI sta {} parathyra'.format(w6))
plt.colorbar()
plt.show()

# Krataw mono osous parousiasan panw apo ta misa parathyra syndeseis
y=(norm_matrix>0.5).astype(int)

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(y,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro diasysxetisis me LAG = 0')
plt.title('Osoi emfanisan syndesi i--j se perissotera apo ta misa parathyra')
plt.colorbar()
plt.show()


####################       LAG = 1     #######################

############## PLOT OLWN TWN PARATHYRWN ME PREWHITE ME LAG = 1

for i in range(0,w6):
    matrix = Pearson61_cor_white[i]  # symmetrikos pinakas 0-1

    fig=plt.figure()

    ax=fig.add_subplot(1,1,1)
    ax.set_aspect('equal')

    plt.imshow(matrix,interpolation='nearest',cmap=plt.cm.Reds)
    plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
    plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
    fig.canvas.set_window_title(str(i+1)+'o Parathyro CROSS CORRELATION ME PREWHITE me LAG = 1')
    plt.title('Starts:'+str(start[i])+'-  Stops:'+str(stop[i]))
    plt.colorbar()
    plt.show()


############## SYGKENTRWTIKO ME PREWHITEN ME LAG=1
x=[[0]*20]*20
for i in range (0,w6):
    x+=np.asarray(Pearson61_cor_white[i])

meg=np.amax(x)
norm_matrix=x.tolist()/meg

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(norm_matrix,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro diasysxetisis me ysterisi 0')
plt.title('Sygkentrwtikes normalized syndeseis i-->j ME YSTERISI me LAG = 1 sta {} parathyra'.format(w6))
plt.colorbar()
plt.show()

# Krataw mono osous parousiasan panw apo ta misa parathyra syndeseis
y=(norm_matrix>0.5).astype(int)

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(y,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro diasysxetisis me LAG = 1')
plt.title('Osoi emfanisan syndesi i-->j se perissotera apo ta misa parathyra')
plt.colorbar()
plt.show()

###########################################################################
###################### AMOIVAIA PLIROFORIA ################################
###########################################################################

############### MUTUAL INFORMATION GIA PREWHITEND TIMESERIES

mut_inf_lag0=[]
mut_inf_lag1=[]

(temp1,temp2,temp3,temp4)=([],[],[],[])

################### ME YSTERISI MIDEN : LAG=0
for s in range(w6):
    temp2=[]
    for j in range(20):
        temp1=[]
        for i in range(20):      
            if i==j : temp1.append(0) #mideniki syndesi sti diagwnio
            else:
                temp1.append(mutual_information_2d(prewhite[s][j],prewhite[s][i]))
        temp2.append(temp1)
    mut_inf_lag0.append(temp2)

################### ME YSTERISI ENA : LAG=1
for s in range(w6):
    temp4=[]
    for j in range(20):
        temp3=[]
        for i in range(20):
            if i==j: temp3.append(0) #mideniki sindesi sti diagwnio
            else:
                temp3.append(mutual_information_2d(prewhite[s][j][:-1],prewhite[s][i][1:]))
            # pernontas tin 1i mexri to prwteleytaio kai ti deyteri apo 1o mexri teleytaio
        temp4.append(temp3)
    mut_inf_lag1.append(temp4)

############################################################################
################  NON PARAMETRIC TEST FOR MUTUAL INFORMATION ###############
    

################### ME YSTERISI MIDEN : LAG=0

MIarray0=[]

for s in range(0,w6): # ana parathyro xronoseiras
    MIarray0_temp1=[]
    for j in range(20):
        MIarray0_temp2=[]
        for i in range(20): 
            m=[]
            m=test_mutual_information_2d(prewhite[s][j],prewhite[s][i])   #######  me prewhite  #######
    ########  Elegxos pou ftiaxnei syndeseis 0  &  1   ########
            if i==j: MIarray0_temp2.append(0) #mideniki sindesi sti diagwnio
            else:
                if m[1]<0.01:
                    MIarray0_temp2.append(1)
                else :
                    MIarray0_temp2.append(0)
        MIarray0_temp1.append(MIarray0_temp2)
    MIarray0.append(MIarray0_temp1)

################### ME YSTERISI MIDEN : LAG=1

MIarray1=[]

for s in range(0,w6): # ana parathyro xronoseiras
    MIarray1_temp1=[]
    for j in range(20):
        MIarray1_temp2=[]
        for i in range(20): 
            m=[]
            m=test_mutual_information_2d(prewhite[s][j][:-1],prewhite[s][i][1:])   #######  me prewhite  #######
    ########  Elegxos pou ftiaxnei syndeseis 0  &  1   ########
            if m[1]<0.01 :
                MIarray1_temp2.append(1)
            else :
                MIarray1_temp2.append(0)
        MIarray1_temp1.append(MIarray1_temp2)
    MIarray1.append(MIarray1_temp1)

########################################################################

##########################################################################
############# PLOT SYGKENTRWTIKWN PARATHYRWN MUTUAL INFO #################

####################       LAG = 0     #######################

############## PLOT OLWN TWN PARATHYRWN MUTUAL INFO 

for i in range(0,w6):
    matrix = MIarray0[i]  # symmetrikos pinakas 0-1

    fig=plt.figure()

    ax=fig.add_subplot(1,1,1)
    ax.set_aspect('equal')

    plt.imshow(matrix,interpolation='nearest',cmap=plt.cm.Reds)
    plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
    plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
    fig.canvas.set_window_title(str(i+1)+'o Parathyro MUTUAL INFORMATION me LAG=0')
    plt.title('Starts:'+str(start[i])+'-  Stops:'+str(stop[i]))
    plt.colorbar()
    plt.show()

############## SYGKENTRWTIKO MUTUAL INFORMATION ME LAG=0

x=[[0]*20]*20
for i in range (0,w6):
    x+=np.asarray(MIarray0[i])

meg=np.amax(x)
norm_matrix=x.tolist()/meg

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(norm_matrix,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro MUTUAL INFORMATION me LAG = 0')
plt.title('Sygkentrwtikes normalized syndeseis i--j MUTUAL INFORMATION xwris YSTERISI sta {} parathyra'.format(w6))
plt.colorbar()
plt.show()

# Krataw mono osous parousiasan panw apo ta misa parathyra syndeseis
y=(norm_matrix>0.5).astype(int)

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(y,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro MUTUAL INFO me LAG = 0')
plt.title('Osoi emfanisan syndesi i--j se perissotera apo ta misa parathyra')
plt.colorbar()
plt.show()

####################       LAG = 1     #######################

############## PLOT OLWN TWN PARATHYRWN MUTUAL INFO ME LAG = 1

for i in range(0,w6):
    matrix = MIarray1[i]  # symmetrikos pinakas 0-1

    fig=plt.figure()

    ax=fig.add_subplot(1,1,1)
    ax.set_aspect('equal')

    plt.imshow(matrix,interpolation='nearest',cmap=plt.cm.Reds)
    plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
    plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
    fig.canvas.set_window_title(str(i+1)+'o Parathyro MUTUAL INFORMATION me LAG = 1')
    plt.title('Starts:'+str(start[i])+'-  Stops:'+str(stop[i]))
    plt.colorbar()
    plt.show()

############## SYGKENTRWTIKO MUTUAL INFORMATION ME LAG = 1

x=[[0]*20]*20
for i in range (0,w6):
    x+=np.asarray(MIarray1[i])

meg=np.amax(x)
norm_matrix=x.tolist()/meg

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(norm_matrix,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro MUTUAL INFORMATION me LAG = 1')
plt.title('Sygkentrwtikes normalized syndeseis i-->j MUTUAL INFORMATION me YSTERISI sta {} parathyra'.format(w6))
plt.colorbar()
plt.show()

# Krataw mono osous parousiasan panw apo ta misa parathyra syndeseis
y=(norm_matrix>0.5).astype(int)

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow(y,interpolation='nearest',cmap=plt.cm.Reds)
plt.xticks(range(20),[str(stocks[x]) for x in range(20)],rotation='vertical')
plt.yticks(range(20),[str(stocks[x]) for x in range(20)],rotation='horizontal')
fig.canvas.set_window_title('Sygkentrwtiko parathyro MUTUAL INFO me LAG = 1')
plt.title('Osoi emfanisan syndesi i-->j se perissotera apo ta misa parathyra')
plt.colorbar()
plt.show()


