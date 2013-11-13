# -*- coding: utf-8 -*-
# pyxsvs.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

r'''
Overview
========
A class for X-ray Speckle Visibility Spectroscopy data analysis.

Classes and functions defined in this file
------------------------------------------

 - :py:class:`pyxsvs() <pyxsvs.pyxsvs>`:
    Main class contiaining all the functions requred to calculate 
    visibility functions and to analyse them.

 - :py:func:`filename() <pyxsvs.filename>`:
    Function returning a list of ID10 style filenames according to the given input 
    parameters (prefix, suffix, extension, first and last file number)

Usage example
-------------

In order to calculate visibility functions the dataset needs to be described by
an input file - a text file with structure similar to MS Windows INI files. 
An example input file is provided. The :Main: section in the file contains 
general settings applicable to the whole XSVS data set. The :Exp_n: sections 
contain information about different exposures: data file suffix and prefix, 
first and last data file number, exposure time.

    >>> import pyxsvs
    >>> calculator = pyxsvs.pyxsvs('./input.txt')
    >>> calculator.calculateVisibility()

Module Documentation
====================
'''

import pylab
from sys import argv,stdout
import fabio
from ConfigParser import RawConfigParser
from os.path import isfile
import os
import numpy
import sys
from time import time
from scipy.stats import nbinom
from scipy.optimize import leastsq
from scipy.special import gamma, gammaln
import string
import pickle
#import lmfit # not used for the moment

# Figure settings
figParams = {
    'legend.fontsize': 10,
    'figure.subplot.bottom': 0.14,
    'figure.subplot.top': 0.95,
    'figure.subplot.left': 0.14,
    'figure.subplot.right': 0.95,
    'figure.subplot.wspace': 0.3,
    'figure.subplot.hspace': 0.3,
    'legend.numpoints': 1,
    'legend.handletextpad': 0.1,
    'font.family': 'serif',
    'axes.formatter.limits': (-2,2),
    }
pylab.rcParams.update(figParams)

swrite = stdout.write
sflush = stdout.flush

class pyxsvs(object):
    r'''Main class containing all functions for speckle visibility
    calculation and analysis.
    '''
    def __init__(self,inputFileName,**kwargs):
        r'''Constructor, initializes the object and reads the input file.
        **Input parameters:**

        :py:attr: inputFileName: string,
                    Path to the input file describing the data set. The input 
                    file is a plain text file with a structure similar to MS Windows
                    INI files.
        :py:attr: \*\*kwargs: key=value pairs,
                    When given, they overwrite the parameters contained in the input 
                    file.

        **Returns:**
        
        :py:class: pyxsvs object
        .. document private Functions
        '''
        # Initialize variables
        config = RawConfigParser()
        self.Parameters = {}
        self.Results = {}
        self.flatField = pylab.array([])
        self.mask = pylab.array([])
        self.initParameters() # initialize Parameters
        self.qVevtor = pylab.array([])
        self.globalPDFArray = []
        self.qBins = []
        self.histStdDev = []
        self.trace = pylab.array([])
        
        # Read input file
        if type(inputFileName) == type(''):
            try:
                config.read(inputFileName)
            except:
                print 'Problems reading %s' % inputFileName
                sys.exit()
        self.parseInput(config) # Parse the input file
        self.setParameters(**kwargs) # Overwrite parameters (if any are given)
        self.flatField = fabio.open(self.Parameters['flatFieldFile']).data # Load flat field
        self.mask = fabio.open(self.Parameters['maskFile']).data # Load mask
        self.constructQVector(self)

    def initParameters(self):
        r'''Acts on the :sefl.Parameters: dictionary. Sets initial parameter values.
        '''
        self.Parameters = {
                'saveDir' : '',
                'dataDir' : '',
                'flatFieldFile' : '',
                'maskFile' : '',
                'q1' : 0.0,
                'q2' : 0.0,
                'qs' : 0.0,
                'dq' : 0.0,
                'outPrefix' : '',
                'wavelength' : 0.0,
                'cenx' : 0.0,
                'ceny' : 0.0,
                'pixSize' : 0.0,
                'sdDist' : 0.0,
                'exposureList' : '',
                'exposureParams' : {},
                }

    def parseInput(self,config):
        r'''Function reading the input file and setting the 
        :self.Parameters: acordingly.
        '''
        self.Parameters['saveDir'] = config.get('Main','save dir')
        self.Parameters['dataDir'] = config.get('Main','data dir')
        self.Parameters['flatFieldFile'] = config.get('Main','flat field')
        self.Parameters['maskFile'] = config.get('Main','mask')
        self.Parameters['q1'] = config.getfloat('Main','q1')
        self.Parameters['q2'] = config.getfloat('Main','q2')
        self.Parameters['qs'] = config.getfloat('Main','qs')
        self.Parameters['dq'] = config.getfloat('Main','dq')
        self.Parameters['outPrefix'] = config.get('Main','sample name')
        self.Parameters['wavelength'] = config.getfloat('Main','wavelength')
        self.Parameters['cenx'] = config.getint('Main','cenx')
        self.Parameters['ceny'] = config.getint('Main','ceny')
        self.Parameters['pixSize'] = config.getfloat('Main','pix')
        self.Parameters['sdDist'] = config.getfloat('Main','sddist')
        exposureList = sort(config.sections().remove('Main'))
        self.Parameters['exposureList'] = exposureList
        for i in xrange(len(exposureList)):
            exposure = exposureList[i]
            currExpParams = {}
            currExpParams['dataSuf'] = config.get(exposure,'data suffix')
            currExpParams['dataPref'] = config.get(exposure,'data prefix')
            currExpParams['n1'] = config.get(exposure,'first data file')
            currExpParams['n2'] = config.get(exposure,'last data file')
            currExpParams['expTime'] = config.getfloat(exposure,'exp time')
            self.Parameters['exposureParams'][exposure] = currExpParams
            self.Results[exposure] = {} # Initialize Results container
            self.Results[exposure]['expTime'] = currExpParams['expTime']

    def setParameters(self,**kwargs):
        r'''Sets the parameters given in keyword - value pairs to known settings 
        keywords. Unknown key - value pairs are skipped.
        '''
        for kw in kwargs:
            if kw in self.Parameters.keys():
                self.Parameters[kw] = kwargs[kw]
            else:
                pass # Ignoring unkwonw settings
        
    def constructQVector(self):
        r'''Generates a q vector acording to settings in :self.Parameters:
        '''
        q1 = self.Parameters['q1']
        q2 = self.Parameters['q2']
        qs = self.Parameters['qs']
        dq = self.Parameters['dq']
        self.qVector = pylab.arange(q1,q2+qs,qs)
        self.qVecLen = len(self.qVector)

    def histogramData(self,fileList,qRings,flatField,bins=arange(10)):
        '''Data reading and processing function. Here's where everything happens.
        '''
        startTime = time()
        # Iterate over files
        n = len(fileList)
        nq = len(qRings)
        trace = zeros((nq,n))
        qHist = list(zeros(nq))
        qBins = list(zeros(nq))
        for i in xrange(n):
            swrite('\r'+str(int(i*100./n))+'%')
            sflush()
            dataFile = fabio.open(fileList[i])
            try:
                rawData = dataFile.GetData(0)
            except ValueError:
                print 'Problems reading file %s' % fileList[i]
                sys.exit()
            if i == 0:
                # Initiate global histogram
                globalPDFArray = list(zeros(nq))
                sqrGlobalPDFArray = list(zeros(nq))
            #rawData /= flatField 
            #rawData = numpy.around(rawData,0)
            # For each file iterate over q rings
            for j in xrange(nq):
                data = rawData[qRings[j]]
                trace[j,i] = numpy.mean(data)
                qHist[j],qBins[j] = numpy.histogram(data,bins=bins[j],normed=True)
                globalPDFArray[j] += qHist[j]
                sqrGlobalPDFArray[j] += numpy.power(qHist[j],2)
            del rawData
        histStddev = range(nq)
        for j in xrange(nq):
            globalPDFArray[j] /= n
            sqrGlobalPDFArray[j] /= n
            histStddev[j] = sqrt(sqrGlobalPDFArray[j] - \
                numpy.power(globalPDFArray[j],2))
        endTime = time()
        print '\nCalculations took %.2f s' % (endTime-startTime)
        return globalPDFArray,qBins,histStddev,trace
    
    def createFastStatic(self,fileList,qRings):
       
        r'''Function creating an averaged 2D SAXS image from data files
        provided in :fileList: and producing bins for photon counting
        histograms based on the numbers of photons found in the averaged image.
        Bins are created independently for each q partition.

        **Input parameters:**

        :py:attr: fileList: list,
                    List of files to use for averaging.
        :py:attr: qRings: list,
                    Coordinates of pixels belonging to different q partitions.
        **Returns:**
        :py:attr: staticFile: numpy.array,
                    Averaged scattering pattern.
        :py:attr: histBins: list of histogram bins for all q partitions.
        '''

        nq = len(qRings)
        histBins = list(zeros(nq))
        res = numpy.zeros(shape(fabio.open(fileList[0]).data))
        for i in xrange(len(fileList)):
           fileName = fileList[i]
           f = fabio.open(fileName)
           try:
               data = f.data
           except ValueError:
               print 'Problems reading file %s' % fileName
               sys.exit()
           res += data
           del data
        staticFile = res/len(fileList)
        for j in xrange(nq):
            data = staticFile[qRings[j]]
            mCnt = numpy.mean(data)
            stddevCnt = numpy.std(data)
            estim = int(mCnt+stddevCnt)
            if estim > 10:
                histBins[j] = arange(int(mCnt+stddevCnt))
            else:
                histBins[j] = arange(10)
        return staticFile,histBins
            
    def calculateVisibility(self):
        r'''Function calculating visibility for each of the exposures
        listed in the input file.
        '''
        exposureList = self.Parameters('exposureList')
        for i in xrange(len(exposureList)):
            exposure = exposureList[i]
            currExpResult = {}
            currExpParams = self.Parameters['exposureParams'][exposure]
            dataSuf = currExpParams['dataSuf']
            dataPref = currExpParams['dataPref']
            n1 = currExpParams['n1']
            n2 = currExpParams['n2']
            expTime = currExpParams['expTime']
            expTimeLabel = 't_exp = %.1e s' % expTime
            fileNames = filename(dataDir+dataPref,dataSuf,n1,n2) # Generate file list
            dim1,dim2 = shape(fabio.open(fileNames[0]).data) # Get the data dimensions
            qArray = ones((dim2,dim1))
            wf = 4*pi/wavelength
            [X,Y] = mgrid[1-ceny:dim2+1-ceny,1-cenx:dim1+1-cenx]
            qArray = wf*sin(arctan(sqrt(X**2+Y**2)*pixSize/sdDist)/2)
            qArray *= mask # Apply mask
            qRings = range(self.self.qVecLen) # Initiate q partition list
            qImg = numpy.ones((dim1,dim2)) # Image to show q partitions
            # Populate q partition list
            for j in xrange(self.qVecLen):
                qRings[j] = where((qArray >= qVector[j] - dq)&(qArray <= qVector[j] + dq))
                qImg[qRings[j]] = 0
            # Create static and bins
            # For the moment, the fast static is created from first 200 files.
            # This should be changed into something smarter.
            fastStatic,histBins = self.createFastStatic(fileNames[:200],qRings)
            qImg = numpy.ma.masked_array(numpy.ones((dim1,dim2)),qImg) # Masked array for plotting
            figQ = figure()
            ax1 = figQ.add_subplot(111)
            ax1.pcolormesh(X,Y,fastStatic*mask)
            ax1.set_aspect(1)
            ax1.pcolormesh(X,Y,qImg,alpha=0.5,cmap=cm.gray)
            savefig(saveDir+outPrefix+dataPref+exposure+'_q_mask.png',dpi=200)
            ###################
            # Start analysis! #
            ###################
            print 'Histograming %s' % exposure
            xsvsRes,hbins,histStddev,trace = self.histogramData(fileNames,qRings,flatField,bins=histBins)
            print 'Done for %s, now plotting and fitting...' % exposure
            # Plot the trace for each q
            figTrace = figure(figsize=(8,7)) # Trace plots
            subplots_adjust(bottom=0.1,top=0.85,wspace=0.4)
            figTrace.suptitle(expTimeLabel,fontsize=12)
            # Plot the histogram for each q
            figHist = figure(figsize=(8,11)) # Histogram plots
            subplots_adjust(bottom=0.1,top=0.9)
            figHist.suptitle(expTimeLabel,fontsize=14)
            nCols = 2
            rest = self.qVecLen % nCols
            if rest > 0:
                nRows = self.qVecLen/nCols + 1
            else:
                nRows = self.qVecLen/nCols
            beta = zeros((self.qVecLen,2))
            MArray = zeros((self.qVecLen,2)) # Number of modes
            KArray = zeros((self.qVecLen,2)) # Mean number of counts
            RArray = zeros((self.qVecLen,2)) # Low-count estimate of 1/M
            for j in xrange(self.qVecLen):
                qKey = 'q%03i' % j
                ax2 = figHist.add_subplot(nRows,nCols,j+1)
                axTrace = figTrace.add_subplot(nRows,nCols,j+1)
                axTrace.plot(trace[j,:])
                initKValue = mean(trace[j,:]) # inital guess for K
                axTrace.set_title('%.2e A^-1'%qVector[j],fontsize=10)
                # Fit a negative binomial distribution pmf to the histograms for each q
                xdata = hbins[j][:-1]
                ydata = xsvsRes[j]
                yerr = histStddev[j]
                logErr = yerr/(ydata*log(10))
                x = linspace(xdata[0],xdata[-1])
                initParams = [5.0] # [M] - try fit with fixed K
                initEval = peval(x,[5.0,initKValue]) # Evaluate model for initial parameters
                #initParams = [5.0,initKValue] # [M,K]
                # Set a roi for the fitting range
                roi = where(ydata>1e-5)
                if len(roi[0]) > len(ydata)-2:
                    roi = (array(roi[0][:-2]),)
                #######
                # Fit #
                #######
                #plsq = leastsq(residuals,initParams,\
                #    args=(log10(ydata[roi]),xdata[roi],logErr[roi]),\
                #    full_output=1)
                #s_sq = (plsq[2]['fvec']**2).sum()/\
                #   (len(ydata[roi])-len(initParams))
                #paramStdDev = sqrt(plsq[1].diagonal()*s_sq)
                #M,K = plsq[0]
                #MArray[j,0],KArray[j,0] = plsq[0]
                #MArray[j,1],KArray[j,1] = paramStdDev
                #resText =  'M = %.2f\nK = %.2f' % (MArray[j,0],KArray[j,0])
                #============================================================
                ################################
                # Alternative fit with fixed K #
                ################################
                plsq = leastsq(residuals2,initParams,\
                    args=(log10(ydata[roi]),xdata[roi],logErr[roi],initKValue),\
                    full_output=1)
                s_sq = (plsq[2]['fvec']**2).sum()/\
                   (len(ydata[roi])-len(initParams))
                paramStdDev = sqrt(plsq[1].diagonal()*s_sq)
                M = plsq[0]
                K = initKValue
                MArray[j,0] = M
                MArray[j,1] = paramStdDev
                KArray[j,0] = K
                KArray[j,1] = std(trace[j,:])
                resText =  'M = %.2f\nK = %.2f' % (MArray[j,0],KArray[j,0])
                #============================================================
                axTrace.axhline(KArray[j,0],ls='-',color='r',lw=2)
                ax2.text(0.75,0.7,resText,transform=ax2.transAxes,fontsize=10)
                fitEval = peval(x,[M,K])
                pEval = poisson(x,K)
                qLabel = '%.2e A^-1' % qVector[j]
                ax2.set_title(qLabel)
                #GammaEval = gammaDist(x,[M,K])
                ax2.errorbar(xdata[roi],ydata[roi],yerr[roi],\
                fmt='o',ms=4,label=qLabel)
                ax2.errorbar(xdata,ydata,yerr,\
                fmt='o',mfc='none',ms=4)
                #ax2.plot(x,initEval,'--g',label='init') # plot initial evaluation
                ax2.plot(x,fitEval,'-r',label='fit') # plot fit result
                ax2.plot(x,pEval,'-b',label='Poisson') # plot Poisson distribution
                #ax2.plot(x,GammaEval,'--m',label='Gamma') # plot Gamma distribution
                ax2.set_yscale('log')
                #ax2.set_ylim(min(where(ydata>0)[0]),max(ydata))
                ax2.set_ylim(1e-10,1)
                beta[j,0] = 1./MArray[j,0]
                beta[j,1] = MArray[j,1]/MArray[j,0]**2
                #legend(loc=3)
                RArray[j,0] = 2.*ydata[2]*(1.-ydata[1])/ydata[1]**2 - 1
                ################
                # Save results #
                ################
                currResults = {}
                currResults[qKey] = {}
                currResults[qKey]['histogramBins'] = xdata
                currResults[qKey]['histogram'] = ydata
                currResults[qKey]['histStdDev'] = yerr
                currResults[qKey]['trace'] = trace[j,:]
                currResults[qKey]['M'] = \
                        {'value': MArray[j,0], 'stddev': MArray[j,1]}
                currResults[qKey]['K'] = \
                        {'value': KArray[j,0], 'stddev': KArray[j,1]}
                currResults[qKey]['beta'] = \
                        {'value': beta[j,0], 'stddev': beta[j,1]}
                currResults[qKey]['R-1'] = \
                        {'value': RArray[j,0], 'stddev': RArray[j,1]}
                self.Results[exposure]['data'] = currResults
            figure(figHist.number)
            savefig(saveDir+outPrefix+dataPref+exposure+'_hist_fits.png',dpi=200)
            figure(figTrace.number)
            savefig(saveDir+outPrefix+dataPref+exposure+'_trace.png',dpi=200)
            # Plot fit results
            figRes = figure(figsize=(6,4)) # Figure for fit results
            figRes.suptitle(expTimeLabel,fontsize=12)
            subplots_adjust(bottom=0.2,right=0.9,top=0.9)
            axB = figRes.add_subplot(111)
            #axB.set_title(outPrefix+' '+dataPref)
            axE = axB.twinx()
            axB.errorbar(qVector,beta[:,0],beta[:,1],\
                fmt='-bo',ms=4,label='contrast')
            axB.plot(qVector,RArray[:,0],\
                '-o',mfc='none',label='R-1')
            axB.set_xlabel('q [A^-1]')
            axB.set_ylabel('contrast')
            axB.legend(loc=2)
            axE.errorbar(qVector,MArray[:,0],MArray[:,1],\
                fmt='-r^',ms=4,label='M')
            axE.errorbar(qVector,KArray[:,0],KArray[:,1],\
                fmt='-gs',ms=4,label='K')
            axE.set_ylabel('M, K')
            axE.legend(loc=1)
            savefig(saveDir+outPrefix+dataPref+exposure+'_fit_params.png',dpi=200)
            pickle.dump(self.Results, open(saveDir+outPrefix+dataPref+exposure+'_results.p', "wb"))

###############################
# Helper function definitions #
###############################

def filename(pref,suf,firstf,lastf):
    '''Function creating file name list 
    '''
    numb=range(firstf,lastf+1)
    fname=numb
    for i in range(len(numb)):
        if numb[i] < 10000:
            fname[i]=pref+string.zfill(str(numb[i]),4)+suf
        else:
            fname[i]=pref+str(numb[i])+suf
    return fname 

def nbinomPMF(x,K,M):
    '''Binomial (Poisson-Gamma) distribution function.
    '''
    coeff = exp(gammaln(x+M)-gammaln(x+1)-gammaln(M))
    Pk = coeff*numpy.power(M/(K+M),M)
    coeff2 = numpy.power(K/(M+K),x)
    Pk *= coeff2
    return Pk

def gammaDist(x,params):
    '''Gamma distribution function
    '''
    M,K = params
    coeff = exp(M*log(M) + (M-1)*log(x) - gammaln(M) - M*log(K))
    Gd = coeff*exp(-M*x/K)
    return Gd

def poisson(x,K):
    '''Poisson distribution function
    '''
    Pk = exp(-K)*numpy.power(K,x)/gamma(x+1)
    return Pk

def residuals(params,y,x,yerr):  
    '''Residuals function used for least squares fitting
    '''
    M,K = params
    pr = M/(K+M)
    result = (y - log10(nbinomPMF(x,K,M)))/yerr
    return result

def residuals2(params,y,x,yerr,K):  
    '''Residuals function used for least squares fitting with 
    *K* parameter fixed.
    '''
    M = params
    pr = M/(K+M)
    result = (y - log10(nbinomPMF(x,K,M)))/yerr
    return result

def peval(x,params):
    '''Function evaluating the binomial distribution for the 
    goven set of input parameters. Redundant - should be removed.
    '''
    M,K = params
    result = nbinomPMF(x,K,M)
    return result

