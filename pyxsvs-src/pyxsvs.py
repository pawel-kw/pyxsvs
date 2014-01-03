#!/usr/bin/python2
# -*- coding: utf-8 -*-
# pyxsvs.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

r'''
.. module:: pyxsvs
    :platform: Unix, Windows

.. moduleauthor:: Pawel Kwasniewski <pawel.kw@gmail.com>

Overview
========
A class for X-ray Speckle Visibility Spectroscopy data analysis.

Classes and functions defined here
----------------------------------

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

PYXSVS Module Documentation
===========================
'''

import pylab
from sys import argv,stdout
import fabio
from ConfigParser import RawConfigParser
from os.path import isfile
import os
import numpy
from numpy import pi,sin,arctan,sqrt,mgrid,where,shape,exp,linspace,std
import sys
from time import time
#from scipy.stats import nbinom
from scipy.optimize import leastsq
from scipy.special import gamma, gammaln
import string
import pickle
from matplotlib.path import Path
import argparse # parsing command line arguments
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
    'axes.formatter.limits': (-4,4),
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
        self.flatField = numpy.array([])
        self.mask = numpy.array([])
        self.static = numpy.array([])
        self.initParameters() # initialize Parameters
        self.qVevtor = numpy.array([])
        self.globalPDFArray = []
        self.qBins = []
        self.histStdDev = []
        self.trace = numpy.array([])
        
        # Read input file
        if type(inputFileName) == type(''):
            try:
                config.read(inputFileName)
            except:
                print 'Problems reading %s' % inputFileName
                sys.exit()
        self.config = config
        self.parseInput(self.config) # Parse the input file
        self.setParameters(**kwargs) # Overwrite parameters (if any are given)
        self.flatField = fabio.open(self.Parameters['flatFieldFile']).data # Load flat field
        self.constructQVector()
        # Create static
        currExpParams = self.Parameters['exposureParams'][self.Parameters['exposureList'][0]]
        dataDir = self.Parameters['dataDir']
        dataSuf = currExpParams['dataSuf']
        n1 = currExpParams['n1']
        n2 = currExpParams['n2']
        dataPref = currExpParams['dataPref']
        fileNames = filename(dataDir+dataPref,dataSuf,n1,n2) # Generate file list
        if len(fileNames) > 200:
            self.static = self.createFastStatic(fileNames[:200])
        else:
            self.static = self.createFastStatic(fileNames)
        try:
            self.mask = fabio.open(self.Parameters['maskFile']).data # Load mask
        except:
            print 'Mask not set, please create one before trying to calculate anything!'

    def initParameters(self):
        r'''Acts on the :sefl.Parameters: dictionary. Sets initial parameter values.
        '''
        self.Parameters = {
                'saveDir' : '',
                'dataDir' : '',
                'flatFieldFile' : '',
                'useFlatField' : False,
                'maskFile' : '',
                'defaultMaskFile' : '',
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
        try:
            self.Parameters['maskFile'] = config.get('Main','mask')
        except:
            print 'Mask not set in the config file'
        self.Parameters['defaultMaskFile'] = config.get('Main','default mask')
        self.Parameters['q1'] = config.getfloat('Main','q1')
        self.Parameters['q2'] = config.getfloat('Main','q2')
        self.Parameters['qs'] = config.getfloat('Main','qs')
        self.Parameters['dq'] = config.getfloat('Main','dq')
        self.Parameters['outPrefix'] = config.get('Main','sample name')
        self.Parameters['wavelength'] = config.getfloat('Main','wavelength')
        self.Parameters['cenx'] = config.getfloat('Main','cenx')
        self.Parameters['ceny'] = config.getfloat('Main','ceny')
        self.Parameters['pixSize'] = config.getfloat('Main','pix')
        self.Parameters['sdDist'] = config.getfloat('Main','sddist')
        secList = config.sections()
        secList.remove('Main')
        exposureList = numpy.sort(secList)
        self.Parameters['exposureList'] = exposureList
        for i in xrange(len(exposureList)):
            exposure = exposureList[i]
            currExpParams = {}
            currExpParams['dataSuf'] = config.get(exposure,'data suffix')
            currExpParams['dataPref'] = config.get(exposure,'data prefix')
            currExpParams['n1'] = config.getint(exposure,'first data file')
            currExpParams['n2'] = config.getint(exposure,'last data file')
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
        self.qVector = numpy.arange(q1,q2+qs,qs)
        self.qVecLen = len(self.qVector)

    def histogramData(self,fileList,qRings,flatField,bins=numpy.arange(10)):
        '''Data reading and processing function. Here's where everything happens.
        By default the data files are histogrammed only after applying the mask.
        When the *useFlatField* parameter is set to *True*, each data file is divided 
        by the flat field and rounded to integers before histogramming.

        *Accepted input:*

        *fileList*: list
            List of data files (full path) to process

        *qRings*: list of pixel indices
            A list containing n lists of pixel indices, corresponding to the chosen q rings

        *flatField*: array
            Flat field to be used to correct each data file

        *bins*: list
            List of bins for the histogram

        *Returns:*
        
        *(globalPDFArray,qBins,histStddev,trace)*: tupile
            Results of histogramming:
            
            *globalPDFArray*: list of histograms for each q ring

            *qBins*: list of histogram bins for each q ring

            *histStddev*: list of standard deviation from the mean (error bars for the histograms)
                for each q ring

            *trace*: list of averaged intensities for each q ring
        '''
        startTime = time()
        # Iterate over files
        n = len(fileList)
        nq = len(qRings)
        trace = numpy.zeros((nq,n))
        qHist = list(numpy.zeros(nq))
        qBins = list(numpy.zeros(nq))
        for i in xrange(n):
            swrite('\r'+str(int(i*100./n))+'%')
            sflush()
            dataFile = fabio.open(fileList[i])
            try:
                rawData = dataFile.data
            except ValueError:
                print 'Problems reading file %s' % fileList[i]
                sys.exit()
            if i == 0:
                # Initiate global histogram
                globalPDFArray = list(numpy.zeros(nq))
                sqrGlobalPDFArray = list(numpy.zeros(nq))
            if self.Parameters['useFlatField']:
                rawData /= flatField 
                rawData = numpy.around(rawData,0) 
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
    
    def createFastStatic(self,fileList,qRings=0):
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
        if qRings != 0:
            nq = len(qRings)
            histBins = list(numpy.zeros(nq))
            for j in xrange(nq):
                data = staticFile[qRings[j]]
                mCnt = numpy.mean(data)
                stddevCnt = numpy.std(data)
                estim = int(mCnt+stddevCnt)
                if estim > 10:
                    histBins[j] = numpy.arange(int(mCnt+stddevCnt))
                else:
                    histBins[j] = numpy.arange(10)
            return staticFile,histBins
        else:
            return staticFile
            
    def createMask(self):
        autoMask = fabio.open(self.Parameters['defaultMaskFile']).data
        saveDir = self.Parameters['saveDir']
        self.mask = maskMaker(self.static,autoMask,saveDir).mask
        pylab.show()

    def calculateVisibility(self):
        r'''Function calculating visibility for each of the exposures
        listed in the input file.
        '''
        exposureList = self.Parameters['exposureList']
        for i in xrange(len(exposureList)):
            exposure = exposureList[i]
            currExpResult = {}
            # Unpack useful variables
            dataDir = self.Parameters['dataDir']
            saveDir = self.Parameters['saveDir']
            outPrefix = self.Parameters['outPrefix']
            wavelength = self.Parameters['wavelength']
            cenx = self.Parameters['cenx']
            ceny = self.Parameters['ceny']
            pixSize = self.Parameters['pixSize']
            sdDist = self.Parameters['sdDist']
            dq = self.Parameters['dq']
            currExpParams = self.Parameters['exposureParams'][exposure]
            dataSuf = currExpParams['dataSuf']
            dataPref = currExpParams['dataPref']
            n1 = currExpParams['n1']
            n2 = currExpParams['n2']
            expTime = currExpParams['expTime']
            expTimeLabel = 't_exp = %.1e s' % expTime
            fileNames = filename(dataDir+dataPref,dataSuf,n1,n2) # Generate file list
            dim1,dim2 = shape(fabio.open(fileNames[0]).data) # Get the data dimensions
            qArray = numpy.ones((dim2,dim1))
            wf = 4*pi/wavelength
            [X,Y] = mgrid[1-ceny:dim2+1-ceny,1-cenx:dim1+1-cenx]
            qArray = wf*sin(arctan(sqrt(X**2+Y**2)*pixSize/sdDist)/2)
            qArray *= (self.mask+1)%2 # Apply mask
            qRings = range(self.qVecLen) # Initiate q partition list
            qImg = numpy.ones((dim1,dim2)) # Image to show q partitions
            # Populate q partition list
            for j in xrange(self.qVecLen):
                qRings[j] = where((qArray >= self.qVector[j] - dq)&(qArray <= self.qVector[j] + dq))
                qImg[qRings[j]] = 0
            # Get static and bins
            if len(fileNames) > 200:
                fastStatic,histBins = self.createFastStatic(fileNames[:200],qRings)
            else:
                fastStatic,histBins = self.createFastStatic(fileNames,qRings)
            qImg = numpy.ma.masked_array(numpy.ones((dim1,dim2)),qImg) # Masked numpy.array for plotting
            figQ = pylab.figure()
            ax1 = figQ.add_subplot(111)
            staticImg = numpy.ma.masked_array(fastStatic,mask=self.mask)
            ax1.pcolormesh(X,Y,staticImg)
            ax1.set_aspect(1)
            ax1.pcolormesh(X,Y,qImg,alpha=0.5,cmap=pylab.cm.gray)
            ax1.set_xlim(numpy.min(X),numpy.max(X))
            ax1.set_ylim(numpy.min(Y),numpy.max(Y))
            pylab.savefig(saveDir+outPrefix+dataPref+exposure+'_q_mask.png',dpi=200)
            ###################
            # Start analysis! #
            ###################
            print 'Histograming %s' % exposure
            xsvsRes,hbins,histStddev,trace = self.histogramData(fileNames,qRings,self.flatField,bins=histBins)
            print 'Done for %s, now plotting and fitting...' % exposure
            # Plot the trace for each q
            figTrace = pylab.figure(figsize=(8,7)) # Trace plots
            pylab.subplots_adjust(bottom=0.1,top=0.85,wspace=0.4)
            figTrace.suptitle(expTimeLabel,fontsize=12)
            # Plot the histogram for each q
            figHist = pylab.figure(figsize=(8,11)) # Histogram plots
            pylab.subplots_adjust(bottom=0.1,top=0.9)
            figHist.suptitle(expTimeLabel,fontsize=14)
            nCols = 2
            rest = self.qVecLen % nCols
            if rest > 0:
                nRows = self.qVecLen/nCols + 1
            else:
                nRows = self.qVecLen/nCols
            beta = numpy.zeros((self.qVecLen,2))
            MArray = numpy.zeros((self.qVecLen,2)) # Number of modes
            KArray = numpy.zeros((self.qVecLen,2)) # Mean number of counts
            RArray = numpy.zeros((self.qVecLen,2)) # Low-count estimate of 1/M
            for j in xrange(self.qVecLen):
                qKey = 'q%03i' % j
                ax2 = figHist.add_subplot(nRows,nCols,j+1)
                axTrace = figTrace.add_subplot(nRows,nCols,j+1)
                axTrace.plot(trace[j,:])
                initKValue = numpy.mean(trace[j,:]) # inital guess for K
                axTrace.set_title('%.2e A^-1' % self.qVector[j],fontsize=10)
                # Fit a negative binomial distribution pmf to the histograms for each q
                xdata = hbins[j][:-1]
                ydata = xsvsRes[j]
                yerr = histStddev[j]
                numpy.logErr = yerr/(ydata*numpy.log(10))
                x = linspace(xdata[0],xdata[-1])
                initParams = [5.0] # [M] - try fit with fixed K
                initEval = peval(x,[5.0,initKValue]) # Evaluate model for initial parameters
                #initParams = [5.0,initKValue] # [M,K]
                # Set a roi for the fitting range
                roi = where(ydata>1e-5)
                if len(roi[0]) > len(ydata)-2:
                    roi = (numpy.array(roi[0][:-2]),)
                #######
                # Fit #
                #######
                #plsq = leastsq(residuals,initParams,\
                #    args=(numpy.log10(ydata[roi]),xdata[roi],numpy.logErr[roi]),\
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
                    args=(numpy.log10(ydata[roi]),xdata[roi],numpy.logErr[roi],initKValue),\
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
                qLabel = '%.2e A^-1' % self.qVector[j]
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
                ax2.set_xlim(xdata[0],xdata[-1])
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
            pylab.figure(figHist.number)
            pylab.savefig(saveDir+outPrefix+dataPref+exposure+'_hist_fits.png',dpi=200)
            pylab.figure(figTrace.number)
            pylab.savefig(saveDir+outPrefix+dataPref+exposure+'_trace.png',dpi=200)
            # Plot fit results
            figRes = pylab.figure(figsize=(6,4)) # Figure for fit results
            figRes.suptitle(expTimeLabel,fontsize=12)
            pylab.subplots_adjust(bottom=0.2,right=0.9,top=0.9)
            axB = figRes.add_subplot(111)
            #axB.set_title(outPrefix+' '+dataPref)
            axE = axB.twinx()
            axB.errorbar(self.qVector,beta[:,0],beta[:,1],\
                fmt='-bo',ms=4,label='contrast')
            #axB.plot(self.qVector,RArray[:,0],\
            #    '-o',mfc='none',label='R-1')
            axB.set_xlabel('q [A^-1]')
            axB.set_ylabel('contrast')
            axB.legend(loc=2)
            axE.errorbar(self.qVector,MArray[:,0],MArray[:,1],\
                fmt='-r^',ms=4,label='M')
            axE.errorbar(self.qVector,KArray[:,0],KArray[:,1],\
                fmt='-gs',ms=4,label='K')
            axE.set_ylabel('M, K')
            axE.legend(loc=1)
            pylab.savefig(saveDir+outPrefix+dataPref+exposure+'_fit_params.png',dpi=200)
        pickle.dump(self.Results, open(saveDir+outPrefix+dataPref+'results.p', "wb"))

class maskMaker:
    '''Interactive mask drawing tool based entirely on matplotlib (not that
    it's a good thing...)
    '''
    def __init__(self, data, auto_mask, savedir):
        self.savedir = savedir
        self.fig = pylab.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title('Select ROI to mask. Press \'m\' to mask, \'u\' to unmask or \'w\' to save and exit ')
        self.canvas = self.ax.figure.canvas
        #self.data = n.log10(data)
        self.data = data
        self.lx, self.ly = shape(self.data)
        self.mask = auto_mask
        self.masked_data = numpy.ma.masked_array(self.data,self.mask)
        self.points = []
        self.key = []
        self.x = 0
        self.y = 0
        self.xy = []
        self.xx = []
        self.yy = []
        self.ind = 0
        self.img = self.ax.imshow(self.masked_data,origin='lower',interpolation='nearest',animated=True)
        self.lc,=self.ax.plot((0,0),(0,0),'-+w',color='black',linewidth=1.5,markersize=8,markeredgewidth=1.5)
        self.lm,=self.ax.plot((0,0),(0,0),'-+w',color='black',linewidth=1.5,markersize=8,markeredgewidth=1.5)
        self.ax.set_xlim(0,self.lx)
        self.ax.set_ylim(0,self.ly)
        for i in range(self.lx):
            for j in range(self.ly):
                self.points.append([i,j])

        cidb = pylab.connect('button_press_event', self.on_click)
        cidk = pylab.connect('key_press_event',self.on_click)
        cidm = pylab.connect('motion_notify_event',self.on_move)

    def on_click(self,event):
        if not event.inaxes:
            self.xy = []
            return
        self.x, self.y = int(event.xdata), int(event.ydata)
        self.key = event.key
        self.xx.append([self.x])
        self.yy.append([self.y])
        self.xy.append([self.y,self.x])
        self.lc.set_data(self.xx,self.yy)
        if self.key == 'm':
            self.xx[-1] = self.xx[0]
            self.yy[-1] = self.yy[0]
            self.xy[-1] = self.xy[0]
            self.ind = Path(self.xy).contains_points(self.points)
            self.mask = self.mask.reshape(self.lx*self.ly,1)
            self.mask[self.ind] = 1
            self.mask = self.mask.reshape(self.lx,self.ly)
            self.update_img()
            self.reset_poly()
        if self.key == 'u':
            self.xx[-1] = self.xx[0]
            self.yy[-1] = self.yy[0]
            self.xy[-1] = self.xy[0]
            self.ind = Path(self.xy).contains_points(self.points)
            self.mask = self.mask.reshape(self.lx*self.ly,1)
            self.mask[self.ind] = 0
            self.mask = self.mask.reshape(self.lx,self.ly)
            self.update_img()
            self.reset_poly()
        draw()
        if self.key == 'w':
            edfImg = fabio.edfimage.edfimage()
            edfImg.data = self.mask
            edfImg.write(self.savedir+'mask.edf')
            print 'Mask saved to %s' % (self.savedir+'mask.edf')
            self.fig.close()

    def on_move(self,event):
        if not event.inaxes: return
        self.xm, self.ym = int(event.xdata), int(event.ydata)
        if self.x != 0:
            self.lm.set_data((self.x,self.xm),(self.y,self.ym))
            draw()

    def update_img(self):
        self.img.set_data(n.ma.masked_array(self.data,self.mask))
        draw()

    def reset_poly(self):
        self.xx = []
        self.yy = []
        self.xy = []
        self.lc.set_data(self.xx,self.yy)
        self.lm.set_data(self.xx,self.yy)
        draw()
        self.x = 0
        self.y = 0

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
    '''Negative Binomial (Poisson-Gamma) distribution function.
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
    coeff = exp(M*numpy.log(M) + (M-1)*numpy.log(x) - gammaln(M) - M*numpy.log(K))
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
    result = (y - numpy.log10(nbinomPMF(x,K,M)))/yerr
    return result

def residuals2(params,y,x,yerr,K):  
    '''Residuals function used for least squares fitting with 
    *K* parameter fixed.
    '''
    M = params
    pr = M/(K+M)
    result = (y - numpy.log10(nbinomPMF(x,K,M)))/yerr
    return result

def peval(x,params):
    '''Function evaluating the binomial distribution for the 
    goven set of input parameters. Redundant - should be removed.
    '''
    M,K = params
    result = nbinomPMF(x,K,M)
    return result

#########################################
# Enable running the module as a script #
#########################################
def main():
    parser = argparse.ArgumentParser(description='Python XSVS data analysis.')
    parser.add_argument('-i',dest='inputFileName', metavar='./input_file.txt', type=str,
                               help='Input file describing the data',required=True)
    args = parser.parse_args()
    calculator = pyxsvs(args.inputFileName)
    calculator.calculateVisibility()

if __name__ == '__main__':
    main()
