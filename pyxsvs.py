from pylab import *
import EdfFile, nfiles
from ConfigParser import RawConfigParser
from os.path import isfile
import os
import numpy
from sys import argv,stdout
import sys
from time import time
from scipy.stats import nbinom
from scipy.optimize import leastsq
from scipy.special import gamma, gammaln
import savedata

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
rcParams.update(figParams)

swrite=stdout.write
sflush=stdout.flush

def nbinomPMF(x,K,M):
    coeff = exp(gammaln(x+M)-gammaln(x+1)-gammaln(M))
    Pk = coeff*numpy.power(M/(K+M),M)
    coeff2 = numpy.power(K/(M+K),x)
    Pk *= coeff2
    return Pk

def gammaDist(x,params):
    M,K = params
    coeff = exp(M*log(M) + (M-1)*log(x) - gammaln(M) - M*log(K))
    Gd = coeff*exp(-M*x/K)
    return Gd

def poisson(x,K):
    Pk = exp(-K)*numpy.power(K,x)/gamma(x+1)
    return Pk

def residuals(params,y,x,yerr):  
    M,K = params
    pr = M/(K+M)
    result = (y - log10(nbinomPMF(x,K,M)))/yerr
    return result

def residuals2(params,y,x,yerr,K):  
    M = params
    pr = M/(K+M)
    result = (y - log10(nbinomPMF(x,K,M)))/yerr
    return result

def peval(x,params):
    M,K = params
    pr = M/(K+M)
    result = nbinomPMF(x,K,M)
    return result

def dataRead(fileList,qRings,flatField,bins=arange(10)):
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
        dataFile = EdfFile.EdfFile(fileList[i])
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

def createFastStatic(fileList,qRings):
    nq = len(qRings)
    histBins = list(zeros(nq))
    res = numpy.zeros(shape(EdfFile.EdfFile(fileList[0]).GetData(0)))
    for i in xrange(len(fileList)):
       fileName = fileList[i]
       f = EdfFile.EdfFile(fileName)
       try:
           data = f.GetData(0)
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
        

# Read input file
inputFile = argv[1]
config = RawConfigParser()
config.read(inputFile)

saveDir = config.get('Main','save dir')
dataDir = config.get('Main','data dir')
flatFieldFile = config.get('Main','flat field')
maskFile = config.get('Main','mask')
q1 = config.getfloat('Main','q1')
q2 = config.getfloat('Main','q2')
qs = config.getfloat('Main','qs')
dq = config.getfloat('Main','dq')
outPrefix = config.get('Main','sample name')
wavelength = config.getfloat('Main','wavelength')
cenx = config.getint('Main','cenx')
ceny = config.getint('Main','ceny')
pixSize = config.getfloat('Main','pix')
sdDist = config.getfloat('Main','sddist')

exposureList = sort(config.sections())
flatField = EdfFile.EdfFile(flatFieldFile).GetData(0)
mask = EdfFile.EdfFile(maskFile).GetData(0)

qVector = arange(q1,q2+qs,qs)
qVecLen = len(qVector)

# Do analysis for each exposure
for i in xrange(len(exposureList)):
    dataToSave = {}
    dataToSave['q'] = qVector
    exposure = exposureList[i]
    if exposure != 'Main':
        dlabel = 'XSVS'
        dataToSave[dlabel] = {} # prepare the data container
        dataSuf = config.get(exposure,'data suffix')
        dataPref = config.get(exposure,'data prefix')
        n1 = config.getint(exposure,'first data file')
        n2 = config.getint(exposure,'last data file')
        fileNames = nfiles.filename(dataDir+dataPref,dataSuf,n1,n2) # Generate file list
        dim1,dim2 = [256,256]
        qArray = ones((dim2,dim1))
        wf = 4*pi/wavelength
        [X,Y] = mgrid[1-ceny:dim2+1-ceny,1-cenx:dim1+1-cenx]
        qArray = wf*sin(arctan(sqrt(X**2+Y**2)*pixSize/sdDist)/2)
        qArray *= mask # Apply mask
        qRings = range(qVecLen)
        qImg = numpy.ones((dim1,dim2))
        for j in xrange(qVecLen):
            qRings[j] = where((qArray >= qVector[j] - dq)&(qArray <= qVector[j] + dq))
            qImg[qRings[j]] = 0
        # Create static and bins
        fastStatic,histBins = createFastStatic(fileNames[:200],qRings) # Create a fast static image for q ring definition
        qImg = numpy.ma.masked_array(numpy.ones((dim1,dim2)),qImg)
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
        #print histBins
        xsvsRes,hbins,histStddev,trace = dataRead(fileNames,qRings,flatField,bins=histBins)
        # Save data to container
        expTime = config.getfloat(exposure,'exp time')
        dataToSave[dlabel]['expTime'] =  expTime # exposure time
        expTimeLabel = 't_exp = %.1e s' % expTime
        # Plot the trace for each q
        figTrace = figure(figsize=(8,7)) # Trace plots
        subplots_adjust(bottom=0.1,top=0.85,wspace=0.4)
        figTrace.suptitle(expTimeLabel,fontsize=12)
        # Plot the histogram for each q
        figHist = figure(figsize=(8,11)) # Histogram plots
        subplots_adjust(bottom=0.1,top=0.9)
        figHist.suptitle(expTimeLabel,fontsize=14)
        nCols = 2
        rest = qVecLen % nCols
        if rest > 0:
            nRows = qVecLen/nCols + 1
        else:
            nRows = qVecLen/nCols
        beta = zeros((qVecLen,2))
        MArray = zeros((qVecLen,2))
        KArray = zeros((qVecLen,2))
        RArray = zeros((qVecLen,2))
        dataToSave[dlabel]['data'] = {}
        for j in xrange(qVecLen):
            qKey = 'q%03i' % j
            dataToSave[dlabel]['data'][qKey] = {}
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
            initEval = peval(x,[M,K])
            fitEval = peval(x,[M,K])
            pEval = poisson(x,K)
            #GammaEval = gammaDist(x,[M,K])
            ax2.errorbar(xdata[roi],ydata[roi],yerr[roi],\
            fmt='o',ms=4,label='%.2e A^-1'%qVector[j])
            ax2.errorbar(xdata,ydata,yerr,\
            fmt='o',mfc='none',ms=4)
            #ax2.plot(x,initEval,'--g',label='init')
            ax2.plot(x,fitEval,'-r',label='fit')
            ax2.plot(x,pEval,'-b',label='Poisson')
            #ax2.plot(x,GammaEval,'--m',label='Gamma')
            ax2.set_yscale('log')
            #ax2.set_ylim(min(where(ydata>0)[0]),max(ydata))
            ax2.set_ylim(1e-10,1)
            beta[j,0] = 1./MArray[j,0]
            beta[j,1] = MArray[j,1]/MArray[j,0]**2
            #legend(loc=3)
            RArray[j,0] = 2.*ydata[2]*(1.-ydata[1])/ydata[1]**2 - 1
            # Save data
            dataToSave[dlabel]['data'][qKey]['histogramBins'] = xdata 
            dataToSave[dlabel]['data'][qKey]['histogram'] = ydata 
            dataToSave[dlabel]['data'][qKey]['histStdDev'] = yerr
            dataToSave[dlabel]['data'][qKey]['trace'] = trace[j,:]
            dataToSave[dlabel]['data'][qKey]['M'] = \
                    {'value': MArray[j,0], 'stddev': MArray[j,1]}
            dataToSave[dlabel]['data'][qKey]['K'] = \
                    {'value': KArray[j,0], 'stddev': KArray[j,1]}
            dataToSave[dlabel]['data'][qKey]['beta'] = \
                    {'value': beta[j,0], 'stddev': beta[j,1]}
            dataToSave[dlabel]['data'][qKey]['R-1'] = \
                    {'value': RArray[j,0], 'stddev': RArray[j,1]}
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
        savedata.savedata(saveDir+outPrefix+dataPref+exposure+'_results.json',dataToSave)        
#show()
