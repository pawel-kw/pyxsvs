import pyxsvs
import sys
import fabio
import pylab
import numpy
from numpy import pi,sin,arctan,sqrt,mgrid,where

defaultMask = fabio.open('/home/kwasniew/Experiments/MAXIPIX/maxipix1_mask_2013.edf').data

inputFile = sys.argv[1]
calculator = pyxsvs.pyxsvs(inputFile)
exposure = calculator.Parameters['exposureList'][0]
expParams = calculator.Parameters['exposureParams'][exposure]
saveDir = calculator.Parameters['saveDir']
outPrefix = calculator.Parameters['outPrefix']
dataDir = calculator.Parameters['dataDir']
dataSuf = expParams['dataSuf']
dataPref = expParams['dataPref']
n1 = expParams['n1']
n2 = expParams['n2']
dq = calculator.Parameters['dq']
wavelength = calculator.Parameters['wavelength']
cenx = calculator.Parameters['cenx']
ceny = calculator.Parameters['ceny']
pixSize = calculator.Parameters['pixSize']
sdDist = calculator.Parameters['sdDist']
fileNames = pyxsvs.filename(dataDir+dataPref,dataSuf,n1,n2)
dim1,dim2 = numpy.shape(fabio.open(fileNames[0]).data) # Get the data dimensions
qArray = numpy.ones((dim2,dim1))
wf = 4*pi/wavelength
[X,Y] = mgrid[1-ceny:dim2+1-ceny,1-cenx:dim1+1-cenx]
qArray = wf*sin(arctan(sqrt(X**2+Y**2)*pixSize/sdDist)/2)
qRings = range(calculator.qVecLen) # Initiate q partition list
qImg = numpy.ones((dim1,dim2)) # Image to show q partitions
# Populate q partition list
for j in xrange(calculator.qVecLen):
    qRings[j] = where((qArray >= calculator.qVector[j] - dq)&(qArray <= calculator.qVector[j] + dq))
    qImg[qRings[j]] = 0
fastStatic,histBins = calculator.createFastStatic(fileNames[:200],qRings)

edfimg = fabio.edfimage.edfimage()
edfimg.data = fastStatic
edfimg.write(saveDir+outPrefix+'_2Dstatic.edf')

fastStatic = numpy.ma.masked_array(fastStatic,mask=defaultMask)
pylab.figure()
pylab.imshow(fastStatic,origin='lower',interpolation='nearest')
pylab.show()
