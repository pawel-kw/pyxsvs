import pyxsvs
import sys
import fabio
import pylab
import numpy
import pyFAI
from numpy import pi,sin,arctan,sqrt,mgrid,where

defaultMask = fabio.open('/home/kwasniew/Experiments/MAXIPIX/maxipix1_mask_2013.edf').data

inputFile = sys.argv[1]
calculator = pyxsvs.pyxsvs(inputFile)
exposure = calculator.Parameters['exposureList'][0]
expParams = calculator.Parameters['exposureParams'][exposure]
maskFile = calculator.Parameters['defaultMaskFile']
mask = fabio.open(maskFile).data
flatFieldFile = calculator.Parameters['flatFieldFile']
flatField = fabio.open(flatFieldFile).data
saveDir = calculator.Parameters['saveDir']
outPrefix = calculator.Parameters['outPrefix']
dataDir = calculator.Parameters['dataDir']
dataSuf = expParams['dataSuf']
dataPref = expParams['dataPref']
n1 = expParams['n1']
n2 = expParams['n2']
q1 = calculator.Parameters['q1']
q2 = calculator.Parameters['q2']
qs = calculator.Parameters['qs']
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
fastStatic,histBins = calculator.createFastStatic(fileNames,qRings)
qImg = numpy.ma.masked_array(numpy.ones((dim1,dim2)),qImg)

# Save the averaged static as edf file
edfimg = fabio.edfimage.edfimage()
edfimg.data = fastStatic
saveFileName = saveDir+outPrefix+'2Dstatic.edf'
edfimg.write(saveFileName)
print 'Static file saved to %s' % saveFileName

# Plot the 2D averaged static + q partitions
#fastStatic = numpy.ma.masked_array(fastStatic,mask=defaultMask)
fastStaticMasked = numpy.ma.masked_array(fastStatic,mask=mask)
fig = pylab.figure()
ax = fig.add_subplot(111)
ax.set_title('Static + q partitions')
plt = ax.pcolormesh(numpy.log10(fastStaticMasked))
ax.pcolormesh(qImg,alpha=0.1,cmap=pylab.cm.gray)
#ax.set_xlim(numpy.min(X),numpy.max(X))
#ax.set_ylim(numpy.min(Y),numpy.max(Y))
ax.set_xlim(0,dim1)
ax.set_ylim(0,dim2)
ax.set_aspect(1)
pylab.colorbar(plt)

# Azimuthally regoup the static image and plot together with q partitions
qv = numpy.arange(q1,q2+qs,qs)

poni1 = ceny*pixSize*1e-3
poni2 = cenx*pixSize*1e-3
integrator = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(poni1=poni1,poni2=poni2,dist=sdDist*1e-3,wavelength=wavelength*1e-10,pixel1=pixSize*1e-3,pixel2=pixSize*1e-3)
staticAzim2d,q,chi = integrator.integrate2d(data=fastStatic/flatField,\
        nbPt_rad=2000,nbPt_azim=360,unit='q_nm^-1',mask=mask)
q /= 10 # to 1/A
fig_azim = pylab.figure()
ax_azim = fig_azim.add_subplot(111)
ax_azim.set_title('Azimuthally regrouped static')
pltA = ax_azim.pcolormesh(q,chi,numpy.log10(staticAzim2d+1))
pylab.colorbar(pltA)
ax_azim.set_xlim(numpy.min(q),numpy.max(q))
ax_azim.set_ylim(numpy.min(chi),numpy.max(chi))

# 1D integration
q,staticAzim1d = integrator.integrate1d(data=fastStatic/flatField,\
        nbPt=2000,unit='q_nm^-1',mask=mask)
q /= 10 # to 1/A
fig1d = pylab.figure()
roi = where(staticAzim1d>0)
ax1d = fig1d.add_subplot(111)
ax1d.set_title('Averaged SAXS + q partitions')
ax1d.semilogy(q[roi],staticAzim1d[roi])
ax1d.set_xlim(numpy.min(q[roi]),numpy.max(q[roi]))

for i in xrange(len(qv)):
    ymin, ymax = pylab.ylim()
    pylab.fill([qv[i]-dq,qv[i]-dq,qv[i]+dq,qv[i]+dq],[ymin,ymax,ymax,ymin],'b',alpha=.2,edgecolor='r')
    pylab.text(qv[i],0.5*ymax,str(i+1),fontsize=7,horizontalalignment='center')


pylab.show()
