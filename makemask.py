#! /usr/bin/env python
"""
Draw a mask on the existing SAXS file to define ROI for calculations.

Usage:
makemask.py </path/to/input.txt>

Requirements:
-mask_medipix.edf must be present in the directory where the script is run.
-SAXS file (output of sumedft.py) must be present.
"""

from pylab import *
from inputparser import *
from utils_xpcs import *
from matplotlib.path import points_in_path as points_inside_poly
import sys
import numpy as n


class Mask:
    def __init__(self, data, auto_mask):
        self.fig = figure()
        title('Select ROI to mask. Press m to mask or w to save and exit ')
        self.ax = self.fig.add_subplot(111)
        self.canvas = self.ax.figure.canvas
        self.data = data
        self.lx, self.ly = shape(self.data)
        self.mask = auto_mask
        self.masked_data = n.ma.masked_array(self.data,self.mask)
        self.points = []
        self.key = []
        self.x = 0
        self.y = 0
        self.xy = []
        self.xx = []
        self.yy = []
        self.ind = 0
        self.img = self.ax.imshow(self.masked_data,vmax=n.max(data)-3,origin='lower',interpolation='nearest',animated=True)
        self.lc,=self.ax.plot((0,0),(0,0),'-+w',color='black',linewidth=1.5,markersize=8,markeredgewidth=1.5)
        self.lm,=self.ax.plot((0,0),(0,0),'-+w',color='black',linewidth=1.5,markersize=8,markeredgewidth=1.5)
        self.ax.set_xlim(0,self.lx)
        self.ax.set_ylim(0,self.ly)
        for i in range(self.lx):
            for j in range(self.ly):
                self.points.append([i,j])

        cidb = connect('button_press_event', self.on_click)
        cidk = connect('key_press_event',self.on_click)
        cidm = connect('motion_notify_event',self.on_move)

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
            self.ind = nonzero(points_inside_poly(self.points,self.xy))
            self.mask = self.mask.reshape(self.lx*self.ly,1)
            self.mask[self.ind] = 1
            self.mask = self.mask.reshape(self.lx,self.ly)
            self.update_img()
            self.reset_poly()
        draw()
        if self.key == 'w':
            # Inverting the mask
            self.mask = (self.mask+1)%2
            saveedf(savdir+'mask_'+sname+'.edf',self.mask)
            close()

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



# Read input parameters
input_file = InputFile()
try: datdir,prefd,sufd,nf1,nf2,darkdir,df1,df2,sname,cx,cy,dt,lambdaw,pix_size,distance,dqv,qv1,qv2,qvs,dis_lev,ttcf_par,savdir=input_file.old_read(sys.argv[1])
except:
    print "Problem reading %s" % sys.argv[1]
    sys.exit()

try:
    saxs_file=getedffile(savdir+'saxs_'+sname+'.edf')
    header=saxs_file.GetStaticHeader(0)
except ValueError:
    print 'There is no SAXS file or something is wrong with it!'
    sys.exit()


auto_maskf=getedffile('/home/kwasniew/Tools/Python/pyxpcs_branch/mask_medipix.edf')
#auto_maskf=getedffile('/buffer/tinc1/PAWEL/pyxpcs_multi_series/mask_medipix.edf')
auto_mask=auto_maskf.GetData(0)
auto_mask[:,143] = 1
init_mask = zeros(shape(auto_mask))
if len(sys.argv)>2:
    init_mask = getedffile(sys.argv[2])
    init_mask = init_mask.GetData(0)
    auto_mask[where(init_mask==0)]=1
    #auto_mask = init_mask
data=log10(saxs_file.GetData(0)+0.0001)

mask = Mask(data,auto_mask)

show()
