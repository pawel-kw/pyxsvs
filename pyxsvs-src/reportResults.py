# reportResults.py

r'''
Overview
========

This is a script which produces a single PDF document, composed of the plots
generated as output from pyxsvs.py script.
'''

import pylab
import pyxsvs 
import pickle
import argparse

#from PIL import Image
from reportlab.platypus.flowables import Image, ParagraphAndImage
from glob import glob

from pdfrw import PdfReader
from pdfrw.buildxobj import pagexobj
from pdfrw.toreportlab import makerl

from reportlab.platypus import Flowable
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.lib import utils

PAGE_HEIGHT=defaultPageSize[1]; PAGE_WIDTH=defaultPageSize[0]
styles = getSampleStyleSheet()

def ReportFirstPage(canvas, doc):
    '''First page of the PDF report
    '''
    Title = 'XSVS Report'
    pageinfo = 'XSVS data analysis'
    canvas.saveState()
    canvas.setFont('Times-Bold',16)
    canvas.drawCentredString(PAGE_WIDTH/2.0, PAGE_HEIGHT-98, Title)
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.75 * inch, "Page %d / %s" % (doc.page, pageinfo))
    canvas.restoreState()

def ReportLaterPages(canvas, doc):
    '''Later pages format
    '''
    pageinfo = 'XSVS data analysis'
    canvas.saveState()
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.75 * inch, "Page %d / %s" % (doc.page, pageinfo))
    canvas.restoreState()

def get_image(path, width=1*inch, caption='Image caption',style=styles['Normal']):
    img = utils.ImageReader(path)
    iw, ih = img.getSize()
    aspect = ih / float(iw)
    resImg = Image(path, width=width, height=(width * aspect))
    caption = Paragraph(caption,style)
    return ParagraphAndImage(caption, resImg, side='bottom')

def main():
    parser = argparse.ArgumentParser(description='Python XSVS data analysis.')
    parser.add_argument('-i',dest='inputFileName', metavar='./input_file.txt', type=str,
                               help='Input file describing the data',required=True)
    args = parser.parse_args()
    calculator = pyxsvs.pyxsvs(args.inputFileName,useFlatField=True)
    expList = sorted(calculator.Parameters['exposureParams'].keys())
    saveDir = calculator.Parameters['saveDir']
    imgFiles = glob(saveDir+'*.png')
    # Create the document
    fileName = calculator.Parameters['outPrefix']+'report.pdf'
    doc = SimpleDocTemplate(fileName)
    Story = [Spacer(1,1*inch)]
    style = styles["Normal"]
    # Add the plots
    for i in xrange(len(expList)):
        exposure = expList[i]
        # Find images for the current exposure:
        currImg = [img for img in sorted(imgFiles) if exposure in img]
        ptext = '<font size=14><b>%s</b></font>' % exposure
        Story.append(Paragraph(ptext, style))
        Story.append(Spacer(1,0.2*inch))
        for j in xrange(len(currImg)):
            bogustext = ("This is Paragraph number %s.  " % i) *20
            Story.append(get_image(currImg[j],width=5*inch,caption=bogustext))
            Story.append(Spacer(1,0.2*inch))
        #cfFig = self.cfFigures.pop(0)
        #cfFigName = cStringIO.StringIO()
        #cfFig.savefig(cfFigName,format='PDF')
        #img = PdfImage(cfFigName)
        #Story.append(img)
        #bogustext = ("This is Paragraph number %s.  " % i) *20
        #p = Paragraph(bogustext, style)
        #Story.append(p)
        Story.append(Spacer(1,0.2*inch))
        Story.append(Spacer(1,0.2*inch))
    doc.build(Story, onFirstPage=ReportFirstPage, onLaterPages=ReportLaterPages)

if __name__ == '__main__':
    main()
