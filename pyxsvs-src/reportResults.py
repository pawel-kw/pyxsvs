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
import os
import time

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

global caption_base,report_title

def ReportFirstPage(canvas, doc):
    '''First page of the PDF report
    '''
    Title = report_title
    pageinfo = report_title
    canvas.saveState()
    canvas.setFont('Times-Bold',16)
    canvas.drawCentredString(PAGE_WIDTH/2.0, PAGE_HEIGHT-98, Title)
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.75 * inch, "Page %d / %s" % (doc.page, pageinfo))
    canvas.restoreState()

def ReportLaterPages(canvas, doc):
    '''Later pages format
    '''
    pageinfo = report_title
    canvas.saveState()
    canvas.setFont('Times-Roman',9)
    canvas.drawString(inch, 0.75 * inch, "Page %d / %s" % (doc.page, pageinfo))
    canvas.restoreState()

def get_image(path, width=1*inch):
    img = utils.ImageReader(path)
    iw, ih = img.getSize()
    aspect = ih / float(iw)
    resImg = Image(path, width=width, height=(width * aspect))
    return resImg

def get_caption(img_file_path):
    '''Generates automated caption based on the image file name
    '''
    file_dir, file_name = os.path.split(img_file_path)
    caption = ''
    for key in caption_base:
        if key in file_name:
            caption = caption_base[key]
    return caption


caption_base = {'fit_params': 'Fit parameters as a function of <i>q</i>. \
                                <font color="blue">Blue</font> - the\
                                value of contrast <i>C = 1/M</i>.\
                                <font color="red">Red</font> - number of modes\
                                <i>M</i>, obtained from fitting the histograms. \
                                <font color="green">Green</font> - \
                                average number of counts <i>K</i>.',
                'hist_fits': 'Photon counts histograms with model fits. \
                                <font color="blue">Blue symbols</font> \
                                mark the histogram values calculated from the data.\
                                The \
                                <font color="blue">blue line</font> \
                                is the Poisson distribution for the \
                                corresponding number of average counts <i>K</i>.\
                                The <font color="red">red</font> \
                                line is the fitted negative-binomial distribution\
                                for <i>K</i> and <i>M</i> given. <i>K</i> is calculated \
                                directly from the data. <i>M</i> is the only free parameter\
                                used in the fit.',
                'q_mask': 'Averaged scattering pattern. The shaded pixels mark <i>q</i> \
                            partitions.',
                'trace': 'Intensity averaged over all pixels in the corresponding <i>q</i>\
                            partition as a function of frame number. The \
                            <font color="red">red line</font> marks\
                            the average number of photons <i>K</i>.'
                }


def main():
    global report_title
    parser = argparse.ArgumentParser(description='Python XSVS data analysis.')
    parser.add_argument('-i',dest='inputFileName', metavar='./input_file.txt', type=str,
                               help='Input file describing the data',required=True)
    args = parser.parse_args()
    calculator = pyxsvs.pyxsvs(args.inputFileName,useFlatField=True)
    expList = sorted(calculator.Parameters['exposureParams'].keys())
    saveDir = calculator.Parameters['saveDir']
    report_title = 'XSVS data analysis: %s' % (calculator.Parameters['figTitle'])
    imgFiles = glob(saveDir+'*.png')
    # Create the document
    fileName = saveDir+calculator.Parameters['outPrefix']+'report.pdf'
    doc = SimpleDocTemplate(fileName)
    Story = [Spacer(1,1*inch)]
    results_file = glob(saveDir+'*.p')
    style = styles["Normal"]
    main_info = 'Analysis done: %s <br />\
                 Partitions of q defined by: <br />\
                 q1 = %.4f <br /> q2 = %.4f <br />\
                 qs = %.4f <br /> dq = %.4f' % (time.ctime(os.path.getctime(results_file[0])),     
                                                calculator.Parameters['q1'],
                                                calculator.Parameters['q2'],
                                                calculator.Parameters['qs'],
                                                calculator.Parameters['dq'])
    Story.append(Paragraph(main_info, style))
    Story.append(Spacer(1,0.4*inch))
    # Add the plots
    for i in xrange(len(expList)):
        exposure = expList[i]
        # Find images for the current exposure:
        currImg = [img for img in sorted(imgFiles) if exposure in img]
        section_caption = 'Exposure %d' % int(exposure.split('_')[-1])
        ptext = '<font size=14><b>%s</b></font>' % section_caption
        Story.append(Paragraph(ptext, style))
        Story.append(Spacer(1,0.4*inch))
        currParams = calculator.Parameters['exposureParams'][exposure]
        exp_info_txt = 'Number of data frames: %d' % (currParams['n2']-currParams['n1']+1,)
        Story.append(Paragraph(exp_info_txt, style))
        Story.append(Spacer(1,0.2*inch))
        data_file_list = pyxsvs.filename(calculator.Parameters['dataDir']+currParams['dataPref'],
                                          currParams['dataSuf'],
                                          currParams['n1'],
                                          currParams['n2'],
                                         )
        first_header = pyxsvs.fabio.open(data_file_list[0]).header
        time_txt = 'Acquisition started on %s' % (first_header['time'])
        Story.append(Paragraph(time_txt, style))
        Story.append(Spacer(1,0.2*inch))
        for j in xrange(len(currImg)):
            bogustext = ("This is Paragraph number %s.  " % i) *20
            Story.append(get_image(currImg[j],width=5*inch))
            Story.append(Spacer(1,0.1*inch))
            ctxt = get_caption(currImg[j])
            caption = Paragraph(r'<b>Fig. %d.%d:</b> %s' % (i+1,j+1,ctxt), style)
            Story.append(caption)
            Story.append(Spacer(1,0.4*inch))
        Story.append(Spacer(1,0.4*inch))
    doc.build(Story, onFirstPage=ReportFirstPage, onLaterPages=ReportLaterPages)

if __name__ == '__main__':
    main()
