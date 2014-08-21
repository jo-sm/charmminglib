#!/usr/bin/env python
"""
This script will produce a histogram plot of a replica exchange (REX)
run, one curve for each replica.  This allows the user to determine
if their REX run has adequate exchange statistics by observing the
overlap between the different histrogram curve.

Histograms are produced using the plaintext data present in the set of
.out files created by the REX script.

Example Usage;
% ./rexhistogram.py -I=run_rex/ld -O=runrex -T=nscale0.85 -F=png
"""


import os
import numpy as np
import matplotlib.pyplot as pyplot
from pychm.scripts.getprop import getProp
from pychm.tools import expandPath


def rexHist(iterable, title, write=False, imgFilePath=None):
    """
    `iterable`  *string         # An iterable of the CHARMM .out files
    `title`     string          # The title of the plot.
    `write`     [`False`,`"svg"`,`"png"`]   # The plot's output format,
                                            # False won't write a file.
    `imgFilePath` string        # The plot's file path, defaults
                                # to `$cwd/default_title`.
    """
    fileList = list(iterable)
    fileList = map(expandPath, fileList)
    tempArray = []
    for fileName in fileList:
        print 'Processing %s...' % fileName
        filePointer = open(fileName)
        rawData = getProp(filePointer, 'averener', 'avertemp')
        tempArray.append(np.array(rawData['avertemp']).mean()) # Calc Temp
        histData = np.histogram(rawData['averener'], 20, new=True)
        pyplot.plot(histData[1][1:],histData[0],label=os.path.basename(fileName))
    pyplot.legend( ['%3dK' % tempArray[i] for i in xrange(len(fileList))]
                ,loc=0, ncol=2)
    pyplot.xlabel(r'$Energy\ (kcal\ mol^{-1})$')
    pyplot.ylabel(r'$Frequency$')
    pyplot.title(r'%s' % title)
    if write:
        if imgFilePath is None:
            imgFilePath = '%s/%s.%s' % (os.path.dirname(fileList[0]), title,
                                        write)
        else:
            imgFilePath = expandPath(imgFilePath)
        pyplot.savefig(imgFilePath, format=write)
    else:
        pyplot.show()


if __name__ == '__main__':


    import sys
    from glob import glob
    from pychm.tools import OptionParser


    # Option Parsing
    useText =\
    """
    DOCME
    """
    optparser = OptionParser(usage=useText, version='%prog 0.1')
    # Optional
    optparser.add_option('-I' ,'--input', default='auto', metavar='DIR',
                        help='DIR where CHARMM .out files live.')
    optparser.add_option('-O', '--output', default='auto', metavar='PATH',
                        help='PATH where figure image is written.')
    optparser.add_option('-T', '--title', default='auto',
                        help='Title of output figure.')
    optparser.add_option('-F', '--format', default=False,
                        choices=[False, 'svg', 'png'],
                        help='Format of output figure, default is not to save.')
    # Parse
    (options, args) = optparser.parse_args(sys.argv)
    # Set defaults
    if options.input == 'auto':
        options.input = os.getcwd()
    if options.output == 'auto':
        options.output = None
    if options.title == 'auto':
        options.title = 'default_title'
    # Do Work
    outList = glob('%s/*.out_*' % options.input)
    if not outList:
        print 'No CHARMM .out files found, please try again!'
        print useText
        sys.exit(1)
    outList = map(expandPath,outList)
    outBase = os.path.basename(outList[0]).split('_')[0]
    outList = [ '%s/%s_%d' % (os.path.dirname(outList[0]), outBase, i)
            for i in xrange(len(outList)) ]
    #
    rexHist(outList, options.title, write=options.format,
            imgFilePath=options.output)
