#!/usr/bin/env python


import os
import sys
from glob import glob
import numpy as np
from pychm.cg.analysis.natq import NatQ
from pychm.scripts.getprop import getProp


#################
# Script Inputs #
#################
taco = NatQ('~/bigbox/chem/cg_paper/1prb/run_rex/1prb_a_pro_7_53.pdb')
nscale = 0.95
pdbcode = '1prb'

# input files
taco.inpFilename = '&/nscale_%4.2f/ld/rex.inp' % nscale
taco.rtfFilename = '&/nscale_%4.2f/ld/merged.rtf' % nscale
taco.prmFilename = '&/nscale_%4.2f/ld/merged.prm' % nscale
taco.psfFilename = '&/cg-vacuum.psf'
taco.crdFilename = '&/cg-vacuum.crd'
taco.dcdPathname = '&/nscale_%4.2f/ld' % nscale

# output path
anlPath = '&/nscale_%4.2f/test/rep_%02d'

# other options
taco.correlArrayLength = 10000
taco.charmmBin = 'cg_mscale'


###############################################################################
###############################################################################
###############################################################################

### Files to extract temperature data from
tempFileList = glob('%s/*out_*' % taco.dcdPathname)
tempFileBase = os.path.basename(tempFileList[0]).split('_')[0]
tempFileList = [ '%s/%s_%d' % (taco.dcdPathname, tempFileBase, i)
                for i in xrange(len(tempFileList)) ]
### Files to extract dynamics data from
dcdFileList = glob('%s/*dcd_*' % taco.dcdPathname)
dcdFileBase = os.path.basename(dcdFileList[0]).split('_')[0]
dcdFileList = [ '%s/%s_%d' % (taco.dcdPathname, dcdFileBase, i)
                for i in xrange(len(dcdFileList)) ]

### Temperature Data
tempList = []
for tempFile in tempFileList:
    tempList.append(np.array(getProp(open(tempFile), 'avertemp',
                                    stopIter=1000)['avertemp']).mean())

#############
# DO Correl #
#############


# Write/Run correl jobs
if 1:
    taco.correlStart = 0
    taco.correlStop = -2
    for i, dcdFile in enumerate(dcdFileList):
        print "Processing Correl jobs on replica %d." % i
        taco.anlPathname = anlPath % (nscale, i)
        taco.inpPathname = taco.anlPathname
        taco.outPathname = taco.anlPathname
        taco.dcdFilename = dcdFile
        comments = ['%s rex natq' % pdbcode,
                    'nscale = %4.2f' % nscale,
                    '%10.5f K' % tempList[i]
                ]
        taco.do_correl(comments=comments)
        taco.data


# Calculate DelG(T)
if 1:


    from pychm.analysis.delg import DelG


    delGList = []
    counts = []
    for i in xrange(len(dcdFileList)):
        taco.anlPathname = anlPath % (nscale, i)
        natQT = taco.get_natQofT()
        # DelG Calculation
        burrito = DelG(natQT, tempList[i])
        burrito.addState('folded', 0.5, 1)
        burrito.addState('unfolded', 0, 0.5)
        burrito.count()
        delGList.append(burrito.get_DelG('unfolded', 'folded'))
    # Filter out infinite energies (for plotting and regression)
    plotTemp = [ temp for i,temp in enumerate(tempList) if delGList[i] != np.inf ]
    plotDelG = [ delg for i,delg in enumerate(delGList) if delGList[i] != np.inf ]
    # Keep data where 250K < temp < 400K
    plotDelG = [ delg for i,delg in enumerate(plotDelG) if 250. <= plotTemp[i] <= 400. ]
    plotTemp = [ temp for i,temp in enumerate(plotTemp) if 250. <= plotTemp[i] <= 400. ]

    # Plotting
    if 1:


        from scipy.optimize import leastsq
        from scipy import linspace
        import matplotlib.pyplot as pyplot


        def residuals(v, x):
            """
            Define the functional form of the curve you want to fit to. 'v' is an array
            of the parameters, and 'x' is the dependant variable.
            """
            return v[0]*(1. - x/v[1]) - v[2]*((v[1]-x) + x*np.log(x/v[1]))


        def error(v, x, y):
            """
            This is my error function. There are many like it, but this one is mine. My
            error function is my best friend. It is my life. I must master it as I must
            master my life. My error function, without me, is useless. Without my error
            function, I am useless.
            """
            return residuals(v,x)-y


        x = np.array(plotTemp)
        y = np.array(plotDelG)

        ## Initial parameter values
        paramArray_0 = np.array([16.8,300,0.28])

        ## Fitting
        paramArray, success = leastsq(error, paramArray_0, args=(x,y), maxfev=10000)

        ## Plot
        print 'Estimater parameters: ', paramArray
        X = linspace(x.min(),x.max(),len(x)*5)
        pyplot.plot(x,y,'g^', X, residuals(paramArray,X),'k-')
        pyplot.xlabel(r'$Temperature\ (K)$')
        #pyplot.ylabel(r'$Radius\ of\ Gyration\ (\AA)$')
        pyplot.ylabel(r'$Gibbs\ Free\ Energy\ (kCal\ mol^-1)$')
        pyplot.title(r'%s Melting Curve' % pdbcode)
        if True:
            Format = 'png'
            pyplot.savefig('%s_delG.%s' % (pdbcode, Format), format=Format)
        else:
            pyplot.show()
