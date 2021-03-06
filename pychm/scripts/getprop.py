#!/usr/bin/env python
"""
A command line CHARMM output file parser.

| Originally implemented by Rick Venable in `csh`
| Rich Pastor Group
| LCB, NHLBI, NIH

:Author: fcp
:Date: 08/05/2010
"""


import sys
from pychm.tools import paragraphs


propPrefixes = ['dyna','aver','fluc','lave','lflc']


propSuffixes = {
    'step':(0,5,14),
    'time':(0,14,27),  'tote':(0,27,40),  'totk':(0,40,53),  'ener':(0,53,66),  'temp':(0,66,79),
    'grms':(1,14,27),  'hfct':(1,27,40),  'hfck':(1,40,53),  'ehfc':(1,53,66),  'virk':(1,66,79),
    'bond':(2,14,27),  'angl':(2,27,40),  'urey':(2,40,53),  'dihe':(2,53,66),  'impr':(2,66,79),
     'vdw':(3,14,27),  'elec':(3,27,40),  'hbon':(3,40,53),   'asp':(3,53,66),  'user':(3,66,79),
    'vire':(4,14,27),  'viri':(4,27,40),'presse':(4,40,53),'pressi':(4,53,66),  'volu':(4,66,79)
    }


helpString = """
Usage %s prop1 prop2 ... propN fileName

Each prop is of the form prefixsuffix, for example avertote.
Valid prefixes are:
    %r\n
Valid suffixes are:
    %r\n
""" % (sys.argv[0], propPrefixes, propSuffixes.keys())


def parse_prop(prop):
    prefix = prop[:4]
    suffix = prop[4:]
    if prefix not in propPrefixes:
        print helpString
        raise TypeError('parse_prop: invalid prefix %s' % prop[:4])
    while 1:
        if suffix:
            if suffix in propSuffixes.keys():
                return (prop, prefix, suffix)
            else:
                suffix = suffix[:-1]
        else:
            print helpString
            raise TypeError('parse_prop: invalid suffix %s' % prop[4:])


def getProp(iterable,*props,**kwargs):
    """
    Returns a dictionary whose keys are the requested properties, and
    whose values are list(prop(t)).  Two key/value pairs are created
    per property: propX and propX_time.

    Usage: getProp(iterable,props*)
    Example:
    >>> getProp(open('charmm.out'),'averener','dynavdw')
    
    Each prop is of the form prefixsuffix, for example 'avertote'

    Valid prefixes are:
        ['dyna','aver','fluc','lave','lflc']

    Valid suffixes are:
        ['ener','hbon','volu','dihe','virk','vire','hfct','tote',
        'totk','presse','pressi','ehfc','asp','hfck','step','user',
        'urey','grms','impr','elec','vdw','temp','angl','time','viri',
        'bond']

    Valid kwargs:
        `stopIter` = int
    """
    props = [ prop.lower() for prop in props ]
    assert len(props) >= 1
    # Make sure properties are valid
    propList = []
    for prop in props:
        propList.append(parse_prop(prop))
    # kwargs
    if 'stopIter' in kwargs.keys():
        stopIter = int(kwargs['stopIter'])
    else:
        stopIter = None
    # Initialize Lists
    outDict = {}
    for prop in propList:
        outDict['%s_time' % prop[0]] = []
        outDict[prop[0]] = []
    stepIterator = paragraphs(iterable, ['DYNA>', ' DYNA>'])
    stepIterator.next() # Discard Header
    stepIterator.next() # Discard Step 0
    # Do Work
    for i,step in enumerate(stepIterator):
        if stopIter:
            if i >= stopIter:
                return outDict
        subSteps = paragraphs(step,[' ----------'])
        for subStep in subSteps:
            # Strip extraneous lines
            tmp = [ line for line in subStep if not line.startswith(' ') ]
            # Ignore empty blocks
            if tmp:
                # Ignore label block
                if tmp[0].startswith('AVER DYN'):
                    continue
                tmpSuffix = tmp[0][0:4].lower()
                for prop in propList:
                    if prop[1] == tmpSuffix:
                        lineNum,start,stop = propSuffixes[prop[2]]
                        outDict['%s_time' % prop[0]].append(i + 1)
                        outDict[prop[0]].append(float(tmp[lineNum][start:stop]))
    return outDict


if __name__ == '__main__':


    import gzip
    import os
    from itertools import izip


    # Parse command line
    try:
        assert len(sys.argv) > 2
    except AssertionError:
        print '\n\nnot enough arguments\n\n'
        print helpString
        sys.exit(0)
    try:
        props = sys.argv[1:-1]
        props = [ prop.lower() for prop in props ]
        fileName = sys.argv[-1]
    except IndexError:
        print '\n\narguments in wrong order\n\n'
        print helpString
        sys.exit(0)

    # Create filepointer
    if fileName.endswith('gz'):
        print '\n\nWARNING - Operating on .gz files can be excruciatingly slow.\n\n'
        filePointer = gzip.open(fileName)
    else:
        filePointer = open(fileName)

    # Main
    outDict = getProp(filePointer, *props)
    for prop in props:
        outFileName = '%s_%s.dat' % (os.path.basename(fileName), prop)
        print 'Writing data to %s' % outFileName
        iterator = izip(outDict['%s_time' % prop], outDict[prop])
        tmp = [ '%12s %15s' % ('time', prop) ]
        for item in iterator:
            tmp.append( '%12d %15.5f' % (item[0], item[1]) )
        write_to = open('%s%s%s' % (os.getcwd(), os.sep, outFileName), 'w')
        write_to.write('\n'.join(tmp))
        write_to.close()
