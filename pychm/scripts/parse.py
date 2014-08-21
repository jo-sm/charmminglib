#!/usr/bin/env python
"""
The standalone command line version of the CHARMMing input parser.

:Author: fcp
:Date: 10/26/2010

:Usage:
    ``parse.py --help`` will give you a help message explaining the various
    options.  "Required options" are marked as such, and defaults
    appear in [brackets].


:Options:
    | ``--version``   show program's version number and exit
    | ``-h, --help``  show this help message and exit
    | ``-I PATH, --input=PATH``
        *Required* PATH of input .pdb file
    | ``-O DIR``, ``--output=DIR``
        DIR where output .pdb files are written [$inputPath]
    | ``-F FORM``, ``--informat=FORM``
        specify the .pdb formatting to expect [auto]
    | ``--outformat=FORM``
        specify the .pdb formatting to output [charmm]
    | ``-M NUM``, ``--model=NUM``
        Specify the model NUM to produce output files from [auto]
    | ``--no_fix_chainid``
        Disable chainid auto fixing
    | ``--no_fix_atom``
        In lieu of quiet auto fixing of atoms and structs, angry
        exceptions are raised.
    | ``-P``, ``--pickle``
        Write `output.chk` files, which are just python pickles of
        PDBFile objects.  Useful for restarting scripts, debugging, etc.
    | ``-V``, ``--verbose``
        Write extra debugging information
    | ``--old_resid``
        Write .pdb files using the canonical `resid` values instead
        of the reindexed CHARMM values.
"""


from cPickle import dump
from pychm.io.pdb import PDBFile
from pychm.tools import expandPath, mkdir, lowerKeys


def parse(pdbFilename, **kwargs):
    """
    Parse a *.pdb* plain text file into its constituent chains and segments, and
    print one CHARMM formatted *.pdb* file per chain/segment combination.

    *kwarg defaults are listed first*

    **kwargs:**
        | ``informat``      ['auto', 'pdborg', 'charmm'] # 'auto' -> detects input formatting
        | ``outformat``     ['charmm', 'pdborg', 'auto'] # 'auto' -> same as input formatting
        | ``outpath``       ['auto', user_specified_path] # 'auto' -> same as pdbFilename path
        | ``fix_chainid``   [True, False]
        | ``autofix``       [True, False]
        | ``modelnum``      ['auto', 0, 1, 2, ...] # 'auto' -> uses the first model found
        | ``pickleme``      [False, True] # pickles the PDBFile object for later use
        | ``verbose``       [False, True]
        | ``old_resid``     [False, True]

    >>> parse('~/pychm/1yjp/1yjp.pdb',outpath='~',pickleme=True)
    """
    # kwargs
    kwargs = lowerKeys(kwargs)
    inFormat = kwargs.get('informat', 'auto')
    outFormat = kwargs.get('outformat', 'charmm')
    outPath = kwargs.get('outpath', 'auto')
    fix_chainid = kwargs.get('fix_chainid', True)
    autoFix = kwargs.get('autofix', True)
    modelNum = kwargs.get('modelnum', 'auto')
    pickleMe = kwargs.get('pickleme', False)
    verbose = kwargs.get('verbose', False)
    old_resid = kwargs.get('old_resid', False)
    # Repackage the PDBFile kwargs, make pdbFile object
    pdbFileArgs = {'informat':inFormat, 'fix_chainid':fix_chainid,
                'autofix':autoFix, 'verbose':verbose}
    pdb = PDBFile(pdbFilename, **pdbFileArgs)
    if verbose:
        print '%s: Output formatting set to `%s`' % (pdb.code, outFormat)
    # Get model number...
    if modelNum == 'auto':
        thisMol = pdb.iter_models().next()
    else:
        thisMol = pdb[modelNum]
    if verbose:
        print '%s: Loading `%s`' % (pdb.code, thisMol.name)
    # Determine explicit output Path
    if outPath == 'auto':
        outPath = pdb.path
    else:
        outPath = expandPath(outPath)
    mkdir(outPath)
    if verbose:
        print '%s: Output path set to `%s`' % (pdb.code, outPath)
    # Do Work
    thisMol.parse()
    if verbose:
        print '%s: Parsing `%s`' % (pdb.code, thisMol.name)
    # Write CHARMMing style output.
    segDict = {'nuc':'nuc', 'pro':'pro', 'good':'goodhet', 'bad':'het',
            'dna':'dna', 'rna':'rna'}
    stdoutList = []
    # Write pdb files
    writeArgs = {'outformat':outFormat, 'ter':True, 'end':False,
                'old_resid':old_resid}
    for seg in thisMol.iter_seg():
        stdoutList.append('%s-%s' % (seg.chainid, segDict[seg.segType]))
        name = '%s/new_%s-%s-%s.pdb' % (outPath, thisMol.code, seg.chainid,
                                        segDict[seg.segType])
        if verbose:
            print '%s: Writing output to file `%s`' % (pdb.code, name)
        seg.write(filename=name, **writeArgs)
    # Write pickle (chk) file
    if pickleMe:
        pickleFilename = '%s/%s.chk' % (outPath, pdb.code)
        pickleFile = open(pickleFilename,'w')
        dump(pdb,pickleFile)
        pickleFile.close()
        if verbose:
            print '%s: Writing pickle to file `%s`' % (pdb.code, pickleFilename)
    if verbose:
        print '\n\nEnd of verbosity\n\n'
    # To STDOUT
    print 'natom=%d' % len(thisMol)
    print 'nwarn=%d' % len(thisMol.warnings)
    if thisMol.warnings:
        print 'warnings=%r' % thisMol.warnings
    print 'seg=%r' % stdoutList


if __name__ == '__main__':


    import sys
    from pychm.tools import OptionParser


    # Option Parsing
    useText =\
    """
    '%prog --help' will give you a help message explaining the various
    options.  "Required options" are marked as such, and defaults
    appear in [brackets].
    """
    optparser = OptionParser(usage=useText, version='%prog 0.1')
    # Required
    optparser.add_option('-I', '--input', required=True, metavar='PATH',
                    help='PATH of input .pdb file')
    # Optional
    optparser.add_option('-O', '--output', default='auto', metavar='DIR',
                    help='DIR where output .pdb files are written [$inputPath]')
    optparser.add_option('-F', '--informat', default='auto',
                        choices=['auto', 'pdborg', 'charmm'], metavar='FORM',
                        help="specify the .pdb formatting to expect [%default]")
    optparser.add_option('--outformat', default='charmm',
                        choices=['charmm', 'pdborg', 'debug', 'xdebug', 'crd',
                                'xcrd'],
                        metavar='FORM',
                        help="specify the .pdb formatting to output [%default]")
    optparser.add_option('-M', '--model', default='auto', metavar='NUM',
                    help='Specify the model NUM to produce output files from [%default]')
    optparser.add_option('--no_fix_chainid', action='store_false', default=True,
                    help='Disable chainid auto fixing')
    optparser.add_option('--no_fix_atom', action='store_false', default=True,
                    help='In lieu of quiet auto fixing of atoms and structs,\
                    angry exceptions are raised.')
    optparser.add_option('-P' ,'--pickle', action='store_true', default=False,
                    help='Write `output.chk` files, which are just python pickles\
                    of PDBFile objects.  Useful for restarting scripts, debugging, etc.')
    optparser.add_option('-V' ,'--verbose' ,action='store_true' ,default=False,
                    help='Write extra debugging information')
    optparser.add_option('--old_resid', action='store_true', default=False,
                    help='Write .pdb files using the canonical `resid` values\
                        instead of the reindexed CHARMM values.')
    # Parse
    (options, args) = optparser.parse_args(sys.argv)
    # Repackage options into kwargs
    kwargs = {
        'outpath':options.output, 'informat':options.informat,
        'outformat':options.outformat, 'modelnum':options.model,
        'fix_chainid':options.no_fix_chainid, 'autofix':options.no_fix_atom,
        'pickleme':options.pickle, 'verbose':options.verbose,
        'old_resid':options.old_resid
        }
    # Do Work
    parse(options.input, **kwargs)
