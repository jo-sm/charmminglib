"""
DOCME
"""
# fcp
# 10/26/2010


import re
import os
import itertools
from charmming.const.etc import alphanum
from charmming.lib.mol import Mol
from charmming.lib.atom import Atom
from charmming.tools import Property, expandPath, cleanStrings, paragraphs,\
        lowerKeys


def get_formatting(input):
    """
    Takes a string representing the location of a .pdb file or an
    iterator of strings as input, and returns a string indicating the
    formatting of the pdb data: "pdborg", "charmm" or "unknown".
    """
    try:
        iterator = ( line.lower() for line in open(input) if
                    line.lower().startswith(('atom', 'hetatm')) )
    except IOError:
        iterator = ( line.lower() for line in input
                    if line.lower().startswith(('atom', 'hetatm')) )
    for line in iterator:
        if line[21:22] == ' ' and line[72:73] in alphanum:
            return 'charmm'
        elif line[21:22] in alphanum and line[12:14].strip() == line[66:].strip():
            return 'pdborg'
    return 'unknown'


class PDBFile(object):
    """
    Class Attributes
        `_autoInFormat`
        `_reCode`
    Public Attributes
        `warnings`
    Properties
        `code`
        `crd`
        `filename`
        `footer`
        `header`
        `inFormat`
        `path`
    Public Methods
        `get_warnings`  TODO :: Sophisticated warning handling
        `iter_all`
        `iter_models`
        `keys`
    Private Methods
        _get_formatting
        _partitions
        _fix_chainids
        _fix_resids     TODO
        _build_models
    Special Methods
    """

    _autoInFormat = 'pdborg'
    """
    A string that defines the default input formatting for all class
    instances.
    """

    _reCode = re.compile(r"(?<![a-z0-9])[0-9][a-z0-9]{3}(?![a-z0-9])")
    """
    A regular expression for deriving the PDB code from the filename of
    a .pdb file named using CHARMMing conventions.
    """

    def __init__(self, filename, **kwargs):
        """
        kwargs:
            `informat`      ['auto','pdborg','charmm']
            `fix_chainid`   [True,False]    # Attempts to fix mangled chainid
            `fix_resid`     [True,False]    # Attempts to fix mangled resid
            `autofix`       [True,False]    # Flag for atom._autoFix
            `verbose`       [False,True]
        """
        super(PDBFile, self).__init__()
        # kwargs
        kwargs = lowerKeys(kwargs)
        inFormat = kwargs.pop('informat', 'auto')
        fix_chainid = kwargs.pop('fix_chainid', True)
        fix_resid = kwargs.pop('fix_resid', True)
        self._autoFix = kwargs.pop('autoFix', True)
        self._verbose = kwargs.pop('verbose', False)
        #
        self.warnings = []
        # Filename
        self.filename = filename
        if self._verbose:
            print 'Parsing data from `%s`\n' % filename
        # Partitioning
        self._partitions()          # self.header / self.crd / self.footer
        # Input Formatting
        if inFormat == 'auto':
            self._inFormat = self._get_formatting(self.filename)
        else:
            self._inFormat = inFormat
        if self._verbose:
            print '%s: Input formatting set to `%s`' % (self.code, self.inFormat)
        # auto fixing
        if fix_chainid:
            if self._verbose:
                print '%s: Fixing `chainid`s' % self.code
            self._fix_chainids()
        #TODO Implement fix_resids
        if fix_resid:
            if self._verbose:
                print '%s: Fixing `resid`s' % self.code
            #self._fix_resids()
        # Processing
        self._models = {}
        self._build_models()

##############
# Properties #
##############

    @Property
    def code():
        doc =\
        """
        The PDB acquisition code of the `PDBFile`.
        """
        def fget(self):
            tmp = [ line for line in self.header if line.startswith('header') ]
            if tmp:
                return tmp[0].split()[-1]
            else:
                searchThis = os.path.basename(self.filename).lower()
                try:
                    return PDBFile._reCode.findall(searchThis)[0]
                except IndexError:
                    return '????'
        return locals()

    @Property
    def crd():
        doc =\
        """
        The atomic coordinate data from the .pdb file, in its
        unadulterated text form.
        """
        def fget(self):
            return self._crd
        return locals()

    @Property
    def filename():
        doc =\
        """
        The name of the .pdb file from which this `PDBFile` instance
        was created.
        """
        def fget(self):
            return self._filename
        def fset(self, value):
            self._filename = expandPath(value)
        return locals()

    @Property
    def footer():
        doc =\
        """
        The metadata in the .pdb file that follows the atomic coordinate
        data.  This property gives access to its unadulterated text form.
        """
        def fget(self):
            return self._footer
        return locals()

    @Property
    def header():
        doc =\
        """
        The metadata in the .pdb file that precedes the atomic
        coordinate data.  This property gives access to its
        unadulterated text form.
        """
        def fget(self):
            return self._header
        return locals()

    @Property
    def inFormat():
        doc =\
        """
        The formatting of the input data: "pdborg" or "charmm".
        """
        def fget(self):
            return self._inFormat
        return locals()

    @Property
    def path():
        doc =\
        """
        The path of the .pdb file from which this `PDBFile` instance
        was created.
        """
        def fget(self):
            return os.path.dirname(self.filename)
        return locals()

##################
# Public Methods #
##################

    def get_warnings(self):
        """
        Return a list of warnings generated by the PDBFile object
        itself, and all of the Mol objects it contains.
        """
        result = self.warnings[:]
        for mol in self.iter_all():
            result.extend(mol.warnings)
        return result

    def iter_all(self):
        """
        Iterate over all of the PDBFile's Mol objects, regardless of
        their origin.
        """
        return itertools.chain(self.iter_models())

    def iter_models(self):
        """
        Iterate over the PDBFile's models. Models are Mol objects
        derived directly from the original pdb file.
        """
        return ( self['model%d' % i] for i in sorted(self._models.keys()) )

    def keys(self):
        """
        Lists the keys to access all of the Mol objects the PDBFile
        object contains.
        """
        models = [ 'model%d' % i for i in sorted(self._models.keys()) ]
        return models

###################
# Private Methods #
###################

    def _build_models(self):
        """
        Parse the crd section, and load the coordinates into `Mol`
        objects, one `Mol` object per model section in the .pdb file.
        """
        models = paragraphs(self.crd, splitter=['model'])
        for model in models:
            if model[0].startswith('model'):
                modelNum = int(model[0].split()[1])
            else:
                modelNum = 0
            iterator = ( Atom(text=line, informat=self.inFormat, index=i,
                        autofix=self._autoFix) for i, line in enumerate(model)
                        if line.startswith(('atom', 'hetatm')) )
            self._models[modelNum] = Mol(iterable=iterator, name='model%d' %
                                        modelNum, code=self.code, autofix=True)

    def _fix_chainids(self):
        """
        Detect and correct chainid mangling.
        """
        # Detect mangling
        def gen():
            for line in self.crd:
                if line.startswith(('atom', 'hetatm')):
                    tmp = Atom(text=line, informat=self.inFormat, autofix=False)
                else:
                    continue
                yield tmp.chainid
        # Correct mangling
        if '' in set(gen()):
#WARN       # Throw a warning if chainids are mangled
            self.warnings.append('chainids mangled')
            print 'One or more chainids appear to be missing, attempting to\
                    autofix.\n'
            print 'Please verify the accuracy of your .pdb file upon\
                    completion.\n'
            chainNum = 0
            for i,line in enumerate(self.crd):
                if line.startswith(('atom', 'hetatm')):
                    tmp = Atom(text=line, informat=self.inFormat, autofix=False)
                    # Blank chainid
                    if not tmp.chainid:
                        tmp.chainid = alphanum[chainNum]
                        self.crd[i] = tmp.Print(outformat=self.inFormat).lower()
                elif line.startswith('ter'):
                    chainNum += 1
                elif line.startswith('model'):
                    chainNum = 0
                else:
                    continue

    def _fix_resids(self):
        """
        Detect and correct resid mangling.
        """
        raise NotImplementedError

    def _get_formatting(self, filename):
        """
        Wrapper method to detect PDB text formatting, defaults to 'pdborg'
        """
        result = get_formatting(filename)
        if result == 'unknown':
#WARN       # Throw a warning if formatting is guessed.
            self.warnings.append('Undetected pdb formatting')
            return self.__class__._autoInFormat
        else:
            return result

    def _partitions(self):
        """
        Partition the file into _header/_crd/_footer sections.
        """
        filePointer = open(self.filename)
        # Populate header
        self._header = []
        tmp = None
        for line in cleanStrings(filePointer):
            if line.startswith(('atom', 'hetatm', 'model')):
                tmp = line
                break
            else:
                self._header.append(line)
        # Populate coordinates
        if tmp is None:
            self._crd = []
        else:
            self._crd = [tmp]
        tmp = None
        for line in cleanStrings(filePointer):
            if not line.startswith(('atom', 'anisou', 'hetatm', 'model', 'ter',
                                    'endmdl')):
                tmp = line
                break
            else:
                self._crd.append(line)
        # Populate footer
        if tmp is None:
            self._footer = []
        else:
            self._footer = [tmp]
        for line in cleanStrings(filePointer):
            self._footer.append(line)
        filePointer.close()

###################
# Special Methods #
###################

    def __getitem__(self, key):
        if key in self.keys():
            try:
                if key.startswith('model'):
                    return self._models[int(key.split('model')[1])]
                else:
                    raise IndexError
            except IndexError:
                raise KeyError('Unknown Key: %s' % key)
        else:
            try:
                dummy = int(key)
                try:
                    return self._models[key]
                except IndexError:
                    raise KeyError('Unknown Key: %s' % key)
            except ValueError:
                raise KeyError('Unknown Key: %s' % key)

    def __len__(self):
        return len(self._models)

    def __contains__(self, key):
        return key in self.keys()

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.keys())
