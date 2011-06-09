"""
:Author: fcp
:Date: 10/26/2010
"""


import re
import os
import itertools
from charmming.const import alphanum
from charmming.tools import Property, expandPath, cleanStrings, paragraphs,\
        lowerKeys
from charmming.lib.atom import Atom
from charmming.lib.mol import Mol, MolError


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


def get_simpleMol(filename, **kwargs):
    """
    A function for creating single ``Mol`` objects from a .pdb text file.
    Useful if you do not need the additional features provided by a
    ``PDBFile`` object.  If the input .pdb text file contains more than one
    model, a ``MolError`` will be raised due to ambiguity.

    **kwargs:**
        *See* ``PDBFile`` *for valid kwargs.*
    """
    tmp = PDBFile(filename, **kwargs)
    if len(tmp) > 1:
        raise MolError("This function is designed to be used on .pdb text \
                    files containing one model.")
    for taco in tmp.iter_models():
        return taco


class PDBFile(object):
    """
    The ``PDBFile`` object serves as a container of ``Mol`` objects which are
    logically united and derive from the same plain text .pdb file.


    Calling the constructor on a plain text .pdb file instantizes one ``PDBFile``
    object and *n* ``Mol`` objects, one for each *model* present in the .pdb file.
    Furthermore, the meta-data contained in the .pdb file is parsed and stored in
    a dictionary for later use.  The ``Mol`` objects created at instantisation
    may be accessed by the ``PDBFile`` object in either a dictionary-like way or
    a list-like way.  Files which contain only one model will store the ``Mol``
    object as *model 0*.

    For example, to access *model 4* from the 1o1o.pdb file one could write:

    >>> taco = PDBFile('1o1o.pdb')[4]

    or...

    >>> taco = PDBFile('1o1o.pdb')['model4']
    
    Furthermore, this object may also be used to couple other ``Mol`` objects to
    each other.  Such as a structure that has been patched, or a minimized
    structure.  The syntax for doing this will be dictionary-like and is TBD.


    **kwargs:**
        | ``informat``      ['auto','pdborg','charmm']
        | ``fix_chainid``   [True,False]    # Attempts to fix mangled chainid
        | **TODO** ``fix_resid``     [True,False]    # Attempts to fix mangled resid
        | ``autofix``       [True,False]    # Flag for atom._autoFix
        | ``verbose``       [False,True]

    :TODO:
        | ``get_warnings`` :: Sophisticated warning handling
        | ``fix_resids`` :: Automagic resid mangling repairs
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
        super(PDBFile, self).__init__()
        # kwargs
        kwargs = lowerKeys(kwargs)
        inFormat = kwargs.get('informat', 'auto')
        fix_chainid = kwargs.get('fix_chainid', True)
        fix_resid = kwargs.get('fix_resid', True)
        self._autoFix = kwargs.get('autoFix', True)
        self._verbose = kwargs.get('verbose', False)
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
        A ``property`` for the PDB accession code of the ``PDBFile``.
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
        A ``property`` for the atomic coordinate data from the .pdb
        file, in its unadulterated text form.
        """
        def fget(self):
            return self._crd
        return locals()

    @Property
    def filename():
        doc =\
        """
        A ``property`` for the name of the .pdb file from which this
        ``PDBFile`` instance was created.
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
        A ``property`` for the metadata in the .pdb file that follows
        the atomic coordinate data.  This property gives access to its
        unadulterated text form.
        """
        def fget(self):
            return self._footer
        return locals()

    @Property
    def header():
        doc =\
        """
        A ``property`` for the metadata in the .pdb file that precedes
        the atomic coordinate data.  This property gives access to its
        unadulterated text form.
        """
        def fget(self):
            return self._header
        return locals()

    @Property
    def inFormat():
        doc =\
        """
        A ``property`` for the formatting of the input .pdb text file
        data, either *"pdborg"* or *"charmm"*.
        """
        def fget(self):
            return self._inFormat
        return locals()

    @Property
    def path():
        doc =\
        """
        A ``property`` for the dirname where the ``PDBFile`` instance's
        parent .pdb text file resides.
        """
        def fget(self):
            return os.path.dirname(self.filename)
        return locals()

##################
# Public Methods #
##################

    def get_warnings(self):
        """
        Return a list of warnings generated by the ``PDBFile`` object
        itself, and all of the ``Mol`` objects it contains.
        """
        result = self.warnings[:]
        for mol in self.iter_all():
            result.extend(mol.warnings)
        return result

    def get_metaData(self):
        """
        Returns a :class:`dict` containing metadata parsed from the
        ``PDBFile`` object's ``_header`` and ``_footer``.  Each key
        corresponds to a section in the header, delimited by the first
        string in each line.
        """
        tmp = {}
        for line in self._header:
            key = line.split()[0]
            value = line.split(key)[1].lstrip()
            if key not in tmp.keys():
                tmp[key] = [value]
            else:
                tmp[key].append(value)
        for line in self._footer:
            key = line.split()[0]
            value = line.split(key)[1].lstrip()
            if key not in tmp.keys():
                tmp[key] = [value]
            else:
                tmp[key].append(value)
        return tmp

    def iter_all(self):
        """
        Iterate over all of the ``Mol`` objects contained in the
        current ``PDBFile`` instance, regardless of origin.
        """
        return itertools.chain(self.iter_models())

    def iter_models(self):
        """
        Iterate over the ``Mol`` objects derived from the .pdb file's
        *models*.
        """
        return ( self['model%d' % i] for i in sorted(self._models.keys()) )

    def keys(self):
        """
        Lists the keys to access all of the ``Mol`` objects the ``PDBFile``
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
                    tmp = Atom(text=line, informat=self.inFormat, autofix=True)
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
                    tmp = Atom(text=line, informat=self.inFormat, autofix=True)
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
        iterator = ( line for line in cleanStrings(open(self.filename) ))
        # Populate header
        self._header = []
        tmp = None
        for line in iterator:
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
        for line in iterator:
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
        for line in iterator:
            self._footer.append(line)

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
