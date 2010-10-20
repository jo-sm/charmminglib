import re
import os
import itertools
from charmming.const.etc import alphanum
from charmming.lib.mol import Mol
from charmming.lib.atom import Atom
from charmming.tools import Property,expandPath,cleanStrings,paragraphs


def get_formatting(input):
    """
    Takes a pdbFile or an iterator of strings, and returns a string
    indicating if it is 'pdborg', 'charmm' or 'unknown' format.
    """
    try:
        iterator = ( line.lower() for line in open(input) if
                    line.lower().startswith(('atom','hetatm')) )
    except IOError:
        iterator = ( line.lower() for line in input if line.lower().startswith(('atom','hetatm')) )
    for line in iterator:
        if line[21:22] == ' ' and line[72:73] in alphanum:
            return 'charmm'
        elif ( line[21:22] in alphanum ) and ( line[12:14].strip() == line[66:].strip() ):
            return 'pdborg'
    return 'unknown'


class PDBFile(object):
    """
    Private Attributes
    Properties
        filename
        path
        code
        format
        header
        crd
        footer
    Public Methods
        keys
        iter_models
        iter_structs
        iter_jobs
        iter_all
        append_struct   Untested
        append_job      Untested
        get_warnings    #TODO# Sophisticated warning handling
    Private Methods
        _get_formatting
        _partitions
        _fix_chainids
        _fix_resids     STUB
        _build_models
    Special Methods
    """

    _reCode = re.compile(r"(?<![a-z0-9])[0-9][a-z0-9]{3}(?![a-z0-9])")
    _defaultFormatting = 'pdborg'

    def __init__(self,filename,**kwargs):
        """
        kwargs:
            format      ['auto','pdborg','charmm']
            fix_chainid [True,False]    # Attempts to fix mangled chainid
            fix_resid   [True,False]    # Attempts to fix mangled resid
            autoFix     [True,False]    # Flag for atom._autoFix
            verbose     [False,True]
        """
        # kwargs
        format      = kwargs.pop('format','auto')
        fix_chainid = kwargs.pop('fix_chainid',True)
        fix_resid   = kwargs.pop('fix_resid',True)
        self._autoFix = kwargs.pop('autoFix',True)
        self._verbose = kwargs.pop('verbose',False)
        if kwargs.keys():
            raise TypeError('Unprocessed kwargs(%r)' % kwargs.keys())
        #
        self.warnings = []
        # Filename
        self.filename = filename
        if self._verbose:
            print 'Parsing data from `%s`\n' % filename
        # Partitioning
        self._partitions()          # self.header / self.crd / self.footer
        # Input Formatting
        if format == 'auto':
            self._format = self._get_formatting(self.filename)
        else:
            self._format = format
        if self._verbose:
            print '%s: Input formatting set to `%s`' % (self.code,self.format)
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
        self._structs = {}
        self._jobs = {}
        self._build_models()

    @Property
    def filename():
        doc = "The filename property."
        def fget(self):
            return self._filename
        def fset(self, value):
            self._filename = expandPath(value)
        return locals()

    @Property
    def path():
        doc = "The path property."
        def fget(self):
            return os.path.dirname(self.filename)
        return locals()

    @Property
    def code():
        doc = "The code property."
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
    def format():
        doc = "The format property."
        def fget(self):
            return self._format
        return locals()

    @Property
    def header():
        doc = "The header property."
        def fget(self):
            return self._header
        return locals()

    @Property
    def crd():
        doc = "The crd property."
        def fget(self):
            return self._crd
        return locals()

    @Property
    def footer():
        doc = "The footer property."
        def fget(self):
            return self._footer
        return locals()

##################
# Public Methods #
##################

    def keys(self):
        """
        Lists the keys to access all of the Mol objects the PDBFile object
        contains.
        """
        models = [ 'model%d' % i for i in sorted(self._models.keys()) ]
        structs = [ 'struct%d' % i for i in sorted(self._structs.keys()) ]
        jobs = [ 'job%d' % i for i in sorted(self._jobs.keys()) ]
        return models + structs + jobs

    def iter_models(self):
        """
        Iterate over the PDBFile's models. Models are Mol objects derived
        directly from the original pdb file.
        """
        return ( self['model%d' % i] for i in sorted(self._models.keys()) )

    def iter_structs(self):
        """
        Iterate over the PDBFile's structs. Structs are user specified Mol
        objects.
        """
        return ( self['struct%d' % i] for i in sorted(self._models.keys()) )

    def iter_jobs(self):
        """
        Iterate over the PDBFile's jobs. Jobs are Mol objects created as the
        result of a CHARMMing process.
        """
        return ( self['job%d' % i] for i in sorted(self._models.keys()) )

    def iter_all(self):
        """
        Iterate over all of the PDBFile's Mol objects, regardless of their
        origin.
        """
        return itertools.chain(self.iter_models(),self.iter_structs(),self.iter_jobs())

    def append_struct(self,input,structid,**kwargs):
        """
        This method should be used for user specified pdb files. 'Input' may be
        a file location or an iterable. This method should only be called by
        the appropriate wrapper method, and not by itself. All valid Mol kwargs
        are acceptable.
        """
        structid = int(structid)
        try:
            iterator = ( line.lower() for line in open(input) if line.lower().startswith(('atom','hetatm')) )
        except IOError:
            iterator = ( line.lower() for line in input if line.lower().startswith(('atom','hetatm')) )
        self._structs[structid] = Mol(iterable=iterator,**kwargs)

    def append_job(self,input,jobid,**kwargs):
        """
        This method should be used for CHARMMing generated pdb files. 'Input'
        may be a file location or an iterable. This method should only be
        called by the appropriate wrapper method, and not by itself. All valid
        Mol kwargs are acceptable.
        """
        jobid = int(jobid)
        try:
            iterator = ( line.lower() for line in open(input) if line.lower().startswith(('atom','hetatm')) )
        except IOError:
            iterator = ( line.lower() for line in input if line.lower().startswith(('atom','hetatm')) )
        self._jobs[jobid] = Mol(iterable=iterator,**kwargs)

    def get_warnings(self):
        """
        Return a list of warnings generated by the PDBFile object itself, and
        all of the Mol objects it contains.
        """
        tmp = self.warnings[:]
        for mol in self.iter_all():
            tmp.extend(mol.warnings)
        return tmp

###################
# Private Methods #
###################

    def _get_formatting(self,fileName):
        """
        Wrapper method to detect PDB text formatting, defaults to 'pdborg'
        """
        tmp = get_formatting(fileName)
        if tmp == 'unknown':
#WARN       # Throw a warning if formatting is guessed.
            self.warnings.append('Undetected pdb formatting')
            return self._defaultFormatting
        else:
            return tmp

    def _partitions(self):
        """
        Partition the file into _header/_crd/_footer sections.
        """
        filePointer = open(self.filename)
        # Populate header
        self._header = []
        tmp = None
        for line in cleanStrings(filePointer):
            if line.startswith(('atom','hetatm','model')):
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
            if not line.startswith(('atom','anisou','hetatm','model','ter','endmdl')):
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

    def _fix_chainids(self):
        """
        Detect and correct chainid mangling.
        """
        # Detect mangling
        def gen():
            for line in self.crd:
                if line.startswith(('atom','hetatm')):
                    tmp = Atom(line,format=self.format,autoFix=False)
                else:
                    continue
                yield tmp.chainid
        # Correct mangling
        if '' in set(gen()):
#WARN       # Throw a warning if chainids are mangled
            self.warnings.append('chainids mangled')
            print 'One or more chainids appear to be missing, attempting to autofix.\n'
            print 'Please verify the accuracy of your .pdb file upon completion.\n'
            chainNum = 0
            for i,line in enumerate(self.crd):
                if line.startswith(('atom','hetatm')):
                    tmp = Atom(line,format=self.format,autoFix=False)
                    # Blank chainid
                    if not tmp.chainid:
                        tmp.chainid = alphanum[chainNum]
                        self.crd[i] = tmp.Print(format=self.format).lower()
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

    def _build_models(self):
        """
        Parse the crd section, and load the coordinates into Mol objects, one
        per model section.
        """
        models = paragraphs(self.crd,splitter=['model'])
        for model in models:
            if model[0].startswith('model'):
                modelNum = int(model[0].split()[1])
            else:
                modelNum = 0
            iterator = ( Atom(line,format=self.format,index=i,autoFix=self._autoFix) for i,line in enumerate(model) if line.startswith(('atom','hetatm')) )
            self._models[modelNum] = Mol(iterable=iterator,name='model%d'%modelNum,code=self.code,autoFix=True)

###################
# Special Methods #
###################

    def __getitem__(self,key):
        if key in self.keys():
            try:
                if key.startswith('model'):
                    return self._models[int(key.split('model')[1])]
                elif key.startswith('struct'):
                    return self._structs[int(key.split('struct')[1])]
                elif key.startswith('job'):
                    return self._jobs[int(key.split('job')[1])]
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
                #try:
                #    return self._structs[key]
                #except IndexError: pass
                #try:
                #    return self._jobs[key]
                #except IndexError: pass
            except ValueError:
                raise KeyError('Unknown Key: %s' % key)

    def __len__(self):
        return len(self._models) + len(self._structs) + len(self._jobs)

    def __contains__(self,key):
        return key in self.keys()

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__,self.keys())
