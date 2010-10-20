from itertools import tee
from charmming.lib.baseatom import BaseAtom
from charmming.tools import Property,expandPath

class StructError(Exception):
    """
    Exception to raise when errors occur involving the BaseStruct class.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class BaseStruct(list):
    """
    This is more or less a dummy class, until we figure out the exact details!
    Candidates for the BaseStruct to inherit from include...
        np.array
        OrderedSet
        OrderedDict
        blist
    """
    def __init__(self,iterable=None,**kwargs):
        """
        BaseStruct
        kwargs:
            name        string
            code        string  # pdbcode
            autoFix     [True,False]
        """
        super(BaseStruct, self).__init__()
        # kwargs
        self._name  = kwargs.pop('name','None')
        self._code  = kwargs.pop('code','None')
        autoFix     = kwargs.pop('autoFix',True)
        if kwargs.keys():
            raise TypeError('Unprocessed kwargs(%r)' % kwargs.keys())
        # Gatekeeper
        if iterable is None:
            list.__init__(self)
        else:
            if autoFix:
                iterable = ( key for key in iterable if isinstance(key,BaseAtom) )
            else:
                iterable,iterchk = tee(iterable,2)  # make a copy of the iterator for checking
                for key in iterchk:
                    if not isinstance(key,BaseAtom):
                        raise StructError('Only objects derived from BaseAtom class may be added to a BaseStruct object')
            list.__init__(self,iterable)

    @Property
    def name():
        doc = "The name property."
        def fget(self):
            return self._name
        def fset(self, value):
            self._name = value
        return locals()

    @Property
    def code():
        doc = "The pdb code of a struct."
        def fget(self):
            return self._code
        def fset(self, value):
            self._code = value
        return locals()

    def del_atoms(self,iterable):
        """
        Deletes, in place, one or more atoms as specified by iterable. The atoms themselves should
        go in the iterable which indicates the atoms to delete.
        >>> BaseStruct.del_atoms(BaseStruct[3:6])
        """
        for key in iterable:
            if not isinstance(key,BaseAtom):
                raise StructError('del_atoms: only objects from BaseAtom class may be deleted')
        del_atoms = set(iterable)
        del_indicies = [ i for i,atom in enumerate(self) if atom in del_atoms ]
        for index in reversed(del_indicies):
            self.pop(index)

    def write(self,fileName,**kwargs):
        """
        Writes a file containing the molecular information contained in the
        BaseStruct object.
        kwargs:
            format          ['charmm','pdborg','debug','xdebug','crd','xcrd']
            old_chainid     [False,True]
            old_segType     [False,True]
            old_resid       [False,True]
            old_atomNum     [False,True]
            ter             [False,True]
            end             True if format in ['pdborg','charmm']
        >>> thisSeg.write('~/charmming/1yjp/1yjp.pdb',Format='charmm',old_resid=True)
        """
        # kwargs
        Format  = kwargs.get('format','charmm')
        ter     = kwargs.pop('ter',None)
        end     = kwargs.pop('end',None)
        #
        fileName = expandPath(fileName)
        writeMe = []
        if Format in ['pdborg','charmm']:
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
                if end is None:
                    end = True
        elif Format in ['debug','xdebug']:
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
        elif Format in ['crd','xcrd']:
            writeMe.append('*\n')
            writeMe.append('   %d\n' % len(self))
            for atom in self:
                writeMe.append(atom.Print(**kwargs))
        # TER/END
        if ter:
            writeMe.append('TER\n')
        if end:
            writeMe.append('END\n')
        # Write file
        writeMe = ''.join(writeMe)
        writeTo = open(fileName,'w')
        writeTo.write(writeMe)
        writeTo.close()
