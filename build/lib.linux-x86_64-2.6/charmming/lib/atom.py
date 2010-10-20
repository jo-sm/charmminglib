from charmming.const.bio import atomMass,backbone,charmm2pdbAtomNames,good,pro,nuc
from charmming.const.etc import alphanum
from charmming.lib.baseatom import BaseAtom,AtomError
from charmming.tools import Property


class Atom(BaseAtom):
    """
    Private Attributes
        _format
        _resIndex
        _hash
    Parse Definition
        parse
    Properties
        addr
        resIndex
    Public Methods
        Print
        isBackbone
        isGood
    Private Methods
        _compliance_resName
        _compliance_atomType
        _sort
    """

    _tagDict = {'nuc':'atom','pro':'atom','good':'atom','bad':'hetatm','dna':'atom','rna':'atom'}

    def __init__(self,text=None,**kwargs):
        """
        kwargs:
            format      ['pdborg','charmm']
            cc          string  # comment character
            index       int
            autoFix     bool
        """
        # kwargs
        self._format = kwargs.pop('format','pdborg')
        # Setup
        super(Atom, self).__init__(text,**kwargs)
        if text is not None:
            # Parse
            self.parse()
            if self._element == 'x':
                self.fix_element()
        # Finally
        self._set_hash()

####################
# Parse Definition #
####################

    def parse(self):
        """
        Parses _text into instance variables.
        """
        if self._format == 'pdborg':
            self.atomNum    = self._text[6:11]
            self.atomType   = self._text[12:16]
            self.resName    = self._text[16:20]
            self.segType    = 'auto'
            self.chainid    = self._text[21]
            self.resid      = self._text[22:26]
            self.cart       = (self._text[30:38],self._text[38:46],self._text[46:54])
            self.weight     = self._text[55:60]
            self.bFactor    = self._text[61:66]
            self.element    = self._text[66:]
            self.resIndex   = self.resid
        elif self._format == 'charmm':
            self.atomNum    = self._text[6:11]
            self.atomType   = self._text[12:16]
            self.resName    = self._text[17:21]
            self.segType    = 'auto'
            self.chainid    = self._text[72:76]
            self.resid      = self._text[22:26]
            self.cart       = (self._text[30:38],self._text[38:46],self._text[46:54])
            self.weight     = self._text[55:60]
            self.bFactor    = self._text[61:66]
            self.element    = 'x'
            self.resIndex   = self.resid
        elif self._format == 'crd':
            pass
        elif self._format == 'amber':
            pass
        else:
            raise AtomError('parse: unknown _format %s' % self._format)
        #
        self._chainid0 = self.chainid
        self._segType0 = self.segType
        self._resid0 = self.resid
        self._atomNum0 = self.atomNum

##############
# Properties #
##############

    @Property
    def addr():
        doc = "The addr property."
        def fget(self):
            return '%s.%s.%d.%d' % (self.chainid,self.segType,self.resid,self.atomNum)
        return locals()

    @Property
    def addr0():
        doc = "The addr0 property."
        def fget(self):
            try:
                return '%s.%s.%d.%d' % (self.chainid0,self.segType0,self.resid0,self.atomNum0)
            except AttributeError:
                return 'Index: %d' % self._index
        return locals()

    @Property
    def atomNum():
        doc = "The atomNum property."
        def fget(self):
            return self._atomNum
        def fset(self, value):
            value = int(value)
            if value < 0:
                if self._autoFix:
                    value = 0
                else:
                    raise AtomError('atomNum %s: %d is less than 0' % (self.addr0,value))
            elif value > 10000:
                if self._autoFix:
                    value = 10000
                else:
                    raise AtomError('atomNum %s: %d is less than 0' % (self.addr0,value))
            self._atomNum = value
        return locals()

    @Property
    def atomNum0():
        doc = "The atomNum0 property."
        def fget(self):
            return self._atomNum0
        return locals()

    @Property
    def atomType():
        doc = "The atomType property."
        def fget(self):
            return self._atomType
        def fset(self, value):
            value = str(value).lower()
            if len(value) > 4:
                if self._autoFix:
                    value = value[:4]
                else:
                    raise AtomError('atomType %s: %s is longer than 4 characters' % (self.addr0,value))
            self._atomType = value
        return locals()

    @Property
    def bFactor():
        doc = "The bFactor property."
        def fget(self):
            return self._bFactor
        def fset(self, value):
            try:
                value = float(value)
            except ValueError:
                if self._autoFix:
                    value = 100.0
                else:
                    raise AtomError('bFactor %s: bFactor value is missing' % self.addr0)
            if value > 100.:
                if self._autoFix:
                    value = 100.
                else:
                    raise AtomError('bFactor %s: %5.2f is greater than 100.0' % (self.addr0,value))
            elif value < 0.:
                if self._autoFix:
                    value = 0.
                else:
                    raise AtomError('bFactor %s: %5.2f is less than 0.0' % (self.addr0,value))
            self._bFactor = value
        return locals()

    @Property
    def chainid():
        doc = "The chainid property."
        def fget(self):
            return self._chainid
        def fset(self, value):
            value = str(value).strip().lower()
            if len(value) > 1:
                if self._autoFix:
                    value = value[0]
                else:
                    raise AtomError('chainid %s: %s is longer than 1 character' % (self.addr0,value))
            self._chainid = value
        return locals()

    @Property
    def chainid0():
        doc = "The chainid0 property."
        def fget(self):
            return self._chainid0
        return locals()

    @Property
    def element():
        doc = "The element property."
        def fget(self):
            return self._element
        def fset(self, value):
            value = str(value).strip().lower()
            if value in atomMass.keys():
                self._element = value
                self._mass = atomMass[value]
                return
            if value in charmm2pdbAtomNames.keys():
                self._element = charmm2pdbAtomNames[value]
                self._mass = atomMass[self._element]
                return
            if self._autoFix:
                self._element = 'x'
                self._mass = 'badMass'
            else:
                raise AtomError('element %s: %s is an invalid element specification' % (self.addr0,value))
        def fdel(self):
            del self._element
        return locals()

    @Property
    def resid():
        doc = "The resid property."
        def fget(self):
            return self._resid
        def fset(self, value):
            value = int(value)
            if value < -1000:
                if self._autoFix:
                    value = -1000
                else:
                    raise AtomError('resid %s: %d is less than -1000' % (self.addr0,value))
            elif value > 10000:
                if self._autoFix:
                    value = 10000
                else:
                    raise AtomError('resid %s: %d is greater than 10000' % (self.addr0,value))
            self._resid = value
        return locals()

    @Property
    def resid0():
        doc = "The resid0 property."
        def fget(self):
            return self._resid0
        return locals()

    @Property
    def resIndex():
        doc = "The resIndex property."
        def fget(self):
            return self._resIndex
        def fset(self, value):
            value = int(value)
            if value < -1000:
                if self._autoFix:
                    value = -1000
                else:
                    raise AtomError('resIndex %s: %d is less than -1000' % (self.addr0,value))
            elif value > 10000:
                if self._autoFix:
                    value = 10000
                else:
                    raise AtomError('resIndex %s: %d is greater than 10000' % (self.addr0,value))
            self._resIndex = value
        return locals()

    @Property
    def resName():
        doc = "The resName property."
        def fget(self):
            return self._resName
        def fset(self, value):
            value = str(value).strip().lower()
            if len(value) > 4:
                if self._autoFix:
                    value = value[:4]
                else:
                   raise AtomError('resName %s: %s is longer than 4 characters' % (self.addr0,value))
            self._resName = value
        return locals()

    @Property
    def segid():
        doc = "The segid property."
        def fget(self):
            return '%s%s' % (self.chainid,self.segType)
        return locals()

    @Property
    def segType():
        doc = "The segType property."
        def fget(self):
            if self._segType == 'auto':
                if self.isPro():
                    return 'pro'
                elif self.isNuc():
                    return 'nuc'
                elif self.isGood():
                    return 'good'
                else:
                    return 'bad'
            else:
                return self._segType
        def fset(self, value):
            self._segType = str(value).strip().lower()
        return locals()

    @Property
    def segType0():
        doc = "The segType0 property."
        def fget(self):
            return self._segType0
        return locals()

    @Property
    def tag():
        doc = "The tag property."
        def fget(self):
            return self._tagDict[self.segType]
        return locals()

    @Property
    def weight():
        doc = "The weight property."
        def fget(self):
            return self._weight
        def fset(self, value):
            try:
                value = float(value)
            except ValueError:
                if self._autoFix:
                    value = 0.
                else:
                    raise AtomError('weight %s: weight value is missing' % self.addr0)
            if value > 1.:
                if self._autoFix:
                    value = 1.
                else:
                    raise AtomError('weight %s: %5.2f is greater than 1.0' % (self.addr0,value))
            elif value < 0.:
                if self._autoFix:
                    value = 0.
                else:
                    raise AtomError('weight %s: %5.2f is less than 0.0' % (self.addr0,value))
            self._weight = value
        return locals()

##################
# Public Methods #
##################

    def fix_element(self):
        if self._segType in set(['pro','nuc','rna','dna']):
            self.element = self._atomType.strip()[0]
        elif self._segType == 'good':
            if self._atomType in charmm2pdbAtomNames.keys():
                self.element = charmm2pdbAtomNames[self._atomType]
            else:
                self.element = self._atomType[:2].strip()
        else:
            self.element = self._atomType[:2].strip()

    def isBackbone(self):
        """
        atomType -> Bool
        """
        return self.atomType in backbone

    def isGood(self):
        """
        resName -> Bool
        """
        return self.resName in good

    def isPro(self):
        """
        resName -> Bool
        """
        return self.resName[-3:] in pro

    def isNuc(self):
        """
        resName -> Bool
        """
        return self.resName[-3:] in nuc

    def Print(self,**kwargs):
        """
        Returns a string representation of the Atom object. Defaults to the
        input formatting.
        kwargs:
            format          ['charmm','pdborg','debug','xdebug','crd','xcrd']
            old_chainid     [False,True]
            old_segType     [False,True]
            old_resid       [False,True]
            old_atomNum     [False,True]
        >>> thisAtom.Print(format='pdborg',old_resid=True)
        """

        # kwargs
        format      = kwargs.pop('format',self._format)
        old_chainid = kwargs.pop('old_chainid',False)
        old_segType = kwargs.pop('old_segType',False)
        old_resid   = kwargs.pop('old_resid',False)
        old_atomNum = kwargs.pop('old_atomNum',False)
        if kwargs.keys():
            raise TypeError('Unprocessed kwargs(%r)' % kwargs.keys())

        # Localize variables
        x,y,z   = self.cart
        if old_chainid:
            chainid = self.chainid0
        else:
            chainid = self.chainid
        if old_segType:
            segType = self.segType0
        else:
            segType = self.segType
        if old_resid:
            resid = self.resid0
        else:
            resid = self.resid
        if old_atomNum:
            atomNum = self.atomNum0
        else:
            atomNum = self.atomNum

        # Set formatting
        if format == 'pdborg':
            taco = '%-6s%5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n' % \
                (self.tag,atomNum,self.atomType,self.resName,chainid,resid,
                x,y,z,self.weight,self.bFactor,self.element)
        elif format == 'charmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n' % \
                (self.tag,atomNum,self.atomType,self.resName,resid,
                x,y,z,self.weight,self.bFactor,chainid)
        elif format == 'debug':
            taco = '%-6s%5i %4s %4s %1s%4i    %8.3f%8.3f%8.3f %3s' % \
                (segType,atomNum,self.atomType,self.resName,chainid,resid,
                x,y,z,self.element)
        elif format == 'xdebug':
            taco = '%15s > %15s: %5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s' % \
                (self.addr0,self.addr,atomNum,self.atomType,self.resName,chainid,resid,
                x,y,z,self.weight,self.bFactor,self.element)
        elif format == 'repr':
            taco = '%15s :: %4s%4s    %8.3f%8.3f%8.3f' % (self.addr,self.atomType,self.resName,x,y,z)
        elif format in ['crd','cor','card']:
            taco = '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f\n' % \
                (atomNum,self.resIndex,self.resName,self.atomType,
                x,y,z,chainid,resid,self.weight)
        elif format in ['xcrd','xcor','xcard']:
            taco = '%10i%10i  %-4s      %-4s    %20.10f%20.10f%20.10f  %-5s    %-4i    %20.10f\n' % \
                (atomNum,self.resIndex,self.resName,self.atomType,
                x,y,z,chainid + segType,resid,self.weight)
        else:
            raise AtomError('Print: unknown format %s' % format)
        return taco.upper()

###################
# Private Methods #
###################

    def _sort(self):
        """
        Scoring method that determines sorting order, for rich comparisons and
        container.sort() methods.

        chainid > segType > resid > atomNum
        """
        chainidScore = dict(((char,i) for i,char in enumerate(alphanum)))
        typeScore = {'pro':1,'dna':2,'rna':3,'nuc':4,'good':5,'bad':6}
        return chainidScore[self.chainid] * 1E10 + typeScore[self.segType] * 1E9 + self.resid * 1E5 + self.atomNum

# CHARMM Compliance
    def _compliance_resName(self):
        """
        Make resName CHARMM compliant.
        """
        if self.resName in ['hoh','tip3']:
           self.resName = 'tip3'
        elif self.resName == 'his':
           self.resName = 'hsd'
        elif self.resName in ['a','da']:
           self.resName = 'ade'
        elif self.resName in ['t','dt']:
           self.resName = 'thy'
        elif self.resName in ['c','dc']:
           self.resName = 'cyt'
        elif self.resName in ['g','dg']:
           self.resName = 'gua'
        elif self.resName in ['u','du']:
           self.resName = 'ura'
        elif (self.atomType == 'zn  ' and self.resName == 'zn'):
           self.resName = 'zn2'
        elif (self.atomType == 'na  ' and self.resName == 'na'):
           self.resName = 'sod'
        elif (self.atomType == 'cs  ' and self.resName == 'cs'):
           self.resName = 'ces'
        elif (self.atomType == 'cl  ' and self.resName == 'cl'):
           self.resName = 'cla'
        elif (self.atomType == 'ca  ' and self.resName == 'ca'):
           self.resName = 'cal'
        elif (self.atomType == ' k  ' and self.resName ==  'k'):
           self.resName = 'pot'

    def _compliance_atomType(self):
        """
        Make atomType CHARMM compliant.
        """
        if self.resName in ['hoh','tip3']:
           self.atomType = ' oh2'
        elif (self.resName == 'ile' and self.atomType == ' cd1'):
           self.atomType = ' cd '
        elif (self.atomType == 'na  ' and self.resName == 'na'):
           self.atomType = 'sod '
        elif (self.atomType == 'cs  ' and self.resName == 'cs'):
           self.atomType = 'ces '
        elif (self.atomType == 'cl  ' and self.resName == 'cl'):
           self.atomType = 'cla '
        elif (self.atomType == 'ca  ' and self.resName == 'ca'):
           self.atomType = 'cal '
        elif (self.atomType == ' k  ' and self.resName ==  'k'):
           self.atomType = 'pot '
