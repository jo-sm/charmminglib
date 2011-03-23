"""
DOCME
"""
# Tim Miller
# 03/22/2011


from commands import getstatusoutput
from numpy import ndarray
from tempfile import NamedTemporaryFile
from charmming.const.bio import aaVDW
from charmming.tools import Property, lowerKeys, modPi
from charmming.lib.pro import NoAlphaCarbonError
from charmming.lib.mol import Mol
from charmming.cg.const import bt_matrix, bt_map, kgs_matrix, kgs_map, \
        mj_matrix, mj_map
from charmming.cg.cgpro import CGPro


class SansomBLN(Mol):
    """
    Class Attributes
        `_parameters`
    Properties
        `nScale`
    Public Methods
        `set_cgStruct`
    Private Methods
        ``
    """

    # Tim Miller: big fat note ... default parameters
    # here are in kcal/mol. Check with Frank to see if
    # this is OK.
    _parameters = {
        'nscale':1.00,
        'domainscale':1,
        'contactrad':4.5,
        'kBondHelix':2.988,
        'kBondSheet':2.988,
        'kBondCoil':2.988,
        'kAngleHelix':8.37,
        'kAngleSheet':8.37,
        'kAngleCoil':5.98,
        }

    def __init__(self, iterable=None, **kwargs):
        self.warnings = []
        super(SansomBLN, self).__init__(iterable=None, **kwargs)
        # kwargs aren't needed here, as this is a pretty simple model
        # 
        if iterable is None:
            raise TypeError('You must initialize a SansomBLN instance with an\
                            iterable of `Atom` like objects.')
        else:
            self.set_allAtoms(iterable, **kwargs)
            self.set_allHeavyAtoms(iterable, **kwargs)
            self.extend(self.gen_cgStruct(**kwargs))
            self.reindex_atomNum()

##############
# Properties #
##############

    @Property
    def contactSet():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return "bln"
        def fset(self, value):
            raise NotImplementedError("Only one BLN contact set is implemented currently")
        return locals()

    @Property
    def domainScale():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._parameters['domainscale']
        def fset(self, value):
            self._parameters['domainscale'] = value
        return locals()

    @Property
    def nScale():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._parameters['nscale']
        def fset(self, value):
            self._parameters['nscale'] = value
            self.calc_contactMatrix()
        return locals()

##################
# Public Methods #
##################

    def calc_contactMatrix(self):
        """
        """
        matrix = self.rawContactMatrix
        self.contactMatrix = self.nScale * matrix

    def gen_cgStruct(self, **kwargs):
        """
        """
        kwargs = lowerKeys(kwargs)
        calcmass = kwargs.get('calcmass', None)
        if calcmass is None:
            iterator = self.allHeavyAtoms.iter_res(restype=CGPro)
        else:
            iterator = self.allAtoms.iter_res(restype=CGPro)
        for res in iterator:
            yield res.get_goBB(**kwargs)
            try:
                yield res.get_goSC(**kwargs)
            except NoAlphaCarbonError:
                pass

    def get_NBFIX(self):
        tmp = []
        return tmp

    def get_nativeBBSC(self):
        bb = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro) ]
        sc = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
            if res.resName != 'gly' ]
        iterator = ( (res_i, res_j) for res_i in bb for res_j in sc
                    if abs(res_j.resid - res_i.resid) > 2 )
        #
        tmp = []
        contactRad = self._parameters['contactrad']
        for res_i, res_j in iterator:
            try:
                for atom_i in res_i.iter_bbAtoms():
                    for atom_j in res_j.iter_scAtoms():
                        if self.get_Rij(atom_i, atom_j) < contactRad:
                            tmp_bb = self.find(chainid=res_i.chainid, resid=res_i.resid)[0]
                            tmp_sc = self.find(chainid=res_j.chainid, resid=res_j.resid)[1]
                            tmp.append( (tmp_bb, tmp_sc) )
                            raise AssertionError
            except AssertionError:
                pass
        return tmp

    def get_nativeSCSC(self):
        sc1 = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
            if res.resName != 'gly' ]
        sc2 = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
            if res.resName != 'gly' ]
        iterator = ( (res_i, res_j) for res_i in sc1 for res_j in sc2
                    if res_j.resid - res_i.resid > 2 )
        tmp = []
        contactRad = self._parameters['contactrad']
        for res_i, res_j in iterator:
            try:
                for atom_i in res_i.iter_scAtoms():
                    for atom_j in res_j.iter_scAtoms():
                        if self.get_Rij(atom_i, atom_j) < contactRad:
                            tmp_sc1 = self.find(chainid=res_i.chainid, resid=res_i.resid)[1]
                            tmp_sc2 = self.find(chainid=res_j.chainid, resid=res_j.resid)[1]
                            tmp.append( (tmp_sc1, tmp_sc2) )
                            raise AssertionError
            except AssertionError:
                pass
        return tmp

    def get_parm(self, resName_i, resName_j):
        """
        """
        try:
            res_i = self.contactMap[resName_i.lower()]
        except KeyError:
            raise KeyError('The %s contact set has no parameters for the residue: %s'
                        % (self.contactSet, resName_i))
        try:
            res_j = self.contactMap[resName_j.lower()]
        except KeyError:
            raise KeyError('The %s contact set has no parameters for the residue: %s'
                        % (self.contactSet, resName_j))
        return self.contactMatrix[res_i][res_j]

    def get_Rij(self, atom_i, atom_j):
        try:
            return self._Rij[(atom_i.addr, atom_j.addr)]
        except KeyError:
            try:
                return self._Rij[(atom_j.addr, atom_i.addr)]
            except KeyError:
                distance = atom_i.calc_length(atom_j)
                self._Rij[(atom_i.addr, atom_j.addr)] = distance
                return distance

    def run_stride(self):
        # tmp file
        tmp = NamedTemporaryFile()
        tmpOut = []
        for atom in self.allHeavyAtoms:
            tmpOut.append(atom.Print(outformat='pdborg'))
        tmpOut = '\n'.join(tmpOut)
        tmp.file.write(tmpOut)
        tmp.file.flush()
        # execute stride
        strideOut = getstatusoutput('%s -h %s' %
                                    (self.strideBin, tmp.name))[1].split('\n')
        if strideOut[0].startswith('Error'):
            raise IOError('run_stride: %s' % strideOut[0])
        tmp.close()
        self.structureOut = [ line for line in strideOut
                            if line.startswith('ASG') ]
        self.hbondOut = [ line for line in strideOut
                        if line.startswith(('ACC', 'DNR')) and int(line[15:20])
                        < int(line[35:40]) ]
        # parse structure
        iterator = ( atom for atom in self if atom.atomType == 'b' )
        for i, atom in enumerate(iterator):
            atom.structure = self.structureOut[i][25:39].strip().lower()

    def set_allAtoms(self, iterable, **kwargs):
        iterable = ( atom for atom in iterable if atom.segType == 'pro' )
        self.allAtoms = Mol(iterable, **kwargs)

    def set_allHeavyAtoms(self, iterable, **kwargs):
        iterable = ( atom for atom in iterable if atom.segType == 'pro' and
                    atom.element != 'h' )
        self.allHeavyAtoms = Mol(iterable, **kwargs)
        self._Rij = {}

    def write_crd(self, filename=None):
        String = []
        String.append('*')
        String.append('%5d' % len(self))
        for atom in self:
            String.append(atom.Print(outformat='crd'))
        if filename is None:
            for line in String:
                print line
        else:
            String = '\n'.join(String)
            writeTo = open(filename, 'w')
            writeTo.write(String)
            writeTo.close()

    def write_pdb(self, filename=None):
        String = []
        for atom in self:
            String.append(atom.Print(outformat='charmm'))
        String.append('TER')
        #
        if filename is None:
            for line in String:
                print line
        else:
            String = '\n'.join(String)
            writeTo = open(filename, 'w')
            writeTo.write(String)
            writeTo.close()

    def write_rtf(self, filename=None):
        String = []
        String.extend(self._rtf_header())
        String.append('')
        String.extend(self._rtf_mass())
        String.append('')
        String.extend(self._rtf_declare())
        String.append('')
        String.extend(self._rtf_residue())
        String.append('')
        String.append('END')
        #
        if filename is None:
            for line in String:
                print line.upper()
        else:
            String = '\n'.join(String).upper()
            writeTo = open(filename, 'w')
            writeTo.write(String)
            writeTo.close()

    def write_prm(self, filename=None):
        # compute some things...
        self.run_stride()
        self.bbAtoms = []
        self.scAtoms = []
        for res in self.iter_res():
            self.bbAtoms.append(res[0])
            try:
                self.scAtoms.append(res[1])
            except IndexError:
                self.scAtoms.append(None)
        #
        String = []
        String.extend(self._prm_header())
        String.append('')
        String.extend(self._prm_bond())
        String.append('')
        String.extend(self._prm_angle())
        String.append('')
        String.extend(self._prm_nonbond())
        String.append('')
        String.extend(self._prm_nbfix())
        String.append('')
        String.append('END')
        #
        del self.bbAtoms
        del self.scAtoms
        #
        if filename is None:
            for line in String:
                print line.upper()
        else:
            String = '\n'.join(String).upper()
            writeTo = open(filename, 'w')
            writeTo.write(String)
            writeTo.close()

###################
# Private Methods #
###################

    def _rtf_header(self):
        String = []
        #
        String.append('* This CHARMM .rtf file describes a Bond-Sansom BLN model of %s' %
                    self.code)
        String.append('*')
        String.append('%5d%5d' % (20,1))
        return String

    def _rtf_mass(self):
        String = []
        #
        for cg in self:
            String.append('MASS %-5d%-8s%10.6f' % (cg.atomNum, cg.prmString,
                                                cg.mass))
        return String

    def _rtf_declare(self):
        String = []
        #
        String.append('DECL +B')
        String.append('DECL -B')
        String.append('DECL #B')
        String.append('DEFAULT FIRST NONE LAST NONE')
        return String

    def _rtf_residue(self):
        String = []
        #
        for res in self.iter_res():
            if len(res) == 1:
                String.append('RESI %-5s       0.0' % res.resName)
                String.append('GROUP')
                String.append('ATOM   B  %-5s  0.0' % res[0].prmString)
                String.append('BOND   B +B')
                String.append('ANGLE -B  B +B')
                String.append('DIHE  -B  B +B #B')
            else:
                String.append('RESI %-5s       0.0' % res.resName)
                String.append('GROUP')
                String.append('ATOM   B  %-5s  0.0' % res[0].prmString)
                String.append('ATOM   S  %-5s  0.0' % res[1].prmString)
                String.append('BOND   B  S  B +B')
                String.append('ANGLE -B  B  S  S  B +B -B  B +B')
                String.append('DIHE  -B  B +B #B')
                String.append('IMPH   B -B +B  S')
            String.append('')
        return String

    def _prm_header(self):
        String = []
        taco = self._parameters
        #
        String.append('* This CHARMM .param file describes a Bond-Sansom BLN model of %s' % self.code)
        String.append('* contactSet             = bln'
        String.append('* nscale                 = %6.3f' % taco['nscale'])
        String.append('* domainscale            = %6.3f' % taco['domainscale'])
        String.append('* kBondHelix             = %6.3f' % taco['kBondHelix'])
        String.append('* kBondSheet             = %6.3f' % taco['kBondSheet'])
        String.append('* kBondCoil              = %6.3f' % taco['kBondCoil'])
        String.append('* kAngleHelix            = %6.3f' % taco['kAngleHelix'])
        String.append('* kAngleSheet            = %6.3f' % taco['kAngleSheet'])
        String.append('* kAngleCoil             = %6.3f' % taco['kAngleCoil'])
        String.append('*')
        return String

    def _prm_bond(self):
        String = []
        kbond = self._parameters['kbond']
        #
        String.append('BOND')
        # bb(i)/bb(i+1)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb_0 = atom
                bb_1 = self.bbAtoms[i+1]
                tmp = '%-8s%-8s%14.6f%12.6f' % \
                        (bb_0.prmString, bb_1.prmString, kbond,
                        bb_0.calc_length(bb_1))
                String.append(tmp)
            except IndexError:
                pass
        # bb(i)/sc(i)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb = atom
                sc = self.scAtoms[i]
                tmp = '%-8s%-8s%14.6f%12.6f' % \
                        (bb.prmString, sc.prmString, kbond,
                        bb.calc_length(sc))
                String.append(tmp)
            except AttributeError:
                    pass
        return String

    def _prm_angle(self):
        String = []
        kangle = self._parameters['kangle']
        #
        String.append('ANGLE')
        # bb(i)/bb(i+1)/bb(i+2)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb_0 = atom
                bb_1 = self.bbAtoms[i+1]
                bb_2 = self.bbAtoms[i+2]
                tmp = '%-8s%-8s%-8s%14.6f%12.6f' % \
                        (bb_0.prmString, bb_1.prmString, bb_2.prmString,
                        kangle, bb_0.calc_angle(bb_1, bb_2))
                String.append(tmp)
            except IndexError:
                pass
        # sc(i)/bb(i)/bb(i+1)
        for i, atom in enumerate(self.scAtoms):
            try:
                sc = atom
                bb_0 = self.bbAtoms[i]
                bb_1 = self.bbAtoms[i+1]
                tmp = '%-8s%-8s%-8s%14.6f%12.6f' % \
                        (sc.prmString, bb_0.prmString, bb_1.prmString,
                        kangle, sc.calc_angle(bb_0, bb_1))
                String.append(tmp)
            except IndexError:
                pass
            except AttributeError:
                pass
        # bb(i)/bb(i+1)/sc(i+1)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb_0 = atom
                bb_1 = self.bbAtoms[i+1]
                sc_1 = self.scAtoms[i+1]
                tmp = '%-8s%-8s%-8s%14.6f%12.6f' % \
                        (bb_0.prmString, bb_1.prmString, sc_1.prmString,
                        kangle, bb_0.calc_angle(bb_1, sc_1))
                String.append(tmp)
            except IndexError:
                pass
            except AttributeError:
                pass
        return String

    def _prm_nonbond(self):
        String = []
        fnn = self._parameters['fnn']
        epsilonnn = self._parameters['epsilonnn']
        #
        String.append('NONBONDED ATOM CDIEL SWITCH VATOM VDISTANCE VSHIFT -')
        String.append('CUTNB 16.0 CTOFNB 12.0 CTONNB 9.0 EPS 15.0 WMIN 1.5 E14FAC 0.7')
        String.append('!atom')
        for atom in self:
            rMinDiv2 = 4.7
            tmp = '%-8s%5.1f%10.2e%12.6f' % (atom.prmString, 0, -epsilonnn, rMinDiv2)
            String.append(tmp)
        return String

    def _prm_nbfix(self):
        String = []
        nscale = self._parameters['nscale']
        domainscale = self._parameters['domainscale']
        bbscinteraction = self._parameters['bbscinteraction']
        #
        String.append('NBFIX')
        # hbonds
        String.append('! backbone hydrogen bonding')
        hbondEnergySum = 0.
        for hbond in self.get_hbond():
            hb_0, hb_1, energy, mult = hbond
            comment = ' '
            if hb_0.domain != hb_1.domain:
                comment += '! Interface between domains %d, %d ' % \
                        (hb_0.domain, hb_1.domain)
                energy *= domainscale
            if mult != 1:
                comment += '! H-Bond multiplicty is %d' % mult
                energy *= mult
            hbondEnergySum += energy
            tmp = '%-8s%-8s%14.6f%12.6f%s' % (hb_0.prmString, hb_1.prmString,
                                            energy, hb_0.calc_length(hb_1), comment)
            String.append(tmp)
        # scsc
        String.append('! native side-chain interactions')
        scscEnergySum = 0.
        for scsc in self.get_nativeSCSC():
            sc_0, sc_1 = scsc
            ljDepth = self.get_parm(sc_0.derivedResName, sc_1.derivedResName)
            comment = ' '
            if sc_0.domain != sc_1.domain:
                comment += '! Interface between domains %d, %d ' % \
                        (sc_0.domain, sc_1.domain)
                ljDepth *= domainscale/nscale
            scscEnergySum += ljDepth
            tmp = '%-8s%-8s%14.6f%12.6f%s' % (sc_0.prmString, sc_1.prmString,
                                            ljDepth, sc_0.calc_length(sc_1), comment)
            String.append(tmp)
        # bbsc
        String.append('! backbone side-chain interactions')
        bbscEnergySum = 0.
        for bbsc in self.get_nativeBBSC():
            bb, sc = bbsc
            ljDepth = bbscinteraction
            comment = ' '
            if bb.domain != sc.domain:
                comment += '! Interface between domains %d, %d ' % (bb.domain, sc.domain)
                ljDepth *= domainscale
            bbscEnergySum += ljDepth
            tmp = '%-8s%-8s%14.6f%12.6f%s' % (bb.prmString, sc.prmString,
                                            ljDepth, bb.calc_length(sc), comment)
            String.append(tmp)
        # czech sum
        String.append('! Czech Sum Info:%8.2f,%8.2f,%8.2f' %
                    (hbondEnergySum, bbscEnergySum, scscEnergySum) )
        #
        return String

################
# Over-written #
################

    def parse(self):
        """
        Not implemented.
        """
        raise NotImplementedError

