"""
DOCME
"""
# fcp
# 02/17/2011


from commands import getstatusoutput
from numpy import ndarray
from tempfile import NamedTemporaryFile
from pychm.const.bio import aaVDW
from pychm.tools import Property, lowerKeys, modPi
from pychm.lib.bond import Bond
from pychm.lib.pro import NoAlphaCarbonError
from pychm.lib.mol import Mol
from pychm.cg.const import bt_matrix, bt_map, kgs_matrix, kgs_map, \
        mj_matrix, mj_map
from pychm.cg.cgpro import CGPro


class KTGo(Mol):
    """
    Class Attributes
        `_parameters`
    Properties
        `nScale`
        `contactSet`
    Public Methods
        `set_cgStruct`
    Private Methods
        ``
    """

    _parameters = {
        'nscale':1.00,
        'domainscale':1,
        'fnn':0.65,
        'contactrad':4.5,
        'kbond':50,
        'kangle':30,
        'kdihedralalphahelix_1':0.30,
        'kdihedralalphahelix_3':0.15,
        'kdihedral310helix_1':0.30,
        'kdihedral310helix_3':0.15,
        'kdihedralnohelix_1':0.55,
        'kdihedralnohelix_3':0.275,
        'hbondenergyalphahelix':-0.25,
        'hbondenergy310helix':-0.25,
        'hbondenergynohelix':-0.50,
        'epsilonnn':1e-12,
        'bbscinteraction':-0.37,
        }

    def __init__(self, iterable=None, **kwargs):
        self.warnings = []
        super(KTGo, self).__init__(iterable=None, **kwargs)
        # kwargs
        self.strideBin = kwargs.get('stridebin', 'stride')
        self.contactSet = kwargs.get('contactset', 'kgs')
        #
        if iterable is None:
            raise TypeError('You must initialize a KTGo instance with an\
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
            return self._contactSet
        def fset(self, value):
            self._contactSet = value
            if value == 'kgs':
                self.rawContactMatrix = kgs_matrix
                self.contactMap = kgs_map
                self.contactOffset = 1.8
            elif value == 'bt':
                self.rawContactMatrix = bt_matrix
                self.contactMap = bt_map
                self.contactOffset = 0.6
            elif value == 'mj':
                self.rawContactMatrix = mj_matrix
                self.contactMap = mj_map
                self.contactOffset = 1.2
            elif value == 'other':
                raise NotImplementedError
            else:
                raise TypeError('Invalid `contactSet` specified specified')
            self.calc_contactMatrix()
        return locals()
    @Property
    def kBond():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._parameters['kbond']
        def fset(self, value):
            self._parameters['kbond'] = value
        return locals()

    @Property
    def kAngle():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._parameters['kangle']
        def fset(self, value):
            self._parameters['kangle'] = value
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
    def fnn():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._parameters['fnn']
        def fset(self, value):
            self._parameters['fnn'] = value
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
        self.contactMatrix = self.nScale * (matrix - self.contactOffset)

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

    #def get_hbond(self):
    #    # multiplicity
    #    hbondMult = {}
    #    iterator = ( (int(line[15:20]), int(line[35:40])) for line in self.hbondOut )
    #    for pair in iterator:
    #        if pair in hbondMult:
    #            hbondMult[pair] += 1
    #        else:
    #            hbondMult[pair] = 1
    #    # ugly hack to map from res index -> cgatom index
    #    keys = xrange(self[-1].resid)
    #    values = ( i for i, cg in enumerate(self) if cg.atomType == 'b' )
    #    bb2cg = dict(zip(keys, values))
    #    # hbond list
    #    energyAlhel = self._parameters['hbondenergyalphahelix']
    #    energy310hel = self._parameters['hbondenergy310helix']
    #    energyHel = (energyAlhel + energy310hel) / 2
    #    energyNohel = self._parameters['hbondenergynohelix']
    #    tmp = []
    #    for i, j in sorted(hbondMult.keys()):
    #        bb_i = self[bb2cg[i]]
    #        bb_j = self[bb2cg[j]]
    #        if 'alphahelix' == bb_i.structure == bb_j.structure:
    #            energy = energyAlhel
    #        elif '310helix' == bb_i.structure == bb_j.structure:
    #            energy = energy310hel
    #        elif ('alphahelix' == bb_i.structure and '310helix' ==
    #            bb_j.structure) or ('310helix' == bb_i.structure and
    #                                'alphahelix' == bb_j.structure):
    #            energy = energyHel
    #        else:
    #            energy = energyNohel
    #        tmp.append( (bb_i, bb_j, energy, hbondMult[(i, j)]) )
    #    return tmp

    def get_hbond(self):
        iterator = ( (int(line[10:15]), int(line[30:35])) for line in self.hbondOut )
        tmp = {}
        for i, j in iterator:
            cg_i = self.find(resid=i, atomtype='b')[0]
            cg_j = self.find(resid=j, atomtype='b')[0]
            bond = Bond(cg_i, cg_j)
            if bond.key in tmp.iterkeys():
                tmp[bond.key].order += 1
            else:
                bond.order = 1
                tmp[bond.key] = bond
        # energy
        energyAlhel = self._parameters['hbondenergyalphahelix']
        energy310hel = self._parameters['hbondenergy310helix']
        energyHel = (energyAlhel + energy310hel) / 2
        energyNohel = self._parameters['hbondenergynohelix']
        for bond in tmp.itervalues():
            if 'alphahelix' == bond.i.structure == bond.j.structure:
                bond.energy = energyAlhel
            elif '310helix' == bond.i.structure == bond.j.structure:
                bond.energy = energy310hel
            elif ('alphahelix' == bond.i.structure and '310helix' ==
                bond.j.structure) or ('310helix' == bond.i.structure and
                                    'alphahelix' == bond.j.structure):
                bond.energy = energyHel
            else:
                bond.energy = energyNohel
        return sorted(tmp.values())

    #def get_nativeBBSC(self):
    #    bb = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro) ]
    #    sc = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
    #        if res.resName != 'gly' ]
    #    iterator = ( (res_i, res_j) for res_i in bb for res_j in sc
    #                if abs(res_j.resid - res_i.resid) > 2 )
    #    #
    #    tmp = []
    #    contactRad = self._parameters['contactrad']
    #    get_Rij = self.get_Rij
    #    for res_i, res_j in iterator:
    #        try:
    #            for atom_i in res_i.iter_bbAtoms():
    #                for atom_j in res_j.iter_scAtoms():
    #                    if get_Rij(atom_i, atom_j) < contactRad:
    #                        tmp_bb = self.find(chainid=res_i.chainid, resid=res_i.resid)[0]
    #                        tmp_sc = self.find(chainid=res_j.chainid, resid=res_j.resid)[1]
    #                        tmp.append( (tmp_bb, tmp_sc) )
    #                        raise AssertionError
    #        except AssertionError:
    #            pass
    #    return tmp

    def get_nativeBBSC(self):
        bb = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro) ]
        sc = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
            if res.resName != 'gly' ]
        iterator = ( (res_i, res_j) for res_i in bb for res_j in sc
                    if abs(res_j.resid - res_i.resid) > 2 )
        #
        tmp = []
        contactRad = self._parameters['contactrad']
        get_Rij = self.get_Rij
        for res_i, res_j in iterator:
            try:
                for atom_i in res_i.iter_bbAtoms():
                    for atom_j in res_j.iter_scAtoms():
                        if get_Rij(atom_i, atom_j) < contactRad:
                            tmp_bb = self.find(chainid=res_i.chainid, resid=res_i.resid)[0]
                            tmp_sc = self.find(chainid=res_j.chainid, resid=res_j.resid)[1]
                            tmp.append( Bond(tmp_bb, tmp_sc) )
                            raise AssertionError
            except AssertionError:
                pass
        return sorted(tmp)

    #def get_nativeSCSC(self):
    #    sc1 = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
    #        if res.resName != 'gly' ]
    #    sc2 = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
    #        if res.resName != 'gly' ]
    #    iterator = ( (res_i, res_j) for res_i in sc1 for res_j in sc2
    #                if res_j.resid - res_i.resid > 2 )
    #    tmp = []
    #    contactRad = self._parameters['contactrad']
    #    get_Rij = self.get_Rij
    #    for res_i, res_j in iterator:
    #        try:
    #            for atom_i in res_i.iter_scAtoms():
    #                for atom_j in res_j.iter_scAtoms():
    #                    if get_Rij(atom_i, atom_j) < contactRad:
    #                        tmp_sc1 = self.find(chainid=res_i.chainid, resid=res_i.resid)[1]
    #                        tmp_sc2 = self.find(chainid=res_j.chainid, resid=res_j.resid)[1]
    #                        tmp.append( (tmp_sc1, tmp_sc2) )
    #                        raise AssertionError
    #        except AssertionError:
    #            pass
    #    return tmp

    def get_nativeSCSC(self):
        sc1 = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
            if res.resName != 'gly' ]
        sc2 = [ res for res in self.allHeavyAtoms.iter_res(restype=CGPro)
            if res.resName != 'gly' ]
        iterator = ( (res_i, res_j) for res_i in sc1 for res_j in sc2
                    if res_j.resid - res_i.resid > 2 )
        tmp = []
        contactRad = self._parameters['contactrad']
        get_Rij = self.get_Rij
        for res_i, res_j in iterator:
            try:
                for atom_i in res_i.iter_scAtoms():
                    for atom_j in res_j.iter_scAtoms():
                        if get_Rij(atom_i, atom_j) < contactRad:
                            tmp_sc1 = self.find(chainid=res_i.chainid, resid=res_i.resid)[1]
                            tmp_sc2 = self.find(chainid=res_j.chainid, resid=res_j.resid)[1]
                            tmp.append( Bond(tmp_sc1, tmp_sc2) )
                            raise AssertionError
            except AssertionError:
                pass
        return sorted(tmp)

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
        String.append('END\n')
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
        String.append('END\n')
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
        String.extend(self._prm_dihedral())
        String.append('')
        String.extend(self._prm_improper())
        String.append('')
        String.extend(self._prm_nonbond())
        String.append('')
        String.extend(self._prm_nbfix())
        String.append('')
        String.append('END\n')
        #
        del self.bbAtoms
        del self.scAtoms
        #
        if filename is None:
            for line in String:
                print line.upper()
        else:
            writeTo = open(filename, 'w')
            writeTo.write('\n'.join(String).upper())
            writeTo.close()

###################
# Private Methods #
###################

    def _rtf_header(self):
        String = []
        #
        String.append('* This CHARMM .rtf file describes a KTGo model of %s' %
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
        String.append('* This CHARMM .param file describes a Go model of %s' % self.code)
        String.append('* contactSet             = %s' % self.contactSet)
        String.append('* nscale                 = %6.3f' % taco['nscale'])
        String.append('* domainscale            = %6.3f' % taco['domainscale'])
        String.append('* fnn                    = %6.3f' % taco['fnn'])
        String.append('* contactrad             = %6.3f' % taco['contactrad'])
        String.append('* kbond                  = %6.3f' % taco['kbond'])
        String.append('* kangle                 = %6.3f' % taco['kangle'])
        String.append('* kdihedralalphahelix_1  = %6.3f' % taco['kdihedralalphahelix_1'])
        String.append('* kdihedralalphahelix_3  = %6.3f' % taco['kdihedralalphahelix_3'])
        String.append('* kdihedral310helix_1    = %6.3f' % taco['kdihedral310helix_1'])
        String.append('* kdihedral310helix_3    = %6.3f' % taco['kdihedral310helix_3'])
        String.append('* kdihedralnohelix_1     = %6.3f' % taco['kdihedralnohelix_1'])
        String.append('* kdihedralnohelix_3     = %6.3f' % taco['kdihedralnohelix_3'])
        String.append('* hBondenergyalphahelix  = %6.3f' % taco['hbondenergyalphahelix'])
        String.append('* hBondenergy310helix    = %6.3f' % taco['hbondenergy310helix'])
        String.append('* hBondenergynohelix     = %6.3f' % taco['hbondenergynohelix'])
        String.append('* epsilonnn              = %6.3e' % taco['epsilonnn'])
        String.append('* bbscinteraction        = %6.3f' % taco['bbscinteraction'])
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
                        kangle, bb_0.calc_angle(bb_1, bb_2, units='deg'))
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
                        kangle, sc.calc_angle(bb_0, bb_1, units='deg'))
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
                        kangle, bb_0.calc_angle(bb_1, sc_1, units='deg'))
                String.append(tmp)
            except IndexError:
                pass
            except AttributeError:
                pass
        return String

    def _prm_dihedral(self):
        String = []
        kdialhel_1 = self._parameters['kdihedralalphahelix_1']
        kdialhel_3 = self._parameters['kdihedralalphahelix_3']
        kdi310hel_1 = self._parameters['kdihedral310helix_1']
        kdi310hel_3 = self._parameters['kdihedral310helix_3']
        kdinohel_1 = self._parameters['kdihedralnohelix_1']
        kdinohel_3 = self._parameters['kdihedralnohelix_3']
        #
        kdihel_1 = (kdialhel_1 + kdi310hel_1) / 2
        kdihel_3 = (kdialhel_3 + kdi310hel_3) / 2
        #
        String.append('DIHEDRAL')
        String.append('! Backbone')
        # bb(i)/bb(i+1)/bb(i+2)/bb(i+3)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb_0 = atom
                bb_1 = self.bbAtoms[i+1]
                bb_2 = self.bbAtoms[i+2]
                bb_3 = self.bbAtoms[i+3]
                #
                dihedral = bb_0.calc_signedDihedral(bb_1, bb_2, bb_3, units='deg')
                #
                if 'alphahelix' == bb_1.structure == bb_2.structure:
                    k = [kdialhel_1, kdialhel_3]
                elif '310helix' == bb_1.structure == bb_2.structure:
                    k = [kdi310hel_1, kdi310hel_3]
                elif ('alphahelix' == bb_1.structure and '310helix' ==
                      bb_2.structure) or ('310helix' == bb_1.structure and
                                        'alphahelix' == bb_2.structure):
                    k = [kdihel_1, kdihel_3]
                else:
                    k = [kdinohel_1, kdinohel_3]
                for i, mult in enumerate([1, 3]):
                    delta = modPi( mult * dihedral - 180. )
                    tmp = '%-8s%-8s%-8s%-8s%14.6f%3d%12.6f' % \
                            (bb_0.prmString, bb_1.prmString, bb_2.prmString,
                            bb_3.prmString, k[i], mult, delta)
                    String.append(tmp)
            except IndexError:
                pass
        return String

    def _prm_improper(self):
        String = []
        #
        String.append('IMPROPER')
        String.append('! Sidechain')
        # bb(i+1)/bb(i)/bb(i+2)/sc(i+1)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb_1 = self.bbAtoms[i+1]
                bb_0 = atom
                bb_2 = self.bbAtoms[i+2]
                sc_1 = self.scAtoms[i+1]
                #
                dihedral = bb_1.calc_signedDihedral(bb_0, bb_2, sc_1, units='deg')
                #
                k = 20. * abs(self.contactMatrix.mean())
                mult = 1
                delta = dihedral + 180.
                tmp = '%-8s%-8s%-8s%-8s%14.6f%3d%12.6f' % \
                        (bb_1.prmString, bb_0.prmString, bb_2.prmString,
                        sc_1.prmString, k, mult, delta)
                String.append(tmp)
            except IndexError:
                pass
            except AttributeError:
                pass
        return String

    def _prm_nonbond(self):
        String = []
        epsilonnn = self._parameters['epsilonnn']
        #
        String.append('NONBONDED NBXMOD 4 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -')
        String.append('CUTNB 23 CTOFNB 20 CTONNB 18 EPS 1.0 WMIN 1.5 E14FAC 0.7')
        String.append('!atom')
        for atom in self:
            if atom.atomType == 'b':
                rMinDiv2 = 20
            elif atom.atomType == 's':
                rMinDiv2 = 10 * aaVDW[atom.derivedResName] * 2**(1/6.) * self.fnn
            else:
                raise AssertionError('How did I get here?')
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
        for hb in self.get_hbond():
            energy = hb.energy
            comment = ' '
            if hb.i.domain != hb.j.domain:
                comment += '! Interface between domains %d, %d ' % \
                        (hb.i.domain, hb.j.domain)
                energy *= domainscale
            if hb.order != 1:
                comment += '! H-Bond multiplicty is %d' % hb.order
                energy *= hb.order
            hbondEnergySum += energy
            tmp = '%-8s%-8s%14.6f%12.6f%s' % (hb.i.prmString, hb.j.prmString,
                                            energy, hb.length, comment)
            String.append(tmp)
        # scsc
        String.append('! native side-chain interactions')
        scscEnergySum = 0.
        for nat in self.get_nativeSCSC():
            ljDepth = self.get_parm(nat.i.derivedResName, nat.j.derivedResName)
            comment = ' '
            if nat.i.domain != nat.j.domain:
                comment += '! Interface between domains %d, %d ' % \
                        (nat.i.domain, nat.j.domain)
                ljDepth *= domainscale/nscale
            scscEnergySum += ljDepth
            tmp = '%-8s%-8s%14.6f%12.6f%s' % (nat.i.prmString, nat.j.prmString,
                                            ljDepth, nat.length, comment)
            String.append(tmp)
        # bbsc
        String.append('! backbone side-chain interactions')
        bbscEnergySum = 0.
        for nat in self.get_nativeBBSC():
            ljDepth = bbscinteraction
            comment = ' '
            if nat.i.domain != nat.j.domain:
                comment += '! Interface between domains %d, %d ' % \
                        (nat.i.domain, nat.j.domain)
                ljDepth *= domainscale
            bbscEnergySum += ljDepth
            tmp = '%-8s%-8s%14.6f%12.6f%s' % (nat.i.prmString, nat.j.prmString,
                                            ljDepth, nat.length, comment)
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

