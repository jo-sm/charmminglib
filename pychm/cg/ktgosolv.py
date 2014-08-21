"""
DOCME
"""
# fcp
# 02/24/2011


from pychm.tools import Property
from pychm.cg.const import blnSolv_map, blnSolv_matrix, blnSolvSC_map
from pychm.cg.ktgo import KTGo


class KTGoSolv(KTGo):
    """
    DOCME
    """

    _parameters = dict(KTGo._parameters)
    _parameters.update({
        'solvscale':1.00
        })

    def __init__(self, iterable=None, **kwargs):
        super(KTGoSolv, self).__init__(iterable, **kwargs)

    @Property
    def solvScale():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._parameters['solvscale']
        def fset(self, value):
            self._parameters['solvscale'] = value
        return locals()

    def set_hbondDonorAcceptor(self):
        acceptors = []
        donors = []
        for line in self.hbondOut:
            if line.startswith('ACC'):
                acceptors.append(line[10:15])
                donors.append(line[30:35])
            elif line.startswith('DNR'):
                donors.append(line[10:15])
                acceptors.append(line[30:35])
        acceptors = set(map(int, acceptors))
        donors = set(map(int, donors))
        #
        for atom in self:
            if atom.resid in acceptors:
                atom.hbondAcceptor = True
            else:
                atom.hbondAcceptor = False
            if atom.resid in donors:
                atom.hbondDonor = True
            else:
                atom.hbondDonor = False

    def get_solvParm(self, resName_i, resName_j):
        """
        """
        try:
            res_i = self.solventMap[resName_i.lower()]
        except KeyError:
            raise KeyError('The solvent contact set has no parameters for the residue: %s'
                        % resName_i)
        try:
            res_j = self.solventMap[resName_j.lower()]
        except KeyError:
            raise KeyError('The solvent contact set has no parameters for the residue: %s'
                        % resName_j)
        return self.solventMatrix[res_i][res_j]

    def write_prm(self, filename=None):
        self.solventMap = blnSolv_map
        self.solventMatrix = self.solvScale * abs(blnSolv_matrix)
        super(KTGoSolv, self).write_prm(filename)

###################
# Private Methods #
###################

    def _rtf_mass(self):
        String = super(KTGoSolv, self)._rtf_mass()
        atomNum = self[-1].atomNum + 1
        String.append('MASS %-5d%-8s%10.6f' % (atomNum, 'p', 72))
        return String

    def _rtf_residue(self):
        String = super(KTGoSolv, self)._rtf_residue()
        String.append('RESI %-5s       0.0' % 'WAT')
        String.append('GROUP')
        String.append('ATOM   P  %-5s  0.0' % 'P')
        return String

    def _prm_header(self):
        String = []
        taco = self._parameters
        #
        String.append('* This CHARMM .param file describes a GoSolv model of %s' % self.code)
        String.append('* contactSet             = %s' % self.contactSet)
        String.append('* nscale                 = %6.3f' % taco['nscale'])
        String.append('* solvscale              = %6.3f' % taco['solvscale'])
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

    def _prm_nonbond(self):
        String = super(KTGoSolv, self)._prm_nonbond()
        fnn = self._parameters['fnn']
        epsilonnn = self._parameters['epsilonnn']
        rMinDiv2 = 4.7 * 2**(1/6.) * fnn
        String.append('! Solvent Parameter')
        String.append('%-8s%5.1f%10.2e%12.6f' % ('p', 0, -epsilonnn, rMinDiv2))
        return String

    def _prm_nbfix(self):
        String = super(KTGoSolv, self)._prm_nbfix()
        self.set_hbondDonorAcceptor()
        fnn = self._parameters['fnn']
        rMinDiv2 = 4.7 * 2**(1/6.) * fnn
        #
        String.append('! Solvent Parameters')
        String.append('%-8s%-8s%14.6f%12.6f' % ('p', 'p',
                        -self.get_solvParm('p', 'p'), rMinDiv2) )
        #
        solvEnergySum = 0
        for atom in self:
            if atom.atomType == 'b':
                key = 'n'
                if atom.hbondDonor:
                    key += 'd'
                if atom.hbondAcceptor:
                    key += 'a'
            elif atom.atomType == 's':
                key = blnSolvSC_map[atom.derivedResName]
            ljDepth = -self.get_solvParm('p', key)
            solvEnergySum += ljDepth
            tmp = '%-8s%-8s%14.6f%12.6f' % ('p', atom.prmString, ljDepth,
                                            rMinDiv2)
            String.append(tmp)
        String.append('! Solvent Czech Sum: %8.2f' % solvEnergySum)
        return String
