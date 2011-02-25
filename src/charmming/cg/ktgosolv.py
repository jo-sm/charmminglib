"""
DOCME
"""
# fcp
# 02/24/2011


from charmming.cg.const import blnSolv_map, blnSolv_matrix, blnSolvSC_map
from charmming.cg.ktgo import KTGo


class KTGoSolv(KTGo):
    """
    DOCME
    """
    def __init__(self, iterable=None, **kwargs):
        super(KTGoSolv, self).__init__(iterable, **kwargs)

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
        self.solventMatrix = self.nScale * abs(blnSolv_matrix)
        super(KTGoSolv, self).write_prm(filename)

    def _rtf_mass(self):
        String = super(KTGoSolv, self)._rtf_mass()
        atomNum = self[-1].atomNum + 1
        String.append('MASS %-5d%-8s%10.6f' % (atomNum, 'p', 72))
        return String

    def _prm_nonbond(self):
        String = super(KTGoSolv, self)._prm_nonbond()
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        String.append('! Solvent Parameter goes here')
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        # PUT WATER PRM HERE
        return String

    def _prm_nbfix(self):
        String = super(KTGoSolv, self)._prm_nbfix()
        self.set_hbondDonorAcceptor()
        #
        String.append('! Solvent Parameters')
        String.append('%-8s%-8s%14.6f%12.6f' % ('p', 'p',
                        -self.get_solvParm('p', 'p'), 2.63779) )
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
                                            2.63779)
            String.append(tmp)
        String.append('! Solvent Czech Sum: %8.2f' % solvEnergySum)
        return String
