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
from charmming/lib/basestruct import BaseStruct
from charmming.cg.const import bt_matrix, bt_map, kgs_matrix, kgs_map, \
        mj_matrix, mj_map
from charmming.cg.cgpro import CGPro
import copy

class SansomBLN(KTGoSolv):

    def __init__(self):
        """
        DOCME
        """
        super(SansomBLN, self).__init__(iterable, **kwargs)

##################
# Public Methods #
##################


    def gen_cgStruct(self, **kwargs):
        """
        """
        # BTM: I don't understand how calcmass
        # works, so I am going to leave it out
        # for now.
        kwargs = lowerKeys(kwargs)
        
        for res in self.allAtoms.iter_res(restype=CGPro):
            tmp = res.get_goBB(**kwargs)
            tmp.atomType = 'n' # we can change this later for secondary struct
            tmp.segType = 'bln'
            yield tmp

            if res.resName = 'gly':
                # one site residue
                continue
            elif res.resName in \
                 ['phe', 'tyr', 'hsd', 'his', 'trp', 'lys', 'arg']:
                # three site residue
                
                # get a BaseStruct of all of the SC atoms
                # named taco in honor of Frank ;-).
                taco = BaseStruct([atom for atom in res if not atom.is_backbone()])

                if res.resName == 'phe':
                    # phenlyalanine gets two "C" type beads, one at the beta carbon
                    # position and one at the ring COM.
                    sc1  = copy.deepcopy(res.find(atomtype=' cb ')[0])
                    sc2l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['cg', 'cd1', 'hd1', 'cd2', 'hd2', 'ce1', \
                                       'he1','ce2', 'he2', 'cz', 'hz']
                    sc1.atomNum = 1
                    sc1.atomType = 'c'
                    sc1.derivedResName = tmp.resName
                    sc1.resName = '%s%d' % (self.chainid, self.resid)
                    sc1.segType = 'bln'

                    sc2 = res.get_alphaCarbon()
                    sc2.atomNum = 2
                    sc2.cart = sc2l.com
                    sc2.atomType = 'c'
                    sc2.derivedResName = tmp.resName
                    sc2.resName = '%s%d' % (self.chainid, self.resid)
                    sc2.segType = 'bln'

                    yield sc1
                    yield sc2

                if res.resName == 'arg':
                    sc2
            else:
                # two site residue
                tmp = res.get_goSC(**kwargs)
                tmp.atomNum = 1
                tmp.segType = 'bln'
                if res.resName in ['ala', 'ile', 'leu', 'pro', 'val']:
                    tmp.atomType = 'c'
                elif res.resName in ['cys', 'met']:
                    tmp.atomType = 'n0'
                elif res.resName in ['asn', 'gln']:
                    tmp.atomType = 'nda'
                elif res.resName in ['ser', 'thr']:
                    tmp.atomType = 'p'
                elif res.resName in ['asp', 'glu']:
                    tmp.atomType = 'qa'
                else:
                    raise AssertionError('What kind of residue is this??!')

                yield tmp

    def write_prm(self, filename=None):
        # compute some things...
        self.run_stride()
        self.bbAtoms = []
        self.sc1Atoms = []
        self.sc2Atoms = []
        for res in self.iter_res():
            self.bbAtoms.append(res[0])
            try:
                self.sc1Atoms.append(res[1])
            except IndexError:
                self.sc1Atoms.append(None)
            try:
                self.sc2Atoms.append(res[2])
            except IndexError:
                self.sc2Atoms.append(None)
        #
        String = []
        String.extend(self._prm_header())
        String.append('')
        String.extend(self._prm_bond())
        String.append('')
        #String.extend(self._prm_angle())
        #String.append('')
        #String.extend(self._prm_nonbond())
        #String.append('')
        #String.extend(self._prm_nbfix())
        #String.append('')
        String.append('END')
        #
        # BTM: why do we do this??!
        #del self.bbAtoms
        #del self.scAtoms
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
        String.append('* mThetaHelix            = %6.3f' % taco['mThetaHelix'])
        String.append('* mThetaSheet            = %6.3f' % taco['mThetaSheet'])
        String.append('* mThetaCoil             = %6.3f' % taco['mThetaCoil'])
        String.append('*')
        return String

    def _prm_bond(self):
        String = []
        #
        String.append('BOND')
        # bb(i)/bb(i+1)
        for i, atom in enumerate(self.bbAtoms):
            try:
                bb_0 = atom
                bb_1 = self.bbAtoms[i+1]

                if 'alphahelix' == bb_0.structure == bb_1.structure:
                    kbond = self._parameters.kBondHelix 
                elif '310helix' == bb_0.structure == bb_1.structure:
                    kbond = self._parameters.kBondHelix 
                elif 'betasheet' == bb_0.structure == bb_1.structure:
                    kbond = self._parameters.kBondSheet
                else:
                    kbond = self._parameters.kBondCoil

                tmp = '%-8s%-8s%14.6f%12.6f' % \
                        (bb_0.prmString, bb_1.prmString, kbond,
                        3.8)
                String.append(tmp)
            except IndexError:
                pass
        # bb(i)/sc(i)
        for i, atom in enumerate(self.bbAtoms):

            bb  = atom
            sc1 = self.sc1Atoms[i]
            sc2 = self.sc2Atoms[i]

            if 'alphahelix' == bb_0.structure:
                kbond = self._parameters.kBondHelix
            elif '310helix' == bb_0.structure:
                kbond = self._parameters.kBondHelix
            elif 'betasheet' == bb_0.structure:
                kbond = self._parameters.kBondSheet
            else:
                kbond = self._parameters.kBondCoil

            if sc1:
                tmp = '%-8s%-8s%14.6f%12.6f' % \
                        (bb.prmString, sc1.prmString, kbond,
                        3.8)
                String.append(tmp)

            if sc1 and sc2:
                tmp = '%-8s%-8s%14.6f%12.6f' % \
                        (sc1.prmString, sc2.prmString, kbond,
                        3.8)
                String.append(tmp)

        return Strin
