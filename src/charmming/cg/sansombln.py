"""
DOCME
"""
# Tim Miller
# 03/22/2011


from commands import getstatusoutput
from numpy import ndarray
from tempfile import NamedTemporaryFile
from charmming.cg.ktgosolv import KTGoSolv
from charmming.const.bio import aaVDW
from charmming.tools import Property, lowerKeys, modPi
from charmming.lib.pro import NoAlphaCarbonError
from charmming.lib.mol import Mol
from charmming.lib.basestruct import BaseStruct
from charmming.cg.const import bt_matrix, bt_map, kgs_matrix, kgs_map, \
        mj_matrix, mj_map
from charmming.cg.cgpro import CGPro
import copy

class SansomBLN(KTGoSolv):

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
        'kBondInternal':2.988,
        'kAngleHelix':8.37,
        'kAngleSheet':8.37,
        'kAngleCoil':5.98,
        'mThetaHelix':-270.0,
        'mThetaSheet':-230.0,
        'mThetaCoil':-240.0
        }


    def __init__(self, iterable=None, **kwargs):
        """
        DOCME
        """
        self.aTypeList = []
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
            tmp.atomType = '   b' # we can change this later for secondary struct
            tmp.resName = 'b' + res.resName
            tmp.segType = 'bln'
            yield tmp

            if res.resName == 'gly':
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
                                       'he1','ce2', 'he2', 'cz', 'hz']])
                    sc1.atomNum = 1
                    sc1.atomType = '   c'
                    sc1.derivedResName = sc1.resName
                    sc1.resName = 'bphe' 
                    sc1.segType = 'bln'

                    sc2 = res.get_alphaCarbon()
                    sc2.atomNum = 2
                    sc2.cart = sc2l.com
                    sc2.atomType = '   c'
                    sc2.derivedResName = tmp.resName
                    sc2.resName = 'bphe'
                    sc2.segType = 'bln'

                    yield sc1
                    yield sc2

                elif res.resName == 'arg' or res.resName == 'lys':
                    # Arginine and Lysine gets a C + Qd, C at the CB, CG, CD COM
                    # and Qd at the COM of the rest of the SC heavy atoms.
                    sc1l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['cb', 'cg', 'cd']])
                    sc2l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['ce', 'nz', 'cz', 'nh1', 'nh2']])

                    sc1 = res.get_alphaCarbon()
                    sc1.cart = sc1l.com
                    sc1.atomNum = 1
                    sc1.atomType = '   s'
                    sc1.derivedResName = sc1.resName
                    sc1.resName = 'b'+ res.resName
                    sc1.segType = 'bln'
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com
                    sc2.atomNum = 2   
                    sc2.atomType = '  s2' 
                    sc2.derivedResName = sc2.resName 
                    sc2.resName = 'b'+ res.resName
                    sc2.segType = 'bln'
                    yield sc2

                elif res.resName == 'tyr':
                    # tyrosine gets C + Nd, C for CB, CG, and half the ring and Nd
                    # for the second half of the ring + the hydroxyl group
                    sc1l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['cb', 'cg', 'cd1', 'cd2']])
                    sc2l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['ce1', 'ce2', 'cz', 'oh']])

                    sc1 = res.get_alphaCarbon()
                    sc1.cart = sc1l.com
                    sc1.atomNum = 1   
                    sc1.atomType = '   s' 
                    sc1.derivedResName = sc1.resName 
                    sc1.resName = 'btyr'
                    sc1.segType = 'bln'
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com
                    sc2.atomNum = 2
                    sc2.atomType = '  s2'
                    sc2.derivedResName = sc2.resName
                    sc2.resName = 'btyr'
                    sc2.segType = 'bln'
                    yield sc2

                elif res.resName == 'his' or res.resName == 'hsd':
                    # histadine is C + Nda, note that we only deal with
                    # neutral histadine here -- charged histadine would
                    # have a Qda bead.
                    sc1l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['cb', 'cg']])
                    sc2l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['nd1', 'cd2', 'ce1', 'ne2']])

                    sc1 = res.get_alphaCarbon()
                    sc1.cart = sc1l.com
                    sc1.atomNum = 1
                    sc1.atomType = '   s'
                    sc1.derivedResName = sc1.resName
                    sc1.resName = 'bhsd'
                    sc1.segType = 'bln' 
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com 
                    sc2.atomNum = 2  
                    sc2.atomType = '  s2'
                    sc2.derivedResName = sc2.resName
                    sc2.resName = 'bhsd'
                    sc2.segType = 'bln' 
                    yield sc2

                elif res.resName == 'trp':
                    # tryptophan is similar to histadine, but with a C
                    # bead on its end
                    sc1l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['cb', 'cg', 'cd1', 'ne1']])
                    sc2l = BaseStruct([a for a in taco if a.atomType.strip() in \
                                      ['cd2', 'ce2', 'cz2', 'ce3', 'cz3', 'ch2']])

                    sc1 = res.get_alphaCarbon()
                    sc1.cart = sc1l.com
                    sc1.atomNum = 1
                    sc1.atomType = '   b'
                    sc1.derivedResName = sc1.resName
                    sc1.resName = 'bhis'
                    sc1.segType = 'bln'
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com
                    sc2.atomNum = 2    
                    sc2.atomType = '  s2'
                    sc2.derivedResName = sc2.resName
                    sc2.resName = 'bhis'
                    sc2.segType = 'bln'
                    yield sc2
                    
            else:
                # two site residue
                tmp = res.get_goSC(**kwargs)
                tmp.atomNum = 1
                tmp.resName = 'b' + res.resName
                tmp.segType = 'bln'
                tmp.atomType = '   s'
                yield tmp

    def write_rtf(self, filename=None):
        String = []
        String.append('* Topology file for three-site BLN model')
        String.extend(['*',''])
        String.extend(self._rtf_masses())

        # declare that we can bond backbone atoms
        String.extend(["", 'DECL +B', 'DECL -B', 'DECL #B', ""])
        String.append('! defaults, note there is no auto-generation')
        String.append('DEFAULT FIRST NONE LAST NONE')
        String.append('')

        String.extend(self._rtf_residues())
        String.append('END')
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

    def _rtf_masses(self):
        String = []

        # Loop over all of the standard atom types
        # NB: the Marrink-style BLN model sets all particle masses to 72,
        # so mass will not be conserved.
        atypes = ['n', 'nd', 'na', 'nda', 'c', 'p', 'q', 'qd', 'qa', 'qda']
        for i, atyp in enumerate(atypes):
            String.append('MASS %-5d%-8s%10.6f' % (i, atyp.upper(), 72.0))
            self.aTypeList.append(atyp)

        # now handle helix and sheet types
        cnt = i + 1
        for sfx in ['h','s']:
            for atyp in atypes:
                prmstr = atyp + sfx
                String.append('MASS %-5d%-8s%10.6f' % (cnt, prmstr.upper(), 72.0))
                self.aTypeList.append(prmstr)
                cnt += 1

        return String

    def _rtf_residues(self):
        String = []

        aares = ['ala','ile','leu','pro','val','phe','cys','met','asn','gln', \
                 'ser','thr','tyr','hsd','trp','asp','glu','lys','arg']

        # prefixes: b = normal BLN (coil), h = helix, s = sheet
        prefix = ['b', 'h', 's']

        for p in prefix:
            for r in aares:
                String.append('RESI %s' % p.upper() + r.upper())
 
                # backbone atom
                if p == 'h':
                   pass
                elif p == 's':
                   pass
                else:
                   pass
                String.append('')

        return String

    def _prm_header(self):
        String = []
        taco = self._parameters
        #
        String.append('* This CHARMM .param file describes a Bond-Sansom BLN model of %s' % self.code)
        String.append('* contactSet             = bln')
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

        for t1 in self.aTypeList:
            for t2 in self.aTypeList:
                if t1.endswith('h') and t2.endswith('h'):
                    kbond = self._parameters['kBondHelix']
                elif t1.endswith('s') and t2.endswith('s'):
                    kbond = self._parameters['kBondSheet']
                else:
                    kbond = self._parameters['kBondCoil']

                tmp = '%-8s%-8s%14.6f%12.6f' % \
                      (t1.upper(), t2.upper(), kbond, 3.8)
        String.append(tmp)

        return String
