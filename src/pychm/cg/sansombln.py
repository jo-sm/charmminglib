"""
DOCME
"""
# Tim Miller
# 03/22/2011


from commands import getstatusoutput
from numpy import ndarray
from tempfile import NamedTemporaryFile
from pychm.cg.ktgosolv import KTGoSolv
from pychm.const.bio import aaVDW
from pychm.tools import Property, lowerKeys, modPi
from pychm.lib.pro import NoAlphaCarbonError
from pychm.lib.mol import Mol
from pychm.lib.basestruct import BaseStruct
from pychm.cg.const import bt_matrix, bt_map, kgs_matrix, kgs_map, \
        mj_matrix, mj_map, bln_matrix, blnSolv_map
from pychm.cg.cgpro import CGPro
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
        'mThetaCoil':-240.0,
        'hbondFC':2.390
        }

    _shortname = {
        'arg': 'r',
        'his': 'h',
        'hsd': 'h',
        'lys': 'k',
        'asp': 'd',
        'glu': 'e',
        'ser': 's',
        'thr': 't',
        'asn': 'n',
        'gln': 'q',
        'cys': 'c',
        'gly': 'g',
        'pro': 'p',
        'ala': 'a',
        'ile': 'i',
        'leu': 'l',
        'met': 'm',
        'phe': 'f',
        'trp': 'w',
        'tyr': 'y',
        'val': 'v'
    }

    def __init__(self, iterable=None, **kwargs):
        """
        DOCME
        """
        self.aTypeList = []
        super(SansomBLN, self).__init__(iterable, **kwargs)

        # run_stride will populate two arrays, self.structureOut and self.hbondOut
        # and atom.structure we need to use these to assign the backbone hydrogen bond types and the
        # helix, sheet, or coild types of the backbones and sidechains. We have to run this
        # AFTER gen_CGstruct (which will be called by the superclass constructor) because we need
        # the backbone and sidechain lists intact...
        self.run_stride()

        # take care of assigning backbone hydrogen bonding, possibly write a stream file with harmonic
        # restraints to hold secondary structural elements
        if kwargs.has_key('hbondstream'):
            hbfp = open(kwargs['hbondstream'], 'w')
            hbfp.write('* hydrogen bond restraints\n')
            hbfp.write('*\n\n')

        dnrlist = set()
        acclist = set()

        for line in self.hbondOut:
            if line.startswith('DNR'):
                # deal solely with acceptors for now since
                # this is the way Frank's code works
                continue

            # I am going to assume that the hydrogen bonds from STRIDE are correct, however
            # the mechanism of Bond & Sansom has slight differences (mostly because they
            # assume we can find the hydrogen coordinates). It might not be a bad idea to
            # investigate this.
            # NB, res2 and res1 are inverted because the acceptor is listed first, then the
            # donor.
            res2 = int(line[10:15])
            res1 = int(line[30:35])
            dist = float(line[42:45])

            dnrlist.add(res1)
            acclist.add(res2)
            if kwargs.has_key('hbondstream'):
                # fixme, figure out how to get the correct segment ID
                hbfp.write('resd kval %10.6f %s %d B %s %d B\n' % (self._parameters['hbondFC'],'a-bln',res1,'a-bln',res2))
        if kwargs.has_key('hbondstream'):
            hbfp.close()

        for atom in self:
            if atom.resid in dnrlist and atom.resid in acclist:
                atom.resName += 'b'
            elif atom.resid in dnrlist:
                atom.resName += 'd'
            elif atom.resid in acclist:
                atom.resName += 'a'

        # take care of assigning Helix and Sheet types
        for line in self.structureOut:
            resnum = int(line[10:15])
            struct = line[24]

            print "Got resid %d struct %s" % (resnum,struct)
            if struct == 'H' or struct == 'G':
                # helix -- Bond and Sansom explicitly state that they only
                # include alpha-helices (H), but I am including 310-helices
                # (G) as well. We need to check this, though...
                for atom in self:
                    if atom.resid == resnum:
                        atom.resName = 'h' + atom.resName[1:]
                
            elif struct == 'E':
                # beta sheet
                for atom in self:
                    if atom.resid == resnum:
                        atom.resName = 's' + atom.resName[1:]

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
            tmp.derivedResName = res.resName
            tmp.resName = 'b' + self._shortname[res.resName]
            tmp.segType = 'bln'
            print "DEBUG: assign %s -> %s" % (tmp.derivedResName,tmp.resName)
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
                    sc1.atomType = '   s'
                    sc1.derivedResName = sc1.resName
                    sc1.resName = 'bf' 
                    sc1.segType = 'bln'

                    sc2 = res.get_alphaCarbon()
                    sc2.atomNum = 2
                    sc2.cart = sc2l.com
                    sc2.atomType = '  s2'
                    sc2.derivedResName = tmp.resName
                    sc2.resName = 'bf'
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
                    sc1.resName = 'b'+ self._shortname[res.resName]
                    sc1.segType = 'bln'
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com
                    sc2.atomNum = 2   
                    sc2.atomType = '  s2' 
                    sc2.derivedResName = sc2.resName 
                    sc2.resName = 'b'+ self._shortname[res.resName]
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
                    sc1.resName = 'by'
                    sc1.segType = 'bln'
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com
                    sc2.atomNum = 2
                    sc2.atomType = '  s2'
                    sc2.derivedResName = sc2.resName
                    sc2.resName = 'by'
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
                    sc1.resName = 'bh'
                    sc1.segType = 'bln' 
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com 
                    sc2.atomNum = 2  
                    sc2.atomType = '  s2'
                    sc2.derivedResName = sc2.resName
                    sc2.resName = 'bh'
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
                    sc1.resName = 'bw'
                    sc1.segType = 'bln'
                    yield sc1

                    sc2 = res.get_alphaCarbon()
                    sc2.cart = sc2l.com
                    sc2.atomNum = 2    
                    sc2.atomType = '  s2'
                    sc2.derivedResName = sc2.resName
                    sc2.resName = 'bw'
                    sc2.segType = 'bln'
                    yield sc2
                    
            else:
                # two site residue
                tmp = res.get_goSC(**kwargs)
                tmp.atomNum = 1
                tmp.resName = 'b' + self._shortname[res.resName]
                tmp.segType = 'bln'
                tmp.atomType = '   s'
                yield tmp

     
    def write_rtf(self, filename=None):
        String = []
        String.append('* Topology file for three-site BLN model')
        String.extend(['*','31  1',''])
        String.extend(self._rtf_masses())

        # declare that we can bond backbone atoms
        String.extend(["", 'DECL +B', 'DECL -B', 'DECL #B', ""])
        String.append('! defaults, note there are no dihedrals in this model')
        String.append('DEFAULT FIRST NONE LAST NONE')
        String.append('AUTO ANGLES')
        String.append('')

        String.extend(self._rtf_residues())
        String.append('END\n')
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
        String.extend(self._prm_angle())
        String.append('')
        String.extend(self._prm_nonbond())
        String.append('')
        String.extend(self._prm_nbfix())
        String.append('')
        String.append('END\n')
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
            String.append('MASS %-5d%-8s%10.6f' % (i+1, atyp.upper(), 72.0))
            self.aTypeList.append(atyp)

        # now handle helix and sheet types
        cnt = i + 2
        for sfx in ['h','s']:
            for atyp in atypes:
                prmstr = atyp + sfx
                String.append('MASS %-5d%-8s%10.6f' % (cnt, prmstr.upper(), 72.0))
                self.aTypeList.append(prmstr)
                cnt += 1

        return String

    def _rtf_residues(self):
        String = []

        # prefixes: b = normal BLN (coil), h = helix, s = sheet
        prefix = ['b', 'h', 's']
        donacc = ['0', 'd', 'a', 'b']

        for p in prefix:
            for r in set(self._shortname.values()):
                for h in donacc:
                    bondstr = ''

                    if r in ['d', 'e']:
                        chrg = -1.0
                    elif r in ['k','r']:
                        chrg = 1.0
                    else:
                        chrg = 0.0
                    if h != '0':
                        String.append('RESI %s       %6.4f' % (p.upper() + r.upper() + h.upper(),chrg))
                    else:
                        String.append('RESI %s       %6.4f' % (p.upper() + r.upper(),chrg))
 
                    # backbone atom
                    # ToDo: check for backbone hydrogen bonding here...
                    if h == '0':
                        bbatom = 'n'
                    elif h == 'd':
                        bbatom = 'nd'
                    elif h == 'a':
                        bbatom = 'na'
                    elif h == 'b':
                        bbatom = 'nda'

                    if p == 'h':
                        bbatom += 'h'
                    elif p == 's':
                        bbatom += 's'

                    String.append('Group')
                    String.append('Atom  B %-4s 0.00' % bbatom)

                    # sidechain 1
                    if r != 'g':
                        String.append('Group')
                        if r in ['a','i','l','p','v','f','y','h','k','r']:
                            chrg = 0.0
                            type = 'c'
                        elif r in ['c', 'm']:
                            chrg = 0.0
                            type = 'n'
                        elif r in ['n', 'q']:
                           chrg = 0.0
                           type = 'nda'
                        elif r in ['s', 't']:
                           chrg = 0.0
                           type = 'p'
                        elif r == 'w':
                           chrg = 0.0
                           type = 'nd'
                        elif r in ['d', 'e']:
                           chrg = -1.0
                           type = 'qa'
                        else:
                           raise AssertionError("WTF")
                        bondstr += 'B S'
                        String.append('Atom  S %-4s %4.2f' % (type,chrg))

                    # sidechain 2
                    hassc2 = False
                    if r in ['r','k']:
                        hassc2 = True
                        chrg = 1.0
                        type = 'qd'
                    elif r == 'y':
                        hassc2 = True
                        chrg = 0.0
                        type = 'nd'
                    elif r == 'h':
                        hassc2 = True
                        chrg = 0.0
                        type = 'nda'
                    elif r == 'f':
                        hassc2 = True
                        chrg = 0.0
                        type = 'c'
                    elif r == 'w':
                        hassc2 = True
                        chrg = 0.0
                        type = 'c'

                    if hassc2:
                        bondstr += ' S S2'
                        String.append('Group')
                        String.append('Atom S2 %-4s %4.2f' % (type,chrg))

                    bondstr += ' B +B'

                    # bonds etc.
                    String.append('Bond %s' % bondstr)
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

        for i, t1 in enumerate(self.aTypeList):
            for t2 in self.aTypeList[i:]:
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

    def _prm_angle(self):
        String = ['THETA']
 
        # OK, this is an ugly hack, but it's late and
        # I'm having trouble thinking of an elegant
        # algorithm to remove combinations that only
        # differ in order
        doneatms = []

        for t1 in self.aTypeList:
            t1h = t1.endswith('h')
            t1s = t1.endswith('s')
            for t2 in self.aTypeList:
                t2h = t2.endswith('h')
                t2s = t2.endswith('s')
                for t3 in self.aTypeList:
                    if t3 in doneatms:
                        continue

                    t3h = t3.endswith('h')
                    t3s = t3.endswith('s')

                    if t1h and t2h and t3h:
                        kangle = self._parameters['kAngleHelix']
                        mtheta = self._parameters['mThetaHelix']
                    elif t1s and t2s and t3s:
                        kangle = self._parameters['kAngleSheet']
                        mtheta = self._parameters['mThetaSheet']
                    else:
                        kangle = self._parameters['kAngleCoil']
                        mtheta = self._parameters['mThetaCoil']

                    tmp = '%-8s%-8s%-8s%14.6f%12.6f' % \
                          (t1.upper(),t2.upper(),t3.upper(), kangle, mtheta)
                    String.append(tmp)

            doneatms.append(t1)
                    
        return String

    def _prm_nonbond(self):
        String = ['NBOND NBXMOD 4 ATOM CDIEL SHIFT VATOM VDISTANCE VSHIFT -']
        String.append(' CUTNB 16 CTOFNB 12 EPS 1.0 WMIN 1.5 E14FAC 0.7')
        for t in self.aTypeList:
            String.append('%-4s    0.0    1.2 4.700 ! Don\'t believe this parameter, it will be overwritten by NBFIX below' % t)
        return String

    def _prm_nbfix(self):
        String = ['NBFIX']
        for t1 in self.aTypeList:
            for t2 in self.aTypeList:

                # strip off helix and sheet designations, as they matter not here
                if t1.endswith('h') or t1.endswith('s'):
                    typ1 = t1[:-1]
                else:
                    typ1 = t1
                if t2.endswith('h') or t2.endswith('s'):
                    typ2 = t2[:-1]
                else:
                    typ2 = t2

                iidx = blnSolv_map[typ1]
                jidx = blnSolv_map[typ2]
                String.append('%-8s%-8s%14.6f%12.6f' % (t1, t2, bln_matrix[iidx][jidx]/1.69, 5.2755))

        return String
