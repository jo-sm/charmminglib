#!/usr/bin/env python

from BaseStruct import BaseStruct
from Mol import Mol

class ktgo(Mol):
    def __init__(self,iterable=None,autoFix=True):
        if iterable:
            iterable = ( atom for atom in iterable if atom._segType == 'pro' )
        Mol.__init__(self,iterable,autoFix)

    def set_backboneAtoms(self):
        self._backboneAtoms = BaseStruct( ( atom for atom in self if atom.isBackbone() ) )

    def set_sidechainAtoms(self):
        self._sidechainAtoms = BaseStruct( ( atom for atom in self if not atom.isBackbone() ) )

    def set_cgStructure(self):
        """
        Reductively map the all-atom structure in the the 2-center cg structure.
        """
        tmp = BaseStruct()
        for res in self.iter_res():
            tmp.add(res.get_alphaCarbon)
            try:
                tmp.add(res.get_sc())
            except ProError: pass
        self._cgStructure = tmp

    def set_domains(self,arg=None):
        """
        Sets _chainid to indicate domain membership
        """
        if type(arg) == NoneType:
            return
        elif type(arg) == str:
            def fixString(String):
                if String.startswith(':'):
                    String = '0' + String
                else:
                    String = '0:' + String
                if String.endswith(':'):
                    String = String + '10000'
                else:
                    String = String + ':10000'
                return String

            def string2list(String):
                String = String.split(':')
                return map(int,String)

            def get_chainid(arg,arg2):
                for i,domain in enumerate(arg2):
                    if arg <= domain:
                        return i
            for cgAtom in self._cgStructure:
                if last_resid < cgAtom._resid <= next_resid:
                    cgAtom._chainid = get_chainid()
                elif
            def get_chainid(iterate=False):
                domainid = iter(Aux.alphanum)
                while 1:
                    if iterate:
                        try:
                            thisDomainid = domainid.next()
                        except StopIteration: pass
                    yield thisDomainid
            domains = arg.split(':')
            domains = iter(map(int,domains))
            nextDomain = 0
            for atom in self:
                if atom._resid <= nextDomain:
                    atom._chainid = get_chainid()
                else:
                    try:
                        nextDomain = domains.next()
                    except StopIteration: pass
                    atom._chainid = get_chainid(iterate=True)
            pass
        elif type(arg) == dict:
            pass
        elif type(arg) == file: # TODO Fix this
            pass
        # Defaults to different seg gets different domain
###

# ( 3 ) Assign each residue to a domain
#    (a) Assign
if domainString:    # When explicitly specified
    domains = domainString.split(':')
    domains = map(int,domains)
    assert len(resDict) == domains[-1]
    thisDomain = 0
    for i,res in enumerate(resList):
        if i+1 <= domains[thisDomain]:
            pass
        else:
            thisDomain += 1
        resDict[res].domainid       = thisDomain
        cgbbDict[res].domainid      = thisDomain
        try:
            cgscDict[res].domainid  = thisDomain
        except KeyError: pass
else:               # Default: each segment gets it's own domain
    alpha2int   = dict( [ (char,i+1) for i,char in enumerate(string.uppercase) ] )
    for res in resList:
        resDict[res].domainid       = alpha2int[res[:1]]
        cgbbDict[res].domainid        = alpha2int[res[:1]]
        try:
            cgscDict[res].domainid  = alpha2int[res[:1]]
        except KeyError: pass
#   (b) Czech that all residues have a domain assignment
for res in resList:
    try:
        dummy = resDict[res].domainid
        dummy = cgbbDict[res].domainid
        try:
            dummy = cgscDict[res].domainid
        except KeyError: pass
    except AttributeError:
        print 'domains: residue %s has no domain assignment' % res
        raise GhastlyDeath

# ( 4 ) Load contact potential
#   (a) Parse Specified Contact File
resNameMap      = {}
contactMatrix   = []
for line in Etc.stripBare(open(contactFileName)):
    if line.startswith('AA'):
        buffer = line.split(':')[1].split()
        resNameMap = dict( zip(buffer,range(len(buffer))) )
    else:
        scaledContacts = [ nScale * abs(float(contact) - contactOffset) for contact in line.split() ]
        contactMatrix.append(scaledContacts)
#   (b) Calculate average interaction
avgEpsilon = sum(Etc.flatten(contactMatrix))/len(Etc.flatten(contactMatrix))
#   (c) Create ResName -> ContactEnergy mapper
def get_contact_energy(resName1,resName2):
    res1 = resNameMap[resName1]
    res2 = resNameMap[resName2]
    try:
        return float(contactMatrix[res1][res2])
    except IndexError:
        return float(contactMatrix[res2][res1])

# ( 5 ) Define more efficient Rij machinery
Rij = {}
def get_Rij(atom_i,atom_j):
    i = atom_i.segid + str(atom_i.atomNumber)
    j = atom_j.segid + str(atom_j.atomNumber)
    try:
        return Rij[(i,j)]
    except KeyError:
        Rij[(i,j)] = atom_i.bond_length(atom_j)
        return Rij[(i,j)]

# ( 6 ) Determine Native Contacts
#   (a) SC/SC Contacts
lowerOffAlmostDiagonal = ( (resCode_i,resCode_j) for i,resCode_i in enumerate(resList) for j,resCode_j in enumerate(resList) if j - i > 2 and resDict[resCode_i].resName != 'GLY' and resDict[resCode_j].resName != 'GLY' )
scContact = []
for resCode_i,resCode_j in lowerOffAlmostDiagonal:
    contact = False
    try:
        for atom_i in scDict[resCode_i]:
            for atom_j in scDict[resCode_j]:
                if get_Rij(atom_i,atom_j) < contactRad:
                    scContact.append( (cgscDict[resCode_i],cgscDict[resCode_j]) )
                    raise AssertionError
    except AssertionError: pass

#   (b) BB/SC Contacts
lowerOffAlmostDiagonal = ( (resCode_i,resCode_j) for i,resCode_i in enumerate(resList) for j,resCode_j in enumerate(resList) if abs(j - i) > 2 and resDict[resCode_j].resName != 'GLY' )
bbscContact = []
for resCode_i,resCode_j in lowerOffAlmostDiagonal:
    contact = False
    try:
        for atom_i in bbDict[resCode_i]:
            for atom_j in scDict[resCode_j]:
                if get_Rij(atom_i,atom_j) < contactRad:
                    bbscContact.append( (cgbbDict[resCode_i],cgscDict[resCode_j]) )
                    raise AssertionError
    except AssertionError:
        pass

# ( 7 ) STRIDE
strideOutput = commands.getstatusoutput('%s -h %s' % (stridePath,inputPath+inputFile) )[1].split('\n')
if strideOutput[0].startswith('Error'): raise IOError('stride: %s' % strideOutput[0])
#   (a) Get helical secondary structure
strideHelixOutput = [ line for line in strideOutput if line.startswith('ASG') ]
for i,line in enumerate(strideHelixOutput):
    if line[34:39] == 'Helix':
        cgbbDict[resList[i]].helix = True
    else:
        cgbbDict[resList[i]].helix = False

#   (b) Get hydrogen bonds
strideHBondOutput = [ line for line in strideOutput if line[:3] in ['ACC','DNR'] and int(line[15:20]) < int(line[35:40]) ]
#       (i) Multiplicity, ie anti-parallel beta sheets will have 'double' hbonds
hBondMult = {}
for line in strideHBondOutput:
    res_i = int(line[15:20])
    res_j = int(line[35:40])
    if (res_i,res_j) in hBondMult:
        hBondMult[(res_i,res_j)] += 1
    else:
        hBondMult[(res_i,res_j)] = 1
#       (ii) Build hbond list, 'double' bonds, get double energy
orderDict = ( (i,j) for i in range(len(resList)) for j in range(len(resList)) if i < j )
bbHBond = []
for (i,j) in orderDict:
    try:
        bb_i = cgbbDict[resList[i]]
        bb_j = cgbbDict[resList[j]]
        if bb_i.helix and bb_j.helix:
            hBondEnergy = hBondEnergyHelix
        else:
            hBondEnergy = hBondEnergyNoHelix
        bbHBond.append( (bb_i,bb_j,hBondEnergy,hBondMult[(i,j)]) )
    except KeyError:
        pass

# ( 8 ) Reindex CG.atomNumber
for i,cg in enumerate(cgList):
    cg.atomNumber = i+1

# ( 9 ) Write .crd file
outputName = 'cg_' + pdbName + '.crd'
writeTo = open(outputPath + outputName,'w')
writeTo.write('*\n')
writeTo.write('%5d\n' % len(cgList) )
for line in cgList:
    writeTo.write(line.Print('cgcrd'))
writeTo.close()

# (10 ) Write a .seq file
outputName = 'cg_' + pdbName + '.seq'
writeTo = open(outputPath + outputName,'w')
writeTo.write(' '.join(resList))
writeTo.close()

# (11 ) Write .pdb file
outputName = 'cg_' + pdbName + '.pdb'
writeTo = open(outputPath + outputName,'w')
for line in cgList:
    writeTo.write(line.Print('cgcharmm'))
writeTo.write('TER\n')
writeTo.close()

# (12 ) Write .rtf file
outputName = 'cg_' + pdbName + '.rtf'
writeTo = open(outputPath + outputName,'w')
writeTo.write('* This CHARMM .rtf file describes a Go model of %s\n' % pdbName)
writeTo.write('*\n')
writeTo.write('%5d%5d\n' % (20,1) )
writeTo.write('\n')

# Mass Section
for i,cg in enumerate(cgList):
    stringBuffer = 'MASS %-5d%-8s%10.6f\n' % (i+1,cg.resCode+cg.atomType.strip(),cg.mass)
    writeTo.write(stringBuffer)
writeTo.write('\n')

# Declare statements & Defaults
writeTo.write('DECL +B\n')
writeTo.write('DECL -B\n')
writeTo.write('DECL #B\n')
writeTo.write('DEFAULT FIRST NONE LAST NONE\n')
writeTo.write('\n')

# Residue Topology Section
for res in resList:
    if resDict[res].resName == 'GLY':
        writeTo.write('RESI %-5s       0.0\n' % res )
        writeTo.write('GROUP\n')
        writeTo.write('ATOM   B  %-5s  0.0\n' % (res+'B') )
        writeTo.write('BOND   B +B\n')
        writeTo.write('ANGLE -B  B +B\n')
        writeTo.write('DIHE  -B  B +B #B\n')
        writeTo.write('\n')
    else:
        writeTo.write('RESI %-5s       0.0\n' % res )
        writeTo.write('GROUP\n')
        writeTo.write('ATOM   B  %-5s  0.0\n' % (res+'B') )
        writeTo.write('ATOM   S  %-5s  0.0\n' % (res+'S') )
        writeTo.write('BOND   B  S  B +B\n')
        writeTo.write('ANGLE -B  B  S  S  B +B -B  B +B\n')
        writeTo.write('DIHE  -B  B +B #B\n')
        writeTo.write('IMPH   B -B +B  S\n')
        writeTo.write('\n')
writeTo.write('\n')
writeTo.write('END\n')
writeTo.close()

# (13 ) Write a .prm file
#   (0) CG Parameter Values
#   (a) Bond
#   (b) Angle
#   (c) Dihedral
#   (d) Improper Dihedral
#   (e) Non-Bonded
#   (f) Backbone hydrogen bonding
#   (g) Native sidechain interaction
#   (h) Backbone sidechain interaction
#   (i) END of file information & Checksum

outputName = 'cg_' + pdbName + '.prm'
writeTo = open(outputPath + outputName,'w')
writeTo.write('* This CHARMM .param file describes a Go model of %s\n' % pdbName)
writeTo.write('* contactParmSet     = %s   \n' % contactParmSet)
writeTo.write('* nScale             = %6.3f\n' % nScale)
writeTo.write('* domainScale        = %6.3f\n' % domainScale)
writeTo.write('* fnn                = %6.3f\n' % fnn)
writeTo.write('* contactRad         = %6.3f\n' % contactRad)
writeTo.write('* kBond              = %6.3f\n' % kBond)
writeTo.write('* kAngle             = %6.3f\n' % kAngle)
writeTo.write('* kDiheHelix_1       = %6.3f\n' % kDiheHelix_1)
writeTo.write('* kDiheHelix_3       = %6.3f\n' % kDiheHelix_3)
writeTo.write('* kDiheNoHelix_1     = %6.3f\n' % kDiheNoHelix_1)
writeTo.write('* kDiheNoHelix_3     = %6.3f\n' % kDiheNoHelix_3)
writeTo.write('* hBondEnergyHelix   = %6.3f\n' % hBondEnergyHelix)
writeTo.write('* hBondEnergyNoHelix = %6.3f\n' % hBondEnergyNoHelix)
writeTo.write('* epsilonNN          = %6.3e\n' % epsilonNN)
writeTo.write('* bbscInteraction    = %6.3f\n' % bbscInteraction)
writeTo.write('*\n')
writeTo.write('\n')

writeTo.write('\n')                                    # (a)
writeTo.write('BOND\n')

#   BB(i)/BB(i+1) length
for i,res in enumerate(resList):
    try:
        bondLength = cgbbDict[res].bond_length(cgbbDict[resList[i+1]])
        stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',resList[i+1]+'B',kBond,bondLength)
        writeTo.write(stringBuffer)
    except IndexError:
        pass

#   BB(i)/SC(i) length
for res in resList:
    try:
        bondLength = cgbbDict[res].bond_length(cgscDict[res])
        stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',res+'S',kBond,bondLength)
        writeTo.write(stringBuffer)
    except KeyError:
        pass

writeTo.write('\n')                                    # (b)
writeTo.write('ANGLE\n')

#   BB(i)/BB(i+1)/BB(i+2) angle
for i,res in enumerate(resList):
    try:
        bondAngle = cgbbDict[res].bond_angle(cgbbDict[resList[i+1]],cgbbDict[resList[i+2]],'deg')
        stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',resList[i+1]+'B',resList[i+2]+'B',kAngle,bondAngle)
        writeTo.write(stringBuffer)
    except IndexError:
        pass

#   SC(i)/BB(i)/BB(i+1) angle
for i,res in enumerate(resList):
    try:
        bondAngle = cgscDict[res].bond_angle(cgbbDict[res],cgbbDict[resList[i+1]],'deg')
        stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'S',res+'B',resList[i+1]+'B',kAngle,bondAngle)
        writeTo.write(stringBuffer)
    except IndexError:
        pass
    except KeyError:
        pass

#   BB(i)/BB(i+1)/SC(i+1) angle
for i,res in enumerate(resList):
    try:
        bondAngle = cgbbDict[res].bond_angle(cgbbDict[resList[i+1]],cgscDict[resList[i+1]],'deg')
        stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',resList[i+1]+'B',resList[i+1]+'S',kAngle,bondAngle)
        writeTo.write(stringBuffer)
    except IndexError:
        pass
    except KeyError:
        pass

writeTo.write('\n')                                    # (c)
writeTo.write('DIHEDRAL\n')
writeTo.write('! Backbone\n')

#   BB(i)/BB(i+1)/BB(i+2)/BB(i+3) dihedral
for i,res in enumerate(resList):
    try:
        dihedral = cgbbDict[res].bond_signed_dihedral(cgbbDict[resList[i+1]],cgbbDict[resList[i+2]],cgbbDict[resList[i+3]],'deg')
        if cgbbDict[resList[i+1]].helix and cgbbDict[resList[i+2]].helix:
            kDihe   = ['',kDiheHelix_1,'',kDiheHelix_3]
        else:
            kDihe   = ['',kDiheNoHelix_1,'',kDiheNoHelix_3]

        for multiplicity in [1,3]:
            delta       = Etc.dihedral_mod( multiplicity * dihedral - 180. )
            kDihedral   = kDihe[multiplicity]

            stringBuffer = '%-8s%-8s%-8s%-8s  %12.6f%3d%12.6f\n' % (\
                    res+'B',resList[i+1]+'B',resList[i+2]+'B',resList[i+3]+'B',kDihedral,multiplicity,delta)
            writeTo.write(stringBuffer)
    except IndexError:
        pass

writeTo.write('\n')                                    # (d)
writeTo.write('IMPROPERS\n')
writeTo.write('! Sidechain\n')

#   BB(i+1)/BB(i)/BB(i+2)/SC(i+1) dihedral
for i,res in enumerate(resList):
    try:
        dihedral = cgbbDict[resList[i+1]].bond_signed_dihedral(cgbbDict[res],cgbbDict[resList[i+2]],cgscDict[resList[i+1]],'deg')
        multiplicity = 1
        kDihedral    = 20. * abs(avgEpsilon)
        delta        = dihedral + 180.
        stringBuffer = '%-8s%-8s%-8s%-8s  %12.6f%3d%12.6f\n' % (\
                resList[i+1]+'B',res+'B',resList[i+2]+'B',resList[i+1]+'S',kDihedral,multiplicity,delta)
    except IndexError:
        pass
    except KeyError:
        pass

writeTo.write('\n')                                    # (e)
writeTo.write('NONBONDED NBXMOD 4 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -\n')
writeTo.write('CUTNB 23 CTOFNB 20 CTONNB 18 EPS 1.0 WMIN 1.5 E14FAC 0.7\n')
writeTo.write('!atom\n')

for res in cgList:
    if res.atomType == '   B':
        rMinDiv2 = 20.
    elif res.atomType == '   S':
        rMinDiv2 = 10. * Etc.get_aa_vdw( res.resName ) * 2.**(1./6.) * fnn
    else:
        raise DeadlyErrorOfDeath
    stringBuffer = '%-8s  %3.1f%10s%12.6f\n' % (\
            res.resCode+res.atomType.strip(),0,str(-1*epsilonNN),rMinDiv2)
    writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (f)
writeTo.write('NBFIX\n')
writeTo.write('! backbone hydrogen bonding\n')

bbEnergySum = 0
for hBond in bbHBond:
    bondLength  = hBond[0].bond_length(hBond[1])
    comment     = ''
    hBondEnergy = hBond[2]
    if hBond[0].domainid != hBond[1].domainid:
        comment     += '! Interface between domains %d, %d ' % (hBond[0].domainid,hBond[1].domainid)
        hBondEnergy *= domainScale
    if hBond[3] != 1:
        comment     += '! hBond multiplcity is %d ' % hBond[3]
        hBondEnergy *= hBond[3]
    bbEnergySum += hBondEnergy
    stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            hBond[0].resCode+'B',hBond[1].resCode+'B',hBondEnergy,bondLength,comment)
    writeTo.write(stringBuffer)

writeTo.write('! native side-chain interactions\n')    # (g)
scEnergySum = 0
for contact in scContact:
    bondLength  = contact[0].bond_length(contact[1])
    comment     = ''
    ljDepth     = -1. * get_contact_energy( contact[0].resName,contact[1].resName )
    if contact[0].domainid != contact[1].domainid:
        comment     += '! Interface between domains %d, %d ' % (contact[0].domainid,contact[1].domainid)
        ljDepth     *= domainScale/nScale
    scEnergySum += ljDepth
    stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            contact[0].resCode+'S',contact[1].resCode+'S',ljDepth,bondLength,comment)
    writeTo.write(stringBuffer)

writeTo.write('! backbone side-chain interactions\n')  # (h)
bbscEnergySum = 0
for contact in bbscContact:
    bondLength  = contact[0].bond_length(contact[1])
    comment     = ''
    ljDepth     = bbscInteraction
    if contact[0].domainid != contact[1].domainid:
        comment     += '! Interface between domains %d, %d ' % (contact[0].domainid,contact[1].domainid)
        ljDepth     *= domainScale
    bbscEnergySum += ljDepth
    stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            contact[0].resCode+'B',contact[1].resCode+'S',ljDepth,bondLength,comment)
    writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (i)
writeTo.write('! Czech Sum Info:%5d,%8.2f,%8.2f\n' % (bbEnergySum,scEnergySum,bbscEnergySum) )
writeTo.write('END\n')
writeTo.close()

