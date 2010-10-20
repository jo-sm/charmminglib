from charmming.const.etc import alphanum
from charmming.lib.basestruct import BaseStruct
from charmming.lib.chain import Chain
from charmming.tools import Property


class Mol(BaseStruct):
    """
    Private Attributes
    Properties
    Public Methods
        iter_chain
        iter_seg
        iter_res
        reindex_atomNum
        reindex_resid
        reindex_resIndex
        parse
        write               STUB
        write_charmming
    Private Methods
        _fix_multi_models
        _diff_nuc
        _compliance_charmm_terminalOxygen
        _compliance_charmm_names
    Special Methods
    """
    def __init__(self,iterable=None,**kwargs):
        super(Mol,self).__init__(iterable,**kwargs)
        self.warnings = []

    def iter_chain(self,**kwargs):
        """
        For looping over chains.
        kwargs:
            chainids
        >>> thisMol.iter_chain(chainids=['a','b'])
        """
        # kwargs
        chainids = kwargs.pop('chainids',None)
        kwargs.pop('segTypes',None)
        if kwargs.keys():
            raise TypeError('Unprocessed kwargs(%r)' % kwargs.keys())
        #
        if chainids is None:
            chainids = list(set(( atom.chainid for atom in self )))
            chainids.sort()
        # Do Work
        for chainid in chainids:
            iterator = ( atom for atom in self if atom.chainid == chainid )
            result = Chain(iterable=iterator,code=self.code,autoFix=False)
            if result:
                yield result

    def iter_seg(self,**kwargs):
        """
        For looping over segments.
        kwargs:
            chainids
            segTypes
        >>> thisMol.iter_seg(chainids=['a'],segTypes=['pro','dna'])
        """
        for chain in self.iter_chain(**kwargs):
            for seg in chain.iter_seg(**kwargs):
                if seg:
                    yield seg

    def iter_res(self,**kwargs):
        """
        For looping over residues.
        kwargs:
            chainids
            segTypes
        >>> thisMol.iter_res(chainids=['a','d'],segTypes=['pro','dna'])
        """
        for chain in self.iter_chain(**kwargs):
            for seg in chain.iter_seg(**kwargs):
                for res in seg.iter_res():
                    yield res

    def reindex_atomNum(self,start=1):
        """
        Reindex the atom numbers of a Mol object, starting at 'start' value.
        """
        for seg in self.iter_seg():
            i = start
            for atom in seg:
                atom.atomNum = i
                i += 1

    def reindex_resid(self,start=1):
        """
        Reindex the resids of a residue, starting at 'start' value.
        """
        for seg in self.iter_seg():
            i = start
            for res in seg.iter_res():
                res.resid = i
                i += 1

    def reindex_resIndex(self,start=1):
        """
        Creates a master residue index on a per Mol basis. The index ordering
        will be alpha ordering, with priority: chainids > segTypes > resid.
        """
        i = start
        for res in self.iter_res():
            res.resIndex = i
            i += 1

    def write(self,fileName,**kwargs):
        """
        Writes a single output file containing all of the Mol object's data.
        kwargs:
            format
        """
        Format      = kwargs.pop('format','charmm')
        raise NotImplementedError

    def write_charmming(self):
        """
        Macro for writing CHARMMing style output.
        """
        segDict = {'nuc':'nuc','pro':'pro','good':'goodhet','bad':'het','dna':'dna','rna':'rna'}
        stdoutList = []
        # Write files
        for seg in self.iter_seg():
            stdoutList.append('%s-%s' % (seg.chainid,segDict[seg.segType]))
            name = 'new_%s-%s-%s.pdb' % (self.pdbCode,seg.chainid,segDict[seg.segType])
            seg.write(filename=name,format='charmm',ter=True,end=False)
        # To STDOUT
        print 'natom=%d' % len(self)
        print 'nwarn=%d' % len(self.warnings)
        if self.warnings:
            print 'warnings=%r' % self.warnings
        print 'seg=%r' % stdoutList

###########
# Parsing #
###########

    def parse(self):
        """
        Performs all of the functions of the CHARMMing parser. Multi-model
        sections are tidied up, meaning the most heavily weighted atom in each
        section is kept, and the rest are deleted. DNA/RNA are differentiated
        based on chemical information and naming conventions of the PDB. Atoms,
        residues are reindexed to be CHARMM compatable. Finally residues and
        atom types are renamed for CHARMM compatability.
        """
        self._fix_multi_models()
        self._diff_nuc()
        self.sort()
        self.reindex_atomNum(start=1)
        self.reindex_resid(start=1)
        self.reindex_resIndex(start=1)
        self._compliance_charmm_terminalOxygen()
        self._compliance_charmm_names()

    def _fix_multi_models(self):
        """
        Searches for multi-model sections, and 'fixes' them. Renames
        residues, and deletes atoms.
        For example...
        15  CE AMET A   1      16.764  26.891  42.696  0.73 60.01           C
        16  CE BMET A   1      18.983  29.323  38.169  0.27 60.30           C
        maps to:
        15  CE  MET A   1      16.764  26.891  42.696  0.73 60.01           C
        """
        delete_these = []
        residues = ( res for res in self.iter_res() if len(res.resName) == 4 )
        for res in residues:
            resLen = len(res)
            if resLen < 2: continue
            tmp = []
            keep = False
            for i,atom in enumerate(res):
                try:
                    # Comparison atom, check ahead if first, check behind for all others
                    if i == 0:
                        comp = res[1]
                    else:
                        comp = res[i-1]
                    # Check last 3 chars are different, and first char is the same.
                    # for example: AMET and BMET
                    if atom.resName[-3:] == comp.resName[-3:] and atom.resName[0] != comp.resName[0]:
                        # If buffer is empty, load it up.
                        if not tmp:
                            tmp.append(atom)
                            continue
                        # If leading char increases, toss atom in buffer
                        if atom.resName[0] > tmp[-1].resName[0]:
                            tmp.append(atom)
                        # If it doesn't increase we're in a new section, cleanup!
                        else:
                            tmp2 = atom
                            raise AssertionError
                        # Last atom in a section, cleanup!
                        if i == resLen - 1:
                            raise AssertionError
                    # ResName has changed, or leading char is gone, cleanup!
                    else:
                        if tmp:
                            tmp2 = ''
                            raise AssertionError
                except AssertionError:
                    maxWeight = max(( atom1.weight for atom1 in tmp ))
                    # Write Notes
                    for atom2 in tmp:
                        if keep == False:
                            if atom2.weight == maxWeight:
                                atom2.note = 'rename'
                                keep = True
                            else:
                                atom2.note = 'delete'
                        else:
                            atom2.note = 'delete'
                    # Process Notes
                    for atom3 in tmp:
                        if atom3.note == 'rename':
                            atom3.resName = atom3.resName[-3:]
                        elif atom3.note == 'delete':
                            delete_these.append(atom3)
                        del atom3.note
                    if tmp2:
                        tmp = [tmp2]
                    else:
                        tmp = []
                    keep = False
        # Delete things
        self.del_atoms(delete_these)

    def _diff_nuc(self):
        """
        Loop over all the segments in the Mol, and differentiate them between
        dna/rna, then set their seg type appropriately. Throw a warning if the
        determination cannot be made. The list of warnings is returned, segType
        is set inline.
        """
        badSegAddr = []
        warnings = []
        # Use Thymine/Uracil to determine dna/rna.
        for seg in self.iter_seg(segTypes=['nuc']):
            thymine = False
            uracil  = False
            thyCheck = set(['t','thy','dt'])
            uraCheck = set(['u','ura','du'])
            for atom in seg:
                if atom.resName in thyCheck:
                    thymine = True
                if atom.resName in uraCheck:
                    uracil  = True
#WARN       # Throw a warning when Uracil and Thymine are found in the same segment
            if thymine and uracil:
                msg = '_diff_nuc: URA & THY in same segment, seg.addr = "%s"' % seg.addr
                badSegAddr.append(seg.addr)
                warnings.append(msg)
            elif thymine:
                seg.segType = 'dna'
            elif uracil:
                seg.segType = 'rna'
        # Use pdb residue names to infer dna/rna.
        for seg in self.iter_seg(segTypes=['nuc']):
            ribo = False
            deoxy = False
            riboCheck = ['a','c','g','u','i',]
            deoxyCheck = ['da','dc','dg','dt','di']
            for atom in seg:
                if atom.resName in riboCheck:
                    ribo = True
                if atom.resName in deoxyCheck:
                    deoxy = True
#WARN       # Throw a warning when ribo- and deoxy- are found in the same segment
            if ribo and deoxy:
                msg = '_diff_nuc: ribo- & deoxy- in same segment, seg.addr = "%s"' % seg.addr
                badSegAddr.append(seg.addr)
                warnings.append(msg)
            elif ribo:
                seg.segType = 'rna'
            elif deoxy:
                seg.segType = 'dna'
        # Final check
        for seg in self.iter_seg(segTypes=['nuc']):
            if seg.addr not in badSegAddr:
                msg = '_diff_nuc: undetermined nucleotide, seg.addr = "%s"' % seg.addr
                badSegAddr.append(seg.addr)
                warnings.append(msg)
        self.warnings = warnings

    def _compliance_charmm_terminalOxygen(self):
        """
        Fix the atom types of terminal oxygen atoms to be CHARMM compliant.
        """
        for seg in self.iter_seg(segTypes=['pro']):
            resList = list(seg.iter_res())
            lastRes = resList[-1]
            for atom in lastRes:
                if atom.atomType == ' o  ':
                    atom.atomType = ' ot1'
                elif atom.atomType == ' oxt':
                    atom.atomType = ' ot2'

    def _compliance_charmm_names(self):
        """
        On a per atom basis, fix residue names and atom types to be CHARMM
        compliant.
        """
        for atom in self:
            atom._compliance_resName()
            atom._compliance_atomType()
