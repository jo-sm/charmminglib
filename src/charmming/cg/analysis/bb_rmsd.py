"""
DOCME
"""
# Filename: bb_rmsd.py
# Owner: fcp
# 10/29/2010


import os
import cPickle as pickle
from numpy import ndarray
from charmming.analysis.cg.basecganalysis import BaseAnalysis, load_correlOutput
from charmming.tools import Property


class BB_RMSD(BaseAnalysis):
    """
    DOCME
    """
    def __init__(self, aaCrdFile):
        """
        DOCME
        """
        super(BB_RMSD, self).__init__(aaCrdFile)

    @Property
    def RMSDofT():
        doc =\
        """
        A numpy array representing the backbone RMSD as a function of time.
        """
        def fget(self):
            filename = '%s/RMSDofT.pickle' % self.anlPath
            try:
                self._RMSDofT = pickle.load(open(filename))
                print 'found pickled RMSDofT data in %s ...' % filename
                return self._RMSDofT
            except IOError:
                print 'processing RMSDofT data in %s ...' % self.anlPath
                self._RMSDofT = load_correlOutput('%s/bb_rmsd.anl' % self.anlPath)
                pickle.dump(self._RMSDofT, open(filename, 'w'))
                return self._RMSDofT
        def fset(self, value):
            assert isinstance(value, ndarray)
            self._RMSDofT = value
        def fdel(self):
            filename = '%s/RMSDofT.pickle' % self.anlPath
            del self._RMSDofT
            try:
                os.remove(filename)
            except IOError:
                pass
        return locals()

# Charmm input creation
    def write_singleReorientInput(self, header=[]):
        """
        This method writes a single charmm .inp file for the purpose
        or reorienting the trajectory to maximally align with the
        reference coordinates.  This facilitates the calculation of
        backbone RMSD wrt to said reference coordinates.
        """
        String = []
        String.append(self.get_correlInputHeader(header))
        String.append('!open files for writing')
        String.append('open unit 100 write unform name reorient.dcd')
        String.append('')
        String.append('! copy coordinates to comparison set')
        String.append('coor copy comp')
        String.append('')
        String.append('traj query unit 10')
        String.append('! reorient')
        String.append('merge firstu 10 nunit 1 begin %d stop %d skip %d output 100 -'
                    % (self.correlStart, self.correlStop, self.correlSkip))
        String.append('  orient mass')
        String.append('close unit 10')
        String.append('')
        String.append('stop')
        String = '\n'.join(String)
        # Write file
        write_to = open('%s/reorient.inp' % self.anlPath, 'w')
        write_to.write(String)
        write_to.close()

    def write_singleCorrelInput(self, header=[], comments=[]):
        """
        This method writes a single CHARMM .inp file for the purpose
        of performing a backbone RMSD calculation.  This RMSD is done
        between the original coarse grained .crd file and the
        reoriented trajectory 'reorient.dcd'.
        """
        self.dcdFile, tmpDcdFile = '%s/reorient.dcd' % self.anlPath, self.dcdFile
        #
        String = []
        String.append(self.get_correlInputHeader(header))
        String.append('! anl :: write')
        String.append('open unit 100 write card name bb_rmsd.anl')
        String.append('')
        String.append('! copy coordinates to comparison set')
        String.append('coor copy comp')
        String.append('')
        String.append('! select backbone atoms')
        String.append('defi bb select type b end')
        String.append('')
        String.append('traj query unit 10')
        String.append('correl maxtimesteps %d maxatom %d maxseries 1' % (self.correlArrayLength,self.correlMaxAtom))
        String.append('enter rmsd rms orient')
        String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select bb end' % (self.correlStart,self.correlStop,self.correlSkip))
        String.append('')
        String.append('write rmsd card unit 100 ')
        String.append('* Backbone RMSD')
        for comment in comments:
            String.append('* %s' % comment)
        String.append('*')
        String.append('')
        String.append('stop')
        String = '\n'.join(String)
        # Write file
        write_to = open('%s/bb_rmsd.inp' % self.anlPath, 'w')
        write_to.write(String)
        write_to.close()
        #
        self.dcdFile = tmpDcdFile

    def run_singleReorientInput(self):
        """
        This method runs a single charmm .inp file for a trajectory
        reorientation.
        """
        lastDir = os.getcwd()
        os.chdir(self.anlPath)
        try:
            os.remove('reorient.out')
        except OSError:
            pass
        os.system('%s < reorient.inp > reorient.out' % self.charmmBin)
        os.chdir(lastDir)

    def run_singleCorrelInput(self):
        """
        This method runs a single charmm .inp file for a backbone RMSD
        calculation.
        """
        lastDir = os.getcwd()
        os.chdir(self.anlPath)
        try:
            os.remove('bb_rmsd.out')
        except OSError:
            pass
        os.system('%s < bb_rmsd.inp > bb_rmsd.out' % self.charmmBin)
        os.chdir(lastDir)

