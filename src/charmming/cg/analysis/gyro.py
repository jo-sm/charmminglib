

import os
import cPickle as pickle
from math import floor,ceil
from numpy import ndarray
import matplotlib.pyplot as pyplot
from charmming.tools import Property
from charmming.analysis.baseanalysis import BaseAnalysis, load_correlOutput


class Gyro(BaseAnalysis):
    """
    """
    def __init__(self,aaCrdFile):
        super(Gyro,self).__init__(aaCrdFile)

    @Property
    def RgOfT():
        doc = "The RgOfT property."
        def fget(self):
            fileName = '%s/RgOfT_%s.pickle' % (self.anlPath,self.correlAtomSelect)
            try:
                self._RgOfT = pickle.load(open(fileName))
                print 'found pickled RgOfT data in %s ...' % fileName
                return self._RgOfT
            except IOError:
                print 'processing RgOfT data in %s ...' % self.anlPath
                self._RgOfT = load_correlOutput('%s/gyro_%s.anl' % (self.anlPath,self.correlAtomSelect))
                pickle.dump(self._RgOfT,open(fileName,'w'))
                return self._RgOfT
        def fset(self, value):
            assert isinstance(value,ndarray)
            self._RgOfT = value
        def fdel(self):
            fileName = '%s/RgOfT_%s.pickle' % (self.anlPath,self.correlAtomSelect)
            del self._RgOfT
            try:
                os.remove(fileName)
            except IOError:
                pass
        return locals()

# Charmm input creation
    def write_singleCorrelInput(self, header=[], comments=[]):
        """
        This method writes a single charmm .inp file for a radius of gyration calculation.
        """
        String = []
        String.append(self.get_correlInputHeader(header))
        String.append('! anl :: write')
        String.append('open unit 100 write card name gyro_%s.anl' % self.correlAtomSelect)
        String.append('')
        if self.correlAtomSelect == 'all':
            String.append('defi taco select all end')
        else:
            String.append('defi taco select ')
            for i,letter in enumerate(self.correlAtomSelect):
                String.append('segid %s ' % letter.upper())
                if i != len(self.correlAtomSelect)-1:
                    String.append('.or. ')
            String.append('end')
        String.append('')
        String.append('traj query unit 10')
        String.append('correl maxtimesteps %d maxatom %d maxseries 1' % (self.correlArrayLength,self.correlMaxAtom))
        String.append('enter gyro gyration')
        String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select taco end' % (self.correlStart,self.correlStop,self.correlSkip))
        String.append('')
        String.append('write gyro card unit 100 ')
        String.append('* Radius of Gyration - Selection: %s' % self.correlAtomSelect)
        for comment in comments:
            String.append('* %s' % comment)
        String.append('*')
        String.append('')
        String.append('stop')
        String = '\n'.join(String)
        # Write file
        write_to = open('%s/gyro_%s.inp' % (self.anlPath,self.correlAtomSelect),'w')
        write_to.write(String)
        write_to.close()

    def run_singleCorrelInput(self):
        """
        This method runs a single charmm .inp file for a radius of gyration calculation.
        """
        lastDir = os.getcwd()
        os.chdir(self.anlPath)
        try:
            os.remove('gyro_%s.out' % self.correlAtomSelect)
        except OSError: pass
        os.system('%s < gyro_%s.inp > gyro_%s.out' % (self.charmmBin,self.correlAtomSelect,self.correlAtomSelect))
        os.chdir(lastDir)

