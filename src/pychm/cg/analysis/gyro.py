

import os
import cPickle as pickle
from numpy import ndarray
import matplotlib.pyplot as pyplot
from pychm.tools import Property, lowerKeys
from pychm.analysis.baseanalysis import BaseAnalysis, load_correlOutput


class Gyro(BaseAnalysis):
    """
    """
    def __init__(self, aaCrdFile, **kwargs):
        super(Gyro,self).__init__(aaCrdFile, **kwargs)
        self.inpFilename = '&/gyro_%s.inp' % self.correlAtomSelection
        self.outFilename = '&/gyro_%s.out' % self.correlAtomSelection
        self.anlFilename = '&/gyro_%s.anl' % self.correlAtomSelection

    @Property
    def anlFilename():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._anlFilename
        def fset(self, value):
            self._anlFilename = self.expandPath(value)
        return locals()

    @Property
    def data():
        doc =\
        """
        The data property.
        """
        def fget(self):
            fileName = '%s.pickle' % self.anlFilename
            try:
                self._data = pickle.load(open(fileName))
                print 'found pickled data data in %s ...' % fileName
                return self._data
            except IOError:
                print 'processing data data in %s ...' % self.anlFilename
                self._data = load_correlOutput('%s' % self.anlFilename)
                pickle.dump(self._data ,open(fileName, 'w'))
                return self._data
        def fset(self, value):
            assert isinstance(value,ndarray)
            self._data = value
        def fdel(self):
            fileName = '%s.pickle' % self.anlFilename
            del self._data
            try:
                os.remove(fileName)
            except IOError:
                pass
        return locals()

# Charmm input creation

    def do_correl(self, **kwargs):
        """
        """
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        comments = kwargs.get('comments', [])
        #
        self.write_correlInput(**kwargs)
        self.run_correlInput(**kwargs)

    def write_correlInput(self, **kwargs):
        """
        This method writes a single charmm .inp file for a radius of gyration calculation.
        """
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        comments = kwargs.get('comments', [])
        #
        String = self.get_correlInputHeader(header)
        String.append('! anl :: write')
        String.append('open unit 100 write card name %s' % self.anlFilename)
        String.append('')
        if self.correlAtomSelection == 'all':
            String.append('defi taco select all end')
        else:
            String.append('defi taco select ')
            for i,letter in enumerate(self.correlAtomSelection):
                String.append('segid %s ' % letter.upper())
                if i != len(self.correlAtomSelection)-1:
                    String.append('.or. ')
            String.append('end')
        String.append('')
        String.append('traj query unit 10')
        String.append('correl maxtimesteps %d maxatom %d maxseries 1' % (self.correlArrayLength, self.maxatom))
        String.append('enter gyro gyration')
        String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select taco end' % (self.correlStart, self.correlStop, self.correlSkip))
        String.append('')
        String.append('write gyro card unit 100 ')
        String.append('* Radius of Gyration - Selection: %s' % self.correlAtomSelection)
        for comment in comments:
            String.append('* %s' % comment)
        String.append('*')
        String.append('')
        String.append('stop')
        String = '\n'.join(String)
        # Write file
        write_to = open('%s' % self.inpFilename, 'w')
        write_to.write(String)
        write_to.close()

    def run_correlInput(self, **kwargs):
        """
        This method runs a single charmm .inp file for a radius of gyration calculation.
        """
        kwargs = lowerKeys(kwargs)
        #
        try:
            os.remove(self.outFilename)
        except OSError:
            pass
        if not os.path.exists(self.inpFilename):
            raise IOError("No such file or directory: '%s'" % self.inpFilename)
        os.system('%s < %s > %s' % (self.charmmBin, self.inpFilename, self.outFilename))
