"""
:Author: fcp
:Date: 08/31/2011
"""


import os
import cPickle as pickle
from numpy import ndarray
from pychm.tools import Property, lowerKeys
from pychm.analysis.baseanalysis import BaseAnalysis, load_correlOutput


class BBRMSD(BaseAnalysis):
    """
    DOCME
    """
    def __init__(self, pdbFilename=None, **kwargs):
        super(BBRMSD, self).__init__(pdbFilename, **kwargs)
        self.inpFilename = '&/bbrmsd_%s.inp' % self.correlAtomSelection
        self.outFilename = '&/bbrmsd_%s.out' % self.correlAtomSelection
        self.anlFilename = '&/bbrmsd_%s.anl' % self.correlAtomSelection

    @Property
    def data():
        doc =\
        """
        A numpy array representing the backbone RMSD as a function of time.
        """
        def fget(self):
            filename = '%s.pickle' % self.anlFilename
            try:
                self._data = pickle.load(open(filename))
                print 'found pickled data in: %s' % filename
                return self._data
            except IOError:
                print 'processing data in %s:' % self.anlFilename
                self._data = load_correlOutput('%s' % self.anlFilename)
                pickle.dump(self._data ,open(filename, 'w'))
                return self._data
        def fset(self, value):
            assert isinstance(value,ndarray)
            self._data = value
        def fdel(self):
            filename = '%s.pickle' % self.anlFilename
            del self._data
            try:
                os.remove(filename)
            except IOError:
                pass
        return locals()

# Charmm input creation
    def do_correl(self, **kwargs):
        """
        """
        self.write_correlInput(**kwargs)
        self.run_correlInput()

    def write_correlInput(self, **kwargs):
        """
        This method writes a single CHARMM .inp file for the purpose
        of performing a backbone RMSD calculation.  This RMSD is done
        between the original coarse grained .crd file and the
        reoriented trajectory 'reorient.dcd'.
        """
        kwargs = lowerKeys(kwargs)
        header = kwargs.get('header', [])
        comments = kwargs.get('comments', [])
        #
        String = self.get_correlInputHeader(header)
        String.append('! anl :: write')
        String.append('open unit 100 write card name %s' % self.anlFilename)
        String.append('')
        String.append('! copy coordinates to comparison set')
        String.append('coor copy comp')
        String.append('')
        String.append('! select backbone atoms')
        String.append('defi bb select type b end')
        String.append('')
        String.append('traj query unit 10')
        String.append('correl maxtimesteps %d maxatom %d maxseries 1' %
                    (self.correlArrayLength, self.maxatom))
        String.append('enter rmsd rms orient')
        String.append('traj firstu 10 nunit 1 begin %d stop %d skip %d select bb end' %
                    (self.correlStart, self.correlStop, self.correlSkip))
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
        write_to = open('%s' % self.inpFilename, 'w')
        write_to.write(String)
        write_to.close()

    def run_correlInput(self):
        """
        This method runs a single charmm .inp file.
        """
        try:
            os.remove(self.outFilename)
        except OSError:
            pass
        if not os.path.exists(self.inpFilename):
            raise IOError("No such file or directory: '%s'" % self.inpFilename)
        os.system('%s < %s > %s' % (self.charmmBin, self.inpFilename, self.outFilename))
