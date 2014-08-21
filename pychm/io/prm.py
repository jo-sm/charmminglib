"""
:Author: fcp
:Date: 02/22/2011
"""


from copy import deepcopy
from pychm.tools import cleanStrings, logicalLines, paragraphs
from pychm.lib.toppar import AnglePRM, BondPRM, DihedralPRM, ImproperPRM, \
        MassPRM, NonBondPRM, NBFixPRM
from pychm.io.basecharmm import BaseCHARMMFile


class PRMFile(BaseCHARMMFile):
    """
    **STUB**
    """

    _sections = ['atom','bond','angl','thet','dihe','phi','impr','imph','cmap',
                'nbon','nonb','nbfi','hbon','end']

    def __init__(self, filename=None, **kwargs):
        super(PRMFile, self).__init__(filename, **kwargs)
        self.cmd = {}
        self.prm = {}
        self.comments = []
        self.warnings = []
        if filename is not None:
            self.parse()

    def parse(self):
        super(PRMFile, self).parse()
        iterator = ( line for line in self.body if not line.startswith('!') )
        iterator = paragraphs(logicalLines(iterator), self._sections)
        for taco in iterator:
            try:
                if taco[0][:4] in self._sections:
                    self.cmd[taco[0][:4]] = taco[0]
                    self.prm[taco[0][:4]] = taco[1:]
                else:
                    self.comments.append(taco)
            except IndexError:
                pass

        # discard empties
        tmp = []
        for key, value in self.prm.iteritems():
            if not value:
                tmp.append(key)
        for key in tmp:
            try:
                del self.cmd[key]
            except KeyError:
                pass
            try:
                del self.prm[key]
            except KeyError:
                pass

        # rename sections
        def rename(old, new):
            if old in self.cmd.keys():
                self.cmd[new] = self.cmd[old]
                del self.cmd[old]
                self.prm[new] = self.prm[old]
                del self.prm[old]
        rename('thet', 'angl')
        rename('phi', 'dihe')
        rename('imph', 'impr')
        rename('nonb', 'nbon')

        # parse sections
        self._parse('atom', MassPRM)
        self._parse('bond', BondPRM)
        self._parse('angl', AnglePRM)
        self._parse('dihe', DihedralPRM)
        self._parse('impr', ImproperPRM)
        self._parse('nbon', NonBondPRM)
        self._parse('nbfi', NBFixPRM)
        # self._parse('hbon', HBondPRM) TODO -- implement HBondPRM objects in lib.toppar
        self._parse_cmap()

    def write(self, filename=None):
        """
        """
        def get_section(name):
            tmp = []
            try:
                tmp.append(self.cmd[name])
                iterator = ( prm.Print() for prm in self.prm[name] )
                tmp.extend(iterator)
                tmp.append('')
            except KeyError:
                pass
            return tmp
        #
        tmp = []
        # header
        iterator = ( '* %s' % line for line in self.header if line )
        tmp.extend(iterator)
        tmp.append('*')
        tmp.append('')
        #
        # atom
        try:
            i = 1
            for prm in self.prm['atom']:
                prm.index = i
                i += 1
        except KeyError:
            pass
        tmp.extend(get_section('atom'))
        tmp.extend(get_section('bond'))
        tmp.extend(get_section('angl'))
        tmp.extend(get_section('dihe'))
        tmp.extend(get_section('impr'))
        tmp.extend(get_section('nbon'))
        tmp.extend(get_section('nbfi'))

        # cmap -- OO cmap not implemented
        try:
            tmp.append(self.cmd['cmap'])
            tmp.extend(self.prm['cmap'])
        except KeyError:
            pass
        # TODO hbon -- hbon not yet implemented

        #
        if filename is None:
            for line in tmp:
                print line.upper()
        else:
            writeTo = open(filename, 'w')
            writeTo.write('\n'.join(tmp).upper())
            writeTo.close()

###################
# Private Methods #
###################

    def _parse_cmap(self):
        try:
            self.prm['cmap'] = self.prm['cmap']
        except KeyError:
            pass

    def _parse(self, section, obj):
        try:
            self.prm[section] = sorted(map(obj, self.prm[section]))
            wild = []
            not_wild = []
            for prm in self.prm[section]:
                if prm.is_wild():
                    wild.append(prm)
                else:
                    not_wild.append(prm)
            self.prm[section] = wild + not_wild
        except KeyError:
            pass

###################
# Special Methods #
###################

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.prm.keys())

    def __add__(self, other):
        def merge_section(obj, name):
            tmp = []
            try:
                tmp.extend(self.prm[name])
            except KeyError:
                pass
            try:
                for prm in other.prm[name]:
                    if prm not in tmp:
                        tmp.append(prm)
                    else:
                        print 'Warning, duplicate parameters found in %s section!\n' % name
                        print '%r\n' % prm
                        print 'Keeping parameter from original (left) prm object\n'
                        obj.warnings.append('Duplicate PRM: %r' % prm)
            except KeyError:
                pass
            if tmp:
                # prm
                obj.prm[name] = tmp
                # cmd
                obj.cmd[name] = ''
                try:
                    obj.cmd[name] = self.cmd[name]
                except KeyError:
                    try:
                        obj.cmd[name] = other.cmd[name]
                    except KeyError:
                        pass
                if not obj.cmd[name]:
                    print 'Warning, no CMD found for %s section!\n' % name
                    obj.warnings.append('No CMD for %s' % name)
            return
        #
        taco = PRMFile()
        taco.header = self.header + other.header
        #
        merge_section(taco, 'atom')
        merge_section(taco, 'bond')
        merge_section(taco, 'angl')
        merge_section(taco, 'dihe')
        merge_section(taco, 'impr')
        merge_section(taco, 'nbon')
        merge_section(taco, 'nbfi')
        merge_section(taco, 'hbon')
        # cmap
        try:
            taco.prm['cmap'] = self.prm['cmap'][:]
        except KeyError:
            try:
                taco.prm['cmap'] = other.prm['cmap'][:]
            except KeyError:
                pass
        if 'cmap' in taco.prm.keys():
            taco.cmd['cmap'] = 'cmap'
        #
        return taco
