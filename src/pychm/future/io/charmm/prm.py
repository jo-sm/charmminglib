"""DOCME
"""

__author__ = ("Frank C. Pickard IV <frank.pickard@nih.gov>")
__all__ = ["open_prm"]

import warnings

from pychm.future.lib import toppar as tp
from pychm.future.io.charmm.base import CharmmCard
from pychm.future.io.charmm import rtf
from pychm.future.tools import _mydict


def open_prm(fname, mode='r', buffering=None):
    """
    The public function responsible for mediating access to PRM file-
    like objects. Opens a file and returns a stream. If the file cannot be
    opened, an IOError is raised. If a file is opened for reading, 


    """
    if not isinstance(fname, basestring):
        raise TypeError("Invalid fname: %r" % fname)
    if not isinstance(mode, basestring):
        raise TypeError("Invalid mode: %r" % mode)
    if buffering is not None and not isinstance(buffering, int):
        raise TypeError("Invalid buffering: %r" % buffering)
    # parse modes
    modes = set(mode)
    if modes - set("arw+x") or len(mode) > len(modes):
        raise ValueError("invalid mode: %r" % mode)
    reading = "r" in modes
    writing = "w" in modes or "x" in modes
    appending = "a" in modes
    updating = "+" in modes
    if reading + writing + appending > 1:
        raise ValueError("Invalid mode: %r" % mode)
    if not (reading or writing or appending):
        raise ValueError("must have exactly one read/write/append mode")
    if "x" in modes and os.path.isfile(fname):
        raise ValueError("you may not set `mode=x` for existing files")
    # instantiate!
    tmp = PRMFile(fname,
                (reading and "r" or "") +
                (writing and "w" or "") +
                (appending and "a" or "") +
                (updating and "+" or ""),
                buffering)
    if 'r' in mode:
        tmp.parse()
    return tmp


def prmobj_from_charmm(arg, cls):
    if not isinstance(arg, basestring):
        raise TypeError("Invalid arg: %r" % arg)
    if not isinstance(cls, tp.PRM):
        raise TypeError("Invalid class: %s" % cls.__class__.__name__)
    return cls(*arg.split())


def prmobj_printer(prm):
    if isinstance(prm, tp.BondPRM):
        args = (prm.atom0, prm.atom1, prm.k, prm.eq)
        return '%-8s%-8s%14.6f%12.6f' % args
    elif isinstance(prm, tp.AnglePRM):
        args = (prm.atom0, prm.atom1, prm.atom2, prm.k, prm.eq, prm.k13, prm.eq13)
        if None in args:
            return '%-8s%-8s%-8s%14.6f%9.3f' % args[:5]
        else:
            return '%-8s%-8s%-8s%14.6f%9.3f%14.6f%9.3f' % args
    elif isinstance(prm, tp.DihedralPRM):
        args = (prm.atom0, prm.atom1, prm.atom2, prm.atom3, prm.k, prm.mult, prm.eq)
        return '%-8s%-8s%-8s%-8s%14.6f%3d%9.3f' % args
    elif isinstance(prm, tp.ImproperPRM):
        args = (prm.atom0, prm.atom1, prm.atom2, prm.atom3, prm.k, prm.mult, prm.eq)
        return '%-8s%-8s%-8s%-8s%14.6f%3d%9.3f' % args
    elif isinstance(prm, tp.CmapPRM):
        args = (prm.text,)
        return '%s' % args
    elif isinstance(prm, tp.NonbondPRM):
        args = (prm.atom, prm.ig, prm.k, prm.eq, prm.ig14, prm.k14, prm.eq14)
        if None in args:
            return '%-8s%5.1f%10s%12.6f' % args[:4]
        else:
            return '%-8s%5.1f%10s%12.6f%12.6f%12.6f%12.6f' % args
    elif isinstance(prm, tp.NBFixPRM):
        args = (prm.atom0, prm.atom1, prm.k, prm.eq, prm.ig14, prm.k14, prm.eq14)
        if None in args:
            return '%-8s%-8s%14.6f%12.6f' % args[:4]
        else:
            return '%-8s%-8s%14.6f%12.6f%12.6f%12.6f' % args
    elif isinstance(prm, tp.HBondPRM):
        args = (prm.atom0, prm.atom1, prm.k, prm.eq)
        return '%-8s%-8s%14.6f%12.6f' % args
    elif isinstance(prm, tp.Mass):
        args = (prm.id, prm.atom, prm.mass)
        return 'mass %-5d%-8s%14.6f' % args
    else:
        raise TypeError("Invalid class: %s" % prm.__class__.__name__)


class PRMFile(CharmmCard):
    """
    :TODO:
        allow kwargs to inject attributes
    """
    _headers = ('atom', 'bond', 'angl', 'thet', 'dihe', 'phi', 'impr',
                    'imph', 'cmap', 'nbon', 'nonb', 'nbfi', 'hbon', 'end')

    sections = ('atom', 'bond', 'angle', 'dihedral', 'improper', 'cmap',
                'nonbond', 'nbfix', 'hbond')

    prm_class_map = {'bond': tp.BondPRM,
                    'angle': tp.AnglePRM,
                    'dihedral': tp.DihedralPRM,
                    'improper': tp.ImproperPRM,
                    'cmap': tp.CmapPRM,
                    'nonbond': tp.NonbondPRM,
                    'nbfix': tp.NBFixPRM,
                    'hbond': tp.HBondPRM}

    def __init__(self, fname, mode='r', buffering=None):
        self.atom = None
        self.bond = None
        self.angle = None
        self.dihedral = None
        self.improper = None
        self.cmap = None
        self.nonbond = None
        self.nbfix = None
        self.hbond = None
        #
        self.commands = _mydict()
        super(PRMFile, self).__init__(fname=fname, mode=mode,
                                    buffering=buffering)

    def parse(self):
        """
        """
        super(PRMFile, self).parse()
        while 1:
            try:
                if self.deque[0].startswith('atom'):
                    self._parse_atom()
                if self.deque[0].startswith('bond'):
                    self._parse_section('bond')
                elif self.deque[0].startswith(('angl', 'thet')):
                    self._parse_section('angle')
                elif self.deque[0].startswith(('dihe', 'phi')):
                    self._parse_section('dihedral')
                elif self.deque[0].startswith(('impr', 'imph')):
                    self._parse_section('improper')
                elif self.deque[0].startswith('cmap'):
                    self._parse_section('cmap')
                elif self.deque[0].startswith(('nbon', 'nonb')):
                    self._parse_section('nonbond')
                elif self.deque[0].startswith('nbfi'):
                    self._parse_section('nbfi')
                elif self.deque[0].startswith('hbon'):
                    self._parse_section('hbond')
                elif self.deque[0].startswith('end'):
                    break
                else:
                    raise AssertionError("Parse: How did I get here?")
            except IndexError:
                warnings.warn("No END statement found.")
                break
        if self.atom is None:
            warnings.warn("No atom section found.")
        if self.bond is None:
            warnings.warn("No bond section found.")
        if self.angle is None:
            warnings.warn("No default angle section found.")
        if self.dihedral is None:
            warnings.warn("No dihedral section found.")
        if self.improper is None:
            warnings.warn("No improper section found.")
        if self.cmap is None:
            warnings.warn("No cmap section found.")
        if self.nonbond is None:
            warnings.warn("No nonbond section found.")
        if self.nbfix is None:
            warnings.warn("No nbfix section found.")
        if self.hbond is None:
            warnings.warn("No hbond section found.")

    def _parse_atom(self):
        self.commands['atom'] = self.deque.popleft()
        tmp = []
        while self.deque[0].startswith('mass'):
            tmp.append(self.deque.popleft())
        tmp = [ ' '.join(line.split()[1:]) for line in tmp ]
        if self.atom is None and tmp:
            self.atom = tmp
        elif self.atom is not None and tmp:
            self.atom.extend(tmp)

    def _parse_section(self, section):
        if self.commands[section] is None:
            self.commands[section] = self.deque.popleft()
        else:
            self.deque.popleft()
        tmp = []
        while 1:
            if self.deque[0].startswith(self._headers):
                break
            else:
                tmp.append(self.deque.popleft())
            if not self.deque:
                break
        if getattr(self, section) is None:
            setattr(self, section, tmp)
        else:
            getattr(self, section).extend(tmp)

    def export_to_toppar(self, toppar):
        for section in self.sections:
            if section == 'cmap':
                self._export_cmap(toppar)
            elif section == 'atom':
                self._export_atom(toppar)
            else:
                self._export_section(section, toppar)

    def _export_cmap(self, toppar):
        toppar.commands['cmap'] = tp._merge_command(self.commands['cmap'],
                                                    toppar.commands['cmap'])
        if self.cmap is not None:
            cmap = [ tp.CmapPRM('\n'.join(self.cmap)) ]
            toppar.cmap = tp._merge_cmap(cmap, toppar.cmap)

    def _export_atom(self, toppar):
        if self.atom is not None:
            masses = [ rtf.massobj_from_charmm(line) for line in self.atom ]
            toppar.mass = tp._merge_mass(masses, toppar.mass)

    def _export_section(self, section, toppar):
        toppar.commands[section] = tp._merge_command(self.commands[section],
                                                    toppar.commands[section])
        self_section = getattr(self, section)
        if self_section is not None:
            cls = self.prm_class_map[section]
            objs = [ prmobj_from_charmm(line, cls) for line in self_section ]
            setattr(toppar, section, tp._merge_section(objs,
                                                    getattr(toppar, section)))

    def import_from_toppar(self, toppar):
        for section in self.sections:
            if section == 'atom':
                self._import_atom(toppar)
            else:
                self._import_section(section, toppar)

    def _import_atom(self, toppar):
        if toppar.mass is not None:
            self.commands['atom'] = 'atom'
            self.atom = [ prmobj_printer(prm) for prm in toppar.mass ]

    def _import_section(self, section, toppar):
        self.commands[section] = toppar.commands[section]
        toppar_section = getattr(toppar, section)
        if toppar_section is not None:
            tmp = [ prmobj_printer(prm) for prm in getattr(toppar, section) ]
            setattr(self, section, tmp)

    def pack(self):
        tmp = []
        tmp.append(self.pack_title())
        tmp.append('')
        tmp.append('')
        for section in self.sections:
            if self.commands[section] is not None:
                tmp.append(self.commands[section])
            if getattr(self, section) is not None:
                tmp.extend(getattr(self, section))
                tmp.append('')
                tmp.append('')
        tmp.append('')
        tmp.append('')
        tmp.append('end')
        tmp.append('')
        tmp.append('')
        return '\n'.join(tmp)

    def write_all(self):
        self.write(self.pack())
