

__author__ = ("Frank C. Pickard IV <frank.pickard@nih.gov>")
__all__ = ["open_rtf"]

import warnings

from pychm.future.lib import toppar as tp
from pychm.future.io.charmm.base import CharmmCard


def open_rtf(fname, mode='r', buffering=None, **kwargs):
    tmp = RTFFile(fname, mode, buffering)
    if 'r' in mode:
        tmp.parse()
    return tmp


def massobj_from_charmm(arg):
    if not isinstance(arg, basestring):
        raise TypeError("Invalid arg: %r" % arg)
    return tp.Mass(*arg.split())


def resiobj_from_charmm(arg):
    x, name, charge = arg[0].split()
    body = arg[1:]
    return tp.Residue(name=name, charge=charge, body=body)


def presobj_from_charmm(arg):
    x, name, charge = arg[0].split()
    body = arg[1:]
    return tp.Patch(name=name, charge=charge, body=body)


def massobj_printer(arg):
    if isinstance(arg, tp.Mass):
        args = (arg.id, arg.atom, arg.mass, arg.element)
        if None in args:
            return 'mass %-5d%-8s%14.6f' % args[:3]
        else:
            return 'mass %-5d%-8s%14.6f%-3s' % args
    else:
        raise TypeError("Invalid class: %s" % arg.__class__.__name__)


def resiobj_printer(arg):
    if isinstance(arg, tp.Residue):
        tmp = []
        tmp.append('resi %6s%8.2f' % (arg.name, arg.charge))
        tmp.extend(arg.body)
        tmp.append('')
        return '\n'.join(tmp)
    else:
        raise TypeError("Invalid class: %s" % arg.__class__.__name__)


def presobj_printer(arg):
    if isinstance(arg, tp.Patch):
        tmp = []
        tmp.append('pres %6s%8.2f' % (arg.name, arg.charge))
        tmp.extend(arg.body)
        tmp.append('')
        return '\n'.join(tmp)
    else:
        raise TypeError("Invalid class: %s" % arg.__class__.__name__)


class RTFFile(CharmmCard):
    """
    """
    def __init__(self, fname, mode='r', buffering=None):
        super(RTFFile, self).__init__(fname, mode, buffering)
        self.mass = None
        self.declare = None
        self.default = None
        self.autogen = None
        self.residue = None
        self.patch = None

    def parse(self):
        """
        """
        super(RTFFile, self).parse()
        while 1:
            try:
                if self.deque[0].startswith('mass'):
                    self._parse_mass()
                elif self.deque[0].startswith('decl'):
                    self._parse_declare()
                elif self.deque[0].startswith('defa'):
                    self._parse_default()
                elif self.deque[0].startswith('auto'):
                    self._parse_autogen()
                elif self.deque[0].startswith(('resi', 'pres')):
                    self._parse_residue()
                elif self.deque[0].startswith('end'):
                    break
                else:
                    raise AssertionError("Parse: How did I get here?")
            except IndexError:
                warnings.warn("No END statement found.")
        if self.mass is None:
            warnings.warn("No mass section found.")
        if self.declare is None:
            warnings.warn("No declaration section found.")
        if self.default is None:
            warnings.warn("No default patching section found.")
        if self.autogen is None:
            warnings.warn("No autogenerate section found.")
        if self.residue is None:
            warnings.warn("No residue information found.")
        if self.patch is None:
            warnings.warn("No patch residue information found.")

    def _parse_mass(self):
        tmp = []
        while self.deque[0].startswith('mass'):
            tmp.append(self.deque.popleft())
        if tmp:
            self.mass = [ ' '.join(line.split()[1:]) for line in tmp ]

    def _parse_declare(self):
        tmp = []
        while self.deque[0].startswith('decl'):
            tmp.append(self.deque.popleft())
        if tmp:
            self.declare = tmp

    def _parse_default(self):
        self.default = self.deque.popleft()

    def _parse_autogen(self):
        self.autogen = self.deque.popleft()

    def _parse_residue(self):
        resi = []
        pres = []
        while 1:
            tmp = []
            try:
                if self.deque[0].startswith('resi'):
                    res_type = 'resi'
                elif self.deque[0].startswith('pres'):
                    res_type = 'pres'
                else:
                    break
                while 1:
                    tmp.append(self.deque.popleft())
                    try:
                        if self.deque[0].startswith(('resi', 'pres', 'end')):
                            break
                    except IndexError:
                        break
                if res_type == 'resi':
                    resi.append(tmp)
                elif res_type == 'pres':
                    pres.append(tmp)
            except IndexError:
                break
        self.residue = resi
        self.patch = pres

    def export_to_toppar(self, toppar):
        toppar.mass = self._export_mass()
        toppar.residue = self._export_residue()
        toppar.patch = self._export_patch()

    def _export_mass(self):
        return [ massobj_from_charmm(line) for line in self.mass ]

    def _export_residue(self):
        return [ resiobj_from_charmm(arg) for arg in self.residue ]

    def _export_patch(self):
        return [ presobj_from_charmm(arg) for arg in self.patch ]
