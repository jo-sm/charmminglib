

__author__ = ("Frank C. Pickard IV <frank.pickard@nih.gov>")
__all__ = ["open_rtf"]

from collections import deque
import warnings

from pychm.future.lib import toppar as tp
from pychm.future.io.charmm.base import CharmmCard
from pychm.future.tools import _mydict, paragraphs
import pdb


def open_rtf(fname, mode='r', buffering=None, **kwargs):
    tmp = RTFFile(fname, mode, buffering)
    if 'r' in mode:
        tmp.parse()
    if 'w' in mode:
        tmp.version = (35, 1)
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
            return 'mass %-5d%-8s%10.4f' % args[:3]
        else:
            return 'mass %-5d%-8s%10.4f%3s' % args
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
        self.mass = None
        self.residue = None
        self.patch = None
        #
        self.commands = _mydict()
        super(RTFFile, self).__init__(fname, mode, buffering)

    def parse(self):
        """
        """
        super(RTFFile, self).parse()
        while 1:
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
        if self.mass is None:
            warnings.warn("No mass section found.")
        if self.residue is None:
            warnings.warn("No residue information found.")
        if self.patch is None:
            warnings.warn("No patch residue information found.")
        if self.commands['declare'] is None:
            warnings.warn("No declaration section found.")
        if self.commands['default'] is None:
            warnings.warn("No default patching section found.")
        if self.commands['autogen'] is None:
            warnings.warn("No autogenerate section found.")

    def _parse_mass(self):
        tmp = []
        while self.deque[0].startswith('mass'):
            tmp.append(self.deque.popleft())
        tmp = [ ' '.join(line.split()[1:]) for line in tmp ]
        if self.mass is None and tmp:
            self.mass = tmp
        elif self.mass is not None and tmp:
            self.mass.extend(tmp)

    def _parse_declare(self):
        tmp = []
        while self.deque[0].startswith('decl'):
            tmp.append(self.deque.popleft())
        if self.commands['declare'] is None and tmp:
            self.commands['declare'] = '\n'.join(tmp)
        elif self.commands['declare'] is not None and tmp:
            self.commands['declare'] += '\n'
            self.commands['declare'] += '\n'.join(tmp)

    def _parse_default(self):
        self.commands['default'] = self.deque.popleft()

    def _parse_autogen(self):
        self.commands['autogen'] = self.deque.popleft()

    def _parse_residue(self):
        resi = []
        pres = []
        other = []
        for block in paragraphs(self.deque, ('resi', 'pres', 'end')):
            if block[0].startswith('resi'):
                resi.append(block)
            elif block[0].startswith('pres'):
                pres.append(block)
            else:
                other.extend(block)
        other.append('end')
        if self.residue is None and resi:
            self.residue = resi
        elif self.residue is not None and resi:
            self.residue.extend(resi)
        if self.patch is None and pres:
            self.patch = pres
        elif self.patch is not None and pres:
            self.patch.extend(pres)
        self.deque = deque(other)

    def export_to_toppar(self, toppar):
        self._export_mass(toppar)
        self._export_residue(toppar)
        self._export_patch(toppar)
        self._export_commands(toppar)

    def _export_mass(self, toppar):
        if self.mass is not None:
            masses = [ massobj_from_charmm(line) for line in self.mass ]
        merged = tp._merge_mass(masses, toppar.mass)
        if merged is not None:
            toppar.mass = merged

    def _export_residue(self, toppar):
        if self.residue is not None:
            resi = [ resiobj_from_charmm(block) for block in self.residue ]
        merged = tp._merge_section(resi, toppar.residue)
        if merged is not None:
            toppar.residue = merged

    def _export_patch(self, toppar):
        if self.patch is not None:
            pres = [ presobj_from_charmm(block) for block in self.patch ]
        merged = tp._merge_section(pres, toppar.patch)
        if merged is not None:
            toppar.patch = merged

    def _export_commands(self, toppar):
        commands = ('declare', 'autogen', 'default')
        for command in commands:
            cmd = self.commands[command]
            if cmd is not None:
                merged = tp._merge_command(cmd, toppar.commands[command])
                if merged is not None:
                    toppar.commands[command] = cmd

    def import_from_toppar(self, toppar):
        self._import_mass(toppar)
        self._import_residue(toppar)
        self._import_patch(toppar)
        self._import_commands(toppar)

    def _import_mass(self, toppar):
        if toppar.mass is not None:
            self.mass = [ massobj_printer(mass) for mass in toppar.mass ]

    def _import_residue(self, toppar):
        if toppar.residue is None:
            self.residue = None
        else:
            self.residue = [ resiobj_printer(resi) for resi in toppar.residue ]

    def _import_patch(self, toppar):
        if toppar.patch is None:
            self.patch = None
        else:
            self.patch = [ presobj_printer(pres) for pres in toppar.patch ]

    def _import_commands(self, toppar):
        self.commands['declare'] = toppar.commands['declare']
        self.commands['autogen'] = toppar.commands['autogen']
        self.commands['default'] = toppar.commands['default']

    def pack(self):
        tmp = []
        tmp.append(self.pack_title())
        tmp.append(self.pack_version())
        tmp.append('')
        tmp.extend(self.mass)
        tmp.append('')
        tmp.append(self.commands['declare'])
        tmp.append('')
        tmp.append(self.commands['default'])
        tmp.append('')
        tmp.append(self.commands['autogen'])
        tmp.append('')
        tmp.extend(self.residue)
        tmp.append('')
        tmp.extend(self.patch)
        tmp.append('')
        tmp.append('end')
        tmp.append('')
        tmp.append('')
        return '\n'.join(tmp)

    def write_all(self):
        self.write(self.pack())
