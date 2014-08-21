

__author__ = ("Frank C. Pickard IV <frank.pickard@nih.gov>")
__all__ = ["open_rtf"]

from collections import deque
import warnings

from pychm.future.lib import toppar as tp
from pychm.future.io.charmm.base import CharmmCard
from pychm.future.tools import mydict, paragraphs
from pychm.future.io.charmm import readwrite as rw


def open_rtf(fname, mode='r', buffering=None, inline_comments=False):
    tmp = RTFFile(fname, mode, buffering)
    if 'r' in mode:
        tmp.parse(inline_comments=inline_comments)
    if 'w' in mode:
        tmp.version = (35, 1)
    return tmp


class RTFFile(CharmmCard):
    """
    """
    def __init__(self, fname, mode='r', buffering=None):
        self.mass = None
        self.residue = None
        self.patch = None
        #
        self.commands = mydict()
        super(RTFFile, self).__init__(fname, mode, buffering)

    def parse(self, inline_comments=False):
        """
        """
        self.deque = deque(self.iter_normalize_card(comments=inline_comments))
        self.title = self._parse_title()
        self.version = self._parse_version()
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
                print self.deque[0]
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

    def export_to_toppar(self, toppar):
        self._export_mass(toppar)
        self._export_residue(toppar)
        self._export_patch(toppar)
        self._export_commands(toppar)

    def import_from_toppar(self, toppar):
        self._import_mass(toppar)
        self._import_residue(toppar)
        self._import_patch(toppar)
        self._import_commands(toppar)

    def pack(self):
        tmp = []
        tmp.append(self.pack_title())
        tmp.append(self.pack_version())
        tmp.append('')
        if self.mass is not None:
            tmp.extend(self.mass)
            tmp.append('')
        if self.commands['declare'] is not None:
            tmp.append(self.commands['declare'])
            tmp.append('')
        if self.commands['default'] is not None:
            tmp.append(self.commands['default'])
            tmp.append('')
        if self.commands['autogen'] is not None:
            tmp.append(self.commands['autogen'])
            tmp.append('')
        if self.residue is not None:
            tmp.extend(self.residue)
            tmp.append('')
        if self.patch is not None:
            tmp.extend(self.patch)
            tmp.append('')
        tmp.append('end')
        tmp.append('')
        tmp.append('')
        return '\n'.join(tmp)

    def write_all(self):
        self.write(self.pack())

# parsing private methods
    def _parse_mass(self):
        tmp = []
        while self.deque[0].startswith('mass'):
            tmp.append(self.deque.popleft())
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

# exporting private methods
    def _export_mass(self, toppar):
        if self.mass is not None:
            masses = [ rw.mass_reader(line) for line in self.mass ]
            merged = tp._merge_mass(masses, toppar.mass)
            if merged is not None:
                toppar.mass = merged

    def _export_residue(self, toppar):
        if self.residue is not None:
            resi = [ rw.residue_reader(block) for block in self.residue ]
            merged = tp._merge_section(resi, toppar.residue)
            if merged is not None:
                toppar.residue = merged

    def _export_patch(self, toppar):
        if self.patch is not None:
            pres = [ rw.patch_reader(block) for block in self.patch ]
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

# importing private methods
    def _import_mass(self, toppar):
        if toppar.mass is not None:
            self.mass = [ rw.mass_writer(mass) for mass in toppar.mass ]

    def _import_residue(self, toppar):
        if toppar.residue is None:
            self.residue = None
        else:
            self.residue = [ rw.residue_writer(resi) for resi in toppar.residue ]

    def _import_patch(self, toppar):
        if toppar.patch is None:
            self.patch = None
        else:
            self.patch = [ rw.patch_writer(pres) for pres in toppar.patch ]

    def _import_commands(self, toppar):
        self.commands['declare'] = toppar.commands['declare']
        self.commands['autogen'] = toppar.commands['autogen']
        self.commands['default'] = toppar.commands['default']

