

from __future__ import division

__author__ = ("Frank C. Pickard IV <frank.pickard@nih.gov>")
__all__ = ["open_prm"]

import warnings

from pychm.future.lib import toppar as tp
from pychm.future.io.charmm.base import CharmmCard


def open_prm(fname, mode='r', buffering=None, **kwargs):
    tmp = PRMFile(fname, mode, buffering)
    if 'r' in mode:
        tmp.parse()
    return tmp


def prmobj_from_charmm(datastring, cls):
    if not isinstance(datastring, basestring):
        raise TypeError("Invalid datastring: %r" % datastring)
    if not isinstance(cls, tp.PRM):
        raise TypeError("Invalid class: %s" % cls.__class__.__name__)
    return cls(*datastring.split())


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
        args = (prm.atom0, prm.ig, prm.k, prm.eq, prm.ig14, prm.k14, prm.eq14)
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
    else:
        raise TypeError("Invalid class: %s" % prm.__class__.__name__)


class PRMFile(CharmmCard):
    """
    """
    _headers = ('atom', 'bond', 'angl', 'thet', 'dihe', 'phi', 'impr',
                    'imph', 'cmap', 'nbon', 'nonb', 'nbfi', 'hbon', 'end')

    sections = ('bond', 'angle', 'dihedral', 'improper', 'cmap', 'nonbond',
                'nbfix', 'hbond')

    prm_class_map = {'bond': tp.BondPRM,
                    'angle': tp.AnglePRM,
                    'dihedral': tp.DihedralPRM,
                    'improper': tp.ImproperPRM,
                    'cmap': tp.CmapPRM,
                    'nonbond': tp.NonbondPRM,
                    'nbfix': tp.NBFixPRM,
                    'hbond': tp.HBondPRM}

    def __init__(self, fname, mode='r', buffering=None):
        self.bond = None
        self.angle = None
        self.dihedral = None
        self.improper = None
        self.cmap = None
        self.nonbond = None
        self.nbfix = None
        self.hbond = None
        super(PRMFile, self).__init__(fname=fname, mode=mode,
                                    buffering=buffering)

    def parse(self):
        """
        """
        super(PRMFile, self).parse()
        while 1:
            if self.deque[0].startswith('bond'):
                self._parse_subsection('bond')
            elif self.deque[0].startswith(('angl', 'thet')):
                self._parse_subsection('angle')
            elif self.deque[0].startswith(('dihe', 'phi')):
                self._parse_subsection('dihedral')
            elif self.deque[0].startswith(('impr', 'imph')):
                self._parse_subsection('improper')
            elif self.deque[0].startswith('cmap'):
                self._parse_subsection('cmap')
            elif self.deque[0].startswith(('nbon', 'nonb')):
                self._parse_subsection('nonbond')
            elif self.deque[0].startswith('nbfi'):
                self._parse_subsection('nbfi')
            elif self.deque[0].startswith('hbon'):
                self._parse_subsection('hbond')
            else:
                break

    def _parse_subsection(self, arg):
        if arg not in self.sections:
            raise ValueError("invalid subsection specified %s" % arg)
        #
        if getattr(self, arg) is None:
            setattr(self, arg, [])
            getattr(self, arg).append(self.deque.popleft())
            while 1:
                if self.deque[0].startswith(self._headers):
                    break
                else:
                    getattr(self, arg).append(self.deque.popleft())
                if not self.deque:
                    break
        else:
            pass
            # This shouldn't happen ... yet!

    def export_to_toppar(self, toppar):
        for section in self.sections:
            if section == 'cmap':
                if self.cmap is not None:
                    toppar.cmap_opts = self.cmap[0]
                    toppar.cmap = [ tp.CmapPRM(line) for line in self.cmap[1:] ]
            else:
                opts, objs = self._export(section)
                if opts is None:
                    continue
                # TODO
                # this needs to be changed to allow multiple read ins, currently
                # second read in will nuke the first
                # TODO
                setattr(toppar, "%s_opts" % section, opts)
                setattr(toppar, section, objs)

    def _export(self, arg):
        if arg not in self.sections:
            raise ValueError("invalid subsection specified %s" % arg)
        #
        self_section = getattr(self, arg)
        if self_section is None:
            return (None, None)
        opts = self_section[0]
        cls = self.prm_class_map[arg]
        objs = [ prmobj_from_charmm(line, cls) for line in self_section[1:] ]
        return (opts, objs)

    def import_from_toppar(self, toppar):
        for section in self.sections:
            setattr(self, section, self._import(toppar, section))

    def _import(self, toppar, arg):
        if arg not in self.sections:
            raise ValueError("invalid subsection specified %s" % arg)
        #
        if getattr(toppar, arg) is None:
            return None
        opts = [ getattr(toppar, "%s_opts" % arg) ]
        prms = [ prmobj_printer(prm) for prm in getattr(toppar, arg) ]
        return opts + prms

    def pack(self):
        tmp = []
        tmp.append(self.pack_title())
        tmp.append('')
        for section in self.sections:
            if getattr(self, section) is not None:
                tmp += getattr(self, section)
                tmp.append('')
        return '\n'.join(tmp)

    def write_all(self):
        self.write(self.pack())
