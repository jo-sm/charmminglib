"""This modedule contains CHARMM specific readers and writers to be used by
.rtf and .prm file objects.
"""

from pychm.future.lib import toppar as tp
from pychm.future.io.charmm import COMMENT_CHAR as CC


# Exceptions
class ReaderError(Exception):
    pass


class WriterError(Exception):
    pass


# Convenience Functions
def _base_reader(arg):
    if not isinstance(arg, basestring):
        raise TypeError("Invalid arg: %r" % arg)
    if CC in arg:
        arg, comment = arg.split(CC, 1)
        arg = arg.strip()
        comment = comment.strip()
    else:
        comment = None
    return (arg.split(), comment)


# Readers
def bond_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom0 = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom0 information")
    try:
        atom1 = str(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom1 information")
    try:
        k = float(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        eq = float(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read eq information")
    return tp.BondPRM(atom0=atom0, atom1=atom1, k=k, eq=eq, comment=comment)


def angle_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom0 = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom0 information")
    try:
        atom1 = str(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom1 information")
    try:
        atom2 = str(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom2 information")
    try:
        k = float(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        eq = float(args[4])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read eq information")
    try:
        k13 = float(args[5])
        eq13 = float(args[6])
    except:
        k13 = None
        eq13 = None
    return tp.AnglePRM(atom0=atom0, atom1=atom1, atom2=atom2, k=k, eq=eq,
                    k13=k13, eq13=eq13, comment=comment)


def dihedral_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom0 = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom0 information")
    try:
        atom1 = str(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom1 information")
    try:
        atom2 = str(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom2 information")
    try:
        atom3 = str(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom3 information")
    try:
        k = float(args[4])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        mult = int(args[5])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read mult information")
    try:
        eq = float(args[6])
    except:
        # add warning here
        eq = 0.0
    return tp.DihedralPRM(atom0=atom0, atom1=atom1, atom2=atom2, atom3=atom3,
                        k=k, mult=mult, eq=eq, comment=comment)


def improper_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom0 = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom0 information")
    try:
        atom1 = str(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom1 information")
    try:
        atom2 = str(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom2 information")
    try:
        atom3 = str(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom3 information")
    try:
        k = float(args[4])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        mult = int(args[5])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read mult information")
    try:
        eq = float(args[6])
    except:
        # add warning here
        eq = 0.0
    return tp.ImproperPRM(atom0=atom0, atom1=atom1, atom2=atom2, atom3=atom3,
                        k=k, mult=mult, eq=eq, comment=comment)


def cmap_reader(arg):
    args, comment = _base_reader(arg)
    tmp = tp.CmapPRM(*args)
    tmp.comment = comment
    return tmp


def nonbond_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom information")
    try:
        ig = float(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read ig information")
    try:
        k = float(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        eq = float(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read eq information")
    try:
        ig14 = float(args[4])
        k14 = float(args[5])
        eq14 = float(args[6])
    except:
        ig14 = None
        k14 = None
        eq14 = None
    return tp.NonbondPRM(atom=atom, ig=ig, k=k, eq=eq, ig14=ig14, k14=k14,
                        eq14=eq14, comment=comment)


def nbfix_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom0 = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom0 information")
    try:
        atom1 = str(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom1 information")
    try:
        k = float(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        eq = float(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read eq information")
    try:
        k14 = float(args[4])
        eq14 = float(args[5])
    except:
        k14 = None
        eq14 = None
    return tp.NBFixPRM(atom0=atom0, atom1=atom1, k=k, eq=eq, k14=k14,
                    eq14=eq14, comment=comment)


def hbond_reader(arg):
    args, comment = _base_reader(arg)
    try:
        atom0 = str(args[0])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom0 information")
    try:
        atom1 = str(args[1])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom1 information")
    try:
        k = float(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read k information")
    try:
        eq = float(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read eq information")
    return tp.HBondPRM(atom0=atom0, atom1=atom1, k=k, eq=eq, comment=comment)


def mass_reader(arg):
    args, comment = _base_reader(arg)
    try:
        id = int(args[1])
    except:
        id = 0
    try:
        atom = str(args[2])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read atom information")
    try:
        mass = float(args[3])
    except:
        print 'On line, "%s"' % arg
        raise ReaderError("Can't read mass information")
    try:
        element = str(args[4])
    except:
        element = None
    return tp.Mass(id=id, atom=atom, mass=mass, element=element, comment=comment)


def residue_reader(arg):
    args, comment = _base_reader(arg[0])
    try:
        name = args[1]
    except:
        print 'On line, "%s"' % arg[0]
        raise ReaderError("Can't read resi name information")
    try:
        charge = float(args[2])
    except:
        print 'On line, "%s"' % arg[0]
        raise ReaderError("Can't read resi charge information")
    return tp.Residue(name=name, charge=charge, body=arg[1:], comment=comment)


def patch_reader(arg):
    args, comment = _base_reader(arg[0])
    try:
        name = args[1]
    except:
        print 'On line, "%s"' % arg[0]
        raise ReaderError("Can't read pres name information")
    try:
        charge = float(args[2])
    except:
        print 'On line, "%s"' % arg[0]
        raise ReaderError("Can't read pres charge information")
    return tp.Patch(name=name, charge=charge, body=arg[1:], comment=comment)


# Writers
def bond_writer(arg):
    if not isinstance(arg, tp.BondPRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom0)
    except:
        raise WriterError("Can't write atom0 information")
    try:
        tmp.append('%-8s' % arg.atom1)
    except:
        raise WriterError("Can't write atom1 information")
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%12.6f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def angle_writer(arg):
    if not isinstance(arg, tp.AnglePRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom0)
    except:
        raise WriterError("Can't write atom0 information")
    try:
        tmp.append('%-8s' % arg.atom1)
    except:
        raise WriterError("Can't write atom1 information")
    try:
        tmp.append('%-8s' % arg.atom2)
    except:
        raise WriterError("Can't write atom2 information")
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%9.3f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.k13 is not None and arg.eq13 is not None:
        tmp.append('%14.6f' % arg.k13)
        tmp.append('%9.3f' % arg.eq13)
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def dihedral_writer(arg):
    if not isinstance(arg, tp.DihedralPRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom0)
    except:
        raise WriterError("Can't write atom0 information")
    try:
        tmp.append('%-8s' % arg.atom1)
    except:
        raise WriterError("Can't write atom1 information")
    try:
        tmp.append('%-8s' % arg.atom2)
    except:
        raise WriterError("Can't write atom2 information")
    try:
        tmp.append('%-8s' % arg.atom3)
    except:
        raise WriterError("Can't write atom3 information")
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%3d' % arg.mult)
    except:
        raise WriterError("Can't write mult information")
    try:
        tmp.append('%9.3f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def improper_writer(arg):
    if not isinstance(arg, tp.ImproperPRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom0)
    except:
        raise WriterError("Can't write atom0 information")
    try:
        tmp.append('%-8s' % arg.atom1)
    except:
        raise WriterError("Can't write atom1 information")
    try:
        tmp.append('%-8s' % arg.atom2)
    except:
        raise WriterError("Can't write atom2 information")
    try:
        tmp.append('%-8s' % arg.atom3)
    except:
        raise WriterError("Can't write atom3 information")
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%3d' % arg.mult)
    except:
        raise WriterError("Can't write mult information")
    try:
        tmp.append('%9.3f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def cmap_writer(arg):
    if not isinstance(arg, tp.CmapPRM):
        raise TypeError("Invalid arg: %r" % arg)
    return arg.text


def nonbond_writer(arg):
    if not isinstance(arg, tp.NonbondPRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom)
    except:
        raise WriterError("Can't write atom information")
    try:
        tmp.append('%5.1f' % arg.ig)
    except:
        tmp.append('  0.0')
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%12.6f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.k14 is not None and arg.eq14 is not None:
        tmp.append('  0.0')
        tmp.append('%14.6f' % arg.k14)
        tmp.append('%12.6f' % arg.eq14)
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def nbfix_writer(arg):
    if not isinstance(arg, tp.NBFixPRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom0)
    except:
        raise WriterError("Can't write atom0 information")
    try:
        tmp.append('%-8s' % arg.atom1)
    except:
        raise WriterError("Can't write atom1 information")
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%12.6f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.k14 is not None and arg.eq14 is not None:
        tmp.append('%14.6f' % arg.k14)
        tmp.append('%12.6f' % arg.eq14)
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def hbond_writer(arg):
    if not isinstance(arg, tp.HBondPRM):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = []
    try:
        tmp.append('%-8s' % arg.atom0)
    except:
        raise WriterError("Can't write atom0 information")
    try:
        tmp.append('%-8s' % arg.atom1)
    except:
        raise WriterError("Can't write atom1 information")
    try:
        tmp.append('%14.6f' % arg.k)
    except:
        raise WriterError("Can't write k information")
    try:
        tmp.append('%12.6f' % arg.eq)
    except:
        raise WriterError("Can't write eq information")
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def mass_writer(arg):
    if not isinstance(arg, tp.Mass):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = ['mass ']
    try:
        tmp.append('%-5d' % arg.id)
    except:
        raise WriterError("Can't write id information")
    try:
        tmp.append('%-8s' % arg.atom)
    except:
        raise WriterError("Can't write atom information")
    try:
        tmp.append('%14.6f' % arg.mass)
    except:
        raise WriterError("Can't write mass information")
    # TODO should we put element writer here?
    if arg.comment is not None:
        tmp.append(' ! %s' % arg.comment)
    return ''.join(tmp)


def residue_writer(arg):
    if not isinstance(arg, tp.Residue):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = 'resi %6s%8.2f' % (arg.name, arg.charge)
    if arg.comment is not None:
        tmp += ' ! %s' % arg.comment
    tmp = [tmp, ]
    tmp.extend(arg.body)
    tmp.append('')
    return '\n'.join(tmp)


def patch_writer(arg):
    if not isinstance(arg, tp.Patch):
        raise TypeError("Invalid arg: %r" % arg)
    tmp = 'pres %6s%8.2f' % (arg.name, arg.charge)
    if arg.comment is not None:
        tmp += ' ! %s' % arg.comment
    tmp = [tmp, ]
    tmp.extend(arg.body)
    tmp.append('')
    return '\n'.join(tmp)


