"""
DOCME
"""
# fcp
# 02/22/2011


from charmming.const.etc import alphanum2num
from charmming.tools import Property, lowerKeys


def prm2int(arg):
    tmp = 0
    for i, char in enumerate(arg):
        tmp += alphanum2num[char]*36**i
    return tmp


class BasePRM(object):
    """
    DOCME
    """

    _formatting = ''
    _altFormatting = ''

    def __init__(self, arg=None):
        """
        """
        super(BasePRM, self).__init__()
        #
        if isinstance(arg, str):
            arg = arg.lower()
            try:
                arg, self.comment = arg.split('!',1)
            except ValueError:
                self.comment = ''
            self.parse(arg)
        elif arg is None:
            self._init_null()
        else:
            raise TypeError('argument should be a string')

    def parse(self, arg):
        self.body = arg.split()
        for i in range(len(self.body)):
            try:
                self.body[i] = float(self.body[i])
            except ValueError:
                pass
        self._set_hash()

    def Print(self):
        if self.body:
            if self.comment:
                tmp = self.body[:] # <- oh so very important!
                tmp.append(self.comment)
                taco = self.__class__._formatting + ' ! %s'
                try:
                    taco = taco % tuple(tmp)
                except TyprError:
                    taco = self.__class__._altFormatting + ' ! %s'
                    taco = taco % tuple(tmp)
            else:
                try:
                    taco = self.__class__._formatting % tuple(self.body)
                except TypeError:
                    taco = self.__class__._altFormatting % tuple(self.body)
            return taco.upper()
        else:
            return ''

    def _set_hash(self):
        raise NotImplementedError

    def _init_null(self):
        self.body = []
        self.comment = ''
        self._hash = 0

###################
# Special Methods #
###################

    def __repr__(self):
        try:
            tmp = self.__class__._formatting % tuple(self.body)
        except TypeError:
            tmp = self.__class__._altFormatting % tuple(self.body)
        return '%s(%s)' % (self.__class__.__name__, tmp)

    def __hash__(self):
        raise NotImplementedError

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.__hash__() == other.__hash__()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self.__hash__() < other.__hash__()

    def __le__(self, other):
        return self.__hash__() <= other.__hash__()

    def __gt__(self, other):
        return self.__hash__() > other.__hash__()

    def __ge__(self, other):
        return self.__hash__() >= other.__hash__()

    def __hash__(self):
        return self._hash


class AnglePRM(BasePRM):
    """
    """

    _formatting = '%-8s%-8s%-8s%14.6f%7.3f'
    _altFormatting = '%-8s%-8s%-8s%14.6f%7.3f%7.3f%7.3f'

    def __init__(self, arg=None):
        super(AnglePRM, self).__init__(arg)

    def _set_hash(self):
        a1, b1, a2 = self.body[:3]
        self._hash = prm2int(a1) + prm2int(a2) + 36**4 * prm2int(b1)


class BondPRM(BasePRM):
    """
    """

    _formatting = '%-8s%-8s%14.6f%12.6f'

    def __init__(self, arg=None):
        super(BondPRM, self).__init__(arg)

    def _set_hash(self):
        a1, a2 = self.body[:2]
        self._hash = prm2int(a1) + prm2int(a2)


class DihedralPRM(BasePRM):
    """
    """

    _formatting = '%-8s%-8s%-8s%-8s%14.6f%3d%12.6f'

    def __init__(self, arg=None):
        super(DihedralPRM, self).__init__(arg)

    def _set_hash(self):
        a1, b1, b2, a2 = self.body[:4]
        c1 = str(int(self.body[5]))
        self._hash = prm2int(a1) + prm2int(a2) + 36**4 * \
                (prm2int(b1) + prm2int(b2)) + 36**8 * prm2int(c1)


class ImproperPRM(DihedralPRM):
    """
    """
    def __init__(self, arg=None):
        super(ImproperPRM, self).__init__(arg)


class MassPRM(BasePRM):
    """
    """

    _formatting = 'mass %-5d%-8s%14.6f'

    def __init__(self, arg=None):
        super(MassPRM, self).__init__(arg)

    def parse(self, arg):
        self.body = arg.split()[2:]
        for i in range(len(self.body)):
            try:
                self.body[i] = float(self.body[i])
            except ValueError:
                pass
        self._set_hash()
        self.index = 0

    def Print(self):
        if self.body:
            tmp = [self.index, self.body[0], self.body[1]]
            if self.comment:
                tmp.append(self.comment)
                taco = self.__class__._formatting + ' ! %s'
                taco = taco % tuple(tmp)
            else:
                taco = self.__class__._formatting % tuple(tmp)
            return taco.upper()
        else:
            return ''

    def _init_null(self):
        super(MassPRM, self)._init_null()
        self.index = 0

    def _set_hash(self):
        a1 = self.body[0]
        self._hash = prm2int(a1)


class NonBondPRM(BasePRM):
    """
    """

    _formatting = '%-8s%5.1f%10s%12.6f'
    _altFormatting = '%-8s%5.1f%10s%12.6f%12.6f%12.6f%12.6f'

    def __init__(self, arg=None):
        super(NonBondPRM, self).__init__(arg)

    def _set_hash(self):
        a1 = self.body[0]
        self._hash = prm2int(a1)


class NBFixPRM(BondPRM):
    """
    """

    _formatting = '%-8s%-8s%14.6f%12.6f'
    _altFormatting = '%-8s%-8s%14.6f%12.6f%12.6f%12.6f'

    def __init__(self, arg=None):
        super(NBFixPRM, self).__init__(arg)


class Residue(object):
    """
    TODO
    """
    def __init__(self, arg):
        super(Residue, self).__init__()


"""
\\For giggles...

def PRMFactory(hashfunc, formatting=''):
    tmp = BasePRM
    tmp._formatting = formatting
    tmp._set_hash = hashfunc
    return tmp


formatting = {
    'angle': '%-8s%-8s%-8s%14.6f%12.6f',
    'bond': '%-8s%-8s%14.6f%12.6f',
    'dihedral': '%-8s%-8s%-8s%-8s%14.6f%3d%12.6f',
    # improper -> dihedral
    'nonbond': '%-8s%5.1f%10s%12.6f'
    # nbfix -> bond
    }


def _hashing():
    def angle(self):
        a1, b1, a2 = self.body[:3]
        self._hash = prm2int(a1) + prm2int(a2) + 36**4 * prm2int(b1)

    def bond(self):
        a1, a2 = self.body[:2]
        self._hash = prm2int(a1) + prm2int(a2)

    def dihedral(self):
        a1, b1, b2, a2 = self.body[:4]
        c1 = str(int(self.body[5]))
        self._hash = prm2int(a1) + prm2int(a2) + 36**4 * \
                (prm2int(b1) + prm2int(b2)) + 36**8 * prm2int(c1)

    # improper -> dihedral

    def nonbond(self):
        a1 = self.body[0]
        self._hash = prm2int(a1)

    # nbfix -> bond

    return locals()
hashing = _hashing()

\\An example:
Angle = PRMFactory(hashing['angle'], formatting['angle'])

\\Should be logically equivalent to AnglePRM
"""
