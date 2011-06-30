"""
:Author: fcp
:Date: 02/22/2011
"""


from pychm.const import alphanum2num
from pychm.tools import Property, lowerKeys, cleanStrings, paragraphs, flatten


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
        super(BasePRM, self).__init__()
        #
        if isinstance(arg, str):
            arg = arg.lower()
            tmp = arg.split('!')
            arg = tmp[0]
            try:
                self.comment = '!'.join(tmp[1:])
            except IndexError:
                self.comment = ''
            self.parse(arg)
        elif arg is None:
            self._init_null()
        else:
            raise TypeError('argument should be a string')

    def is_wild(self):
        return 'x' in self.body

    def parse(self, arg):
        self.body = arg.split()
        for i in range(len(self.body)):
            try:
                self.body[i] = float(self.body[i])
            except ValueError:
                pass

    def Print(self):
        if self.body:
            if self.comment:
                tmp = self.body[:] # <- oh so very important!
                tmp.append(self.comment)
                taco = self.__class__._formatting + ' ! %s'
                try:
                    taco = taco % tuple(tmp)
                except TypeError:
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

    def _set_sort(self):
        raise NotImplementedError

    def _init_null(self):
        self.body = []
        self.comment = ''
        self._sort = 0

###################
# Special Methods #
###################

    def __repr__(self):
        try:
            tmp = self.__class__._formatting % tuple(self.body)
        except TypeError:
            tmp = self.__class__._altFormatting % tuple(self.body)
        return '%s(%s)' % (self.__class__.__name__, tmp)

    # STUB
    def __eq__(self, other):
        raise NotImplementedError

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return self._sort < other._sort

    def __le__(self, other):
        return self._sort <= other._sort

    def __gt__(self, other):
        return self._sort > other._sort

    def __ge__(self, other):
        return self._sort >= other._sort

    def __hash__(self):
        return hash(self._sort)


class AnglePRM(BasePRM):
    """
    """

    _formatting = '%-8s%-8s%-8s%14.6f%7.3f'
    _altFormatting = '%-8s%-8s%-8s%14.6f%7.3f%7.3f%7.3f'

    def __init__(self, arg=None):
        super(AnglePRM, self).__init__(arg)
        if not self.body[0] < self.body[2]:
            self.body[0], self.body[2] = self.body[2], self.body[0]
        self._set_sort()

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.body[:3] == other.body[:3]

    def _set_sort(self):
        self._sort = '%-8s%-8s%-8s' % tuple(self.body[:3])


class BondPRM(BasePRM):
    """
    """

    _formatting = '%-8s%-8s%14.6f%12.6f'

    def __init__(self, arg=None):
        super(BondPRM, self).__init__(arg)
        if not self.body[0] < self.body[1]:
            self.body[0], self.body[1] = self.body[1], self.body[0]
        self._set_sort()

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.body[:2] == other.body[:2]

    def _set_sort(self):
        self._sort = '%-8s%-8s' % tuple(self.body[:2])


class DihedralPRM(BasePRM):
    """
    """

    _formatting = '%-8s%-8s%-8s%-8s%14.6f%3d%12.6f'

    def __init__(self, arg=None):
        super(DihedralPRM, self).__init__(arg)
        if not self.body[0] < self.body[3]:
            self.body[0], self.body[3] = self.body[3], self.body[0]
            self.body[1], self.body[2] = self.body[2], self.body[1]
        self._set_sort()

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.body[:4] == other.body[:4] and self.body[5] == other.body[5]

    def _set_sort(self):
        tmp = self.body[:6]
        tmp.pop(4)
        self._sort = '%-8s%-8s%-8s%-8s%2s' % tuple(map(str,tmp[:5]))


class HBondPRM(BasePRM):
    """
    """
    def __init__(self, arg=None):
        super(HBondPRM, self).__init__(arg)
        raise NotImplementedError
        self._set_sort()

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
        self._set_sort()

    def parse(self, arg):
        self.body = arg.split()[2:]
        for i in range(len(self.body)):
            try:
                self.body[i] = float(self.body[i])
            except ValueError:
                pass
        self._set_sort()
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
            return taco
        else:
            return ''

    def _init_null(self):
        super(MassPRM, self)._init_null()
        self.index = 0

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.body[0] == other.body[0]

    def __repr__(self):
        tmp = self.__class__._formatting % tuple([0] + self.body)
        return '%s(%s)' % (self.__class__.__name__, tmp)

    def _set_sort(self):
        self._sort = '%-8s' % self.body[0]


class NonBondPRM(BasePRM):
    """
    """

    _formatting = '%-8s%5.1f%10s%12.6f'
    _altFormatting = '%-8s%5.1f%10s%12.6f%12.6f%12.6f%12.6f'

    def __init__(self, arg=None):
        super(NonBondPRM, self).__init__(arg)
        self._set_sort()

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        return self.body[0] == other.body[0]

    def _set_sort(self):
        self._sort = '%-8s' % self.body[0]


class NBFixPRM(BondPRM):
    """
    """

    _formatting = '%-8s%-8s%14.6f%12.6f'
    _altFormatting = '%-8s%-8s%14.6f%12.6f%12.6f%12.6f'

    def __init__(self, arg=None):
        super(NBFixPRM, self).__init__(arg)


class TopRes(object):
    """
    **TODO**
    """
    def __init__(self, arg):
        super(TopRes, self).__init__()
        tmp1 = []
        tmp2 = []
        for line in cleanStrings(arg, CC='!'):
            if line.startswith(('group', 'atom')):
                tmp1.append(line) # atom data
            else:
                tmp2.append(line) # everything else
        # turn tmp1 into a list
        self.groups = []
        for group in paragraphs(tmp1, ['grou']):
            self.groups.append([ tuple(line.split()[1:]) for line in group[1:] ])
        self.atoms = flatten(self.groups, ltypes=list)
        # turn tmp2 into a dict
        meta = {}
        for line in tmp2:
            key = line.split()[0]
            value = line.split(key)[1].lstrip()
            if key not in meta.keys():
                meta[key] = [value]
            else:
                meta[key].append(value)
        self.meta = meta
        # dicts
        self.chemDict = dict(( line[:2] for line in self.atoms ))
        self.chargeDict = dict(( [line[0], float(line[2])] for line in self.atoms ))



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
