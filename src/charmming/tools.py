"""
This module contains a 'random' assortment of useful utility functions
and classes.  If a utility function is closely associated with a single
particular module, it should be located in said module; otherwise, it
should be put here for general use.

TODO:
    The optparse module is depricated starting with python 2.6, and so
    we should probably move to the newer argparse.
"""
# fcp
# 10/22/2010


import optparse
import os


def paragraphs(iterable,splitter):
    """
    Cut a text stream up into 'paragraphs,' where partitions are
    determined by a list named splitter.
    """
    assert isinstance(splitter, (tuple, list))
    splitter = tuple(splitter)
    paragraph = []
    for line in iterable:
        if line.startswith(splitter):
            if paragraph:
                yield paragraph
            paragraph = [line]
        else:
            paragraph.append(line)
    if paragraph:
        yield paragraph


def expandPath(String):
    """
    A better version of the default path expander.
    """
    if '~' in String:
        return os.path.expanduser(String)
    else:
        return os.path.abspath(String)


def cleanStrings(iterable,CC=None):
    """
    Takes an iterable of strings, strips out comments, blank lines,
    forces lower case and left justifies.
    """
    if CC is not None:
        iterable = ( line.split(CC)[0] for line in iterable )
    iterable = ( line.strip() for line in iterable if line )
    iterable = ( line.lower() for line in iterable )
    return iterable


def flatten(l,ltypes=(list,tuple)):
    """
    Takes a nested list of arbitrary depth, and returns a flattened
    one. This function will only flatten types specified by `ltypes`.
    """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)


def Property(func):
    """
    Decorator function for property. Prevents namespace pollution in
    classes.

    >>> @Property
    >>> def foo():
    >>>     doc = "The foo property."
    >>>     def fget(self):
    >>>         return self._foo
    >>>     def fset(self, value):
    >>>         self._foo = value
    >>>     return locals()
    """
    return property(**func())


def mkdir(value):
    """
    Recursively make a directory at specified path. The following
    example will first try to make '~', then '~/python', etc.

    >>> mkdir('~/python/projects/taco/sauce')
    """
    value = expandPath(value).split(os.sep)[1:]
    tmp = os.sep
    for entry in value:
        tmp += '%s%s' % (entry,os.sep)
        try:
            os.mkdir(tmp)
        except OSError: pass


def out2inp(iterable,lookFor='CHARMM>',CC='!',maxDrought=1000):
    """
    Parse through a CHARMM .out file, and reconstruct the corresponding
    .inp file.

    WARNING -- This can be potentially (and quietly) very bad if you
    use it on .out files that have been concatenated!
    """
    n = 0
    for line in iterable:
        line = line.strip()
        if line.startswith(lookFor):
            line = line.split(lookFor)[1]
            line = line.split(CC)[0]
            if line:
                yield line
            n = 0
        else:
            n += 1
        if n > maxDrought:
            return


def logicalLines(iterable,continueChar='-'):
    """
    Convert an iterable of physical lines with `continuationChar` into
    an iterator of logical lines.
    """
    tmp = ''
    for line in iterable:
        line = line.strip()
        if line.endswith(continueChar):
            line = line[:-len(continueChar)]
            tmp += line
        else:
            if tmp:
                yield tmp
                tmp = ''
            yield line


def get_inpProp(prop,iterable):
    """
    Parse through the iterable, looking for string 'prop'.  If it is
    found, return the word directly following prop in the same line.
    """
    for line in cleanStrings(iterable,CC='!'):
        llist = line.split()
        try:
            propIndex = index(prop,llist) + 1
        except ValueError:
            pass
        else:
            try:
                return float(llist[propIndex])
            except:
                raise ValueError


def index(predicate,iterable):
    """
    Like the builtin string method 'index' only instead of looking for
    an exact match, this only looks for the first instance of
    'startswith'.  If it fails, then a ValueError is raised.
    """
    for i,value in enumerate(iterable):
        if value.startswith(predicate):
            return i
    raise ValueError


class walk(object):
    """
    A forward iterator that traverses a directory tree.
    """
    def __init__(self, directory):
        self.stack = [directory]
        self.files = []
        self.index = 0

    def __getitem__(self, index):
        while 1:
            try:
                file = self.files[self.index]
                self.index = self.index + 1
            except IndexError:
                # pop next directory from stack
                self.directory = self.stack.pop()
                self.files = os.listdir(self.directory)
                self.index = 0
            else:
                # got a filename
                fullname = os.path.join(self.directory, file)
                if os.path.isdir(fullname) and not os.path.islink(fullname):
                    self.stack.append(fullname)
                return fullname


class OptionWithDefault(optparse.Option):
    """
    Exension of optparse to include 'required options.'  Based on code found at
    http://code.activestate.com/recipes/573441/
    """
    ATTRS = optparse.Option.ATTRS + ['required']

    def __init__(self, *opts, **attrs):
        if attrs.get('required', False):
            attrs['help'] = '(Required) ' + attrs.get('help', "")
        optparse.Option.__init__(self, *opts, **attrs)


class OptionParser(optparse.OptionParser):
    """
    Exension of optparse to include 'required options.'  Based on code found at
    http://code.activestate.com/recipes/573441/

    >>> parser = OptionParser( ... )
    >>> parser.add_option('-I','--input',required=True,metavar='PATH',
                        help='PATH of input .pdb file')
    """
    def __init__(self, **kwargs):
        kwargs['option_class'] = OptionWithDefault
        optparse.OptionParser.__init__(self, **kwargs)

    def check_values(self, values, args):
        for option in self.option_list:
            if hasattr(option, 'required') and option.required:
                if not getattr(values, option.dest):
                    self.error("option %s is required" % (str(option)))
        return optparse.OptionParser.check_values(self, values, args)

def lowerKeys(dictionary):
    """
    Modifies a dictionary by applying `string.lower` to each of its keys.

    This is helpful for making kwargs case insensitive.
    """
    return dict(((key.lower(), value) for key, value in dictionary.iteritems()))
