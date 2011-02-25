"""
DOCME
"""
# fcp
# 02/22/2011


from charmming.tools import Property, expandPath, lowerKeys, logicalLines


class BaseCHARMMFile(object):
    """
    DOCME
    """

    _properties = {
        'header': [],
        'body': []
    }

    def __init__(self, arg=None, **kwargs):
        """
        DOCME
        """
        super(BaseCHARMMFile, self).__init__()
        # kwargs
        kwargs = lowerKeys(kwargs)
        #
        if arg is None:
            self._init_null()
        elif isinstance(arg, str):
            self.read(open(expandPath(arg)))
        else:
            self.read(arg)

##############
# Properties #
##############

    @Property
    def header():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._header
        def fset(self, value):
            self._header = value
        return locals()

    @Property
    def body():
        doc =\
        """
        DOCME
        """
        def fget(self):
            return self._body
        def fset(self, value):
            self._body = value
        return locals()

##################
# Public Methods #
##################

    def read(self, iterable):
        """
        DOCME
        """
        self.header = []
        self.body = []
        iterable = ( line.strip() for line in iterable )
        iterable = ( line for line in iterable if line )
        iterable = logicalLines(iterable)
        #
        start = False
        for line in iterable:
            if start:
                if line.startswith('*'):
                    self.header.append(line)
                else:
                    self.body.append(line.lower())
            else:
                if line.startswith('*'):
                    self.header.append(line)
                    start = True

    def write(self, filename):
        """
        DOCME
        """
        NotImplementedError

###################
# Private Methods #
###################

    def _init_null(self):
        """
        """
        for key, value in self.__class__._properties.iteritems():
            setattr(self, key, value)

