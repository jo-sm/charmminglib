"""
:Author: fcp
:Date: 06/14/2011
"""


from itertools import chain
from pychm.tools import chomp, expandPath


class BaseCHARMMFile(object):
    """
    **STUB**
    """
    def __init__(self, filename, **kwargs):
        super(BaseCHARMMFile, self).__init__()
        self.filename = expandPath(filename)
        self.header = []
        self.body = []

    def parse(self):
        """
        """
        iterator = ( line.strip() for line in open(self.filename) )
        iterator = ( line.lower() for line in iterator if line )
        for line in iterator:
            if line.startswith('*'):
                self.header.append(line[1:])
            else:
                iterator = chain([line], iterator) # rewind by one
                break
        for line in iterator:
            self.body.append(line)

