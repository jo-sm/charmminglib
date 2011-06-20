"""
:Author: fcp
:Date: 06/14/2011
"""


from itertools import chain
from pychm.tools import chomp


class BaseCHARMMFile(object):
    """
    **STUB**
    """
    def __init__(self, filename, **kwargs):
        super(BaseCHARMMFile, self).__init__()
        self.filename = filename

    def parse(self):
        """
        """
        self.header = []
        self.body = []
        iterator = ( line.strip() for line in open(self.filename) )
        iterator = ( line for line in iterator if line )
        for line in iterator:
            if line.startswith('*'):
                self.header.append(line[1:].lower())
            else:
                iterator = chain([line], iterator) # rewind by one
                break
        for line in iterator:
            self.body.append(line.lower())

