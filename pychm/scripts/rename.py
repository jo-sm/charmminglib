#!/usr/bin/env python

import os

class walk(object):
    """
    A forward iterator that traverses a directory tree, and returns strings,
    one for each file found.
    """
    def __init__(self, directory):
        self.stack = [directory]
        self.files = []
        self.index = 0

    def __getitem__(self, index):
        while 1:
            try:
                file = self.files[self.index]
                self.index += 1
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

iterator = ( filename for filename in walk('.') if filename.endswith('.py') )
for filename in iterator:
    os.popen("sed -i 's/pychm/pychm/g' %s" % filename)
