

from charmming.tools import Property


class Contact(object):
    """
    docstring for Contact
    """
    __index = 0
    nativeRad = 7.10
    def __init__(self, scI, scJ):
        if scI.addr < scJ.addr:
            self._scI = scI
            self._scJ = scJ
        else:
            self._scI = scJ
            self._scJ = scI

        self._index = Contact.__index
        Contact.__index += 1

    @Property
    def index():
        doc = "The index property."
        def fget(self):
            return self._index
        return locals()

    @Property
    def residI():
        doc = "The residI property."
        def fget(self):
            return self._scI.resid
    #       return self._scI.resIndex
        return locals()

    @Property
    def residJ():
        doc = "The residJ property."
        def fget(self):
            return self._scJ.resid
    #       return self._scJ.resIndex
        return locals()

    @Property
    def chainidI():
        doc = "The chainidI property."
        def fget(self):
            return self._scI.chainid
        return locals()

    @Property
    def chainidJ():
        doc = "The chainidJ property."
        def fget(self):
            return self._scJ.chainid
        return locals()

    @Property
    def isNative():
        doc = "The isNative property."
        def fget(self):
            if Contact.nativeRad < 0:
                raise ValueError('isNative: nativeRad %5.2f, should be greater than zero.')
            else:
                return self._scI.calc_length(self.scJ) < Contact.nativeRad
        return locals()


#class ContactSet(OrderedSet):
#    """
#    docstring for ContactSet
#    """
#    def __init__(self,aaCrdFile,contactRad):
#        pass
#        # parse(aaCrdFile)
#        # pdb = aaCrdFile.get_model(0)
#        # ktModel = ktgo(pdb)
