.. _atoms:

Atom Modules
============

Collectively, the atom modules comprise one-half of the workhorse modules of
charmminglib (the other half being their associated container objects).  The
:class:`Atom` class is the only operative class of the three.  The :class:`BaseAtom`
and poorly named :class:`MetaAtom` classes both serve to abstract away common
elements to other :class:`Atom`-like objects such as the :class:`cgAtom` class
and any potential future QM implementations.

:class:`Atom`-like instances are typically instantized using a string from either
CHARMM output or from the PDB.  However null :class:`Atom` instances may also
be created for debugging purposes, or using any object that descents from
:class:`MetaAtom` provided the appropriate class variables are defined.  For more
information, please see the :mod:`charmming.lib.metaatom` documentation.

Currently, :class:`Atom`-like objects are hashable, despite the fact they are
mutable, a big python no-no.  This was originally done so that different container
objects like :class:`OrderedDict` and :class:`OrderedSet` could be considered
as container classes (we currently use a :class:`list` based object), however
because hashes are created at instanization, and are used to determine uniqueness,
putting :class:`Atom`-like objects into :class:`dict` and :class:`set` objects
can have unintended consequences, as seemingly different objects collide.

:class:`Atom` objects support rich comparisons, which can be customized by
tweaking the :meth:`_sort` method, for more information see the
:mod:`charmming.lib.baseatom` documentation.

Finally, the major purpose of the :class:`Atom` class is to serve as an easy store
of numerical coordinate data, and biochemical metadata.  Numerical coordinate data
is stored using :class:`numpy.array` objects.  This allows a great speed up
for doing algebraic manipulations versus using pure python.  Metadata is
gleaned from text strings using the :meth:`parse` method, and written to strings
using the :meth:`Print` method, for more information see the
:mod:`charmming.lib.atom` documentation.

**Contents:**

.. toctree::
  :maxdepth: 2
  
  atom
  baseatom
  metaatom
