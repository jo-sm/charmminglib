#!/usr/bin/env python


import pychm
import emap


taco = pychm.io.pdb.PDBFile('1A7N.pdb')[0]
molA = pychm.lib.mol.Mol(taco.iter_seg(chainids=['l'], segtypes=['pro']).next())
molA.parse()
rtfA = pychm.io.rtf.RTFFile('~/bigbox/toppar/top_all27_prot_na.rtf')
molA.populate_charges(rtfA)
lastRes = pychm.lib.basestruct.BaseStruct(list(molA.iter_res())[-1])
lastRes.find(atomtype=' c  ')[0].charge =  0.34
lastRes.find(atomtype=' ot1')[0].charge = -0.67
lastRes.find(atomtype=' ot2')[0].charge = -0.67

mapA = emap.EMap(mol=molA)
