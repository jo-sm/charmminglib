"""
A collection of useful biological and chemical constants.

:Attributes:
    | ``aaAlphabet`` - A :class:`dict` from the 3-letter residue
        abbreviation to the 1-letter residue abbreviation.
    | ``aaMass`` - A :class:`dict` from the 3-letter residue abbreviation
        to its mass in AMU.
    | ``aaVDV`` - A :class:`dict` from the 3-letter residue abbreviation to
        the residue's van der Waals radius in Angstroms.
    | ``atomMass`` - A :class:`dict` from the element's abbreviation to its
        mass in AMU.
    | ``backBone`` - A :class:`set` of strings representing the
        atomType of atoms present in a protein back bone.
    | ``charmm2pdbAtomNames`` - A :class:`dict` from CHARMM element abbreviations
        to standard element abbreviations.
    | ``dna`` - A :class:`set` of strings representing the resName of
        residues present in DNA.
    | ``good`` - A :class:`set` of strings represenitng the atomType of
        hetero atoms recognized internally by CHARMM.
    | ``nuc`` - A :class:`set` of strings representing the resName of
        residues present in generic nucleic acids.
    | ``pro`` - A :class:`set` of the 3-letter residue abbreviations
        for residues present in proteins.
    | ``rna`` - A :class:`set` of strings representing the resName of
        residues present in RNA.
"""

# Atom.resName
aaAlphabet = {
    'ala':'a',
    'arg':'r',
    'asn':'n',
    'asp':'d',
    'cys':'c',
    'gln':'q',
    'glu':'e',
    'gly':'g',
    'his':'h',
    'hsd':'h',
    'ile':'i',
    'leu':'l',
    'lys':'k',
    'met':'m',
    'phe':'f',
    'pro':'p',
    'sec':'u',
    'ser':'s',
    'thr':'t',
    'trp':'w',
    'tyr':'y',
    'val':'v'
    }


# Atom.resName
aaMass = {
    'ala': 71.079,
    'arg':156.188,
    'asn':114.104,
    'asp':115.089,
    'cys':103.145,
    'gln':128.131,
    'glu':129.116,
    'gly': 57.052,
    'his':137.141,
    'hsd':137.141,
    'hse':137.141,
    'hsp':137.141,
    'ile':113.160,
    'leu':113.160,
    'lys':128.17 ,
    'met':131.199,
    'phe':147.117,
    'pro': 97.117,
    'ser': 87.078,
    'thr':101.105,
    'trp':186.213,
    'tyr':163.176,
    'val': 99.133
    }


# Atom.resName
aaVDW = {
    'ala':2.51958406732374,
    'arg':3.28138980397453,
    'asn':2.84049696898525,
    'asp':2.79030096923572,
    'cys':2.73823091624513,
    'gln':3.00796101305807,
    'glu':2.96332591119925,
    'gly':2.25450393833984,
    'his':3.04273820988499,
    'hsd':3.04273820988499,
    'hse':3.04273820988499,
    'hsp':3.04273820988499,
    'ile':3.09345983013354,
    'leu':3.09345983013354,
    'lys':3.18235414984794,
    'met':3.09345983013354,
    'phe':3.18235414984794,
    'pro':2.78004241717965,
    'ser':2.59265585208464,
    'thr':2.81059478021734,
    'trp':3.38869998431408,
    'tyr':3.22881842919248,
    'val':2.92662460060742
    }


# Atom.element
atomMass = {
     'h':  1.0079,'he':  4.0026,'li':  6.941 ,'be':  9.0122, 'b': 10.811 ,
     'c': 12.0107, 'n': 14.0067, 'o': 15.9994, 'f': 18.9984,'ne': 20.1797,
    'na': 22.9897,'mg': 24.3050,'al': 26.9815,'si': 28.0855, 'p': 30.9738,
     's': 32.065 ,'cl': 35.453 ,'ar': 39.948 , 'k': 39.098 ,'ca': 40.078 ,
    'sc': 44.9559,'ti': 47.867 , 'v': 50.9415,'cr': 51.9961,'mn': 54.9380,
    'fe': 55.845 ,'co': 58.9332,'ni': 58.6934,'cu': 63.546 ,'zn': 65.38  ,
    'ga': 69.723 ,'ge': 72.64  ,'as': 74.9216,'se': 78.96  ,'br': 79.904 ,
    'kr': 83.798 ,'rb': 85.4678,'sr': 87.62  , 'y': 88.9059,'zr': 91.224 ,
    'nb': 92.9064,'mo': 95.94  ,'tc': 98.    ,'ru':101.07  ,'rh':102.9055,
    'pd':106.42  ,'ag':107.8682,'cd':112.411 ,'in':114.818 ,'sn':118.710 ,
    'sb':121.760 ,'te':127.60  , 'i':126.9045,'xe':131.293 ,'cs':132.905 ,
    'ba':137.327 ,'la':138.9055,'ce':140.116 ,'pr':140.9077,'nd':144.24  ,
    'pm':145.    ,'sm':150.36  ,'eu':151.964 ,'gd':157.25  ,'tb':158.9253,
    'dy':162.500 ,'ho':164.9303,'er':167.259 ,'tm':168.9342,'yb':173.04  ,
    'lu':174.967 ,'hf':178.49  ,'ta':180.9479, 'w':183.84  ,'re':186.207 ,
    'os':190.23  ,'ir':192.217 ,'pt':195.078 ,'au':196.9666,'hg':200.59  ,
    'tl':204.3833,'pb':207.2   ,'bi':208.9804,'po':209.    ,'at':210.    ,
    'rn':222.    ,'fr':223.    ,'ra':226.    ,'ac':227.    ,'th':232.0381,
    'pa':231.0359, 'u':238.0289,'np':237.    ,'pu':244.    ,'am':243.    ,
    'cm':247.    ,'bk':247.    ,'cf':251.
    }


# Atom.atomType
backbone = set([' n  ',' ca ',' c  ',' o  ',' ot1'])


# Atom.atomType
charmm2pdbAtomNames = {
    'sod' :'na','sod ':'na',' sod':'na',
    'ces' :'cs','ces ':'cs',' ces':'cs',
    'cla' :'cl','cla ':'cl',' cla':'cl',
    'pot' :'k' ,'pot ':'k' ,' pot':'k'
    }


# Atom.resName
dna = set(['da','dg','dt','dc'])


# Atom.resName
good = set(['hoh','tip3','ip3''zn2','sod','ces','cla','cal','pot','zn','fe',
            'na','ca','mg','cs','k','cl'])


# Atom.resName
nuc = set([
    'a','g','t','c','u','ade','thy','gua','cyt','ura','da','dg','dt','dc','du'
    ])


# Atom.resName
pro = set(aaVDW.keys())


# Atom.resName
rna = set(['a','g','c','u'])
