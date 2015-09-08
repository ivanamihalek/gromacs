#!/usr/bin/python

from sys import argv
from getopt import getopt
from string import *

opts, args = getopt(argv[1:], 'f:o:')

fnIN = 'in.pdb'
fnOUT = 'index.ndx'

for o, a in opts:
    if o == '-f':
        fnIN = a
    if o == '-o':
        fnOUT = a

if (fnIN[-4:] == '.gro'):
    ID = (15, 20)
    RESN = (5, 10)
    ATMN = (10, 15)
else:
    ID = (6, 11)
    RESN = (17, 20) 
    ATMN = (12, 16)

aminoacids = ['ABU', 'AIB', 'ALA', 'ARG', 'ARGN', 'ASN', 'ASN1', 'ASP', 
              'ASP1', 'ASPH', 'CYS', 'CYS1', 'CYS2', 'CYSH', 'CSZ', 'CS-', 
              'GLN', 'GLU', 'GLUH', 'GLY', 'HIS', 'HIS1', 'HISA', 'HISB', 
              'HISH', 'HYP', 'ILE', 'LEU', 'LYS', 'LYSH', 'MET', 'PHE', 
              'PRO', 'PRN', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'DALA', 
              'MELEU', 'MEVAL', 'PHL', 'PHEU', 'TRPU', 'TYRU', 'PHEH', 
              'TRPH', 'TYRH', 'ACE', 'NAC', 'NH2', 'LY+']

NDX = {'CA': [], 
       'C': [], 
       'N': [], 
       'H': [], 
       'O': [], 
       'CB': [], 
       'Hsc': [], 
       'SCo': [], 
       'NonP': [],
       'SOL': [],
       'NA+': [],
       'CL-': [],
       'CA2+': [],
       'HOH': []}

fIN = open(fnIN, 'r')
line = fIN.readline()
if fnIN[-4:] == '.gro':
    line = fIN.readline()
    natoms = int(line)
    line = fIN.readline()
else:
    while not line[0:6] in ['ATOM  ', 'HETATM']:
        line = fIN.readline()
    natoms = 1000000
linenum = 1

oldid = 0
add = 0
while line:
    if fnIN[-4:] == '.pdb' and not line[0:6] in ['ATOM  ', 'HETATM']:
        pass
    elif linenum > natoms:
        break
    else:
        try:
            id = int(line[ID[0]:ID[1]]) + add
            if id == oldid - 99999:
                add = add + 100000
                id = id + 100000
            resn = strip(line[RESN[0]:RESN[1]])
            atmn = strip(line[ATMN[0]:ATMN[1]])
            if fnIN[-4:] == '.pdb':
                elem = strip(line[12:14])
            else:
                elem = atmn[0]
                # Dirty, but only necessary for hydrogens...
            if resn in aminoacids:
                if atmn == 'CA':
                    NDX['CA'].append(id)
                elif atmn == 'C':
                    NDX['C'].append(id)
                elif atmn == 'N':
                    NDX['N'].append(id)
                elif atmn == 'H' or ( fnIN[-4:] == '.gro' and resn in ['PRO', 'PRN', 'NH2'] and atmn in ['H1', 'H2']):
                    NDX['H'].append(id)
                elif atmn == 'O':
                    NDX['O'].append(id)
                elif atmn == 'CB':
                    NDX['CB'].append(id)
                else:
                    if elem == 'H':
                        NDX['Hsc'].append(id)
                    else:
                        NDX['SCo'].append(id)
            elif resn == 'SOL':
                NDX['SOL'].append(id)
            elif resn in ['NA', 'NA+', 'Na', 'Na+']:
                NDX['NA+'].append(id)
            elif resn in ['CL', 'CL-', 'Cl', 'Cl-']:
                NDX['CL-'].append(id)
	    elif resn == 'CA2+':
		NDX['CA2+'].append(id)
            elif resn == 'HOH':
                NDX['HOH'].append(id)
            else:
                NDX['NonP'].append(id)
        except ValueError:
            # Probably stumbled upon the box in the gro file...
            pass
    line = fIN.readline()
    linenum = linenum + 1
    
## System = Protein + NonP
## Protein = MainChain+H + SideChain
## Protein-H = MainChain + SideChain-H
## Prot-Masses = MainChain+H + SideChain
## MainChain+H = MainChain + H
## MainChain = Backbone + O
## Backbone = CA + C + N
## SideChain = CB + SCo + Hsc

calpha = NDX['CA']
calpha.sort()

backbone = NDX['CA'] + NDX['C'] + NDX['N']
backbone.sort()

mainchain = backbone + NDX['O']
mainchain.sort()

mainchaincb = mainchain + NDX['CB']
mainchaincb.sort()

mainchainh = mainchain + NDX['H']
mainchainh.sort()

sidechainnoh = NDX['CB'] + NDX['SCo']
sidechainnoh.sort()

sidechain = sidechainnoh + NDX['Hsc']
sidechain.sort()

protein = mainchainh + sidechain
protein.sort()

proteinnoh = mainchain + sidechainnoh
proteinnoh.sort()

protmasses = protein

naplus = NDX['NA+']
naplus.sort()

clmin = NDX['CL-']
clmin.sort()

ca2plus = NDX['CA2+']
ca2plus.sort()

ions = naplus + clmin
ions.sort()

nonprotein = NDX['NonP'] + ions
nonprotein.sort()

basesystem = protein + NDX['NonP']
basesystem.sort()

sol = NDX['SOL']
sol.sort()

solvent = ions + sol + NDX['HOH'] + NDX['CA2+']
solvent.sort()

system = basesystem + solvent
system.sort()

fOUT = open(fnOUT, 'w')

def writendx(list, fh, name):
    if not list:
        return
    
    fh.write('[ %s ]\n' % name)
    teller = 0
    for i in list:
        fh.write('%6d' % i)
        teller = teller + 1
        if not teller % 13:
            fh.write('\n')
    if teller % 13:
        fh.write('\n')

writendx(system, fOUT, 'System')
writendx(protein, fOUT, 'Protein')
writendx(proteinnoh, fOUT, 'Protein-H')
writendx(calpha, fOUT, 'C-alpha')
writendx(backbone, fOUT, 'Backbone')
writendx(mainchain, fOUT, 'MainChain')
writendx(mainchaincb, fOUT, 'MainChain+Cb')
writendx(mainchainh, fOUT, 'MainChain+H')
writendx(sidechain, fOUT, 'SideChain')
writendx(sidechainnoh, fOUT, 'SideChain-H')
writendx(protmasses, fOUT, 'Prot-Masses')
writendx(nonprotein, fOUT, 'Non-Protein')
writendx(basesystem, fOUT, 'Base')
writendx(sol, fOUT, 'SOL')
writendx(solvent, fOUT, 'Solvent+')
writendx(ions, fOUT, 'Ions')
writendx(naplus, fOUT, 'NA+')
writendx(clmin, fOUT, 'CL-')
writendx(ca2plus, fOUT, 'CA2+')
