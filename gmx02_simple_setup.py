#!/usr/bin/python -u

import os, subprocess
import Bio.PDB as biopython
#from Bio.PDB.PDBParser import PDBParser

from argparse import Namespace
from gmx_utils import run_setup
from gmx_utils.run_setup  import WorkdirStructure
from gmx_utils.gmx_engine import GmxEngine
from gmx_utils.gmx_params import GmxParameters

from gmx01_core import core


#########################################
def make_directories(params):
	os.chdir(params.run_options.workdir)
	for dir in WorkdirStructure.workdir_names:
		if not os.path.exists(dir): os.mkdir(dir)

#########################################
def fill_input_dir(params):
	in_dir = "/".join([params.run_options.workdir, WorkdirStructure.workdir_names[0]])
	os.chdir(in_dir)
	pdbname = params.run_options.pdb
	if os.path.exists("../%s.pdb"%pdbname): os.rename("../%s.pdb"%pdbname,"./%s.pdb"%pdbname)
	subprocess.call(["bash", "-c", "cp -f %s/* %s"%(params.run_options.mdp_template_home, in_dir)])


#########################################
def check_pdb(params):
	currdir = params.rundirs.in_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb
	parser = biopython.PDBParser(PERMISSIVE=1)
	structure = parser.get_structure(pdbname, pdbname+".pdb")
	for chain in structure.get_chains():
		sequential_id = 0
		for residue in chain:
			# residue id is a triple; [1] is the id assigned by PDB
			# here can check for the gap in the residue_id sequence (missing residues)
			print chain.id, residue.get_resname(), residue.id[1], sequential_id


			atom_list = biopython.Selection.unfold_entities(residue, 'A')
			# what should be the number of atoms in each residue
			# check for missing sidechains
			for atom in atom_list:
				print "\t",atom.get_name()
			sequential_id += 1
	exit()
# simple setup: single peptide, single run
# create dir structure
# clean up the structure (fix gaps or missing residues)
# make sure  that the itp and gro files for the ligands are provided
#########################################
def main():

	params = Namespace()
	params.run_options  = run_setup.parse_commandline()
	params.physical     = GmxParameters(params.run_options)
	params.gmx_engine   = GmxEngine("/usr/local/gromacs/bin/GMXRC.bash")

	make_directories(params)
	fill_input_dir  (params)
	params.rundirs      = WorkdirStructure(params.run_options)
	params.command_log  = open(params.run_options.workdir+"/commands.log","w")
	######################
	# check pdb file for missing residues and sidechains
	######################
	check_pdb(params)
	exit()
	core(params)
	params.command_log.close()
	return True


#########################################
if __name__ == '__main__':
	main()

