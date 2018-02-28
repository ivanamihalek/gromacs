#!/usr/bin/python -u

import os, subprocess, shutil
import Bio.PDB as biopython
#from Bio.PDB.PDBParser import PDBParser

from argparse import Namespace
from gmx_lib import run_setup
from gmx_lib.run_setup  import WorkdirStructure
from gmx_lib.gmx_engine import GmxEngine
from gmx_lib.gmx_params import GmxParameters
from modeling_lib.scwrl_engine import ScwrlEngine

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
def get_biopython_structure(params):
	currdir = params.rundirs.in_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb
	parser = biopython.PDBParser(PERMISSIVE=1)
	structure = parser.get_structure(pdbname, pdbname+".pdb")
	# TODO: error state?
	return structure


#########################################
def check_pdb_for_missing_atoms(structure):

	number_of_atoms = { "ASN":8, "GLU":9, "CYS":6, "THR":7, "TYR":12, "PRO":7, "HIS":10,
						"SER":6, "LEU":8, "TRP":14, "PHE":11, "MET":8, "ALA":5, "GLY":4,
						"VAL":7, "GLN":9, "LYS":9, "ILE":8, "ASP":8, "ARG":11}
	# TODO: atypical residues, terminal residues
	chain_breaks     = {}
	missing_sidechains = {}
	for chain in structure.get_chains():
		chain_breaks[chain.id] = []
		missing_sidechains[chain.id] = []
		sequential_id = 0
		prev_res_id = None
		for residue in chain:
			# residue id is a triple; [1] is the id assigned by PDB
			# here can check for the gap in the residue_id sequence (missing residues)
			atom_list = biopython.Selection.unfold_entities(residue, 'A')

			# what should be the number of atoms in each residue
			# check for missing sidechains
			if len(atom_list)!= number_of_atoms[residue.get_resname()]:
				print "missing sidechain:"
				print chain.id, residue.get_resname(), residue.id[1], sequential_id
				print "\t", len(atom_list), number_of_atoms[residue.get_resname()]
				missing_sidechains[chain.id].append(sequential_id)
			resid = int(residue.id[1])
			if prev_res_id and (resid-prev_res_id)>1:
				print 'chain break:'
				print chain.id, residue.get_resname(), residue.id[1], sequential_id
				chain_breaks[chain.id].append(sequential_id)
			prev_res_id = resid
			sequential_id += 1
	return chain_breaks, missing_sidechains

#########################################
def fix_sidechains(params, hash_of_residue_positions):
	# TODO: make sure this actually works for multiple chains
	currdir = params.rundirs.in_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname  = params.run_options.pdb
	pdb_cleanup_dir = "pdb_cleanup"
	if not os.path.exists(pdb_cleanup_dir): os.mkdir(pdb_cleanup_dir)
	shutil.copyfile(pdbname+".pdb", "%s.orig.pdb"%(pdbname))
	shutil.copyfile(pdbname+".pdb", "%s/%s.pdb"%(pdb_cleanup_dir, pdbname))
	os.chdir(pdb_cleanup_dir)

	for chain, positions in hash_of_residue_positions.iteritems():
		# get sequence
		pdb2seq = "%s/pdb2seq.pl"%params.run_options.perl_utils
		stdout  = subprocess.Popen([pdb2seq, pdbname+".pdb", chain], stdout=subprocess.PIPE)
		seqlist = list(stdout.communicate()[0].replace("\n","").lower())
		for pos in positions: seqlist[pos] = seqlist[pos].upper()
		seq = ''.join(seqlist)
		print seq
		# send sequence and position to SCWRL
		pass



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
	params.ScwrlEngine  = ScwrlEngine("/usr/local/bin/scwrl4/Scwrl4")

	make_directories(params)
	fill_input_dir  (params)
	params.rundirs      = WorkdirStructure(params.run_options)
	params.command_log  = open(params.run_options.workdir+"/commands.log","w")

	######################
	# check pdb file for missing residues and sidechains
	######################
	structure = get_biopython_structure(params)
	chain_breaks, missing_sidechains = check_pdb_for_missing_atoms(structure)

	fix_sidechains(params, missing_sidechains)

	#  TODO: fix chain breaks
	exit()
	core(params)
	params.command_log.close()
	return True


#########################################
if __name__ == '__main__':
	main()

