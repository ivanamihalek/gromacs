#!/usr/bin/python -u

import os, subprocess


from argparse import Namespace
from gmx_lib import run_setup
from gmx_lib.run_setup  import WorkdirStructure
from gmx_lib.gmx_engine import GmxEngine
from gmx_lib.gmx_params import GmxParameters
from modeling_lib.scwrl_engine import ScwrlEngine
from modeling_lib.pdb_doctor import check_pdb_for_missing_atoms, fix_sidechains

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
	params.scwrl_engine = ScwrlEngine("/usr/local/bin/scwrl4/Scwrl4")

	make_directories(params)
	fill_input_dir  (params)
	params.rundirs      = WorkdirStructure(params.run_options)
	params.command_log  = open(params.run_options.workdir+"/commands.log","w")

	######################
	# check pdb file for missing residues and sidechains
	######################
	chain_breaks, missing_sidechains = check_pdb_for_missing_atoms(params)
	if missing_sidechains: fix_sidechains(params, missing_sidechains)
	#  TODO: fix chain breaks (if they exist)
	# if chain_breaks: fix_chain_breaks(params, missing_sidechains)

	######################
	# adjust the run length
	######################
	# the frequency of output adjusted to 100 frames, no matter what the length of simulation
	params.physical.set_run_lengths(params, em_steep=10000, em_lbfgs=1000,
											pr_nvt="10ps", pr_npt="10ps", md="20ps")

	######################
	# check pdb file for missing residues and side chains
	######################
	core(params)

	######################
	# cleanup and exit
	######################
	params.command_log.close()
	return True


#########################################
if __name__ == '__main__':
	main()

