#!/usr/bin/python -u

import os, subprocess

from argparse import Namespace
from gmx_lib import run_setup, gro_and_top, box, water
from gmx_lib import counterions, local_energy_minimum
from gmx_lib import solvent_equilibration, production
from gmx_lib import postproduction
from gmx_lib.run_setup  import WorkdirStructure
from gmx_lib.gmx_engine import GmxEngine
from gmx_lib.gmx_params import GmxParameters


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
def rna_restraints(params):
	# change to  directory for this stage
	currdir = params.rundirs.production_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	proc = subprocess.Popen(["bash", "-c", "grep include ../%s/*top | grep -v %s | grep RNA "%(params.rundirs.top_dir,params.physical.forcefield)],
					stdout=subprocess.PIPE, stderr=None)

	for line in proc.stdout.readlines():
		itp = line.rstrip().replace("#","").replace('"','').replace("include","").replace(" ","")
		subprocess.call(["bash", "-c", "ln -sf ../%s/%s %s"%(params.rundirs.top_dir,itp,itp)],
					stdout=None, stderr=None)

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
	# note: pdb is not checked for missing residues and sidechains
	######################

	######################
	# adjust the run length
	######################
	# the frequency of output adjusted to 100 frames, no matter what the length of simulation
	params.physical.set_run_lengths(params,  em_steep=12000, em_lbfgs=100,
											pr_nvt="50ps", pr_npt="10ps", md="20ns")

	if False: params.physical.set_annealing_schedule(params, annealing_type='single',
								annealing_npoints='9', # these are picoseconds
								annealing_time= "0    200  400  600  700  800  1000  1200  1400 ",
								annealing_temp= "300  350  400  450  470  450   400   350   300 ")


	###############
	# process pdb into gro and topology files
	###############
	gro_and_top.generate(params)

	###############
	# place the system in a box (define box boundaries, move to center)
	###############
	box.generate(params)

	###############
	# solvate
	###############
	water.add(params)

	###############
	# add counterions if needed
	###############
	counterions.add(params)

	###############
	# geometry "minimization"  - round 1; method: steepest decent
	###############
	local_energy_minimum.find(params, 'em1')

	###############
	# geometry "minimization"  - round 2; method: LBFGS
	###############
	local_energy_minimum.find(params, 'em2')

	###############
	# position restrained md    -  NVT
	###############
	solvent_equilibration.run(params, 'pr1')

	###############
	# position restrained md    -  NPT
	###############
	solvent_equilibration.run(params, 'pr2')

	###############
	# ... and finally ... the production run!
	###############
	rna_restraints(params) # link to restrains from production dir
	params.physical.request_restraints(params,"md.mdp")
	production.run(params, position_restrained_run=True)

	###############
	#  compress the trajectory, strip hydrogens (and gzip)
	###############
	postproduction.produce_viewable_trajectory(params)




	######################
	# cleanup and exit
	######################
	params.command_log.close()
	return True


#########################################
if __name__ == '__main__':
	main()

