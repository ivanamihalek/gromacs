#!/usr/bin/python -u

from argparse import Namespace
from gmx_utils import run_setup, gro_and_top, box, water
from gmx_utils import counterions, local_energy_minimum
from gmx_utils import solvent_equilibration, production
from gmx_utils import postproduction
from gmx_utils.run_setup  import WorkdirStructure
from gmx_utils.gmx_engine import GmxEngine
from gmx_utils.gmx_params import GmxParameters

#########################################
def echo_options(run_options):
	print "pdb:", run_options.pdb
	print "ligands:", run_options.ligands
	print "minimization only?", run_options.minimization

# the core of the pipeline; it assumes that the
# structure is cleaned up (no gaps or missing residues)
# and that the itp and gro files for the ligands are provided


#########################################
def core(params):
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
	production.run(params)

	###############
	#  compress the trajectory, strip hydrogens (and gzip)
	###############
	postproduction.produce_viewable_trajectory(params)


	return True

#########################################
def main():

	###############
	# set up the input and the parameters for the run
	###############
	params = Namespace()
	params.run_options  = run_setup.parse_commandline()
	params.physical     = GmxParameters(params.run_options)
	params.gmx_engine   = GmxEngine("/usr/local/gromacs/bin/GMXRC.bash")
	params.rundirs      = WorkdirStructure(params.run_options)
	params.command_log  = open(params.run_options.workdir+"/commands.log","w")

	core(params)
	params.command_log.close()


#########################################
if __name__ == '__main__':
	main()
