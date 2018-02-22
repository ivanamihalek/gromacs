#!/usr/bin/python -u

from argparse import Namespace
from gmx_utils.python import run_setup, gro_and_top, box, water, grompp
from gmx_utils.python import counterions
from gmx_utils.python.run_setup  import WorkdirStructure
from gmx_utils.python.gmx_engine import GmxEngine
from gmx_utils.python.gmx_params import GmxParameters


#########################################
def echo_options(run_options):
	print "pdb:", run_options.pdb
	print "ligands:", run_options.ligands
	print "minimization only?", run_options.minimization


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

	###############
	# TODO: check pdb - missing residues and sidechains

	###############
	# TODO: check if itp and gro given for all ligands

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
	# grompp = generate parametrized topology file
	###############
	grompp.generate(params)

	###############
	# add counterions if the grompp warned us about the charged system
	###############
	counterions.add(params)

	###############
	# geometry "minimization"  - round 1
	###############

	###############
	# geometry "minimization"  - round 2
	###############

	###############
	# position restained md    -  NVT
	###############

	###############
	# position restained md    -  NPT
	###############

	###############
	# ... and finally ... the production run!
	###############

	###############
	#  compress the trajectory, strip hydrogens and gzip
	###############

	params.command_log.close()

	return True


#########################################
if __name__ == '__main__':
	main()
