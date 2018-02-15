#!/usr/bin/python -u

from gmx_utils.python import run_setup
from gmx_utils.python.run_setup import GmxEngine

#########################################
def echo_options(run_options):
	print "pdb:", run_options.pdb
	print "ligands:", run_options.ligands
	print "minimization only?", run_options.minimization


#########################################
def main():

	run_options = run_setup.parse_commandline()
	echo_options(run_options)
	print "#########################"
	gmx_engine = GmxEngine("/usr/local/gromacs/bin/GMXRC.bash")


	return True


#########################################
if __name__ == '__main__':
	main()
