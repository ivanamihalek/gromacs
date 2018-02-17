#!/usr/bin/python -u

from gmx_utils.python import run_setup
from gmx_utils.python.run_setup import GmxEngine, Parameters, WorkdirStructure

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
	params     = Parameters(run_options)
	gmx_engine = GmxEngine("/usr/local/gromacs/bin/GMXRC.bash")
	workdir_names = ["00_input", "01_topology", "02_em_steepest", "03_em_lbfgs",
					"04_nvt_eq", "05_mpt_eq", "06_production"]
	rundirs    = WorkdirStructure(run_options, workdir_names)

	return True


#########################################
if __name__ == '__main__':
	main()
