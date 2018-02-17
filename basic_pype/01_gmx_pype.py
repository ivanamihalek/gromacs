#!/usr/bin/python -u

import os, subprocess
from argparse import Namespace
from gmx_utils.python import run_setup
from gmx_utils.python.run_setup import GmxEngine, GmxParameters, WorkdirStructure

#########################################
def echo_options(run_options):
	print "pdb:", run_options.pdb
	print "ligands:", run_options.ligands
	print "minimization only?", run_options.minimization

#########################################
def generate_gro_and_top(params):
	# change to topology directory
	os.chdir("/".join([params.run_options.workdir, params.rundirs.top_dir]))
	pdbname = params.run_options.pdb
	pdbfile = pdbname+".pdb"
	topfile = pdbname+".top"
	grofile = pdbname+".gro"
	if os.path.exists(topfile) and os.path.exists(grofile):
		print "\t %s and %s found" % (topfile, grofile)
		return
	program = "pdb2gmx"
	log = "pdb2gmx.log"
	# create topology for the protein
	if pdbname!="none":
		print "\t running", program, ", creating topology for the peptide"
		cmd  = "%s %s " % (params.gmx_engine.gmx, program)
		cmd += "-ignh  -ff  %s  -water %s " % (params.physical.forcefield, params.physical.water)
		cmd += "-f ../%s/%s -o %s -p %s" % (params.rundirs.in_dir, pdbfile,  grofile, topfile)
		# if there are salt bridges or termini modifications,
		# create file with input options for pdb2gmx, and indicate so on the command line
		if os.path.exists("pdb2gmx_in"): os.remove("pdb2gmx_in")
		if os.path.exists("ssbridges"):
			subprocess.call(["bash","-c", "echo ssbridges > pdb2gmx_in"], stdout=subprocess.PIPE)
			cmd  += " -ss"
		if os.path.exists("termini"):
			subprocess.call(["bash","-c", "echo termini >> pdb2gmx_in"], stdout=subprocess.PIPE)
			cmd  += " -ter"
		if os.path.exists("pdb2gmx_in"):
			cmd  += " < pdb2gmx_in"
		params.command_log.write(cmd+"\n")
		proc = subprocess.Popen(["bash","-c", cmd], stdout=subprocess.PIPE)

	# create top for ligands
	# merge

	return

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
	# TODO: check if itp and gro given for all ligands

	###############
	# process pdb into gro and topology files
	###############
	generate_gro_and_top(params)

	###############
	# place the sytem in a box
	###############

	###############
	# add water
	###############

	###############
	# add counterions
	###############

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
