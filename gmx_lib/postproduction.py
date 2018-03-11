#!/usr/bin/python -u
import os, subprocess

from argparse import Namespace
import run_setup
from run_setup  import WorkdirStructure
from gmx_engine import GmxEngine
from gmx_params import GmxParameters

#########################################
def produce_viewable_trajectory(params):

	# change to production directory
	currdir = params.rundirs.post_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))

	pdbname     = params.run_options.pdb
	xtc_file    = pdbname + ".md.xtc"
	trrfile_in  = "../%s/%s.md_out.trr"%(params.rundirs.production_dir, pdbname)
	tprfile_in  = "../%s/%s.md_input.tpr"%(params.rundirs.production_dir, pdbname)
	for infile in [trrfile_in, tprfile_in]:
		if not os.path.exists(infile):
			print "\t in postproduction.produce_viewable_trajectory(): %s not found (?)" % (infile)
			exit(1)

	#########################
	# xtc, a compressed version of trajectory
	program = "trjconv"
	cmdline_args  = "-pbc nojump -f %s -o %s " % (trrfile_in, xtc_file)
	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "creating xtc file", params.command_log)
	params.gmx_engine.check_logs_for_error(program)

	#########################
	# pdb file
	program = "trjconv"
	cmdline_args  = "-s %s -f %s -fit progressive -o %s.trj.pdb " % (tprfile_in, xtc_file, pdbname)
	# pipe in the group to fit on and to output (2 = Protein - Hydrogen in both cases)
	outf = open("trjconv.pipein", "w");  outf.write("2\n2\n"); outf.close()
	params.gmx_engine.run(program, cmdline_args, "producing pdb file", params.command_log, pipein="cat  trjconv.pipein")
	params.gmx_engine.check_logs_for_error(program)


	#########################
	# number of hydrogen bonds within the protein
	program = "hbond"
	cmdline_args  = "-s %s -f %s  -xvg none -num %s.hbonds " % (tprfile_in, xtc_file, pdbname)
	outf = open("hbonds.pipein", "w");  outf.write("1\n1\n"); outf.close()
	params.gmx_engine.run(program, cmdline_args, "counting hbonds", params.command_log, pipein="cat  hbonds.pipein")
	params.gmx_engine.check_logs_for_error(program)

	#########################
	subprocess.Popen(["bash", "-c", "rm -f \#* *.pipein"])
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

	produce_viewable_trajectory(params)
	params.command_log.close()



#########################################
if __name__ == '__main__':
	main()



