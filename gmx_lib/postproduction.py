#!/usr/bin/python -u
import os, subprocess

from argparse import Namespace
import run_setup
from run_setup  import WorkdirStructure
from gmx_engine import GmxEngine
from gmx_params import GmxParameters

#########################################
def make_xtc(params):

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
	params.command_log.write("in %s:\n"%(params.rundirs.post_dir))
	params.gmx_engine.run(program, cmdline_args, "creating xtc file", params.command_log)
	params.gmx_engine.check_logs_for_error(program)

#########################################
# see counterions.py to see explanation for the similar case of genion
def hack_group_numbers_out_of_trjconv(params, tprfile_in, xtc_file):
	protein_group_number = ""
	non_water_group_number = ""

	program = "trjconv"
	cmdline_args  = "-s %s -f %s -fit progressive -o grouphack " % (tprfile_in, xtc_file)
	params.gmx_engine.run(program, cmdline_args, pipein="echo 1")
	[logfile, errlogfile]  = params.gmx_engine.lognames(program)
	infile = open(errlogfile,"r")
	for line in infile:
		field = line.split()
		if len(field)<4: continue
		if field[0]=='Group':
			if  field[3]=='Protein-H)':
				protein_group_number = field[1]
			elif field[3]=='non-Water)':
				non_water_group_number = field[1]
	infile.close()
	for group_number in ["protein_group_number","non_water_group_number"]:
		if eval(group_number) == "":
			print group_number, "not found in %s" % errlogfile
			exit(1)
	subprocess.call(["bash", "-c", "rm -f %s %s grouphack.xtc"%(logfile, errlogfile)])
	return [protein_group_number, non_water_group_number]

#########################################
def make_pdb(params):

	program = "trjconv"
	pdbname     = params.run_options.pdb
	xtc_file    = pdbname + ".md.xtc"
	tprfile_in  = "../%s/%s.md_input.tpr"%(params.rundirs.production_dir, pdbname)

	# not sure how the group numbers are assigned, so hack them out
	[protein_group_number, non_water_group_number] = hack_group_numbers_out_of_trjconv(params, tprfile_in, xtc_file)
	# pipe in the groups to fit on and to output
	outf = open("trjconv.pipein", "w");  outf.write("%s\n%s\n"%(protein_group_number, non_water_group_number)); outf.close()

	pdb_w_hydrogens = "%s.trj.wH.pdb"%pdbname
	cmdline_args  = "-s %s -f %s -fit progressive -o %s" % (tprfile_in, xtc_file, pdb_w_hydrogens)
	params.gmx_engine.run(program, cmdline_args, "producing pdb file", params.command_log, pipein="cat  trjconv.pipein")
	params.gmx_engine.check_logs_for_error(program)

	# strip hydrogen atoms from the pdb
	inf  = open(pdb_w_hydrogens,"r")
	outf = open("%s.trj.pdb"%pdbname,"w")
	for line in inf:
		if line[:4]=="ATOM" and line[-2:-1]=="H": continue
		outf.write(line)
	inf.close()
	outf.close()
	#cleanup
	subprocess.call(["bash", "-c", "rm -f %s "% pdb_w_hydrogens])

#########################################
def count_hbonds(params):

	# number of hydrogen bonds within the protein
	program = "hbond"
	pdbname     = params.run_options.pdb
	xtc_file    = pdbname + ".md.xtc"
	tprfile_in  = "../%s/%s.md_input.tpr"%(params.rundirs.production_dir, pdbname)
	cmdline_args  = "-s %s -f %s  -xvg none -num %s.hbonds " % (tprfile_in, xtc_file, pdbname)
	outf = open("hbonds.pipein", "w");  outf.write("1\n1\n"); outf.close()
	params.gmx_engine.run(program, cmdline_args, "counting hbonds", params.command_log, pipein="cat  hbonds.pipein")
	params.gmx_engine.check_logs_for_error(program)

def cleanup(params):
	#########################
	subprocess.Popen(["bash", "-c", "rm -f \#* *.pipein"])
	return

#########################################
def produce_viewable_trajectory(params):
	# change to production directory
	os.chdir("/".join([params.run_options.workdir, params.rundirs.post_dir]))

	make_xtc(params)
	make_pdb(params)
	cleanup(params)
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
	params.command_log  = open(params.run_options.workdir+"/postproc_commands.log","w")

	# change to production directory
	currdir = params.rundirs.post_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))

	make_xtc(params)
	make_pdb(params)
	count_hbonds(params)
	cleanup(params)
	params.command_log.close()



#########################################
if __name__ == '__main__':
	main()



