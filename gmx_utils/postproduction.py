import os, subprocess


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
	program = "trjconv"
	cmdline_args  = "-pbc nojump -f %s -o %s " % (trrfile_in, xtc_file)
	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "creating xtc file", params.command_log)
	params.gmx_engine.check_logs_for_error(program)

	#########################
	program = "trjconv"
	cmdline_args  = "-s %s -f %s -o tmp.pdb " % (tprfile_in, xtc_file)
	params.gmx_engine.run(program, cmdline_args, "producing pdb file", params.command_log, pipein="echo 1")
	params.gmx_engine.check_logs_for_error(program)
	cmd = "awk -F ''  '$14 != \"H\"' tmp.pdb >  %s.trj.pdb"%pdbname
	params.command_log.write(cmd+"\n")
	subprocess.call(["bash","-c",cmd])
	subprocess.call(["bash","-c","rm -f tmp.pdb"])

	return





