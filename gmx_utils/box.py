import os

#########################################
def generate(params):
	# change to topology directory
	currdir = params.rundirs.top_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname = params.run_options.pdb
	grofile_in  = pdbname + ".gro"
	grofile_out = pdbname + ".box.gro"
	if os.path.exists(grofile_out):
		print "\t %s found" % (grofile_out)
		return
	if not os.path.exists(grofile_in):
		print "\t in box.generate(): %s not found in %s (?)" % (grofile_in, os.getcwd())
		exit(1)

	program = "editconf"
	cmdline_args  = "-f %s -o %s " % (grofile_in, grofile_out)
	# -c is the centering command:
	cmdline_args += "-bt %s -d %s -c" % (params.physical.box_type, params.physical.box_edge)

	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "placing system in a box", params.command_log)
	false_alarms = ["masses will be determined based on residue and atom names", "No boxtype specified"]
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	return



