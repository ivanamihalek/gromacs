import os, subprocess

#########################################
def generate(params):
	# change to topology directory
	currdir = params.rundirs.em1_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb
	topfile      = "../%s/%s.top"%(params.rundirs.top_dir,pdbname)
	grofile_in   = "../%s/%s.water.gro"%(params.rundirs.top_dir,pdbname)
	paramsfile   = "../%s/em_steep.mdp"%(params.rundirs.in_dir)
	tprfile_out  = pdbname+".em_input.tpr"
	if os.path.exists(tprfile_out):
		print "\t %s found" % (tprfile_out)
		return
	for infile in [topfile, grofile_in, paramsfile]:
		if not os.path.exists(infile):
			print "\t in grompp.generate(): %s not found (?)" % (infile)
			exit(1)

	# earlier program  genbox was split itno solvate and insert-molecules
	# the latter inserts -nmol copies of the system specified in the -ci input file - I guess I don'e need this
	# (perhaps for lipid generation in the simulation of a membrane bound system?)

	# a layer of water can be added - should this be tried?
	program = "grompp"
	cmdline_args  = "-c %s -p %s " % (grofile_in, topfile)
	cmdline_args += "-f %s -o %s " % (paramsfile, tprfile_out)
	# I don't want to get killed on a single warning (that I have charge, most likely here)
	cmdline_args += "-maxwarn 3 "

	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "creating parametrization file", params.command_log)
	false_alarms = []
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	# TODO: fix number of waters (if some were present in the initial structure)

	return



