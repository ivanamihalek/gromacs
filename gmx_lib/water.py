import os


#########################################
def add(params):
	# change to topology directory
	currdir = params.rundirs.top_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb
	topfile      = pdbname+".top"
	grofile_in   = pdbname+".box.gro"
	water_coords = params.physical.water_coords_gro # gmx will fnd this in its lib paths
	grofile_out  = pdbname+".water.gro"

	if os.path.exists(grofile_out):
		print "\t %s found" % (grofile_out)
		return
	for infile in [topfile, grofile_in]:
		if not os.path.exists(infile):
			print "\t in water.add(): %s not found (?)" % (infile)
			exit(1)

	# earlier program  genbox was split itno solvate and insert-molecules
	# the latter inserts -nmol copies of the system specified in the -ci input file - I guess I don'e need this
	# (perhaps for lipid generation in the simulation of a membrane bound system?)

	# a layer of water can be added - should this be tried?
	program = "solvate"
	cmdline_args  = "-cp %s -p %s " % (grofile_in, topfile)
	cmdline_args += "-cs %s -o %s " % (water_coords, grofile_out)

	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "adding water", params.command_log)
	false_alarms = ["masses will be determined based on residue and atom names", "radii will be determined"]
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	# TODO: fix number of waters (if some were present in the initial structure)
	# TODO: in that case, when fixing the number of water molecules, check that
	# -p *top option for genion still results in the correct number of SOL or water molecules

	return



