import os, subprocess

#########################################
def minimize(params):
	# change to topology directory
	currdir = params.rundirs.em1_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb

	tprfile_in   = pdbname+".em_input.tpr"
	grofile_out  = pdbname+".em_out.gro"
	traj_out     = pdbname+".em_out.trr"
	edrfile_out  = pdbname+".em_out.edr"
	if os.path.exists(grofile_out):
		print "\t %s found" % (grofile_out)
		return
	for infile in [tprfile_in]:
		if not os.path.exists(infile):
			print "\t in grompp.generate(): %s not found (?)" % (infile)
			exit(1)

	program = "mdrun" # nt 1; run multiple trajectories instead
	# mdrun will produce trajectory and edr (energy) files,  whether we ask for it or not,
	# so we might just as well name them so we can remove them later
	cmdline_args  = "-s %s -c %s -nt 1 -o %s -e %s" % (tprfile_in, grofile_out, traj_out, edrfile_out)

	params.command_log.write("in %s:\n" % (currdir))
	params.gmx_engine.run(program, cmdline_args, "looking for the local energy minimum", params.command_log)
	false_alarms = ["masses will be determined based on residue and atom names"]
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	print "\t ", params.gmx_engine.convergence_line(program)
	os.remove(traj_out)
	os.remove(edrfile_out)
	
	return



