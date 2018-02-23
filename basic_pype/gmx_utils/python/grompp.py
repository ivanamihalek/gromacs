import os, subprocess


#########################################
def dir_and_files(params, previous_stage):

	[currdir, topfile, grofile_in, paramsfile, tprfile_out] = ['']*5
	pdbname      = params.run_options.pdb

	if previous_stage in ['water', 'ions']:
		currdir = params.rundirs.em1_dir
		topfile      = "../%s/%s.top"       % (params.rundirs.top_dir, pdbname)
		grofile_in   = "../%s/%s.%s.gro"    % (params.rundirs.top_dir, pdbname, previous_stage)
		paramsfile   = "../%s/em_steep.mdp" % (params.rundirs.in_dir)
		tprfile_out  = pdbname+".em_input.tpr"

	elif previous_stage== 'em1':
		currdir = params.rundirs.em2_dir
		topfile      = "../%s/%s.top"        % (params.rundirs.em1_dir, pdbname)
		grofile_in   = "../%s/%s.em_out.gro" % (params.rundirs.em1_dir, pdbname)
		paramsfile   = "../%s/em_lbfgs.mdp"  % (params.rundirs.in_dir)
		tprfile_out  = pdbname+".em_input.tpr"

	elif previous_stage== 'em2':
		currdir = params.rundirs.em2_dir
		topfile      = "../%s/%s.top"        % (params.rundirs.em2_dir, pdbname)
		grofile_in   = "../%s/%s.em_out.gro" % (params.rundirs.em2_dir, pdbname)
		paramsfile   = "../%s/em_lbfgs.mdp"  % (params.rundirs.in_dir)
		tprfile_out  = pdbname+".pr_input.tpr"


	return [currdir, topfile, grofile_in, paramsfile, tprfile_out]


#########################################
def generate(params, stage):


	[currdir, topfile, grofile_in, paramsfile, tprfile_out] = dir_and_files(params, stage)
	# change to topology directory
	os.chdir("/".join([params.run_options.workdir, currdir]))
	if os.path.exists(tprfile_out):
		print "\t %s found" % (tprfile_out)
		return
	for infile in [topfile, grofile_in, paramsfile]:
		if not os.path.exists(infile):
			print "\t in grompp.generate(): %s not found (?)" % (infile)
			exit(1)

	program = "grompp"
	cmdline_args  = "-c %s -p %s " % (grofile_in, topfile)
	cmdline_args += "-f %s -o %s " % (paramsfile, tprfile_out)
	# I don't want to get killed on a single warning (that I have charge, most likely here)
	cmdline_args += "-maxwarn 3 "

	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "creating parametrization file", params.command_log)
	false_alarms = []
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	# TODO group.ndx  - when and where do I need that?

	return



