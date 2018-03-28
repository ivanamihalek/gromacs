import os, subprocess


#########################################
def dir_and_files(params, stage):

	[currdir, topfile, grofile_in, paramsfile, tprfile_out] = ['']*5
	pdbname      = params.run_options.pdb
	topfile      = "../%s/%s.top"       % (params.rundirs.top_dir, pdbname)

	if stage == 'ions':
		currdir = params.rundirs.em1_dir
		grofile_in   = "../%s/%s.%s.gro"    % (params.rundirs.top_dir, pdbname, 'water')
		paramsfile   = "../%s/em_steep.mdp" % (params.rundirs.in_dir)
		tprfile_out  = pdbname + ".em_input.tpr"

	elif stage== 'em1':
		currdir = params.rundirs.em1_dir
		grofile_in   = "../%s/%s.%s.gro"    % (params.rundirs.top_dir, pdbname, 'ions')
		paramsfile   = "../%s/em_steep.mdp" % (params.rundirs.in_dir)
		tprfile_out  = pdbname + ".em_input.tpr"

	elif stage== 'em2':
		currdir = params.rundirs.em2_dir
		grofile_in   = "../%s/%s.em_out.gro" % (params.rundirs.em1_dir, pdbname)
		paramsfile   = "../%s/em_lbfgs.mdp"  % (params.rundirs.in_dir)
		tprfile_out  = pdbname + ".em_input.tpr"

	elif stage== 'pr1':
		currdir = params.rundirs.pr1_dir
		grofile_in   = "../%s/%s.em_out.gro" % (params.rundirs.em2_dir, pdbname)
		paramsfile   = "../%s/pr_nvt.mdp"  % (params.rundirs.in_dir)
		tprfile_out  = pdbname + ".pr_input.tpr"

	elif stage== 'pr2':
		currdir = params.rundirs.pr2_dir
		grofile_in   = "../%s/%s.pr_out.gro" % (params.rundirs.pr1_dir, pdbname)
		paramsfile   = "../%s/pr_npt.mdp"  % (params.rundirs.in_dir)
		tprfile_out  = pdbname + ".pr_input.tpr"

	elif stage== 'production':
		currdir = params.rundirs.production_dir
		grofile_in   = "../%s/%s.pr_out.gro" % (params.rundirs.pr2_dir, pdbname)
		paramsfile   = "../%s/md.mdp"  % (params.rundirs.in_dir)
		tprfile_out  = pdbname + ".md_input.tpr"

	return [currdir, topfile, grofile_in, paramsfile, tprfile_out]


#########################################
def generate(params, stage, position_restrained_run=False):


	[currdir, topfile, grofile_in, paramsfile, tprfile_out] = dir_and_files(params, stage)
	# change to topology directory
	os.chdir("/".join([params.run_options.workdir, currdir]))
	if os.path.exists(tprfile_out):
		print "\t in grompp.generate(%s): %s found" % (stage, tprfile_out)
		params.gmx_engine.delete_junk()
		return
	for infile in [topfile, grofile_in, paramsfile]:
		if not os.path.exists(infile):
			print "\t in grompp.generate(%s): %s not found (?)" % (stage, infile)
			exit(1)

	program = "grompp"
	cmdline_args  = "-c %s -p %s " % (grofile_in, topfile)
	# Note:
	# "From GROMACS-2018, you need to specify (option -r) the position restraint coordinate files
	# explicitly to avoid mistakes, although you can still use the same file as you
	# specify for the -c option."
	# ergo repeat the same argument twice with two different flags (whoever came up with this 'solution')
	if position_restrained_run: cmdline_args += "-r %s " % (grofile_in)
	cmdline_args += "-f %s -o %s " % (paramsfile, tprfile_out)
	# I don't want to get killed on a single warning (that I have charge, most likely here)
	cmdline_args += "-maxwarn 3 "

	params.command_log.write("in %s:\n"%(currdir))
	msg = "creating parametrization file (%s)" % stage
	params.gmx_engine.run(program, cmdline_args, msg, params.command_log)
	false_alarms = []
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	# TODO group.ndx  - when and where do I need that?

	return



