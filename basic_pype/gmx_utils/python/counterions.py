import os, subprocess

#########################################
def add(params):

	# change to topology directory
	currdir = params.rundirs.em1_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb
	# previously generated tpr file - we will add counterions
	# if the system complained there was non-zero total charge
	[log_previous, logerr_previous]  = params.gmx_engine.lognames('grompp')
	if not os.path.exists(logerr_previous):
		print "\t in counterions.add(): %s not found " % logerr_previous
		print "\t run grompp once to establish the total charge in the system "
		exit(1)
	# grep -s = suppress errmsg about nonexisting files
	# TODO what do I do if they change the error message?
	command = ['bash', '-c', "grep  -s \'System has non-zero total charge\' %s" % logerr_previous]
	ret = subprocess.Popen(command, stdout=subprocess.PIPE).communicate()
	# communicate() returns a tuple (stdoutdata, stderrdata)
	# and closes the input pype for the subprocess
	zero_charge = not (ret and len(ret[0])>0 and 'charge' in ret[0])

	if zero_charge:
		print "\t the total charge in the system in zero"
		return

	print "\t there is nonzero charge in the system - adding counterions"
	# parse ret to get  the charge
	charge = int(round(float(ret[0].rstrip().split()[-1])))
	if charge>0:
		ion_request = " -nname %s -nn %d " % (params.physical.neg_ion, charge)
	else:
		# notice we turn the charge variable into its absolute value
		ion_request = " -pname %s -np %d " % (params.physical.pos_ion, -charge)

	print ion_request
	exit()

	topfile      = "../%s/%s.top"%(params.rundirs.top_dir,pdbname)

	program = "genion"
	tprfile_previous  = pdbname + ".em_input.tpr"
	ion_gro           = pdbname + ".ion.gro"
	cmdline_args  = "-s %s -o %s" % (tprfile_previous, ion_gro)


	return

'''
	topfile      = "../%s/%s.top"%(params.rundirs.top_dir,pdbname)
	grofile_in   = "../%s/%s.water.gro"%(params.rundirs.top_dir,pdbname)
	paramsfile   = "../%s/em_steep.mdp"%(params.rundirs.in_dir)

	# check whether there is the warning about the charged system in the err.log file ffrom grompp
	# if

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
	program = "counterions"
	cmdline_args  = "-c %s -p %s " % (grofile_in, topfile)
	cmdline_args += "-f %s -o %s " % (paramsfile, tprfile_out)
	# I don't want to get killed on a single warning (that I have charge, most likely here)
	cmdline_args += "-maxwarn 3 "
	
	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "creating parametrization file", params.command_log)
	false_alarms = []
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	# TODO: fix number of waters (if some were present in the initial structure)
'''




