import os, subprocess

from ..python import grompp

#########################################
def add(params):

	# try running the preprocessor, and see if it complains about the charge in the system
	grompp.generate(params, 'ions')

	# change to topology directory
	currdir = params.rundirs.em1_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))

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

	#############################################
	if zero_charge:
		print "\t the total charge in the system is zero"
		# our 'ions' file is the same as the 'water' file
		pdbname      = params.run_options.pdb
		os.chdir("/".join([params.run_options.workdir, params.rundirs.top_dir]))
		# f = force; remove xisting destination
		cmd          = "ln -sf %s.water.gro  %s.ions.gro " % (pdbname,pdbname)
		subprocess.call(["bash", "-c", cmd])
		return

	#############################################
	print "\t there is nonzero charge in the system"
	pdbname      = params.run_options.pdb
	topfile      = "../%s/%s.top"%(params.rundirs.top_dir,pdbname)
	tprfile_previous  = pdbname + ".em_input.tpr"
	for infile in [topfile, tprfile_previous, topfile]:
		if not os.path.exists(infile):
			print "\t in counterions.add(): %s not found (?)" % (infile)
			exit(1)
	# parse ret to get  the charge
	ion_request = parse_non_zero_charge_msg(params, ret)
	# find the number assigned to water by genion (genion input is interactive, so we need to pipe it in)
	water_group_number = hack_water_group_number_out_of_genion(params, tprfile_previous)

	program = "genion"
	# say, we keep all the geometry files in top_dir
	ion_gro = "../%s/%s.ions.gro"%(params.rundirs.top_dir,pdbname)
	# when -p *top is present the original topology file is modified in place:
	# the number of SOL molecules is reduced by the number that was replaced by ions
	cmdline_args  = "-s %s -o %s -p %s" % (tprfile_previous, ion_gro, topfile)
	cmdline_args += ion_request

	params.command_log.write("in %s:\n"%(currdir))
	params.gmx_engine.run(program, cmdline_args, "adding couterions",params.command_log,"echo %s"%water_group_number)
	false_alarms = ["turning of free energy, will use lambda=0"]
	params.gmx_engine.check_logs_for_error(program, false_alarms)
	os.remove(tprfile_previous)


	return

####################################
def parse_non_zero_charge_msg(params, ret):
	charge = int(round(float(ret[0].rstrip().split()[-1])))
	if charge>0:
		ion_request = " -nname %s -nn %d " % (params.physical.neg_ion, charge)
	else:
		# notice we turn the charge variable into its absolute value
		ion_request = " -pname %s -np %d " % (params.physical.pos_ion, -charge)

	return ion_request


####################################
'''
Interactive prompt by genion looks something like this:
"Select a continuous group of solvent molecules
Group     0 (         System) has  5201 elements
Group     1 (        Protein) has   596 elements
Group     2 (      Protein-H) has   295 elements
Group     3 (        C-alpha) has    36 elements
Group     4 (       Backbone) has   108 elements
Group     5 (      MainChain) has   145 elements
Group     6 (   MainChain+Cb) has   179 elements
Group     7 (    MainChain+H) has   182 elements
Group     8 (      SideChain) has   414 elements
Group     9 (    SideChain-H) has   150 elements
Group    10 (    Prot-Masses) has   596 elements
Group    11 (    non-Protein) has  4605 elements
Group    12 (          Water) has  4605 elements
Group    13 (            SOL) has  4605 elements
Group    14 (      non-Water) has   596 elements
Select a group:"

at which point we are supposed to type in 12 for water.
This number may change according to the content of tpr,
which is a binary. The following function is a hack to get the number
by feeding  genion some dummy input.Rather than parsing some 
other file, this way I know exactly what is genion going to ask.'
'''
def hack_water_group_number_out_of_genion(params, tprfile):
	water_group_number = ""
	cmdline_args = "-s %s -pname Na -np 13" % (tprfile)
	program = "genion"
	params.gmx_engine.run(program, cmdline_args, pipein="echo -1")
	[logfile, errlogfile]  = params.gmx_engine.lognames(program)
	infile = open(errlogfile,"r")
	for line in infile:
		field = line.split()
		if len(field)<4: continue
		if field[0]=='Group' and field[3]=='SOL)':
			water_group_number = field[1]
			break
	infile.close()
	if water_group_number == "":
		print "water_group_number not found in %s" % errlogfile
		exit(1)
	subprocess.call(["bash", "-c", "rm -f %s %s"%(logfile, errlogfile)])
	return water_group_number



