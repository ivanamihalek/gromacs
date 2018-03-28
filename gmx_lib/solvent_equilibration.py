import os, subprocess

import grompp

#########################################
# this is a hack to find all restraint files
# we have available
def itp_files(params):
	proc = subprocess.Popen(["bash", "-c", "grep include ../%s/*top | grep -v %s"%(params.rundirs.top_dir,params.physical.forcefield)],
					stdout=subprocess.PIPE, stderr=None)
	itps = []
	for line in proc.stdout.readlines():
		itps.append(line.rstrip().replace("#","").replace('"','').replace("include","").replace(" ",""))
	return itps

#########################################
def run(params, stage):

	#grompp = generate parametrized topology file (tpr; binary; complete input compiled)
	grompp.generate(params, stage, position_restrained_run=True) # use the previous step to construct the tpr file

	# change to  directory for this stage
	currdir = params.rundirs.name[stage]
	os.chdir("/".join([params.run_options.workdir, currdir]))

	# apparently there is no other way to pass yhe restraints to mdrun:
	# (alternatively I could have the *top file each time,
	# it contains the line #include "posre.itp")
	for itp in itp_files(params):
		subprocess.call(["bash", "-c", "ln -sf ../%s/%s %s"%(params.rundirs.top_dir,itp,itp)],
					stdout=None, stderr=None)

	pdbname      = params.run_options.pdb
	tprfile_in   = pdbname+".pr_input.tpr"
	grofile_out  = pdbname+".pr_out.gro"
	traj_out     = pdbname+".pr_out.trr"
	edrfile_out  = pdbname+".pr_out.edr"
	native_log   = pdbname+".pr_native.log"
	if os.path.exists(grofile_out):
		print "\t %s found" % (grofile_out)
		return
	for infile in [tprfile_in]:
		if not os.path.exists(infile):
			print "\t in local_energy_minimum.find(%s): %s not found (?)" % (stage, infile)
			exit(1)

	program = "mdrun" # nt 1; run multiple trajectories instead
	# mdrun will produce trajectory and edr (energy) files,  whether we ask for it or not,
	# so we might just as well name them so we can remove them later
	cmdline_args  = " -s %s -c %s  -nt 4 -o %s -e %s -g %s" % \
	 				(tprfile_in, grofile_out, traj_out, edrfile_out, native_log)
	params.command_log.write("in %s:\n" % (currdir))
	msg = "solvent equilibration with peptide restrained"
	params.gmx_engine.run(program, cmdline_args, msg, params.command_log)

	# check for errors
	false_alarms = ["masses will be determined based on residue and atom names"]
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	# check convergence
	print "\t ", params.gmx_engine.convergence_line(program)

	if os.path.exists(traj_out): os.remove(traj_out)
	if os.path.exists(edrfile_out): os.remove(edrfile_out)


