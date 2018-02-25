import os, subprocess

from ..python import grompp

#########################################
def run(params):

	stage = 'production'

	#grompp = generate parametrized topology file (tpr; binary; complete input compiled)
	grompp.generate(params, stage) # use the previous step to construct the tpr file

	# change to  directory for this stage
	currdir = params.rundirs.name[stage]
	os.chdir("/".join([params.run_options.workdir, currdir]))

	pdbname      = params.run_options.pdb
	tprfile_in   = pdbname+".md_input.tpr"
	grofile_out  = pdbname+".md_out.gro"
	traj_out     = pdbname+".md_out.trr"
	edrfile_out  = pdbname+".md_out.edr"
	native_log   = pdbname+".md_native.log"
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
	cmdline_args  = " -s %s -c %s -nt 1 -o %s -e %s -g %s" % \
					(tprfile_in, grofile_out, traj_out, edrfile_out, native_log)

	params.command_log.write("in %s:\n" % (currdir))
	params.gmx_engine.run(program, cmdline_args, "production run", params.command_log)
	false_alarms = ["defaults to zero instead of generating an error"]
	params.gmx_engine.check_logs_for_error(program, false_alarms)

	print "\t ", params.gmx_engine.convergence_line(program)


	# if we safely got to here, we'll assume that we do not need the checkpoint files
	# (the error checks exit with non-zero code, and thus the cpt files will stay)
	subprocess.call(["bash","-c","rm -f *.cpt"])
