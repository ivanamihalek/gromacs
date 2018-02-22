import os, subprocess

#########################################
def generate(params):
	# change to topology directory
	currdir = params.rundirs.top_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname = params.run_options.pdb
	pdbfile = pdbname+".pdb"
	topfile = pdbname+".top"
	grofile = pdbname+".gro"
	if os.path.exists(topfile) and os.path.exists(grofile):
		print "\t %s and %s found" % (topfile, grofile)
		return
	program = "pdb2gmx"

	#################################
	# create topology for the protein
	if pdbname!="none":
		cmdline_args = "-ignh  -ff  %s  -water %s " % (params.physical.forcefield, params.physical.water)
		cmdline_args += "-f ../%s/%s -o %s -p %s" % (params.rundirs.in_dir, pdbfile,  grofile, topfile)
		# if there are salt bridges or termini modifications,
		# create file with input options for pdb2gmx, and indicate so on the command line
		# -ss is interactive but candidates must exist
		if os.path.exists("pdb2gmx_in"): os.remove("pdb2gmx_in")
		if os.path.exists("ssbridges"):
			subprocess.call(["bash","-c", "echo ssbridges > pdb2gmx_in"], stdout=subprocess.PIPE)
			cmdline_args  += " -ss"
		# notes on terminal mods:
		# the  -ter flag will actually prompt for manual input  - we are funneling in the answers
		# the default ins NH3+ and COO-, which might be the best choice enabling the hydrogen
		# bonds in the secondary structure
		if os.path.exists("termini"): # amber forcefield is incompatible with ter
			subprocess.call(["bash", "-c", "echo termini >> pdb2gmx_in"], stdout=subprocess.PIPE)
			cmdline_args += " -ter"

		if os.path.exists("pdb2gmx_in"): cmdline_args += " < pdb2gmx_in"
		params.command_log.write("in %s:\n"%(currdir))

		##########################################
		# here is where the run actually happens:
		params.gmx_engine.run(program, cmdline_args,"generating top and gro files for the peptide", params.command_log)
		# the things in bracket are some false alarm msgs gmx sometimes issues
		params.gmx_engine.check_logs_for_error(program,["masses will be determined based on residue and atom names"])
		##########################################

		if os.path.exists("pdb2gmx_in"): subprocess.Popen(["bash", "-c", "rm -f pdb2gmx_in"])

	# TODO: create top for ligands

	# merge

	return

