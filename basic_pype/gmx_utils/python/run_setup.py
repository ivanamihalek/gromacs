import argparse
import sys, os
import pprint
import subprocess
#########################################
def parse_commandline():
	if len(sys.argv)<=1:
		print "Usage:"
		print "%s -h" % sys.argv[0]
		print "or"
		print "%s -p <pdb_root> [-l <ligand list>] [-w <workdir>] [--remd] [--min] [--posres]" % sys.argv[0]
		exit(1)

	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pdb', metavar='pdbfile', type=str,  default="none",
						dest="pdb", help='input pdb file (peptide, ions, and ligands)')
	parser.add_argument('-l','--ligands', metavar='ligand', type=str, nargs='+', default=[],
						dest="ligands", help='list of ligand names')
	parser.add_argument('-w','--workdir', metavar='workdir', type=str, default=".",
						dest="workdir", help='work directory - check WorkdirStructure for the expected structure')
	parser.add_argument('-m','--min', action='store_true', default=False,
						dest="minimization", help='energy minimization run (defaults to False)')
	parser.add_argument('-r','--remd', action='store_true', default=False,
						dest="remd", help='REMD run (defaults to False)')
	parser.add_argument('-s','--posres', action='store_true', default=False,
						dest="posres", help='position restrained run (defaults to False)')

	run_options = parser.parse_args(sys.argv[1:])
	return run_options


#########################################
class GmxEngine:

	# we rely on everything being accesible from the paths set in GMXRC.bash
	# that is, on what is currently the standard way of installing gromacs
	def __init__(self, gmx_bash):
		if not os.path.exists(gmx_bash):
			print gmx_bash, "not found"
			exit(2)
		# check if the remaining gmx executables can be found in the paths specified in gmx_bash
		# (TODO: gmx_bash will also set paths for libraries which we are not checking)
		self.gmx_bash = gmx_bash
		self.gmx_executables = ['gmx', 'blub']
		for executable in self.gmx_executables:
			command = ['bash', '-c', "(source %s && which %s) | wc -l " % (self.gmx_bash, executable)]
			proc = subprocess.Popen(command, stdout=subprocess.PIPE)
			# communicate() returns a tuple (stdoutdata, stderrdata)
			# and closes the input pype for the subprocess
			number_of_returned_lines =  int(proc.communicate()[0])
			if number_of_returned_lines==0:
				print executable, "not found in the path specified in", self.gmx_bash


#########################################
class Parameters:

	def __init__(self, run_options):  # maybe one day I can have some parameters passed for the most typical runs
		self.forcefield  = "amber99sb"
		self.water = "tip3p"
		# self.forcefield = "oplsaa"
		# self.forcefield = "gmx"
		# self.box_type   = "cubic"
		self.box_type    = "triclinic"
		self.box_edge    =  1.2 # distance of the box edge from the molecule in nm
		self.neg_ion     = "Cl" # theses names depend on the choice of the forcefield (check "ions.itp")
		self.pos_ion     = "Na"
		self.genion_solvent_code = 12

		if self.forcefield=="oplsaa" or  self.forcefield=="amber99sb":
			self.neg_ion  = "CL"
			self.pos_ion  = "NA"

		if run_options.pdb=="none":
			if run_options.posres:
				self.box_edge = 1.2
			else:
				self.box_edge = 1.6
		else:
			self.box_edge = 0.7

#########################################
class WorkdirStructure:

	def __init__(self, run_options, workdir_names):
		self.check_workdir_existence(run_options.workdir, workdir_names)
		self.in_dir  = workdir_names[0]
		self.top_dir = workdir_names[1]
		self.em1_dir = workdir_names[2]
		self.em2_dir = workdir_names[3]
		self.pr1_dir = workdir_names[4]
		self.pr2_dir = workdir_names[5]
		self.production_dir =  workdir_names[6]
		## check run-specific setup requirements
		wp = run_options.workdir
		if not run_options.pdb=="none":
			self.exist_or_die("/".join([wp, self.in_dir, run_options.pdb+".pdb"]))
		for mdp in ["em_steep.mdp", "em_lbfgs.mdp"]:
			self.exist_or_die("/".join([wp, self.in_dir, mdp]))
		if not run_options.minimization:
			self.exist_or_die("/".join([wp, self.in_dir, "pr_nvt.mdp"]))
		if not run_options.minimization and not run_options.posres:
			for mdp in ["pr_npt.mdp", "md.mdp"]:
				self.exist_or_die("/".join([wp, self.in_dir, mdp]))

	def exist_or_die(self, fullpath):
		if not os.path.exists(fullpath):
			print fullpath, "not found"
			exit(1)
	def bedir_or_die(self, fullpath):
		if not os.path.isdir(fullpath):
			print fullpath, "not found"
			exit(1)
	def check_workdir_existence (self, wkdir_path, workdir_names):
		for wd in workdir_names:
			self.exist_or_die(wkdir_path+"/"+wd)
			self.bedir_or_die(wkdir_path+"/"+wd)