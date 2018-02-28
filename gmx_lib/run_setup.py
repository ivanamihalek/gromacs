import argparse
import sys, os
import subprocess
#########################################
def parse_commandline():
	if len(sys.argv)<=1:
		print "Usage:"
		print "%s -h" % sys.argv[0]
		print "or"
		print "%s -p <pdb_root> [-l <ligand list>] [-w <workdir>] [--remd] [--min] [--posres]" % sys.argv[0]
		exit(1)
	pwd = os.getcwd()
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pdb', metavar='pdbfile', type=str,  default="none",
						dest="pdb", help='input pdb file (peptide, ions, and ligands)')
	parser.add_argument('-l','--ligands', metavar='ligand', type=str, nargs='+', default=[],
						dest="ligands", help='list of ligand names')
	parser.add_argument('-w','--workdir', metavar='workdir', type=str, default=pwd,
						dest="workdir", help='work directory - check WorkdirStructure for the expected structure')
	parser.add_argument('-m','--min', action='store_true', default=False,
						dest="minimization", help='energy minimization run (defaults to False)')
	parser.add_argument('-r','--remd', action='store_true', default=False,
						dest="remd", help='REMD run (defaults to False)')
	parser.add_argument('-s','--posres', action='store_true', default=False,
						dest="posres", help='position restrained run (defaults to False)')

	run_options = parser.parse_args(sys.argv[1:])

	# basic sanity check:
	if run_options.pdb=="none" and run_options.ligands.length()==0:
		print "No protein, no small molecule ... what are we doing here?"
		exit(1)

	# home directory of the gromac pipeline # NOTE THE HARDCODED FILENAME HERE
	run_options.gromacs_pype_home = "/".join(sys.argv[0].split("/")[:-1])
	run_options.mdp_template_home = "/".join([run_options.gromacs_pype_home,"gmx00_templates"])
	run_options.perl_utils =  "/".join([run_options.gromacs_pype_home,"perl_utils"])
	return run_options


#########################################
class WorkdirStructure:

	workdir_names = ["00_input", "01_topology", "02_em_steepest", "03_em_lbfgs",
					"04_nvt_eq", "05_mpt_eq", "06_production", "07_post-production"]

	def __init__(self, run_options):
		self.check_workdir_existence(run_options.workdir, self.workdir_names)
		self.in_dir  = self.workdir_names[0]
		self.top_dir = self.workdir_names[1]
		self.em1_dir = self.workdir_names[2]
		self.em2_dir = self.workdir_names[3]
		self.pr1_dir = self.workdir_names[4]
		self.pr2_dir = self.workdir_names[5]
		self.production_dir =  self.workdir_names[6]
		self.post_dir=  self.workdir_names[7]

		self.name = {'water':self.top_dir, 'ions':self.top_dir,
					'em1':self.em1_dir, 'em2':self.em2_dir,
					'pr1':self.pr1_dir, 'pr2':self.pr2_dir,
					'production':self.production_dir}
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