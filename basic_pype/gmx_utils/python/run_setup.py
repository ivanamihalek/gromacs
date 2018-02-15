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
		print "%s -p <pdb_root> [-l <ligand list>] [--remd] [--min] [--posres]" % sys.argv[0]
		exit(1)

	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pdb', metavar='pdbfile', type=str,  default="none",
						dest="pdb", help='input pdb file (peptide, ions, and ligands)')
	parser.add_argument('-l','--ligands', metavar='ligand', type=str, nargs='+', default=[],
						dest="ligands", help='list of ligand names')
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
			command = ['bash', '-c', "(source %s && which %s) | wc -l  " % (self.gmx_bash, executable)]
			proc = subprocess.Popen(command, stdout=subprocess.PIPE)
			# communicate() returns a tuple (stdoutdata, stderrdata)
			# and closes the input pype for the subprocess
			number_of_returned_lines =  int(proc.communicate()[0])
			if number_of_returned_lines==0:
				print executable, "not found in the path specified in", self.gmx_bash




