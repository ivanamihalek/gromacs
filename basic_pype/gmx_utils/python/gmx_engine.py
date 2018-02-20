import subprocess, os

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
		command = ['bash', '-c', "(source %s && which %s) | wc -l " % (self.gmx_bash, 'gmx')]
		proc = subprocess.Popen(command, stdout=subprocess.PIPE)
		# communicate() returns a tuple (stdoutdata, stderrdata)
		# and closes the input pype for the subprocess
		number_of_returned_lines =  int(proc.communicate()[0])
		if number_of_returned_lines==0:
			print "gmx not found in the path specified in", self.gmx_bash
			exit(1)
		if number_of_returned_lines>1:
			print "multiple gmx found in the path specified in", self.gmx_bash
			exit(1)
		self.gmx = "source %s && gmx" % self.gmx_bash # making sure that gmx sees the expected environment
		return

	def run(self, program, cmdline_args, message=None, higher_level_log=None):
		if message: print "\t running %s: %s" % (program, message)
		log    = open("%s.log"%program, "w")
		errlog = open("%s.err.log"%program, "w")
		cmd  = "%s %s " % (self.gmx, program)
		cmd += cmdline_args
		if higher_level_log: higher_level_log.write(cmd+"\n")
		failed = False
		try:
			# call blocks, Popen does not (we want to wait, because the output is the input for the next step)
			subprocess.call(["bash", "-c", cmd], stdout=log, stderr=errlog)
		except subprocess.CalledProcessError as errc:
			print "error error running", cmd
			print errc.output
			failed = True
		log.close()
		errlog.close()
		# delete gromacs junk
		subprocess.Popen(["bash", "-c", "rm -f \#*"])
		if failed: exit(2)
		return
