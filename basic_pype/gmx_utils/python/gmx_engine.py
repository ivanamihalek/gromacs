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
		return


	###########################
	@staticmethod
	def lognames(program): # static method does not need 'self'
		return ["%s.log"%program,"%s.err.log"%program]

	###########################
	def check_logs_for_error(self, program, tolerated_error_msgs):
		# it should be in logerr, but I do not trust these kids to be systematic
		[logname,errlogname] = self.lognames(program)
		for logf in [logname,errlogname]:
			# gcq is the funny quote, which funnily enough can contain the word error in it
			cmd = "grep -i error %s | grep -v gcq"%logf
			proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE)
			for line in proc.stdout.readlines():
				if not line in tolerated_error_msgs:
					print
					print "%s ended in error state. Please check the file  %s.\n" % (program, logf)
					print line
					exit(3)
		return

	###########################
	def run(self, program, cmdline_args, message=None, higher_level_log=None, pipein=None):
		if message: print "\t running %s: %s" % (program, message)
		[logname,errlogname] = self.lognames(program)
		log    = open(logname, "w")
		errlog = open(errlogname, "w")
		# sourcing the bash script ensures that  gmx sees the expected environment
		if pipein:
			cmd    = "source %s && (%s | gmx %s %s)" % (self.gmx_bash, pipein, program, cmdline_args)
		else:
			cmd    = "source %s && (gmx %s %s)" % (self.gmx_bash, program, cmdline_args)

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

	###########################
	def convergence_line(self, program):
		[logname,errlogname] = self.lognames(program)
		cmd  = "grep converge %s | tail -n1" % errlogname
		return subprocess.Popen(["bash", "-c", cmd]).communicate()[0]
