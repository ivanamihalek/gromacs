import subprocess, os

#########################################
class  ScwrlEngine:

	def __init__(self, scwrl_path):
		if not os.path.exists(scwrl_path):
			print scwrl_path, "not found"
			exit(2)
		if not os.path.isfile(scwrl_path):
			print scwrl_path, "not a file"
			exit(3)
		if not os.access(scwrl_path, os.X_OK):
			print scwrl_path, "not executable"
			exit(3)

		self.scwrl_path = scwrl_path

	###########################
	@staticmethod
	def lognames(): # static method does not need 'self'
		return ["scwrl.log","scwrl.err.log"]

	###########################
	# seq_file specifies the residues that should be kept fixed
	def run(self, in_pdb, out_pdb, seq_file, message=None, higher_level_log=None):
		if message: print "\t running scwrl: %s" % (message)
		[logname,errlogname] = self.lognames()
		log    = open(logname, "w")
		errlog = open(errlogname, "w")
		# -h suppresses the output of hydrogen atoms
		cmd    = "%s  -i %s  -o %s  -s %s -h" % (self.scwrl_path, in_pdb, out_pdb,  seq_file)
		if higher_level_log: higher_level_log.write(cmd+"\n")
		failed = False
		try:
			# call blocks, Popen does not (we want to wait,
			# because the output is the input for the next step)
			subprocess.call(["bash", "-c", cmd], stdout=log, stderr=errlog)
		except subprocess.CalledProcessError as errc:
			print "error error running", cmd
			print errc.output
			failed = True

		log.close()
		errlog.close()
		if failed: exit(2)
		return


