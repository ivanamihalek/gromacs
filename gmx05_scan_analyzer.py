#!/usr/bin/python -u

import os, subprocess,time
import matplotlib.pyplot as plt
import colorsys, random

from argparse import Namespace
from gmx_lib import run_setup
from gmx_lib.gmx_engine import GmxEngine


# mutational scan: simple run for each mutated structure
#########################################
def main():

	graph_only = False

	params = Namespace()
	params.run_options  = run_setup.parse_commandline()
	params.gmx_engine   = GmxEngine("/usr/local/gromacs/bin/GMXRC.bash")
	params.command_log  = open(params.run_options.workdir+"/analyzer_commands.log","w")

	######################
	# take that the mutation dir has the name pof the format <pos>_<from>_<to>
	os.chdir(params.run_options.workdir)
	subdirs = [name for name in os.listdir(params.run_options.workdir) if os.path.isdir(name)]
	subdirs.remove('super-input')
	subprocesses = {}
	for subdir in subdirs:
		field = subdir.split("_")
		if len(field)!=3: continue
		subdir_path = "/".join([params.run_options.workdir, subdir])
		os.chdir(subdir_path)
		cmd = "ls 06_production"
		proc = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE)
		trjfile = None
		for line in proc.stdout.readlines():
			if "md_out.trr" in line:
				trjfile = line.strip()
		if not trjfile:
			print "\t trajectory not found"
			continue
		if graph_only: continue
		postprocessor = "/".join([params.run_options.gromacs_pype_home,"gmx_lib", "postproduction.py"])
		pdbname = params.run_options.pdb+"."+subdir
		cmd = "nice %s -p %s -w %s  > /dev/null 2>&1" % (postprocessor, pdbname, subdir_path)
		params.command_log.write(cmd+"\n")
		subprocesses[subdir] = subprocess.Popen(["bash", "-c", cmd], stdout=None, stderr=None,  close_fds=True)

	# check all processes done
	done = graph_only
	while not done:
		done = True
		time.sleep(3)
		print "checking ..."
		for subdir,p in subprocesses.iteritems():
			if p.poll() == None:
				# the subprocess is still alive
				done = False
	print "all done"

	# plot
	os.chdir(params.run_options.workdir)
	# see https://en.wikipedia.org/wiki/Hue
	hue_step = 240.0/len(subdirs);
	hue = 0
	for subdir in subdirs:
		xvg_file = params.run_options.pdb+"."+subdir+".hbonds.xvg"
		with open("/".join([subdir,"07_post-production", xvg_file])) as f:
			lines = f.readlines()
			# col 1 are hydrogen bonds, and col 2 are apirs within 0.35nm
			# col 0 is time in ps
			x = [int(line.split()[0])/1000 for line in lines]
			y = [(float(int(line.split()[1])) + 2*random.uniform(0, 1)-1) for line in lines]
			hue += hue_step
			markerfacecolor = colorsys.hsv_to_rgb(hue/360.0, 1.0, 1.0)
			linecolor = colorsys.hsv_to_rgb(random.uniform(0, 1), 1.0, 1.0)
			markeredgecolor  = colorsys.hsv_to_rgb(random.uniform(0, 1), 1.0, 1.0)
			plt.plot(x,y, marker='o',
			         markerfacecolor= markerfacecolor,
			         markeredgecolor= markeredgecolor,
			         linestyle = '-', color=linecolor,
			         label=subdir)
	plt.legend()
	plt.ylim(115, 160)
	plt.xlim(-1,8)
	plt.show()




#########################################
if __name__ == '__main__':
	main()
