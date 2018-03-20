#!/usr/bin/python -u

import os, subprocess
import collections

from argparse import Namespace
from gmx_lib import run_setup
from gmx_lib.run_setup  import WorkdirStructure
from gmx_lib.gmx_engine import GmxEngine
from gmx_lib.gmx_params import GmxParameters
from modeling_lib.scwrl_engine import ScwrlEngine
from modeling_lib.pdb_doctor import get_biopython_structure, replace_sidechains
from modeling_lib.pdb_doctor import check_pdb_for_missing_atoms, fix_sidechains


#########################################
def read_mutations(params):

	os.chdir("/".join([params.run_options.workdir,params.rundirs.in_dir]))
	mutfile = open(params.run_options.xtra, "r")
	mutation = {}
	for line in mutfile:
		line=line.strip()
		aa_from=line[:1]
		res_id= int(line[1:-1])
		aa_to = line[-1:]
		mutation[res_id] = [aa_from, aa_to]
	if len(mutation)==0:
		print "no mutations found in %s (?)" % params.run_options.xtra
		exit()
	sorted_keys = mutation.keys()
	sorted_keys.sort()
	mutation_sorted = collections.OrderedDict()
	for key in sorted_keys:
		mutation_sorted[key] = mutation[key]

	return mutation_sorted


#########################################
def fill_in_dir(params):
	os.chdir(params.run_options.workdir)
	in_dir = params.rundirs.in_dir
	if not os.path.exists(in_dir): os.mkdir(in_dir)

	# move the list of mutations to super-input directory
	os.rename(params.run_options.xtra, "/".join([in_dir, params.run_options.xtra]))

	os.chdir(in_dir)
	pdbname = params.run_options.pdb
	if os.path.exists("../%s.pdb"%pdbname): os.rename("../%s.pdb"%pdbname,"./%s.pdb"%pdbname)
	if not os.path.exists("%s.pdb"%pdbname):
		print "%s.pdb  not found"%pdbname
		exit()
	subprocess.call(["bash", "-c", "cp -f %s/* . " % (params.run_options.mdp_template_home)])


#########################################
def remove_nonexistent_positions(params, mutations):

	structure = get_biopython_structure(params)
	# TODO what if I have multiple chains
	for chain in structure.get_chains():
		residue_ids_present = []
		for residue in chain:
			residue_ids_present.append(int(residue.id[1]))

	for key in mutations.keys():
		if not key in residue_ids_present:
			del  mutations[key]
	if len(mutations) == 0:
		print "no positions from the mutation list seem to be present on the structure"
		exit()
	return


#########################################
def res_id_to_seq_id(params):
	biopython_structure = get_biopython_structure(params)
	chains = [c for c in biopython_structure.get_chains()]
	if len(chains)==1:
		chain=chains[0]
	elif params.run_options.chain=="none":
		print "when there are multiple chains in the input pdb, ",
		print "a chain must be specified to do mutation scan"
		exit(1)
	else:
		chain=params.run_options.chain
	# 0 here is mor the first (and only, in this case, model)
	chain_object = biopython_structure[0][chain]

	seq_no = 0
	seq_id = {}
	for residue in chain_object:
		seq_id[int(residue.id[1])] = seq_no
		seq_no += 1
	return seq_id

#########################################
def make_mutant(params,newdir, seq_id, new_aa):
	os.chdir(params.run_options.workdir)
	if not os.path.exists(newdir): os.mkdir(newdir)
	os.chdir(newdir)
	# mutation
	pdbname = params.run_options.pdb
	original_pdb = "../%s/%s.pdb"%(params.rundirs.in_dir, pdbname)
	new_pdb = "%s.%s.pdb"%(pdbname, newdir)
	# this function will call scwrl
	replace_sidechains(params, {seq_id: new_aa}, original_pdb, params.run_options.chain, new_pdb)


# mutational scan: simple run for each mutated structure
#########################################
def main():

	params = Namespace()
	params.run_options  = run_setup.parse_commandline()
	params.physical     = GmxParameters(params.run_options)
	params.gmx_engine   = GmxEngine("/usr/local/gromacs/bin/GMXRC.bash")
	params.scwrl_engine = ScwrlEngine("/usr/local/bin/scwrl4/Scwrl4")
	params.command_log  = open(params.run_options.workdir+"/commands.log","w")

	######################
	# special requirement for mutation scan:
	if (params.run_options.xtra == "none"):
		print "this pipe expects xtra file: list of mutations"
		exit()

	######################
	# create a super-input dir with fixed pdb and mdp files
	params.rundirs = Namespace()
	params.rundirs.in_dir = "super-input"
	fill_in_dir(params)

	######################
	# pdb is assumed to be ok
	# TODO: add pdb cleanup (not that here there might be multiple, non-peptide chains

	######################
	# adjust the run length
	params.physical.set_run_lengths(params,  em_steep=12000, em_lbfgs=100,
											pr_nvt="50ps", pr_npt="10ps", md="2000ps")

	if True: params.physical.set_annealing_schedule(params, annealing_type='single',
								annealing_npoints='8',
								annealing_time= "0    100  200 300 400 500 600 700",
								annealing_temp= "300  320  340 360 380 400 420 440")


	######################
	# read in the mutation list
	mutations = read_mutations(params)
	# remove mutations in positions outside of the given structure
	remove_nonexistent_positions(params, mutations)

	######################
	# for each mutation create directory and the mutant structure
	seq_id = res_id_to_seq_id(params)
	newdirs = []
	for resid, mut_pair in mutations.iteritems():
		newdir = "{}_{}_{}".format(str(resid).zfill(3), mut_pair[0],  mut_pair[1])
		newdirs.append(newdir)
		make_mutant(params, newdir, seq_id[resid], mut_pair[1])

	######################
	# create dir tree, fill input for each new directory
	for newdir in newdirs:
		os.chdir(params.run_options.workdir)
		if not os.path.exists(newdir): os.mkdir(newdir)
		os.chdir(newdir)
		for dir in WorkdirStructure.workdir_names:
			if not os.path.exists(dir): os.mkdir(dir)
		in_dir = WorkdirStructure.workdir_names[0]
		# move pdb to in_dir
		subprocess.call(["bash", "-c", "mv -f *.pdb %s"%in_dir], stdout=None, stderr=None)
		# move mdp files to in_dir
		super_in_dir = "/".join([params.run_options.workdir, params.rundirs.in_dir])
		subprocess.call(["bash", "-c", "cp -f %s/*.mdp %s"%(super_in_dir, in_dir)], stdout=None, stderr=None)
		# cleanup
		subprocess.call(["bash", "-c", "rm -f chain*seq scwrl*"], stdout=None, stderr=None)

	######################
	# start the simulation in each sub-dir
	for newdir in newdirs:
		workdir = "/".join([params.run_options.workdir, newdir])
		core_pipe = "/".join([params.run_options.gromacs_pype_home,"gmx01_core.py"])
		pdbname = params.run_options.pdb+"."+newdir
		cmd = "nice %s -p %s -w %s" % (core_pipe, pdbname, workdir)
		params.command_log.write(cmd+"\n")
		subprocess.Popen( ["bash", "-c", cmd], stdout=None, stderr=None,  close_fds=True)

	######################
	# cleanup and exit
	######################
	params.command_log.close()
	return True


#########################################
if __name__ == '__main__':
	main()
