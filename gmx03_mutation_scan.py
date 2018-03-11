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

from gmx01_core import core

#########################################
def read_mutations(params):
	os.chdir(params.run_options.workdir)
	if (params.run_options.xtra == "none"):
		print "this pipe expects xtra file: list of mutations"
		exit()
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
	structure = get_biopython_structure(params)
	# TODO what if I have multiple chains
	for chain in structure.get_chains():
		seq_no = 0
		seq_id = {}
		for residue in chain:
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
	replace_sidechains(params, {seq_id: new_aa}, original_pdb, '', new_pdb)

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
	# create a super-input dir with fixed pdb and mdp files
	params.rundirs = Namespace()
	params.rundirs.in_dir = "super-input"
	fill_in_dir(params)

	######################
	# put fixed pdb in super-input
	chain_breaks, missing_sidechains = check_pdb_for_missing_atoms(params)
	if missing_sidechains: fix_sidechains(params, missing_sidechains)
	#  TODO: fix chain breaks (if they exist)
	# if chain_breaks: fix_chain_breaks(params, missing_sidechains)

	######################
	# adjust the run length
	params.physical.set_run_lengths(params, em_steep=10000, em_lbfgs=100,
											pr_nvt="50ps", pr_npt="10ps", md="100ns")

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
	#
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
