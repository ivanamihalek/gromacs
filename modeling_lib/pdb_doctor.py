
import Bio.PDB as biopython
import os, subprocess, shutil

#########################################
def get_biopython_structure(params):
	currdir = params.rundirs.in_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname      = params.run_options.pdb
	parser = biopython.PDBParser(QUIET=True) # silence the crappy warnings
	structure = parser.get_structure(pdbname, pdbname+".pdb")
	# TODO: error state?
	return structure

#########################################
def init_chain_hash (chain_hash, chain_id):
	if not chain_hash: chain_hash = {}
	if not chain_hash.has_key(chain_id) or chain_hash[chain_id]==None:
		chain_hash[chain_id] = []
	return chain_hash

#########################################
def check_pdb_for_missing_atoms(params):

	structure = get_biopython_structure(params)
	number_of_atoms = { "ASN":8, "GLU":9, "CYS":6, "THR":7, "TYR":12, "PRO":7, "HIS":10,
						"SER":6, "LEU":8, "TRP":14, "PHE":11, "MET":8, "ALA":5, "GLY":4,
						"VAL":7, "GLN":9, "LYS":9, "ILE":8, "ASP":8, "ARG":11}
	# TODO: atypical residues, terminal residues
	chain_breaks       = None
	missing_sidechains = None
	for chain in structure.get_chains():
		sequential_id = 0
		prev_res_id = None
		for residue in chain:
			# residue id is a triple; [1] is the id assigned by PDB
			# here can check for the gap in the residue_id sequence (missing residues)
			atom_list = biopython.Selection.unfold_entities(residue, 'A')

			# what should be the number of atoms in each residue
			# check for missing sidechains
			if len(atom_list)!= number_of_atoms[residue.get_resname()]:
				print "missing sidechain:"
				print chain.id, residue.get_resname(), residue.id[1], sequential_id
				print "\t", len(atom_list), number_of_atoms[residue.get_resname()]
				missing_sidechains = init_chain_hash (missing_sidechains, chain.id)
				missing_sidechains[chain.id].append(sequential_id)
			resid = int(residue.id[1])
			if prev_res_id and (resid-prev_res_id)>1:
				print 'chain break:'
				print chain.id, residue.get_resname(), residue.id[1], sequential_id
				chain_breaks = init_chain_hash (chain_breaks, chain.id)
				chain_breaks[chain.id].append(sequential_id)
			prev_res_id = resid
			sequential_id += 1
	return chain_breaks, missing_sidechains


#########################################
def make_scwrl_input_sequence(params, positions, original_pdb, chain):
	# get sequence
	pdb2seq = "%s/pdb2seq.pl"%params.run_options.perl_utils
	stdout  = subprocess.Popen([pdb2seq, original_pdb, chain], stdout=subprocess.PIPE)
	# lowercase positions that stay the same, uppercase for those that change
	seqlist = list(stdout.communicate()[0].replace("\n","").lower())
	if type(positions)==list:
		for pos in positions: seqlist[pos] = seqlist[pos].upper()
	elif type(positions)==dict:
		for pos in positions.keys():  seqlist[pos] = positions[pos].upper()
	else:
		print "type of 'positions' variable not recognized"
		exit()
	seq  = ''.join(seqlist)
	# write it to a file (uppercase are residues with sidechains requiring fixing)
	scwrl_input_seq_name = "chain%s.seq"%chain
	outf = open(scwrl_input_seq_name, "w")
	outf.write(seq)
	outf.close()
	return scwrl_input_seq_name


#########################################
def replace_sidechains(params, positions, original_pdb, chain, new_pdb):
	# split the input file into one that contains the chain, and the other one
	# that contains the rest of the structure
	pdb_extract = "%s/pdb_extract_chain.pl"%params.run_options.perl_utils
	tmp_chain_pdb    = "chain%s.pdb" % chain
	tmp_notchain_pdb = "not_chain%s.pdb" % chain
	subprocess.call([pdb_extract, original_pdb, "-c%s"%chain, "-o%s"%tmp_chain_pdb], stdout=None, stderr=None)
	subprocess.call([pdb_extract, original_pdb, "-c%s"%chain, "-o%s"%tmp_notchain_pdb, "-i"], stdout=None, stderr=None)

	# get the input file indicating which positions need to be mutated
	scwrl_input_seq_name = make_scwrl_input_sequence(params, positions, tmp_chain_pdb, chain)
	# here: scwrl run
	mutated_chain_pdb = "mutated_chain%s.pdb" % chain
	params.scwrl_engine.run(tmp_chain_pdb, mutated_chain_pdb, scwrl_input_seq_name,
							message="fixing sidechains", higher_level_log=params.command_log)
	# concatenate mutated peptide and the rest of the pdb into new file
	command = ['bash', '-c', "cat %s %s > %s" % (mutated_chain_pdb, tmp_notchain_pdb, new_pdb)]
	subprocess.call(command, stdout=None, stderr=None)
	# cleanup
	command = ['bash', '-c', "rm %s %s %s" % (mutated_chain_pdb, tmp_notchain_pdb, tmp_chain_pdb)]
	subprocess.call(command, stdout=None, stderr=None)

#########################################
def fix_sidechains(params, hash_of_residue_positions):
	# TODO: make sure this actually works for multiple chains
	currdir = params.rundirs.in_dir
	os.chdir("/".join([params.run_options.workdir, currdir]))
	pdbname  = params.run_options.pdb
	pdb_cleanup_dir = "pdb_cleanup"
	original_pdb = "%s.orig.pdb"%pdbname
	if not os.path.exists(pdb_cleanup_dir): os.mkdir(pdb_cleanup_dir)
	shutil.copyfile(pdbname+".pdb", original_pdb)
	os.chdir(pdb_cleanup_dir)

	for chain, positions in hash_of_residue_positions.iteritems():
		replace_sidechains(params, positions,  "../"+original_pdb, chain, "%s.pdb"%pdbname)
		os.rename("%s.pdb"%pdbname, "../%s.pdb"%pdbname)

