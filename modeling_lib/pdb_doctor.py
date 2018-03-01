
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
def replace_sidechains(params, positions, original_pdb, chain, new_pdb):
	# get sequence
	pdb2seq = "%s/pdb2seq.pl"%params.run_options.perl_utils
	stdout  = subprocess.Popen([pdb2seq, original_pdb, chain], stdout=subprocess.PIPE)
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
	outf = open("chain%s.seq"%chain, "w")
	outf.write(seq)
	outf.close()
	#  here: scwrl run ; scwrl loses the b-factor - not sure whether that could ever matter
	params.scwrl_engine.run(original_pdb, new_pdb, "chain%s.seq"%chain,
							message="fixing sidechains", higher_level_log=params.command_log)


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

