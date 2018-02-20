
#########################################
class GmxParameters:

	def __init__(self, run_options):  # maybe one day I can have some parameters passed for the most typical runs
		self.forcefield  = "amber99sb" # for simulations involving ligands
		# self.forcefield = "oplsaa"
		# self.forcefield = "gmx"
		self.water = "tip3p"
		# self.box_type   = "cubic"
		self.box_type    = "triclinic"
		self.box_edge    = "1.2" # distance of the box edge from the molecule in nm
		self.neg_ion     = "Cl" # theses names depend on the choice of the forcefield (check "ions.itp")
		self.pos_ion     = "Na"
		self.genion_solvent_code = "12"

		if self.forcefield=="oplsaa" or self.forcefield=="amber99sb":
			self.neg_ion  = "CL"
			self.pos_ion  = "NA"

		if run_options.pdb=="none":
			if run_options.posres:
				self.box_edge = "1.2"
			else:
				self.box_edge = "1.6"
		else:
			self.box_edge = "0.7"

		# The default solvent is Simple Point Charge water(SPC), with coordinates from $GMXLIB / spc216.gro.
		# These coordinates can also be used for other 3-site water models,
		# since a short equibilibration will remove the small differences between the models
		self.water_coords_gro = "spc216.gro"
