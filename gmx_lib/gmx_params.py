import os

#########################################
class GmxParameters:

	def __init__(self, run_options):  # maybe one day I can have some parameters passed for the most typical runs
		#self.forcefield  = "amber99sb" # for simulations involving ligands
		self.forcefield  = "amber14sb" # for simulations involving RNA, perhaps also ligands
		# self.forcefield = "oplsaa"
		# self.forcefield = "gmx"
		self.water = "tip3p"
		# self.box_type   = "cubic"
		self.box_type    = "triclinic"
		self.box_edge    = "1.2" # distance of the box edge from the molecule in nm
		self.neg_ion     = "Cl" # these names depend on the choice of the forcefield (check "ions.itp")
		self.pos_ion     = "Na"
		self.genion_solvent_code = "12"

		if self.forcefield=="oplsaa" or "amber" in self.forcefield:
			self.neg_ion  = "CL"
			self.pos_ion  = "NA"

		if run_options.pdb=="none":
			if run_options.posres:
				self.box_edge = "1.2"
			else:
				self.box_edge = "1.6"
		else:
			self.box_edge = "0.7"

		self.timestep_in_fs = 2;

		# The default solvent is Simple Point Charge water(SPC), with coordinates from $GMXLIB / spc216.gro.
		# These coordinates can also be used for other 3-site water models,
		# since a short equibilibration will remove the small differences between the models
		self.water_coords_gro = "spc216.gro"


	#################################
	@staticmethod # static method does not need 'self'
	def to_number_of_steps(provided_value, timestep_in_fs):
		if type(provided_value)==int: return provided_value
		if type(provided_value)!=str:
			print "The value provided as number of steps must be either"
			print "and integer, or a string of the format \"%dunit\", "
			print "with unit one of fs, ps, ns or ms."
			exit(1)
		
		number_of_steps = 0
		time_unit = provided_value[-2:].lower()
		number_of_time_units = int(provided_value[:-2])
		if time_unit=="fs":
			number_of_steps = int(round(number_of_time_units/timestep_in_fs))
		elif time_unit=="ps":
			number_of_steps = int(round(number_of_time_units/timestep_in_fs)*1.0e3)
		elif time_unit=="ns":
			number_of_steps = int(round(number_of_time_units/timestep_in_fs)*1.0e6)
		elif time_unit=="ms":
			number_of_steps = int(round(number_of_time_units/timestep_in_fs)*1.0e9)
		else:
			print "in gmx_params: unrecognized time unit:", time_unit
			exit()

		return number_of_steps

	#################################
	def set_run_lengths(self, run_params, em_steep=100, em_lbfgs=100, pr_nvt=100, pr_npt=100, md=100):
		os.chdir("/".join([run_params.run_options.workdir, run_params.rundirs.in_dir]))
		dt = self.timestep_in_fs
		for name in ['em_steep', 'em_lbfgs']:
			filename = name+".mdp"
			outf = open ("tmp.mdp","w")
			inf =  open (filename,"r")
			number_of_steps = self.to_number_of_steps(eval(name), dt)
			for line in inf:
				parameter = line.lstrip()[:6]
				if parameter=="nsteps":
					outf.write("nsteps   =  %d \n"  % number_of_steps)
				else:
					outf.write(line)
			inf.close()
			outf.close()
			os.rename("tmp.mdp", filename)

		for name in ['pr_nvt', 'pr_npt']:
			filename = name+".mdp"
			outf = open("tmp.mdp","w")
			inf  = open(filename,"r")
			number_of_steps = self.to_number_of_steps(eval(name), dt)
			for line in inf:
				parameter = line.lstrip()[:6]
				if parameter[:2]=="dt":
					outf.write("dt   =  %5.3f ; %d fs \n"  % (self.timestep_in_fs*1.e-3, self.timestep_in_fs))
				elif parameter=="nsteps":
					outf.write("nsteps   =  %d ; %s \n"  % (number_of_steps, eval(name)))
				elif parameter in ["nstxou", "nstvou", "nstfou"]:
					# write coords and velocities only in the last step
					outf.write("%st =  %d \n" % (parameter,number_of_steps))
				elif parameter == "nstlog":
					outf.write("nstlog    =  %d \n" % (number_of_steps/10))
				elif parameter == "nstene":
					outf.write("nstenergy =  %d \n" % (number_of_steps/10))
				else:
					outf.write(line)
			inf.close()
			outf.close()
			os.rename("tmp.mdp", filename)

		for name in ['md']:
			filename = name+".mdp"
			outf = open ("tmp.mdp","w")
			inf =  open (filename,"r")
			number_of_steps = self.to_number_of_steps(eval(name),dt)
			for line in inf:
				parameter = line.lstrip()[:6]
				if parameter[:2]=="dt":
					outf.write("dt   =  %5.3f ; %d fs \n"  % (self.timestep_in_fs*1.e-3, self.timestep_in_fs))
				elif parameter=="nsteps":
					outf.write("nsteps   =  %d ; %s \n"  % (number_of_steps, eval(name)))
				elif parameter in ["nstxou", "nstvou", "nstfou"]:
					outf.write("%st   =  %d \n" % (parameter,number_of_steps/100))
				elif parameter == "nstlog":
					outf.write("nstlog    =  %d \n" % (number_of_steps/1000))
				elif parameter == "nstene":
					outf.write("nstenergy =  %d \n" % (number_of_steps/1000))
				else:
					outf.write(line)
			inf.close()
			outf.close()
			os.rename("tmp.mdp", filename)

		return


	#################################
	# annealing_type: cycle through annealing thing with 'periodic' (alt: 'single')
	# annealing-npoints: A list with the number of annealing reference/control points used for each temperature group.
	#     Use 0 for groups that are not annealed. The number of entries should equal the number of temperature groups.
	# annealing-time: List of times at the annealing reference/control points for each group (ps
	# annealing-temp: List of temperatures at the annealing reference/control points for each group.

	def set_annealing_schedule(self, run_params,
								annealing_type='periodic',
								annealing_npoints='5',
								annealing_time= "0     5    7.5     10     15   20",
								annealing_temp= "300   320   340    360   340   320"):

		os.chdir("/".join([run_params.run_options.workdir, run_params.rundirs.in_dir]))
		filename = "md.mdp"
		outf = open ("tmp.mdp","w")
		inf =  open (filename,"r")
		for line in inf:
			if 'anneal' in line: continue
			outf.write(line)
		inf.close()
		outf.write("\n")
		outf.write("annealing  = %s\n" % annealing_type)
		outf.write("annealing_npoints = %s\n" % annealing_npoints)
		outf.write("annealing_time = %s\n" % annealing_time)
		outf.write("annealing_temp = %s\n" % annealing_temp)
		outf.close()
		os.rename("tmp.mdp", filename)
		return
