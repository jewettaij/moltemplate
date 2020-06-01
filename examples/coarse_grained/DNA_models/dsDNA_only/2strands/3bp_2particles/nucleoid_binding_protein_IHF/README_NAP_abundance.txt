Scratchwork used for figuring out how many DNA-binding proteins to add
to the simulation.  (Upper bounds.)

# How many base pairs of DNA are in the cell,
# and how many are in this simulation?
# Later, we will scale the number of proteins in the simulation by this ratio.
#
Nbp_per_monomer = 3.0
Npolylength = 59105
Nbp_sim = Nbp_per_monomer*Npolylength # = 177315 (bp)
Nbp_per_genome = 4639221.0
N_genome = 2.0
Nbp_cell = N_genome * Nbp_per_genome

# -----  IHF -----

N_IFS_cell = 56000.0  # Azam++Ishihama_JBacteriol1999, Figure 2J (peak at ~4hrs)

N_IFS_sim = N_IFS_cell*(Nbp_sim/Nbp_cell)
   # = 1070.184  # = number of IFS particles to put in the sim

# -----  H-NS -----

N_HNS_cell = 23833.0  # Azam++Ishihama_JBacteriol1999, Figure 2G (peak at ~3hrs)

N_HNS_sim = N_HNS_cell*(Nbp_sim/Nbp_cell)
   # = 455.4589  # = number of H-NS particles to put in the sim
   # = 910.9177  #   (if there is only 1 copy of the genome per cell, N_genome=1)
