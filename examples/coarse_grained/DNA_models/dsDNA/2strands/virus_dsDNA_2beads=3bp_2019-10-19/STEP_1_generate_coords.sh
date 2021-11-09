# This is a simulation of a polymer inside a container (a virus).
# I will approximate the shape of the container as a sphere of radius R.
# We need to generate a random curve of the correct length with a given
# number of monomers which can fit within the container.

# Step 1a)

# Since we will be approximating the shape of the curve as a lattice polymer,
# we first need to figure out how big to make the lattice.

R=24.0   # Assume the container can be approximated by a sphere of this radius.
NMONOMERS=3058  # The number of monomers in the polymer.  (In this example
                # the polymer is a coarse-grained model of DNA.  Each "monomer"
                # represents 3 base-pairs of DNA and has length 0.996nm.)
B=0.996  # The physical distance between consecutive "monomers" along the chain

# L = the physical length of the polymer = B * NMONOMERS.  Calculate using awk:
L=`awk -v n="$NMONOMERS" -v b="$B" 'BEGIN{print n*b}'`

# Let us imagine an Nx*Ny*Nz lattice fitting inside the container.  Let "delta"
# refer to the physical separation between the points in this lattice. The
# self-avoiding polymer must visit every site in this Nx*Ny*Nz lattice.  The
# length of the polymer, L (which is NMONOMERS*B), satisfies this constraint:
# 1) delta * Nx * Ny * Nz = L = NMONOMERS*B
# Suppose it is a cubical lattice (Nx=Ny=Nz)
# 2)  Suppose we want it to fit within the sphere of radius R.  In that case:
#    --> R=delta*sqrt((Nx/2)**2+(Ny/2)**2+(Nz/2)**2) = delta * (Nx/2) * sqrt(3)
#    Solving for Nx and delta in terms of R and L yeilds:
#    --> delta = sqrt(R**3/L) * 2**(3/2) / (3)**(3/4)
#    --> Nx = ceil( (L / delta)**(1/3) )
#
# In this case, this yields:
# Nx=11
# Ny=11
# Nz=11
# More generally, we can use awk to calculate the formula above:

Nx=`awk -v R="$R" -v L="$L" 'BEGIN{delta=sqrt(R^3.0/L)*2^(3.0/2)/(3)^(3.0/4); Nx=(L/delta)^(1.0/3); if (Nx != int(Nx)) Nx=int(Nx)+1; print Nx}'`
Ny=$Nx
Nz=$Ny


# Step 1b)
# Generate a random self-avoiding polymer in a lattice of size Nx * Ny * Nz.
# If the "ndmansfield" program we need has not been downloaded and compiled yet,
# do it now.

if [ ! -f ndmansfield ]; then
  # Check the prerequisites:
  if [ ! -f ndmansfield ] && ! which git; then
    # If "ndmansfield" was not installed then we will need git to install it.
    # ("ndmansfield" is a random curve generator useful for polymer simulations)
    echo 'Error: You must install git before running this script.' >&2
    exit 1
  fi
  # Download the source code for "ndmansfield"
  git clone https://github.com/jewettaij/ndmansfield ndmansfield_src
  # Compile it
  cd ndmansfield_src/src/
  source setup_gcc.sh
  make clean
  make
  # Optional: Move the "ndmansfield" binary to the place where it will be used.
  mv -f ndmansfield ../../
  make clean
  cd ../../
fi


# Step 1c)
# Now run "ndmansfield" to generate a random compact self-avoiding curve        

if [ -z "$RANDOM_SEED" ]; then
    RANDOM_SEED=0  # If RANDOM_SEED was not defined, pick a default value
fi

./ndmansfield -box $Nx $Ny $Nz -cyclic no -seed $RANDOM_SEED  \
              -tsave 1000000 -tstop 3000000    \
              > ndmansfield_traj.raw

# We just need to run the lattice simulation long enough to get a reasonably
# random polymer conformation.   (You can check if we ran it long enough
# by watching the messages from ndmansfield sent to stderr.  When it's random,
# the number of bonds in each direction: x,y,z should be approximately equal,
# especially for large lattices with Nx > 20.

# The trajectory file created by ndmansfield is a 3-column ASCII file c         
# containing the coordinates of the polymer at different moments in "time".     
#     Note: The number of lattice sites is Ntot = Nx*Ny*Nz                      
# Each frame of this file has Ntot lines containing coordinates and 1 blank line
# (separator). We want to extract the last frame of this trajectory using:      

Ntot=$(( Nx*Ny*Nz ))
tail -n ${Ntot} ndmansfield_traj.raw > coords_lattice.raw


# Step 1d)
# Interpolate the coordinates to obtain a smoother curve.
# Also rescale these coordinates (multiply them by a constant).  What is that
# constant?  The physical length of the polymer should be NMONOMERS*B
# Hence, the scaling factor should be (NMONOMERS*B)/(Nx*Ny*Nz)

SCALE=`awk -v Nx="$Nx" -v Ny="$Ny" -v Nz="$Nz" -v L="$L" -v B="$B" 'BEGIN{print L/(Nx*Ny*Nz-1)}'`

# Optional: Center the coordinates around the origin (0,0,0):
recenter_coords.py 0 0 0 < coords_lattice.raw > coords_lattice_centered.raw

# Then interpolate and rescale them:
interpolate_curve.py $NMONOMERS $SCALE \
          < coords_lattice_centered.raw \
          > moltemplate_files/init_crds_polymer_backbone.raw

# The final curve (moltemplate_files/init_crds_polymer_backbone.raw)
# should contain one point at every location where each monomer should go.
# Documentation for "interpolate_curve.py" can be found here:
# https://github.com/jewettaij/moltemplate/blob/master/doc/doc_interpolate_curve.md
# NOTE: This is an extremely small, tight lattice.  This is not a realistic
#       conformation of DNA, but we can relax the shape of the DNA later.

# discard temporary files:
rm -f ndmansfield_traj.raw  coords_lattice.raw coords_lattice_centered.raw


# Step 1e)
# This is a simulation of DNA inside a container.
# The container is made from a collection of immobile particles through which
# the polymer cannot pass (in this case, the capsid wall of a virus).
# It is convenient to copy the coordinates of those immobile particles
# into the same directory with the other moltemplate files.  Later we will
# create another file (in moltemplate format) describing those particles.

cp capsid_locations.txt moltemplate_files/coords_wall.raw

# Note: If you are simulating the behavior of a bulk polymer melt
#       (without a container), this step is not relevant to you.

# Note: Estimate the final polymer density and spacing:
#       After the polymer is allowed to expand to fill up the full container
#       the distance between nearest-neighbor polymer strands should
#           a = sqrt(volume / (sin(pi/3)*L))   <-- assumes hexagonal packing
#       where L is the physical length of the polymer in nm, L=NMONOMERS*B
#       and volume is the volume of the container where it will be confined
#       (or the volume of the whole simulation box, if there is no container).
#       In this example, initially, the polymer is only occupying a
#       cube-shaped subvolume of the sphere (which itself is only occupying
#       a subvolume of the actual container).  Consequently we expect the
#       initial density to be somewhat higher.
