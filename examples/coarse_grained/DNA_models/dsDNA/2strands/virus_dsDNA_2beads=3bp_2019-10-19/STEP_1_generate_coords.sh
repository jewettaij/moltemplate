# Check the prerequisites:

if [ ! -f ndmansfield ] && ! which git; then
    # If "ndmansfield" was not installed then we will need git to install it.
    # ("ndmansfield" is a random curve generator useful for polymer simulations)
    echo 'Error: You must install git before running this script.' >&2
    exit 1
fi


# Step 0)
# Copy the coordinates of the capsid to moltemplate_files/coords_wall.raw

cp capsid_locations.txt moltemplate_files/coords_wall.raw

# Step 1)
# I will approximate the mature HIV capsid as a sphere of radius R
# Generate a random self-avoiding polymer within this spherical region.

R=24.0   # *the interior of the mature HIV capsid is at least 2*R in diameter)
NMONOMERS=3058
B=0.996  # ~= distance between consecutive monomers (3bp) along the chain

# Let "delta" be the physical spacing of the points in the lattice.
# The self-avoiding polymer visits every lattice site on an Nx * Ny * Nz lttice.
# Hence the length of the polymer, L (which is NMONOMERS*B), satisfies
# this constraint:
# 1) delta * Nx * Ny * Nz = L = NMONOMERS*B
#     Suppose it is a cubical lattice (Nx=Ny=Nz)
#     Suppose we want it to fit within the sphere of radius R.  In that case:
# 2) delta*sqrt((Nx/2)**2+(Ny/2)**2+(Nz/2)**2) = delta * (Nx/2) * sqrt(3) = R
#    Solving for Nx and delta in terms of R and L yeilds:
#    --> delta = sqrt(R**3 / L) / ((sqrt(3)/2)**(3/2))
#    --> Nx = ceil( (L / delta)**(1/3) )
#           = 11, in this case
# NOTE: This is an extremely small, tight lattice.  It will take a great deal
#       of careful energy minimization before the polymer is relaxed enough
#       to be numerically stable (much less have a realistic conformation).

Nx=11
Ny=11
Nz=11

# Note: The polymer is only occupying a cube-shaped subvolume of the sphere
#       (which itself is only occupying a subvolume of the actual container).
#       After the polymer is allowed to expand to fill up the full container
#       it should have a spacing equal to:
#           a = sqrt(volume / (sin(pi/3)*L))
#       where L is the physical length of the polymer in nm, L=NMONOMERS*B
#       and volume is the volume of the container where it will be confined.
#       For a sphere, volume = (4/3)*pi*R**3.)


# If the program we need has not been compiled yet, do it now:
# If you don't have the "ndmansfield_src" source code, you can download it using
# git clone https://github.com/jewettaij/ndmansfield ndmansfield_src

if [ ! -f ndmansfield ]; then
  git clone https://github.com/jewettaij/ndmansfield ndmansfield_src
  cd ndmansfield_src/src/
  source setup_gcc.sh  # or "source setup_gcc_mac.sh" if you use a mac    
  make clean
  make
  mv -f ndmansfield ../../
  make clean
  cd ../../
fi

# Now run "ndmansfield" to generate a random compact self-avoiding curve        

if [ -z "$RANDOM_SEED" ]; then
    RANDOM_SEED=0  # In general
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


# We will rescale these coordinates and use the rescaled coordinates to         
# interpolate the coordinates of the physical polymer of physical length
# NMONOMERS*B through the lattice.
# The scaling factor should be (NMONOMERS*B)/(Nx*Ny*Nz)
SCALE=`awk -v Nx="$Nx" -v Ny="$Ny" -v Nz="$Nz" -v NMONOMERS="$NMONOMERS" -v B="$B" 'BEGIN{print NMONOMERS*B/(Nx*Ny*Nz-1)}'`

recenter_coords.py 0 0 0 < coords_lattice.raw \
          | interpolate_curve.py $NMONOMERS $SCALE \
          > moltemplate_files/init_crds_polymer_backbone.raw

# discard temporary files:
rm -f ndmansfield_traj.raw  coords_lattice.raw

