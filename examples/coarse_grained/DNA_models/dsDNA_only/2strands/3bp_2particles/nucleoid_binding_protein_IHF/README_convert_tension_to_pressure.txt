---- Warning: ----
First make sure that you are using the dimensionless ("lj") units in LAMMPS.
(Equivalently, make sure the "system.in.init" file contains "units lj". If you
are using some other unit system, additional conversions are needed. For my own
benefit, I wrote those down below.)  For more details about units in LAMMPS, see
http://lammps.sandia.gov/doc/units.html
------------------

To simulate a chain at constant tension, I use periodic boundary conditions and connect the two ends of the chain at opposite ends of the simulation box together.  The length of the simulation box controls the distance between the two ends of the chain.  The length of the simulation box is controlled by a barostat.

To achieve the desired tension in the chain, a constant negative pressure 
("pxx") is applied in the x direction.

pxx = -tension / (Ly*Lz)
   (in units of energy/distance) 

UNITS:
  In this simulation I use the LAMMPS command "units lj"
  http://lammps.sandia.gov/doc/units.html
  In these simulations, my distance units are nm, 
      (1nm approximately equals the width of the DNA, and a length of 3bp)
  my energy units are kCal/mole, hence pressure is in units of (kCal/mole)/nm.

  If you prefer to measure energy in units of kB*T, 
  then kB = 0.001987207 (kcal/mole)/degreeK
  --> kB*T = 0.5961621, assuming T=300Kelvin

Example

Ly=Lz=1000.0                      # (in nm)
tension = 0.8 * 0.001987207*300   # ie 0.8 kB*T per nm  (assuming T=300Kelvin)
pxx = -tension / (Ly*Lz)

==>     =-4.769e-7  (in units of (kcal/mole)/nm^3)

So, in LAMMPS (with "units lj"), use this command to maintain this tension
#                   pxx_start   pxx_stop   pdamp(time-units, usually 1000 steps)
fix fxnph all nph x -4.769e-7  -4.769e-7   1000.0
(The chain will not collapse at this tension.)

# --------- scratchwork ----------





Comments/scratchwork
1) The size of the box must not be allowed to vary in the Y and Z directions,
   (so be sure not so use "shrink-wrapped" boundary conditions in Y and Z)
    In LAMMPS use "boundary p f f" instead of "boundary p s s")
    http://lammps.sandia.gov/doc/boundary.html
    http://lammps.sandia.gov/doc/change_box.html
    http://lammps.sandia.gov/doc/fix_deform.html
    http://lammps.sandia.gov/doc/fix_nh.html
   (Note: I usually make the y and z directions non-periodic, although sometimes
        I use a cylindrical barrier to confine the polymer in those directions.)

-----------------------  ALTERNATE PRESSURE UNITS ----------------------

2) This discussion can be ignored because we are using the "units lj" command
   LAMMPS in the "settings.in.init" file (<==> moltemplate_files/forcefield.lt)

   However if we had not used the "units lj" command, then LAMMPS would revert
to the the default units ("units real").  In this case pressure is specified in atmospheres, where 1 atmosphere ~= 1 Bar
   For completeness, here is the conversion:

http://lammps.sandia.gov/doc/units.html
   Conversion:
1 Bar = 10000 Pascals = 10000 Newton/m^2
 1 Newton = 1 Joule / m
 1 Joule = (1/4.184) * thermochemical_calorie    (exactly 4.184)
 1 thermochemical_calorie = 0.001 kCal
 1 kCal = Na * (kcal/mole)     ( Na=6.022140857e23 )
 1 Bar = 10000.0  _____________Pascal______________
		  _____1_Joule____  / ___meter^3___
                   1J_in_kcal/mole  / (_m_in_nm_)^3
       = 10000.0 * (Na*0.001/4.184) / ((1e9)**3)
       = 0.00051815743  (kcal/mole) / nm

..so, if you are using "units real" in LAMMPS, then calculate pressure this way:

pxx  = -(tension / (Ly*Lz))   /   0.00051815743



