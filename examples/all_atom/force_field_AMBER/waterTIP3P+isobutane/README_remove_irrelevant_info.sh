
# Note: By default, the system.data and system.in.settings files contain
# extra information for atoms defined in AMBER/GAFF2 ("gaff2.lt") which you
# are not using in this simulation.
# This is harmless, but if you to delete this information from your
# system.in.settings and system.in.data files, run this script:

cleanup_moltemplate.sh

# (Note: Removing unecessary atom types will make it easier to visualize the
#        simulation in VMD.)
