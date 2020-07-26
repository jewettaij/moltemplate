This file contains some suggestions for how to improve the appearance of long
polymers using VMD.  (Other useful programs like OVITO are not covered here.)
These instructions were tested with VMD version 1.9.4a38.
(Later versions of VMD may have different menu options, but hopefully are
 similar enough that these instructions remain helpful.)

-------------------------------
Follow the instructions in the "README_visualize_VMD.txt" file in the local
directory to display the LAMMPS data file and load the DUMP file(if applicable).
(Sometimes the file is named "README_visualize.txt" instead.)

1) In the "Display"->"Render Mode"->"GSL" (this increases rendering speed
later when you change the drawing method to a non-wirefram representation).

2) OPTIONAL: Select the "Graphics"->"Representations" menu option.
Under the "Coloring Method" pull-down menu select "Index".

This will assign color to the atoms according to their "Index" 
(which equals their LAMMPS atom-ID - 1).  It uses color map which by default
fades from blue to white to red. VMD comes with several colormaps built in.
You can select from them using the "Graphics"->"Colors" global menu option
from the main window, clicking on the "Color Scale" tab, and selecting one
of the options from the "Method" pull-down menu.

I personally do not like any of these colormaps, so I created several
additional ones you can load manually.  To do that, select the
"Extensions"->"Tk Console" menu option from the main window.  In the new window
that is created, select the "File"->"Load File" menu option and select one of
the .tcl files in this directory.  For example:

  vmd_colorscale_jet_linear.tcl
  vmd_colorscale_jet_circular.tcl

Then (if you haven't already) select "Index" from the "Coloring Method"
pull-down menu.

3) Optional: Select the "Display"->"Display Settings" main menu option and
increase the magnitude of the "Screen Dist" parameter from the default value
(-2.0) to -10.0 or -30.0 or larger.  (Scroll the mouse wheel to compensate
for the increased magnification that this causes.)

4) Optional: In the main graphics window, hit the "T" key, and hold down the
right mouse button while moving the mouse back and forth.  This will bring b
he object closer to the camera reducing the effect of "Cue" (fog).  This
should increase its brightness (contrast).  Stop before portions of the object
get too close to the viewer and begin to clip.  (Scroll the mousewheel to
compensate for the increased magnification that this causes.)  You can also
increase or decrease the thickness of the fog by changing the "Cue Density".

5) Select the "Graphics"->"Representations" menu option.
Under the "Drawing Method" pull-down menu select "Licorice".
(In my experience, "CPK", "VDW", and "Dotted" and "Quicksurf" are also very
 useful.)
(This will slow down the frame rate during rotation considerably compared
 to a wireframe representation.  You might want to wait until the objects
 are positioned  and rotated correctly beforehand.)

6) Optional: You can exclude certain atoms from being visualized by clicking
on the "Selected Atoms" text field (in the "Graphical Representations" window),
and replacing "all" with something like "type <= 2" (which only displays
atoms of type 1 or 2).  You can represent different atom types in different
ways by clicking on the "Create Rep" button, customizing the atom selection,
and changing the drawing and coloring method for that atom selection.

7) Select the "Display"->"Display Settings" menu option from the main window
and select "Shadows: On", "Amb. Occl.: On".  Then select the "File"->"Render"
menu option and select "Tachyon (internal memory rendering)" from the pull-down
menu, and click on the "Start Rendering" button.  (Depending on your window size
and number of processors, a decent quality rendering of a large DNA structure,
such as a bacterial chromsome, can take an hour or longer.)

