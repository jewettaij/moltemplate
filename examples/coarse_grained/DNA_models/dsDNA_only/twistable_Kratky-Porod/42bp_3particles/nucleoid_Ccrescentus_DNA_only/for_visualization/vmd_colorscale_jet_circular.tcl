
proc lerpcolor { col1 col2 alpha } {
  set dc [vecsub $col2 $col1]
  set nc [vecadd $col1 [vecscale $dc $alpha]]
  return $nc  
}

proc coltogs { col } {
  foreach {r g b} $col {}
  set gray [expr ($r + $g + $b) / 3.0]
  return [list $gray $gray $gray]
}

proc jet_tricolor_scale {} {
  display update off
  set mincolorid [expr [colorinfo num] - 1]
  set maxcolorid [expr [colorinfo max] - 1]
  set colrange [expr $maxcolorid - $mincolorid]
  set colhalf [expr $colrange / 2]
  for {set i $mincolorid} {$i <= $maxcolorid} {incr i} {
    set colpcnt [expr ($i - $mincolorid) / double($colrange)]

    # The following color definitions for "jet" sort-of came from:
    #http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale
    # but it was missing "green", so I inserted a green somewhere in the middle.

    # magenta
    set color0 { 1.0 0.0 1.0 }

    set color1 { 0.08 0.0 0.77 }

    # blue 
    set color2 { 0.0 0.0 1.0 }

    # cyan 
    set color3 { 0.0 1.0 1.0 }

    # turquoisegreen 
    set color4 { 0.0 1.0 0.2 }

    # green 
    set color5 { 0.0 1.0 0.0 }

    # chartreuse 
    set color6 { 0.1 1.0 0.0 }

    # yellow 
    set color7 { 1.0 1.0 0.0 }

    # orange 
    set color8 { 1.0 0.25 0.0 }

    # red 
    set color9 { 1.0 0.0 0.0 }

    # darkred 
    set color10 { 0.93 0.0 0.0 }

    if { $colpcnt < 0.121 } {
      set nc [lerpcolor $color0 $color1 [expr $colpcnt/(0.121-0.0)]]
    } elseif { $colpcnt < 0.244 } {
      set nc [lerpcolor $color1 $color2 [expr ($colpcnt-0.121)/(0.244-0.121)]]
    } elseif { $colpcnt < 0.371 } {
      set nc [lerpcolor $color2 $color3 [expr ($colpcnt-0.244)/(0.371-0.244)]]
    } elseif { $colpcnt < 0.424 } {
      set nc [lerpcolor $color3 $color4 [expr ($colpcnt-0.371)/(0.424-0.371)]]
    } elseif { $colpcnt < 0.471 } {
      set nc [lerpcolor $color4 $color5 [expr ($colpcnt-0.424)/(0.471-0.424)]]
    } elseif { $colpcnt < 0.489 } {
      set nc [lerpcolor $color5 $color6 [expr ($colpcnt-0.471)/(0.489-0.471)]]
    } elseif { $colpcnt < 0.630 } {
      set nc [lerpcolor $color6 $color7 [expr ($colpcnt-0.489)/(0.630-0.489)]]
    } elseif { $colpcnt < 0.714 } {
      set nc [lerpcolor $color7 $color8 [expr ($colpcnt-0.630)/(0.714-0.630)]]
    } elseif { $colpcnt < 0.815 } {
      set nc [lerpcolor $color8 $color9 [expr ($colpcnt-0.714)/(0.815-0.714)]]
    } elseif { $colpcnt < 0.888 } {
      set nc [lerpcolor $color9 $color10 [expr ($colpcnt-0.815)/(0.888-0.815)]]
    } else {
      set nc [lerpcolor $color10 $color0 [expr ($colpcnt-0.888)/(1.0-0.888)]]
    }

    #    set nc [coltogs $nc]
    foreach {r g b} $nc {}
    puts "index: $i $r $g $b   -- $colpcnt"
    display update ui
    color change rgb $i $r $g $b 
  }
  display update on
}

jet_tricolor_scale

