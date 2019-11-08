#arguments are passed from the setup.py script
#these can be modified as user input variables in that script
set args $argv
set input [lindex $args 0]
set d1 [lindex $args 1]
set d2 [lindex $args 2]

mol new "../initial_structures/aligned.cent_$input.pocket.pdb"

set sel [atomselect 0 all]
set r [transaxis z $d1]
set t [transaxis x $d2]

for {set i 0} {$i <= 1} {incr i} {
	for {set j 0} {$j < 24} {incr j} {
		$sel move $r
		set deg [expr 15 * $j]
		puts $deg

		if {$i == 0} {
			set tilt "notilt"
		} else {
			set tilt "tilted"
		}

		$sel writepdb "../initial_structures/$input.$deg.$tilt.pdb"
	}

	$sel move $t
}

exit
