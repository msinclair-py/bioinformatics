###this script takes the input of pocket file name and then
###the normal vector of its longest axis and aligns it with z
set arg [split $argv " "]
set name [lindex $arg 0]
set xvec [lindex $arg 1]
set yvec [lindex $arg 2]
set zvec [lindex $arg 3]
set which [lindex $arg 4]

#if the name is the original protein then use the correct filepath
#otherwise use the filepath for where all the other structures are
#stored
if {$which == 1} {
	set path "../initial_structures/"
} else {
	set path "../PDB/"
}
mol new "$path$name"

#do the vector inversion and output the newly aligned structure
set vector [list $xvec $yvec $zvec ]
set sel [atomselect top all]
set M [transvecinv $vector]
$sel move $M
set M [transaxis y -90]
$sel move $M
$sel writepdb "$path$name.aligned"
exit
