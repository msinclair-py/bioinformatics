###this script measure the solvent-accessible surface area
###of the residues of the pocket. the input is the name
###of the pocket output pbd from fpocket, the pocket number,
###the volume of the original pocket and the volume filter
set input [split $argv " "]
set inp [lindex $input 0]
set pock [lindex $input 1]
set volume [lindex $input 2]
set filter [lindex $input 3]

proc listFromFile {filename} {
	set f [open $filename r]
	set data [split [string trim [read $f]]]
	close $f
	return $data
}


proc average L {
    expr ([join $L +])/[llength $L].
}

if { $pock == 0 } {
	mol new "../initial_structures/$inp.pdb"
	set output [open "../initial_structures/$inp.sasa.txt" w]
	set res [listFromFile "../initial_structures/$inp.pos"]
} else {
	mol new "../PDB/pdb$inp.ent"
	set output [open "../PDB/$inp.$pock.sasa.txt" w]
	set res [listFromFile "../rosetta/$inp\_$pock.pos"]
}

set protein [atomselect top all]

set count 0
set sasa []
foreach r $res {
	set sel [atomselect top "residue $r"]
	set rsasa [measure sasa 1.4 $protein -restrict $sel]
	$sel set user $rsasa
	$sel delete
	lappend sasa $rsasa
	puts $output "residue $r, sasa: $rsasa"
	if { $rsasa > $filter } {
		incr count
	}
}

set a [average $sasa]
set avg [expr $a / $volume]
puts $output "average, volume normalized sasa: $avg"
puts $output "number of residue over filter ($filter): $count"
close $output

exit
