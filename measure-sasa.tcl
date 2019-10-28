###this script measure the solvent-accessible surface area
###of the residues of the pocket. the input is the name
###of the pocket output pbd from fpocket
set input [split $argv " "]
set inp [lindex $input 0]
set pock [lindex $input 1]

mol new "../PDB/pdb$inp.ent"

set input [open "../rosetta/$inp.pose" r]
set output [open "../PDB/$inp.sasa.txt" w]

set protein [atomselect 0 all]
set pocket [atomselect top all]

set res $input
set res [split $res " "]
puts $res

#set sasa [measure sasa 1.4 $protein -restrict $pocket]
#puts $output "$sasa"
#close $output
