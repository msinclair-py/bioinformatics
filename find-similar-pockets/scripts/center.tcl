#this script selects all the pdb's in the PBD
#directory and centers them all by center of mass

set arg [split $argv " "]
set mode [lindex $argv 0]
set name [lindex $argv 1]

if {$mode == 1} {
cd ../initial_structures
set pockets "$name.pocket.pdb"
} elseif {$mode == 2} {
cd ../PDB/
set pockets [glob "*.pdb"]
}

foreach pocket $pockets {
	mol new $pocket

	set pock [atomselect top all]
	set cent [measure center $pock]
	$pock moveby [vecscale -1.0 $cent]

	$pock writepdb cent_$pocket
}
exit
