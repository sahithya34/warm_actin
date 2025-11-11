mol new ../../step4.4_equilibration_protein.pdb
#mol addfile ../npt_unb_protein.xtc waitfor all
mol addfile ../../md_1_protein.nojump.xtc waitfor all

set outfile [open "dihedral_Q137.dat" w]
for {set unit 0} {$unit < 13} {incr unit} {
	set n [molinfo top get numframes ]
	for {set i 0} {$i < $n} {incr i} {
		set sd2 [atomselect top "pfrag $unit and resid 137 and name CB" frame $i]
		set sd1 [atomselect top "pfrag $unit and resid 137 and name CA" frame $i]
		set sd3 [atomselect top "pfrag $unit and resid 137 and name CG" frame $i]
		set sd4 [atomselect top "pfrag $unit and resid 137 and name CD" frame $i]
		set com1 [measure center $sd1 weight mass]
		set com2 [measure center $sd2 weight mass]
		set com3 [measure center $sd3 weight mass]
		set com4 [measure center $sd4 weight mass]
		set v1 [vecsub $com2 $com1 ]
		set mod_v1 [veclength $v1]
		set v2 [vecsub $com3 $com2 ]
		set mod_v2 [veclength $v2]
		set v3 [vecsub $com4 $com3 ]
		set mod_v3 [veclength $v3]
		set n1 [vecnorm [veccross $v1 $v2 ]]
		set n2 [vecnorm [veccross $v2 $v3 ]]
		set dot_prod [vecdot $n1 $n2 ]
		set cos_theta [expr $dot_prod  ]
		set theta [expr acos($cos_theta)*180.0/$M_PI ]
		set n3  [vecnorm [ veccross $n1 $n2 ]]
		set sign [vecdot  $v2 $n3 ]
		if {$sign < 0} {
			set theta [expr {0-$theta}]
		} elseif {$sign == 0} {
			set theta 180.00
		} else {
			set theta $theta }
		puts  $outfile [format "%s %9.5f" $unit  $theta ]
	}
}	
close $outfile
