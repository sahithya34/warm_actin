set outfile [open "dihedral_1micro_with_Dloop.dat" w]
for {set unit 0} {$unit < 27} {incr unit} {
	set n [molinfo top get numframes ]
	for {set i 0} {$i < $n} {incr i} {
		set sd2 [atomselect top "pfrag $unit and resid 33 to 69 and name CA" frame $i]
		set sd1 [atomselect top "pfrag $unit and resid 1 to 32 70 to 144 338 to 375 and name CA" frame $i]
		set sd3 [atomselect top "pfrag $unit and resid 145 to 180 270 to 337 and name CA" frame $i]
		set sd4 [atomselect top "pfrag $unit and resid 181 to 269 and name CA" frame $i]
		set com1 [measure center $sd2 weight mass]
		set com2 [measure center $sd1 weight mass]
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
