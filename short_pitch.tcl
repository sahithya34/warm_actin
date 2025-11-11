set outfile [open "twist_full.align.dat" w]
for {set unit 27} {$unit >= 1} {incr unit -1} {
	set n [molinfo top get numframes ]
	for {set i 0} {$i < $n} {incr i} {
		set next [expr $unit - 1]
		set P1 [atomselect top "pfrag $next and name CA and resid 5 to 374" frame $i]
		set P2 [atomselect top "pfrag $unit and name CA and resid 5 to 374" frame $i]
		set M [ measure fit $P1 $P2 ] 
		set cos_theta [lindex [lindex $M 0] 0]
		set theta [expr acos($cos_theta)*180.0/$M_PI ]
		puts $theta
		puts  $outfile [format "%s %9.5f" $unit  $theta ]
}
}
close $outfile
