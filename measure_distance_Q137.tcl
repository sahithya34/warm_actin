set outfile [open "distance_Q137-Pi.dat" w]
set n [molinfo top get numframes ]
for {set unit 0} {$unit < 13} {incr unit} { 
	for {set i 0} {$i < $n} {incr i} {
	set Q137 [atomselect top "pfrag $unit and resid 137 and name OE1" frame $i]
	set H161 [atomselect top "resname ATP and name PG and within 10 of (pfrag $unit and resid 137)" frame $i]
	set com2 [measure center $Q137 weight mass]
	set com3 [measure center $H161 weight mass]
	set v1 [vecsub $com2 $com3 ]
	set mod_v1 [veclength $v1]
	puts  $outfile [format "%-5d %7.3f" $unit  $mod_v1 ]
	}
}
close $outfile
