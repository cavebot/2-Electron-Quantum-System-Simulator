#!/usr/bin/perl


#
#   run as:  ./scr/run_Nph.pl $exe $system $nof_photons $gauge $r_initial $r_final $dr 
#   help
#


if ($#ARGV != 6 ) {
   print "usage: run_Nph.pl exe system nof_photons gauge r_init r_fin dr \n";
   print "\n";
   print "              exe = ' directory location of the executables' \n";
   print "           system =  hn,he0,li1,mg,ca,.. \n";
   print "      nof_photons =  1,2... \n";
   print "            gauge =  l, v   \n";
   print "           r_init =  initial box radius in a.u.   \n";
   print "           r_fin  =    final box radius in a.u.   \n";
   print "              dr  =    box radius step  in a.u.   \n";
   print "\n";
   exit;
}



$exe=$ARGV[0] ;
$system=$ARGV[1] ;
$nph=$ARGV[2] ;
$gauge=$ARGV[3] ;
$r1= $ARGV[4] ;
$r2=$ARGV[5] ;
$dr=$ARGV[6] ;


system("mkdir -p $nph-ph") ;
system("mkdir -p $nph-ph/$gauge") ;
system("mkdir -p $nph-ph/$gauge/dat") ;
system("mkdir -p $nph-ph/$gauge/inp") ;
system("mkdir -p $nph-ph/$gauge/out") ;


#
# create and write inp/csNph_merge_r.inp
#


$csNph_merge_inp = "inp/csNph_merge_r.inp";
unlink $csNph_merge_inp ;
open merge_file, "+>>$csNph_merge_inp" or die "$!";

$nr = 0;
for ($r= $r1 ; $r <= $r2 ; $r += $dr ){
   $nr = $nr + 1;
}
print merge_file "$nr \n";

for ($r= $r1 ; $r <= $r2 ; $r += $dr ){
    print merge_file "$r\t"; 
#print merge_file "    ";
}
print merge_file "\n";
for ($r= $r1 ; $r <= $r2 ; $r += $dr ){
    print merge_file "1\t" ; 
}
close merge_file;


system( "cp $csNph_merge_inp   $nph-ph/$gauge/inp");

#
#  run multiphoton cross sections N=1-4
#

for ($r= $r1 ; $r <= $r2 ; $r += $dr ){

    $dir="r$r";                        # director for box radius r = .. 
    print " r = $r  dir = $dir \n" ;   # r = ...

    chdir($dir) or die "$!";
    system( "cp ../inp/csNph-$nph.inp inp/csNph.inp");


    system( "$exe/RcsNph_i $nph $gauge $system");
    system( "$exe/RcsNph_l $nph $gauge");
	if($nph == 1) {
	    system( "cp dat/csNph-1.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l1.r$r.dat");
	}
	elsif($nph == 2){
	    system( "cp dat/csNph-0.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l0.r$r.dat");
	    system( "cp dat/csNph-2.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l2.r$r.dat");
	}
	elsif($nph == 3){
	    system( "cp dat/csNph-1.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l1.r$r.dat");
	    system( "cp dat/csNph-3.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l3.r$r.dat");
	}
	elsif($nph == 4){
	    system( "cp dat/csNph-0.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l0.r$r.dat");
	    system( "cp dat/csNph-2.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l2.r$r.dat");
	    system( "cp dat/csNph-4.$gauge.dat   ../$nph-ph/$gauge/dat/csNph-l4.r$r.dat");
	}
	else{
	print " no implementation for nof photons above 4.";
	die "$!";
	}	

#
#
#
    system( "rm dat/dmx.dat");
    system( "rm dat/en.dat");
    chdir('../');

    $nr = $nr + 1
}

#
#  merge the cross sections for various radii
#


$dir="$nph-ph/$gauge"; 
print " dir = $dir";
chdir($dir) or die "$!";
system( "$exe/RcsNph_merge_r $nph");
print "merge cross sections done";
#
#
#

print " nr = ", $nr, " \n";
print " end of run script ";
system("date");

