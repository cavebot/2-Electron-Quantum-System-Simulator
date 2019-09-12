#!/usr/bin/perl


# l1e_min           can be added as well
# l2e_min           can be added as well
# n2e ==            can be added as well
# gauge             can be added as well? 
# help


if ($#ARGV != 10 ) {
   print "usage: bsci2eb.pl exe system znuc l1e n1e l2e n2e rb_min rb_max dr_b\n";
   print "\n";
   print "0              exe =  directory location of the executables   \n";
   print "1           system =  hn,he0,li1,mg,ca,..                     \n";
   print "2                z =  1,2..                                   \n";
   print "3               v0 =  potential depth  (0=atoms)              \n";
   print "4              l1e =  max 1e angular momenta                  \n";
   print "5              n1e =  nof B-splines                           \n";
   print "6              l2e =  max 2e angular momenta                  \n";
   print "7              n2e =  max nof states for each angular momenta \n";
   print "8               r1 =  initial box radius in a.u.              \n";
   print "9               r2 =  final box radius in a.u.                \n";
   print "10               dr =  box radius step  in a.u.                \n";
   print "\n";
   exit;
}


$exe=$ARGV[0] ;
$system=$ARGV[1] ;
$z=$ARGV[2] ;
$z=$ARGV[3] ;
$l1e=$ARGV[4] ;
$n1e=$ARGV[5] ;
$l2e=$ARGV[6] ;
$n2e=$ARGV[7] ;
$r1=$ARGV[8] ;
$r2=$ARGV[9] ;
$dr=$ARGV[10] ;



print "$exe/run-bsci2e.pl  $exe $system $z $v0 $l1e $n1e $l2e $n2e $r1 $r2 $dr\n" ;   
#exit;



#
# assumes inp/h1e.inp file where rmax = $r1 
#

for ($r= $r1 ; $r <= $r2 ; $r += $dr ){

    print " r = $r,   store in  dir = $system/$r \n" ;   

    system("mkdir -p log") ;
    system("mkdir -p dat") ;
    system("mkdir -p out") ;
    system("mkdir -p inp") ;
    system("$exe/Rbspci2e -zn $z -rb $r -nb $n1e -v $v0 ");  # write input file
    system("$exe/run-h1e.scr 0 $l1e $exe"       );  # h1e  (structure + transitions)
 
   for ($l= 0 ; $l <= $l2e ; $l += 1 ){             # h2e  (structure)
	system("$exe/run-h2e.scr  $l $exe");
    }
    for ($li= 0 ; $li < $l2e ; $li += 1 ){      # h2e (transitions)
	$lf = $li + 1;
	system("$exe/run-d2e.scr $li $lf $n2e v $system $exe");
	system("$exe/run-d2e.scr $li $lf $n2e  l $system $exe");
    }
    system("mkdir -p  $system") ;
    system("mkdir -p  $system/r$r") ;
    system("mv    dat $system/r$r");
    system("mv    out $system/r$r");
    system("mv    log $system/r$r");
    system("cp -r inp $system/r$r");

#    $rr = $r1 + $dr;
#    perl -i.bak -pe "s/$r1/$rr/g" h1e.inp or die "$!";

}
print "end of run script ";
system("date");

