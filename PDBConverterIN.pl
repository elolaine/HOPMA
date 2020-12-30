#!/usr/bin/perl

# Elodie Laine april 16th 2010
# convert a PDB file to a csv file
# usage ./PDBConverterIN.pl PDBFile_noext 

$fileIN = $ARGV[0].".pdb";
$fileOUT = $ARGV[0].".csv";

open(OUT,">".$fileOUT)  || die "cannot open OUT file";
open(IN, $fileIN) || die "cannot open IN file";
  while (<IN>) {
    @ff = split(' ',$_);
    if ((index($ff[0],"ATOM")>(-1))||(index($ff[0],"HETATM")>(-1))) {
	chomp($_);
	$record = substr $_, 0, 6;
	$atnum = substr $_, 6, 5;
	$name = substr $_, 12, 4;
	$resname = substr $_, 17, 3;
	$chain = substr $_, 21, 1;
	$resid = substr $_, 22, 4;
	#$resid = substr $_, 22, 5;
	$alter = substr $_, 26, 1;
        $x = substr $_, 30, 8;
	$y = substr $_, 38, 8;
	$z = substr $_, 46, 8;
	$occ = substr $_, 54, 6;
	$bfact = substr $_, 60, 6;
	$id = substr $_, 72, 3;
	$elt = substr $_, 76, 2;
	$chg = substr $_, 78, 2;
	printf OUT $record.",".$atnum.",".$name.",".$resname.",".$chain.",".$resid.",".$alter.",".$x.",".$y.",".$z.",".$occ.",".$bfact.",".$id.",".$elt.",".$chg."\n";
    };
  };

close(IN);
close(OUT);







      






