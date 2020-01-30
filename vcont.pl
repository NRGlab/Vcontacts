#!/usr/bin/perl -w

$pid = $$;
$usage = "USAGE: vcont_f.pl -f <PDB> -r RESNUMC -c -s -m LIST\n";
$usage .= "\t  -f  PDB filename or code. Obligatory\n";
$usage .= "\t  -r  output of Residue RES number NUM chain  C.\n";
$usage .= "\t  -c  output only contacts.\n";
$usage .= "\t  -s  output only solvent accessible area (SAS), doesn't work with -c.\n";
$usage .= "\t  -p  output unique contacts in PDB format (works with -r).\n";
$usage .= "\t  -ca include CA atoms in output (works with -p).\n";
$usage .= "\t  -x exclude RESNUMC in output.\n";
$usage .= "\t  -m matrix of contact areas, diagonals are SAS and LIST restrict output to those in LIST. Required arg, use ALL otherwise\n";

$pdb = "";
$sas_only=0;
$query_rid = "";
$just_contact_list = 0;
$sas_contact_list = 0;
$pdb_format=0;
$add_calpha = 0;
$xclude_qrid = 0;
$contact_matrix=0;
$list = "ALL";

while(@ARGV){
  $arg = shift @ARGV;
  if($arg eq "-f"){$pdb = shift @ARGV;}
  if($arg eq "-r"){$query_rid = shift @ARGV;}
  if($arg eq "-s"){$sas_only = 1;}
  if($arg eq "-c"){$just_contact_list = 1;}
  if($arg eq "-C"){$just_contact_list = 1;$sas_contact_list = 1;}
  if($arg eq "-p"){$pdb_format = 1;}
  if($arg eq "-ca"){$add_calpha = 1;}
  if($arg eq "-x"){$xclude_qrid = 1;}
  if($arg eq "-m"){$contact_matrix = 1;$list=shift @ARGV;}
  if($arg eq "-h"){die $usage;}  
}

$error="ERROR: ";
if($pdb eq ""){
  $error .= "A PDB filename or code is needed\n";
  print $usage;
  die $error;
}

if($query_rid eq ""){$query_rid="ALL";}
if($query_rid =~ /ALL/){$xclude_qrid = 0;}

$delete_pdb_file = 0; 
if($pdb =~ /\.pdb/){
  $pdb_file = $pdb;
}else{
  $pdb_file = $pdb."_".$pid.".pdb";
  
  $com  = "ssh bcb.med.usherbrooke.ca cat /databases/pdb/data/structures/all/pdb/pdb";
  $com .= $pdb.".ent.gz | gunzip ";
  $mand = " >> ".$pdb_file;
  system($com.$mand);
  if(!-e $pdb_file){die "ERROR $pdb doesn't exist\n";}
  $delete_pdb_file = 1;
}

open(VCONT, "~/bin/Vcontacts -radii ~/Vcalcs-v1.2/radii.dat $pdb_file |");
while($line = <VCONT>){
  if($line =~ /^\#/){next;}
  if($line =~ /SAS/){
    ($anum = substr($line,0,6))  =~ s/^\s+//gi;
    ($rnum = substr($line,11,5)) =~ s/^\s+//gi;
    ($rnam = substr($line,18,3)) =~ s/\s/\-/gi;
    ($chn  = substr($line,22,1)) =~ s/\s/\-/gi;
    ($sas  = substr($line,33,18)) =~ s/\s+//gi;
    $rid = $rnam.$rnum.$chn;
    #print $line;
    #print "SAS=",$sas;<STDIN>;
    if(!exists($sas{$rid})){
      $sas{$rid}=$sas; # SAS of the residue
    }else{
      $sas{$rid}+=$sas; # SAS of the residue
    }
    #print $rid," ",$sas{$rid},"\n";
    $sas_line{$rid}{$anum}=$line;
    if($#resids == -1 || $rid ne $resids[$#resids]){
      push(@resids,$rid);
    }
  }else{
    $con_line{$rid}{$anum} .= $line; #
    if($line !~ /^$/){
      ($rnum = substr($line,32,5)) =~ s/^\s+//gi;
      ($rnam = substr($line,39,3)) =~ s/\s/\-/gi;
      ($chn  = substr($line,43,1)) =~ s/\s/\-/gi;
      ($area = substr($line,44,8)) =~ s/^\s+//gi;
      #print "area= ",$area,"\n";
      $c_rid=$rnam.$rnum.$chn;
      #if($c_rid eq $rid){print "$c_rid eq to $rid\n";<STDIN>;}
      if(!exists($contact{$rid}{$c_rid})){ #
	$contact{$rid}{$c_rid}=$area; # contact area between rid and c_rid
	push(@{$con_list{$rid}},$c_rid);
      }else{
	$contact{$rid}{$c_rid} += $area; #
      }
    }
  }
}
close(VCONT);

if($contact_matrix==1){
  if($list eq "ALL"){
    #output matrix for all residues
    foreach $rid (@resids){$matrix_list{$rid}=1;}
  }else{
    #print $list;<STDIN>;
    @tmp=split(/\,/,$list);
    foreach $rid (@tmp){
      $sas_str = sprintf "%4.2f",$sas{$rid};
      print "$rid $sas_str\n";      
      $matrix_list{$rid}=1;
    }
    
  }
  
  foreach $c_rid (@resids){
    if(exists($matrix_list{$c_rid})){print "\t",$c_rid;}
  }
  print "\n";
  foreach $rid (@resids){
    if(exists($matrix_list{$rid})){
      print $rid;
      foreach $c_rid (@resids){
	if(exists($matrix_list{$c_rid})){
	  if(!$contact{$rid}{$c_rid}){
	    print "\t0.0";
	  }else{
	    print "\t",$contact{$rid}{$c_rid};
	  }
	}
      }
      print "\n";
    }
  }
}else{
  if($just_contact_list == 0){
    if($query_rid eq "ALL"){
      foreach $rid (@resids){
	@anums=sort {$a <=>$b} keys %{$sas_line{$rid}};
	foreach $anum (@anums){
	  print $sas_line{$rid}{$anum};
	  if($sas_only==0){print $con_line{$rid}{$anum};}
	}
      }
    }else{
      @anums=sort {$a <=>$b} keys %{$sas_line{$query_rid}};
      %c_anums=();
      if($pdb_format == 1){%record=&pdb_lines_by_anum($pdb_file);}
      foreach $anum (@anums){
	if($pdb_format == 0){
	  print $sas_line{$query_rid}{$anum};
	  if($sas_only==0){print $con_line{$query_rid}{$anum};}
	}else{
	  @tmp=split(/\n/,$con_line{$query_rid}{$anum});
	  foreach $t (@tmp){
	    @tmp1=split(/\s+/,$t);
	    $c_anums{$tmp1[1]}=1;
	  }
	}
      }
      if($pdb_format == 1){
	if($add_calpha==1){%ca_rid=&ca_rid_anum(%record);}
	%pdbout_anums=();
	foreach $c (keys %c_anums){
	  $rnam = uc(substr($record{$c},17,3));
	  $rnum=substr($record{$c},22,4);
	  $rnum =~ s/^\s+//gi;
	  $rchn=substr($record{$c},21,1);
	  $rid=$rnam.$rnum.$rchn;
	  if($xclude_qrid == 1 && $rid eq $query_rid){next;}
	  %pdbout_anums=(%pdbout_anums,$c => "1");
	  if($add_calpha==1 && exists($ca_rid{$rid})){
	    %pdbout_anums=(%pdbout_anums,$ca_rid{$rid} => "1");
	  }
	}
	foreach $c (sort {$a<=>$b} keys %pdbout_anums){print $record{$c};}
      }
    }
  }else{
    if($query_rid eq "ALL"){
      foreach $rid (@resids){
	$sas_str="";
	if($sas_contact_list == 1){
	  $sas=0;
	  @anums=keys %{$sas_line{$rid}};
	  foreach $anum (@anums){
	    @tmp=split(/\s+/,$sas_line{$rid}{$anum});
	    $sas += $tmp[$#tmp];
	  }
	  $sas_str = sprintf "%4.2f",$sas;
	}
	
	printf "%-8s %5s",$rid,$sas_str;
	foreach $c_rid (@{$con_list{$rid}}){
	  if($xclude_qrid == 1 && $c_rid eq $query_rid){next;}
	  print " $c_rid";
	}
	print "\n";
      }
    }else{
      $sas_str="";
      if($sas_contact_list == 1){
	$sas=0;
	@anums=keys %{$sas_line{$query_rid}};
	foreach $anum (@anums){
	  @tmp=split(/\s+/,$sas_line{$query_rid}{$anum});
	  $sas += $tmp[$#tmp];
	}
	$sas_str = sprintf "%4.2f",$sas;
      }
      
      printf "%-8s %5s",$query_rid,$sas_str;
      foreach $c_rid (@{$con_list{$query_rid}}){
	if($xclude_qrid == 1 && $c_rid eq $query_rid){next;}
	print " $c_rid";
      }
      print "\n";
    }
  }
}
if($delete_pdb_file == 1){unlink($pdb_file);}
  
#1234567890123456789012345678901234567890123456789012345678901234567890123456#
#         1         2         3         4         5         6         7      #
#1234567890123456789012345678901234567890123456789012345678901234567890123456#
sub pdb_lines_by_anum(){
    use strict;
    my $pdb_file = shift;
    my ($line,$anum);
    my %records=();
    local (*PDB);

    open(PDB,$pdb_file) || die "cannot open pdb file $pdb_file for read\n";
    while($line=<PDB>){
	if($line =~ /^ATOM/ || $line =~ /^HETATM/){
	    #print $line;
	    $anum = substr($line,6,5);
	    $anum =~ s/^\s+//gi;
	    #print "<$anum>\n";<STDIN>;
	    $records{$anum}=$line;
	}
    }
    close(PDB);
    return %records;
}
#------------------------------------------------------------------------------
sub ca_rid_anum(){
    use strict;
    my %records=@_;
    my %map=();
    my ($r,$anam,$anum,$rnam,$rnum,$rchn,$rcod);
    
    foreach $r (keys %records){
	$anam=substr($records{$r},12,4);
	if($anam eq " CA "){
	    $rnam = uc(substr($records{$r},17,3));
	    $rcod=Three2One($rnam);
	    if($rcod ne "X"){
		$anum=substr($records{$r},6,5);
		$anum =~ s/^\s+//gi;
		$rnum=substr($records{$r},22,4);
		$rnum =~ s/^\s+//gi;
		$rchn=substr($records{$r},21,1);
		#print "$anam $anum $rnam $rnum $rchn\n";
		$map{$rnam.$rnum.$rchn}=$anum;
	    }
	}
	
    }
    return %map;
}
#----------------------------------------------------------------------------
sub Three2One(){
    use strict;
    my $three=lc(shift);
    my $one="x";
    
    if($three eq "gly"){$one="g";}
    if($three eq "ala"){$one="a";}
    if($three eq "val"){$one="v";}
    if($three eq "leu"){$one="l";}
    if($three eq "ile"){$one="i";}
    if($three eq "met"){$one="m";}
    if($three eq "pro"){$one="p";}
    if($three eq "phe"){$one="f";}
    if($three eq "trp"){$one="w";}
    if($three eq "ser"){$one="s";}
    if($three eq "thr"){$one="t";}
    if($three eq "asn"){$one="n";}
    if($three eq "gln"){$one="q";}
    if($three eq "tyr"){$one="y";}
    if($three eq "cys"){$one="c";}
    if($three eq "lys"){$one="k";}
    if($three eq "arg"){$one="r";}
    if($three eq "his"){$one="h";}
    if($three eq "asp"){$one="d";}
    if($three eq "glu"){$one="e";}

    return uc($one);
}
