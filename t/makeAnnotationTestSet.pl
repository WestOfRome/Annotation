#!/usr/bin/perl 

unless ( @ARGV ) {
    print <<USAGE;
    
    $0 gff fasta --fragnum=200 --gapeln=100
	
USAGE
exit;
}

use Getopt::Long;

$fragnum=200;
$gaplen=100;

GetOptions(
    'fragnum=i' => \$fragnum,
    'gaplen=i' => \$gaplen,
    'restrict=s' => \$restrict
    );

####################################
# read in gff 
####################################

$, = "\t";
$\ = "\n";
$/ = "\n";

my $file = $ARGV[0];
open(GFF, shift @ARGV);
while ( <GFF> ) {
    chomp;
    next unless /\w/;
    next if /^\#/;    
    my @r=split/\t/;
    next unless $r[2] eq 'CDS';
    next unless $r[0] =~ /chr[XVI]/;
    push @{$hash{ $r[0] }{GFF}}, \@r;        
}
close GFF;

# sort by start.

foreach my $chr ( keys %hash ) {
    my @sort = sort { $a->[3] <=> $b->[3] } @{$hash{$chr}{GFF}};
    $hash{$chr}{GFF} = \@sort;
}

####################################
# read in fasta
####################################

$/ = ">";

while ( <> ) {
    chomp;
    next unless /\w/;
    my ($id,@r) = split/\n/;
    next unless $id =~ /chr[XVI]/;
    $hash{$id}{SEQ} = join('', @r);
    $total += length($hash{$id}{SEQ});
}

my $meanlen = int($total/$fragnum);

print STDERR $meanlen;

####################################
# split chromomoes 
####################################

$/ = ">";

my %new;
foreach my $chr ( keys %hash ) {
    if ( $restrict ) {
	next unless $chr eq $restrict;
    }
    
    # 

    my $seq = $hash{$chr}{SEQ};
    my $residlen=length($seq);
    my $init=$newid+1;
    print STDERR ">$chr $residlen [$init]";

    unless ( $fragnum==1 ) {
	until ( $residlen < ($meanlen+10000) ) {
	    my $stop = $residlen;
	    my $start = $residlen - $meanlen;
	    my $gap = int(10**rand(4));
	    $residlen = $start - $gap;
	    print STDERR $start, $stop, $gap, $residlen;
	    
	    $newid++;
	    
	    %{$new{$newid}} = (
		START => $start,
		STOP => $stop,
		OLD => $chr,
		NEW => $newid,
		SEQ => substr($seq, $start, $meanlen),
		GFF => []
		);
	}
    }
    
    $newid++;    
    %{$new{$newid}} = (
	START => 1,
	STOP => $residlen,
	OLD => $chr,
	NEW => $newid,
	SEQ => substr($seq, $start, $meanlen),
	GFF => []
	);

    # 

    #print STDERR $lo,$hi;

    my $old;
  GFF:foreach my $gff ( @{$hash{$chr}{GFF} } ) {	
      
      # 
      
      if ( $old ) {
	  if ( $gff->[3] >= $old->{START} && $gff->[3] <= $old->{STOP} ) {
	      push @{$new{$old->{NEW}}{GFF}}, $gff;
	      next;
	  } elsif ( $gff->[4] >= $old->{START} && $gff->[4] <= $old->{STOP} ) {
	      push @{$new{$old->{NEW}}{GFF}}, $gff;
	      next;
	  }
      }
      
      # 
      
      foreach my $bl ( $init..$newid ) {
	  die($bl) unless exists $new{$bl};
	  my $seg = $new{$bl}; 
	  print $seg->{START}, $seg->{STOP}, $gff->[3];
	  if ( $gff->[3] >= $seg->{START} && $gff->[3] <= $seg->{STOP} ) {
	      push @{$new{$bl}{GFF}}, $gff;
	      $old=$seg;
	      next GFF;
	  } elsif ( $gff->[4] >= $seg->{START} && $gff->[4] <= $seg->{STOP} ) {
	      push @{$new{$bl}{GFF}}, $gff;
	      next GFF;
	  }	    
      }	
      # gaps! 
      # die(join("\t",@{$gff}));
  }
}    
 
####################################
# output 
####################################

open(GFF, ">".($restrict ? $restrict."_" : undef).$fragnum."_".$gaplen."_$file");
open(FSA, ">".($restrict ? $restrict."_" : undef).$fragnum."_".$gaplen."_$file".'.fsa');
foreach my $id (sort {$a <=> $b} keys %new) {
    print STDERR $id,$new{$id}{OLD},$#{$new{$id}{GFF}};
    
    foreach my $gff ( @{ $new{$id}{GFF} } ) {
	$gff->[0] = $id;
	$gff->[3] = ($gff->[3] < $new{$id}{START} ? $new{$id}{START}-$gff->[3]+1 : 1);
	$gff->[4] = ($gff->[4] < $new{$id}{STOP} ? $new{$id}{START}-$gff->[4]+1 : $new{$id}{STOP});
	print GFF join("\t", @{$gff});
    }
    
    print FSA ">$id\n$new{$id}{SEQ}";
}
