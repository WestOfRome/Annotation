package YGOB;

use Annotation;
use GlobalVars;
use Storable qw(nstore retrieve);

sub new {
    my $class = shift;
    my $file = shift || $ENV{YGOB_ORTHO_DB};
    die unless -e $file;
    my $sto = $file.'.sto';

    if ( -e $sto && -B $sto ) {
	return( retrieve( $sto ) );
    }

    # 

    my %hash;
    open(YGOB, $file) || die;
    while ( my $line = <YGOB> ) {
	chomp($line);
	next if $line =~ /^\#/;
	next if $line =~ /_CEN\d+\s+/;
	next if $line =~ /_RUF\d+\s+/;
	next if $line =~ /Scer_|rrna|trna|_RDN/; # ncRNAs

	my @parts = split(/\s+/, $line);	
	my $id = shift @parts;
	die unless $id =~ /^YGOB_\d+_DS$/;
	
	$hash{$id}{'ID'} = $id;
	$hash{$id}{'SD'} = pop @parts;
	$hash{$id}{'MEAN'} = pop @parts;
	$hash{$id}{'N'} = pop @parts;
	
	# 
	
	foreach my $gene ( grep { !/_YGOB_/ } grep {!/\-{3}/} @parts ) {
	    next unless my $sp = _species_key( $gene );
	    push @{$hash{$id}{$sp}}, $gene;
	    $hash{ $gene } = $hash{ $id };
	}
    }
    close YGOB;

    my $obj = bless \%hash, $class;
    nstore($obj, $sto) || die;
    print STDERR "Stored: $sto";
    return $obj;
}

sub access {
    my $self = shift;
    my $pillar = shift;
    my $species = shift;
    
    if ( $pillar =~ /YGOB_\d+_DS/ ) {
	die unless exists $self->{$pillar};
    } else {
	return undef unless exists $self->{$pillar};
    }
    return( defined $species ? @{$self->{$pillar}->{$species}} : $self->{$pillar} );
}

sub _species_key {
    die unless my $gene = shift;

    if ($gene =~ /([YA])[A-P][LR]\d{3}.+/) { # SGD 
    } elsif ($gene =~ /([A-Z]{4})\d+[A-Z]\d+[grst]?/) { # Genolevures 
    } elsif ($gene =~ /(\w{3,4})_\d+\.\d+[a-z]?/) { # Standard
    } else { die($gene); } 

    my $sp = uc($1);

    return( exists $HOMOLOGY{$sp} ? $sp :  $SPECIES{$sp} );
}

1;
