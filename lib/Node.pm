package Node;

#########################################
# inheritance 
#########################################

use Carp;
use GlobalVars;

@ISA = qw(GlobalVars);

#########################################
# 
#########################################

# autoloaded methods 

our $AUTOLOAD;

my %AUTOMETH = (
    # floss 
    name => undef,
    description => undef,
    # Required for valid tree structure 
    root => 0,
    ancestor => undef,
    descendents => undef,
    # Numerical 
    branch_length => undef,
    bootstrap => undef
    # Methodology.
    # Removing -- these are tree attributes not node. 
    # data_type => undef, # codons, dN, dS etc
    # data_size => undef, # number of characters
    # reconstruction => undef, # e.g. NJ
    # rooting => undef # e.g. midpoint 
    ); 

sub AUTOLOAD {
    my $self = shift;    
    my $name = uc((split(/\:/, $AUTOLOAD))[-1]); # remove fully qualified 
    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

local $\ = "\n";
local $, = "\t";

#########################################
# Methods 
#########################################

=head2 new()

=cut 

sub new {
    my $class = shift;
    my $args = {@_};
    my $self = bless {}, $class;

    # 

    foreach my $k ( map {uc($_)} keys %AUTOMETH) {
	if (exists $args->{$k}) {
	    $self->{$k} = $args->{$k};
	} elsif ( ref($AUTOMETH{$k}) eq 'ARRAY' ) {
	    $self->{$k} = [];  
	} else {
	    $self->{$k} = $AUTOMETH{$k};	    
	}
	#print $k, $self->{ $k }, $self->lc($k);
    }

    # define root or ancestor 

    if ( $self->root ) {
	#$self->throw unless scalar($self->descendents) == $self->root;
	$self->throw unless $self->root == 2 || $self->root == 3;
	$self->throw if defined $self->branch_length;
	$self->throw if defined $self->ancestor;
    } else {
	#$self->throw unless $self->ancestor;
	#$self->throw unless $self->branch_length;
    }

    # should not have descendents at start. Nevertheless: 

    map { $self->throw unless $self->isa(ref($_)) } $self->descendents;

    return $self;
}   

=head2 stream
=cut

sub stream {  
    my $self = shift;
    return ( (map { $_->stream } $self->descendents), $self );
}

=head2 descendents()
=cut 

sub descendents { 
    my $self = shift; 
    $self->throw if @_;
    return @{$self->{'DESCENDENTS'}};
}

=head2 ancestor()
=cut 

sub ancestor { 
    my $self = shift; 
    $self->throw if @_;
    return $self->{'ANCESTOR'};
}

=head2 add( -oject => $node )
=cut 

sub add {
    my $self = shift;
    my $args = {@_};

    # 

    my $obj = $args->{'-object'};
    $self->throw unless $self->isa(ref($obj));

    # 
    
    $self->throw if $obj->ancestor;
    map { $self->throw if $_ eq $obj } $self->descendents;

    # 
 
    push @{$self->{'DESCENDENTS'}}, $obj;
    $obj->{'ANCESTOR'} = $self ;    
    map { $self->throw unless $_->ancestor eq $self } $self->descendents;

    return $self;
}

=head2 remove( -object => $node )
=cut 

sub remove {
    my $self = shift;
    my $args = {@_};

    # 

    my $obj = $args->{'-object'};
    $self->throw unless $self->isa(ref($obj));

    # 
    
    $self->throw unless $obj->ancestor eq $self;
    $self->throw unless scalar(grep { $_ eq $obj } $self->descendents)==1; 

    # 

    for my $i ( 0 .. @{$self->{'DESCENDENTS'}} ) {
	if ( $self->{'DESCENDENTS'}->[$i] eq $obj ) {
	    splice( @{$self->{'DESCENDENTS'}}, $i, 1);
	    last;
	}
    }
    $obj->{'ANCESTOR'} = undef;
    
    # basic QC

    map { $self->throw unless $_->ancestor eq $self } $self->descendents;

    return $self;
}

#########################################
# IO
#########################################

=head2 newick( -screen => [0|1], -file => [undef|name], -process => 0|1 )

    Returns tree (or subtree) as a newick format text string or 
    an array ref of branchlengths. The latter makes wild assumptions
    and is not generalized for now. 

    With options will optionally also write screen and file. 

    e.g. (Skud_1.2:0.0136875,Sbay_1.5:0.0269125,((Scer_1.22:0.005,
    Spar_1.6:0.0054):0.008,Smik_1.4:0.0152):0.0035375);

=cut 

sub newick {
    my $self = shift;
    my $args = {@_};

    my $fh = STDOUT;

    my @data = $self->_write_tree_helper;
    my $string = join( ',', @data ).";" ;
    
    print {$fh} $string if $args->{'-screen'}==1 || $args->{'-print'}==1;

    if ( defined $args->{'-file'} ) {
	my $file = (! -e $args->{'-file'} && $args->{'-file'} =~ /\.nwk$/ 
		    ? $args->{'-file'}.'.nwk' 
		    : $args->{'-file'});
	
	open(my $fh2, $file);
	print {$fh2} $string;
	close $fh2
    } 

    ########################################
    # (Skud_1.2:0.0136875,Sbay_1.5:0.0269125,(Smik_1.4:0.016,
    # (Scer_1.22:0.005,Spar_1.6:0.0054):0.0072):0.0035375);
    if ( $args->{'-process'} ) {
	my @r;
	push @r, $1 while ($string =~ /\:([\d\.\-]+)/g); 
	#print join(', ', @r);
	return \@r;
    }
    ########################################

    return $string;
}

sub _write_tree_helper {
  my $self = shift;
  my @data;
  
  foreach my $node ( $self->descendents ) {
      push @data, $node->_write_tree_helper();
  }
  
  my $label = 
      ($self->descendents ? $self->bootstrap : $self->name ).
      ($self->root ? undef : ':'.$self->branch_length);
  
  if ( scalar(@data) >= 1) {
      $data[0] = "(" . $data[0];
      $data[-1] .= ")";
      $data[-1] .= $label;
  } else {
      push @data, $label;
  }

  return @data;
}

#########################################
# Other
#########################################

sub throw {
    my $self = shift;
    my $string = shift;
    confess($string);
}

#########################################
# End Methods 
#########################################

1;
