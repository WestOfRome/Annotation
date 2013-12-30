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
    bootstraps => undef
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

sub print {
    my $self = shift;
    my $args = {@_};

    my $fh = $args->{'-fh'} || STDOUT;
    
    my @data = $self->_write_tree_helper;
    my $string = join( ',', @data ).";" ;

    print {$fh} $string;
    exit;
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
