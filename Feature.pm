package Annotation::Feature;

# inheritance 

use Annotation;
@ISA = qw(Annotation);

#########################################
#########################################

# autoloaded methods 

our $AUTOLOAD;

my %AUTOMETH = (
	coord => undef,
	score => undef,
	feature => undef,
	strand => undef
);

sub AUTOLOAD {
    my $self = shift;    
    my $name = uc((split(/\:/, $AUTOLOAD))[-1]); # much faster 
    #$name =~ s/.*://;   # strip fully-qualified portion

    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

sub DESTROY {
    my $self=shift;
    
    if (my $link = $self->link) {           
	$link->link(undef) if $link->link eq $self;
	$self->link(undef);
    }

    return $self->SUPER::DESTROY;
}

#########################################
#########################################

sub new {
    my $class = shift;
    my $args = {@_};
    my $self  = $class->SUPER::new(@_);  

   foreach my $k (keys %AUTOMETH) {
		$k = uc($k);
		if ($self->id == 0) {
		    $self->{$k} = $AUTOMETH{$k};
		} else {
		    $self->throw("$k undefined -- $args->{$k}") unless 
				exists $args->{$k} && defined $args->{$k};
		    $self->{$k} = $args->{$k};
		}
    }

    bless $self, $class;
}

#########################################
#########################################

sub output {
  my $self = shift;
  my $string;
  foreach my $k (sort {$a cmp $b} keys %AUTOMETH) {
  	$string .= $self->$k."\t";
  }
  #print $string;
  return $self;
}

#########################################
#########################################

sub _invert_feature {
    my $self = shift;
    my $feature = $self->feature;    
	if ($feature =~ /_1/) {
		$feature =~ s/_1/_2/;
	} else {$feature =~ s/_2/_1/;}
	$self->feature($feature);
    return $self;
}

sub _invert_coords {
    my $self = shift;
    my $length = $self->up->length;    
    $self->coord($length - $self->coord + 1);
    return $self;
}

1;

