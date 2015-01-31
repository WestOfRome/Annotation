package Annotation;

# inheritance 

use ExApp;
use GlobalVars;
use Storable qw(nstore);
# use Bio::Root::Root;
use Carp;
use Data::Dumper;
use LWP::Simple;

@ISA = qw(ExApp GlobalVars);

# uses closure to get unique ids

$unique_id;# 
BEGIN {
    sub _id_generator {
	my $init = shift;
	my $iterator = sub {$init++;};
	return $iterator;
    }
    $unique_id = &_id_generator( time**2 );
}

my $STORE = 0;

#########################################
#########################################

# autoloaded methods 

our $AUTOLOAD;

my %AUTOMETH = (
		id => undef,
		left => undef,
		right => undef,
		up => undef,
		down  => undef # there is actually a method for this
		); 

sub AUTOLOAD {
    my $self = shift;    
    my $name = $AUTOLOAD;
    my $name = uc((split(/\:/, $AUTOLOAD))[-1]); # remove fully qualified 
    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

sub _method { (split(/::/, (caller(1))[3]))[-1] } 

#########################################
#########################################

=head2 new()

    Make new Annotaion object. 
    Called by Instance object (Genome or ORF etc).

    UP => undef,
    DOWN => undef, 
    LEFT => undef,
    RIGHT => undef,
    ID => undef,

=cut 

sub new {
    my $class = shift;
    my $args = {@_};
    my $self = bless {}, $class;

    foreach my $k ( map {uc($_)} keys %AUTOMETH) {
	if (exists $args->{$k}) {
	    $self->{$k} = $args->{$k};
	} else {
	    $self->{$k} = $AUTOMETH{$k};
	}
    }

	# ID is initially set to _ID unless 
	# supplied by user (as when making a base ID => 0) 

    $self->{_ID} = &{$unique_id};
    $self->id($self->_internal_id) unless defined $self->id;

    return $self;
}   

=head2 clone()

    Makes a shallow copy of calling object.
    Copies attributes which point to HASH and ARRAY refs 
    but NOT objetcs. All objetcs created this way are 
    not in Genome object (though they can access it). 

=cut 

sub clone {
    my $self = shift;
    
    # clones attributes other than references
    # -- a superficial clone that is not part of $genome.
    # can be manipulated without affecting other objects etc
    # does not do deep cloning. takes one level. 
    # will need to be rewritten if more complex objects used. 
    
    my @r;
    foreach my $atr (keys %{$self}) {
	if ($atr eq 'UP') {
	    push @r, ($atr => $self->{$atr});
	} elsif (ref($self->{$atr}) =~ /::/) {
	    next;  # no references
	} elsif (ref($self->{$atr}) eq 'ARRAY') {
	    # copy values not reference
	    push @r, ($atr => [@{$self->{$atr}}]);
	} elsif (ref($self->{$atr}) eq 'HASH') {		
	    my $hash = $self->{$atr};
	    push @r, (
		$atr => 
		{ map { $_ => $hash->{$_} } grep {$hash->{$_} !~ /::/} keys %{$hash} }
	    );
	} else {push @r, ($atr => $self->{$atr});}
    }
    
    return ref($self)->new(@r);
}

=head2 DESTROY()

    Custom DESTROY routine. Call for all objects. 
    Recurses through daughters and extracts from parents.

=cut

sub DESTROY {
	my $self = shift;

	# recurse through daughter objects

	map { $_->DESTROY } grep {defined} 
	($self->stream, $self->_down); # base must be last
	$self->throw if $self->_down || $self->down;

	# what do I need to deal with ?
	
	my $up = $self->up;
	my $left = $self->left;
	my $right = $self->right;	

	# remove from parent object.
	
	if ( $up && $up->_down eq $self ) { # base objects (object 0)
	    # we can only remove base from an object as 
	    # part of ultimate destruction. all other 
	    # dependent objects must be removed first
	    $self->throw if $left || $right; 
	    $up->{'DOWN'}=undef;
	    
	} elsif ( $up && ! ($self->left || $self->right) &&
		  $up->_down && $up->_down->left eq $self ) { # object 1
	    $self->throw unless $up->_down->right eq $self; 
	    # remove from parent .. 
	    $up->remove(-object => $self, -warn => 0);

	} elsif ($up && ($left || $right)) { # objects 2 ... N 
	    # will be tested in remove whether 
	    # left/right have same parent 
	    $up->remove(-object => $self, -warn => 0);

	} elsif ($up) { # other odd objects -- 
	    # pseudos/features may look like this
	    # know their parent (maybe to get access to sequence)
	    # but have no siblings. 
	    
	} elsif ($left || $right) {
	    $self->throw("Objects MUST know parent");			
	} else {}  # no up, no left, no right == genome
	
	# destroy
	
	$self->up(undef);
	$self->left(undef);
	$self->right(undef);
	$self->{'DOWN'} = undef; # no setter for 'DOWN'
	
	# shoudl be clear ....

	undef $self;
	return 1;
}

#########################################
#########################################

# the 'base' system is meant to be fully transparent

=head2 down(-direction => 'left|right|1|-1')

    Returns objects below caller in specified direction. 
    Remember, all object levels are closed circles so 
    the object to the left (-1) of the first gene on a contig 
    is the last gene. 

    Defaults to 1 (right). 

=cut

sub down {
    my $self = shift;	
    my $args = {@_};    
    
    return undef unless $self->_down; # Features || Exons
    
    if ($args->{'-direction'} eq 'left' || $args->{'-direction'} == -1) {
	return $self->_down->left;
    } else {
	return $self->_down->right;
    }
}

=head2 traverse(-direction => 'left|right', -distance => infinity)

    Return all objetcs along datastructure in specified dir from calling object.
    
    @r = $obj->traverse(-dir => 'left', -dist => 20) 

    Gathers up to 20 objects (not incl. caller) in specificed direction and
    returns those objects + the caller. Never retuns base object even if it
    is caller. 

=cut

sub traverse {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-distance'} = $INFINITY unless exists $args->{'-distance'};
    $args->{'-array'} = [] unless exists $args->{'-array'};
    $args->{'-direction'} = 'right' unless exists $args->{'-direction'};
	
    push @{$args->{'-array'}}, $self 
	if ($self && $self->id != 0);
    my $method = $args->{'-direction'}; 		
    my $next = $self->$method();
		
    # scalar(@{$args->{'-array'}}) == $args->{'-distance'}+1 
    unless (($#{$args->{'-array'}} == $args->{'-distance'}) 
	   || ! defined $next)  {
	$next->traverse(
	    -direction => $args->{'-direction'}, 
	    -distance => $args->{'-distance'},
	    -array => $args->{'-array'}
	    );			
    }
    
    return @{$args->{'-array'}};
}

=head2 stream()

    Return all daughter objetcs to caller. 
    e.g. if called on an orf object, returns all exons

=cut

sub stream {
    my $self = shift;
	return () unless $self->_down;     # Features || Exons
	return($self->_down->traverse(@_));
}

sub _stream {
    my $self = shift;
    return () unless $self->_down;     # Features || Exons
    return($self->_down, $self->_down->traverse(@_));
}

=head2 index(-method => 'start')

    Reindexes all daughter objects on calling object.
    e.g. if 5 genes are added to a contig, can renumber all 
    genes (by start coord or any method) so that all genes 
    are sequentially labelled.

=cut

sub index {
    my $self = shift;             
    my $args = {@_};
    return $self unless my $base = $self->_down;
    
    my $index = (exists $args->{'-start'} ? $args->{'-start'} : 1); # -start for super-contigs 
    my $method = (exists $args->{'-method'} ? $args->{'-method'} : 'start'); 
	
    my @r = (
	($self->strand && $self->strand == -1) ? 
	sort {$b->$method <=> $a->$method} $self->stream : 
	sort {$a->$method <=> $b->$method} $self->stream # return 0 unless $a && $b;
	);
	
    my $left;	
    foreach my $orf (@r) {		
	$orf->id($index); 
	
	if ($left) {
	    $orf->left($left);
	    $left->right($orf);
	} else {
	    $orf->left(undef);
	    $base->right($orf);
	}
	$left = $orf;
	$index++;
    }    
    $left->right(undef) if $left;
    $base->left($left);

    return $index-1; # -1 is adjustment for base object 
}

=head2 add(-object => obj, -left obj_23, -right => obj_24)

    The specified daughter object will be added to the caller
    between the objects currently indexed at 23 and 24. 

    Specifying left/right is optional. 

=cut 

sub add {
    my $self = shift;
    my $args = {@_};
    my $new = $args->{'-object'};
    
    $self->throw("Must add object: -object => $new")
    	unless $new && $new =~ /::/;
    #$self->throw if $self->contains($new);
    $self->throw("Cannot add a base to an object -- supplied on creation only")
    	if $new->id == 0; # if this is changed must change code below 	
    
    # if we are moving the last orf from one contig to another 
    # the base in the old object will still point to orf -- 
    # need to do a small intervention to prevent this.
    # only a problem of cource if user has not called ->remove
    
    if ($new->up && ($new->left || $new->right)) {
	$new->up->remove(-object => $new, -force => 1)
	    unless $new->up eq $self;
    } elsif ($new->up) {
	my $base = $new->up->_down;
	if ($base->left eq $new || $base->right eq $new) {
	    $new->up->remove(-object => $new, -force => 1)
		unless $new->up eq $self;
	}
    }
    
    # can add to end of list of genes 
    # to a specified object/location on contig
    
    if (exists $args->{'-left'} && exists $args->{'-right'}) {
	$self->warn("$args->{'-left'} and $args->{'-right'} are not neighbours\n")
	    unless $args->{'-left'}->right eq $args->{'-right'};
    } elsif (exists $args->{'-left'}) {
	$args->{'-right'} = $args->{'-left'}->right;
    } elsif (exists $args->{'-right'}) {
	$args->{'-left'} = $args->{'-right'}->left;
    } else {		
	if (my $left = $self->down(-direction => 'left')) {
	    $args->{'-left'} = $left;
	} else {
	    $args->{'-left'} = $self->_down;
	}
	$args->{'-right'} = $self->_down;
    }
    
    $self->throw("Neighbours not defined: $args->{'-left'}, $args->{'-right'}")
    	unless defined $args->{'-right'} && defined $args->{'-left'};
    
	# make updates 
    
    $new->left($args->{'-left'}) unless $args->{'-left'}->id == 0;
    $new->right($args->{'-right'}) unless $args->{'-right'}->id == 0;
    $args->{'-left'}->right($new); 
    $args->{'-right'}->left($new); 
    
    $new->up($self);
    return $self;
}

=head2 remove(-object => obj, -force => 0|1)

    This is the ONLY sanctioned way to remove a daughter 
    object from the caller. NB: The object is undamaged
    with the sole exception that it is safely excised 
    from parent object. The parent (caller) is healed and returned.

    e.g. $genome = $genome->remove(-object => contig_obj)

    Use -force to suppress a warnign if removing last 
    object from parent. 
    
=cut

sub remove {
    my $self = shift;
    my $args = {@_};
    my $obj = $args->{'-object'};

    # 

    $self->throw("Must supply an object!: $obj ($self)")
	unless ref($obj) eq ref($self->_down) || ref($obj) =~ /Feature/;
    $self->throw("$obj not attached to $self")
    	unless $obj->up eq $self;			
    $self->throw("Not attached to self") 
	unless $self->_down eq $obj || $self->contains($obj);
    $self->throw("Cannot remove base from an object")
    	if $obj->id == 0;

    # defaults 
    
    $args->{'-warn'} = $args->{'-force'} if exists $args->{'-force'};
    $args->{'-warn'} = 1 unless exists $args->{'-warn'};
    
    # get base and  neighbours 

    $self->throw unless my $base = $self->_down;

    my $left = $obj->left;
    my $right = $obj->right;
    map {  $self->throw("$_ (L/R) not attached to $self") unless $_->up eq $self } 
    grep { defined }  ($left,$right,$base);

    # rewire neighbours 
    
    if ($left && $right) {
	$right->left($left);
	$left->right($right);
    } elsif ($left) {
	$base->left($left);
	$left->right(undef);
    } elsif ($right) {
	$base->right($right);
	$right->left(undef);
    } else {
	$base->left(undef);
	$base->right(undef); 
	$self->warn("Removing only ".ref($obj)." from ".ref($self)." ".$self->id)
	    unless $args->{'-warn'}==0;
    }
    
    # destroy pointer to parent and neighbours. 
    # descendants, orthologs, ohnologs unaffected 
    # -> we might still want to do somethign with object 

    $obj->left(undef);
    $obj->right(undef);
    $obj->up(undef);
    $obj = undef;

    return $self;
}

sub delete {
    my $self = shift;
    return $self->remove(@_);
}

=head2 contains(obj)

    Test if object is owned by caller by querying all 
    left/right relationships (does not use UP).
    Does not test multiple presence pathology. Returns
    TRUE (obj) on first detection in the data structure.

=cut 

sub contains {
    my $self = shift;
    my $obj = shift;
    $self->throw unless $self->_down->isa( ref($obj) );

    foreach my $o ( $self->stream ) {
	return $obj if $o eq $obj; 
    }
    
    return();
}

=head2 swap( -object => obj )

    Swap locations of two objects with respect to neighbours. 
    Useful for randomizations. 

    Reversed by calling index().

=cut 

sub swap {
    my $self = shift;
    my $args = {@_};
    my $other = $args->{'-object'};

    $self->throw if $self->left eq $other || $self->right eq $other;
    my ($l1,$r1,$u1) = ($self->left,$self->right,$self->up);
    my ($l2,$r2,$u2) = ($other->left,$other->right,$other->up);
    
    # 
    
    $self->left($l2);
    $l2->right($self);
    $self->right($r2);
    $r2->left($self);
    $self->up($u2);

    # 

    $other->left($l1);
    $l1->right($other);
    $other->right($r1);
    $r1->left($other);
    $other->up($u1);

    return $self;
}

=head2 transfer(-from => , -to => , -warn => 1|0, -index => 1|0, -fast => undef)

    Transfer an orf object from one contig to another. 
    Uses add() to move the caller to the new contig so the 
    binding is reciprocal. Orf knows contigs and vice versa. 
    
    Returns self.

  NB: User is responsible for ensuring coordinates are 
    meaningful in new context (eg $orf->translatable()).  

=cut 

sub transfer {
    my $self = shift;
    my $args = {@_};    

    $args->{'-warn'} = 1 unless exists $args->{'-warn'};
    $args->{'-index'} = 1 unless exists $args->{'-index'};
    $args->{'-fast'} = undef unless exists $args->{'-fast'};

    $args->{'-from'} = $self->up unless exists $args->{'-from'};
    $args->{'-to'} = shift if $#_ == 0; # 
    my ($from, $to) = ($args->{'-from'},$args->{'-to'});

    $self->throw unless ref($from) eq ref($to);
    $self->throw unless $self->isa(ref($from->_down));
    # $self->throw if $to->contains($self) unless $args->{'-fast'};
    # this is redundant with use of contains() below. 

    # 

    if ( $args->{'-fast'} || $from->contains($self) ) {
	$from->remove(-object => $self, -warn => $args->{'-warn'});
    } elsif ( $from eq $self->up ) {
	$self->up(undef);
    } else { $self->throw; }

    # 

    $to->add(-object => $self);
    $to->index unless $args->{'-index'} == 0;

    # 

    return $self;
}

#########################################
#########################################

=head2 output()

    Output caller object. Recurses though as many levels of 
    objects as specified and dispatches work to appropriate 
    object ar each level.

=cut 

sub output {
    my $self = shift;
	foreach my $o ($self->stream) {
		$o->output(@_);
	}
	return $self;
}

=head2 debug()

    Dump object innards. 

=cut 

sub debug {
    my $self = shift;
    my $out = "$self -> ";
    foreach my $k (keys %{$self}) { 
	next if $k =~ /sequence/i;
	$out .= "$k:$self->{$k}, ";
    }       
    print $out;
    return $out;
}

=head2 dump

    Dump structure and die.

=cut 

sub dump {
    my $self = shift;
    my $args = {@_};
    print Dumper($self);
    return $self;
}

=head2 warn

    Custom warn with stack trace.

=cut 

sub warn {
    my $self = shift;   
    my ($string) = grep {!/\:\:/} @_;
    my @objs = grep { ref($_) =~ /\:\:/ } @_;
 
    my $go = $self;
    $go = $go->up until (ref($go) =~ /Genome/i || ! $go->up ); # DEVIN : THIS IS BAD
    $self->throw unless $go;
    my $fh = $go->{_LOG};

    unless ( $fh ) {
	my $stamp = $self->history(0)->{STAMP};
        open($fh, '>'.$go->organism.'.'.$stamp.'.log') || $self->throw;
        $go->{'_LOG'} = $fh;
    }
    
    $string = join(
	"\n",
	'>##################################',
	$string,
	join("\t", (caller(1))[0..4] ),
	$self->output(-string => 1, -creator => 1),
	#( map {$_->output(-string => 1, -creator => 1)} grep {defined} @objs ),     
	'###################################'
	)."\n";

    print $fh $string;
    print ($string);
    return $self;
}

=head2 die 

    Custom warn with stack trace.

=cut 

sub throw {
    my $self = shift;
    my $args = {@_};
    my $string = ($#_ == 0 ? shift(@_) : $args->{'-string'});

    if ( $args->{'-output'} && (caller(1))[3] !~ /output/i ) {
	$self->output( -fh => \*STDERR );
    }

    confess($string);
}

#########################################
#########################################

# private methods 

sub _down {
    my $self = shift;
    return $self->{'DOWN'};
}

sub _internal_id {
    my $self = shift;
    return $self->{'_ID'};
}
        

1;
