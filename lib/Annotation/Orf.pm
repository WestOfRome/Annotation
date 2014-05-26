#!/usr/bin/perl

package Annotation::Orf;

use Annotation;
use GlobalVars;
use YGOB;
use Node;

use File::Copy;
use GD::SVG;

@ISA = qw(Annotation);

use Annotation::Exon;                   

# Memory leakage ...  
# use constant HAS_LEAKTRACE => eval{ require Test::LeakTrace };
# use Test::More HAS_LEAKTRACE ? (tests => 1) : (skip_all => 'require Test::LeakTrace');
# use Test::LeakTrace;
         
#########################################
#########################################

# autoloaded methods orthogroup

my $YGOB;

our $AUTOLOAD;

# this needs to be fixed-- 
# defaults are never used. 

my %AUTOMETH = (
    START => {req => 1, def => undef},
    STOP => {req => 1, def => undef},
    STRAND => {req => 1, def => undef},

    EXONS => {req => 0, def => []},
    DATA => {req => 0, def => {}}, # converted to an hash ref by _init_data 
    EVIDENCE => {req => 0, def => 'NONE'},

    OGID => {req => 0, def => undef},
    FAMILY => {req => 0, def => undef},
    OHNOLOG => {req => 0, def => undef},
    DESCRIPTION => {req => 0, def => undef}, # this is set ONLY by impute()
    
    _CREATOR => {req => 0, def => undef},
    _DEBUG => {req => 0, def => undef}
    ); 

sub AUTOLOAD {
    my $self = shift;    
    my $name = (split(/\:/, $AUTOLOAD))[-1];
    #$name =~ s/.*://;   # strip fully-qualified portion
    $name = uc($name);

    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

sub DESTROY {
    my $self=shift;

    # 0. debugging.. 

    if ( my $debug = $self->_debug ) {
	#$self->_linked_pair_generic_get_set(-attribute => '_DEBUG', -object => undef);
	delete $self->{'_DEBUG'};
	delete $debug->{'_DEBUG'};
	$debug->DESTROY;
    } 
    
    # 1. ohnologs 
    
    if ( my $ohno = $self->{'OHNOLOG'} ) {
	delete $self->{'OHNOLOG'};
	delete $ohno->{'OHNOLOG'};
    } 
    
    # 2. orthogroups 

    $self->throw("Must call _dissolve_orthogroup first.") if $self->ogid;

    # 3. Neighbours / Up / Down 
    
    return $self->SUPER::DESTROY;
}

#########################################
#########################################

=head2 new()

    START => coord, 
    STOP => coord,
    STRAND => ±1,

    # not required 

    EXONS => [Exon objects], 
    UP => CONTIG :: this is usually set after creation by calling $contig->add(orf)

    # not recommended 
    
    _DEBUG => ORF
    _CREATOR => method
    OGID => i
    OHNOLOG => ORF
    DATA => {} # this is set by clone( -data => 1 ) 

    # not allowed

    LEFT => ORF 
    RIGHT => ORF 
    DOWN => EXON

    .. orthologs .. => ORF 
    DESCRIPTION 
    ID
    _INTERNAL_ID 

=cut 

sub new {
    my $class = shift;
    my $args = {@_};

    $args->{'_CREATOR'} = (caller(1))[3] unless exists $args->{'_CREATOR'};

    # 

    my $base = Annotation::Exon->new(ID => 0);
    my $self = $class->SUPER::new(@_, DOWN => $base);
    $base->up($self);

    # "reverse sort keys" is critical below. 
    # we need to set start/stop/strand before we can call _init_exons.
    # _init_data overrides EVIDENCE 
    
    foreach my $k ( map {uc($_)} reverse sort keys %AUTOMETH) {
	if ($self->id == 0) {
	    $self->{$k} = $AUTOMETH{$k}->{'def'};

	} elsif ( $AUTOMETH{$k}->{'req'} == 1) {
	    $self->throw("$k : $args->{$k} ") 
		unless exists $args->{$k} && defined $args->{$k};
	    $self->{$k} = $args->{$k};
	    
	} elsif ($k eq 'EXONS') {
	    $self->_init_exon_structure( @{$args->{$k}} );

	} elsif ($k eq 'DATA') {
	    $self->_init_data_structure( $args->{$k} || () );
	    
	} elsif ($k eq 'EVIDENCE') {
	    $self->{$k} = 'NONE'; # set by _init_data_structure 
	    
	} else {
	    $self->{$k} = $args->{$k}; 
	}
    }

    # coords must be supplied to consruct an Orf object 
    # but responsibility is devolved instantly to Exons

    delete $self->{START}; 
    delete $self->{STOP};  
    # delete $self->{STRAND};  
    
    # 

    return $self;
}

=head2 clone

    Over-ride method that clones Exons as well as Orf. 
    Necessary for sequence comparinsons between Orf 
    and a modified version of itself. 

=cut 

sub clone {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-data'} = undef unless exists $args->{'-data'};

    my %new = map { $_ => $self->data($_) } keys %{$self->_data()} if $args->{'-data'};

    my $new = ref($self)
	->new(
	UP => $self->up,
	START => $self->start,
	STOP => $self->stop,
	STRAND => $self->strand,
	EXONS => [ map { $_->clone } $self->stream ], # exons reattached to clone 
	(%new ? (DATA => \%new) : ()),
	_CREATOR => (caller(1))[3]
	);

    $new->evaluate( -structure => 0, -validate => 0 );
    return $new;
}

=head2 update(-blastdb => YGOB_AA_DB, -hmmerdb => YGOB_HMMER3_DB
    -reference => 'SGD', -full => 1, -evalute => 1)
    
    Re-gather evidence for the callling Orf after a coordinate: 
    BLAST, HMMER, exonerate, synteny, structure etc.

    Runs only on protein-coding genes. 

=cut

sub update {
    my $self=shift;
    my $args = {@_};
    
    return $self if $self->assign =~ /GAP|FEATURE|RNA/;

    # RUN by default 
    $args->{'-blast'} = $ENV{'YGOB_AA_DB'} unless exists $args->{'-blast'};
    $args->{'-hmmer'} = $ENV{'YGOB_HMMER3_DB'} unless exists $args->{'-hmmer'};
    $args->{'-exonerate'} = 1 unless exists $args->{'-exonerate'};
    $args->{'-hsp'} = 1 unless exists $args->{'-hsp'};
    $args->{'-synteny'} = 1 unless exists $args->{'-synteny'};
    # NOT run by default 
    $args->{'-repeats'} = undef unless exists $args->{'-repeats'};
    $args->{'-ncbi'} = undef unless exists $args->{'-ncbi'};
    $args->{'-kaks'} = undef unless exists $args->{'-kaks'};
    # moderate behaviour -- by default we EVAL but do not INIT 
    $args->{'-reference'} = 'SGD' unless exists $args->{'-reference'};   # BLAST -> impute 
    $args->{'-evaluate'} = 1 unless exists $args->{'-evaluate'};         # last action 
    $args->{'-initialize'} = undef unless exists $args->{'-initialize'}; # first action 
    # a hack 
    $args->{'-store'} = undef unless exists $args->{'-store'};
    # 
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};

    my $key = '_STORE_BLAST';
 
    #########################################   
    # 
    #########################################

    if ( $args->{'-initialize'} ) {
	$self->_init_data_structure( 
	    ref($args->{'-initialize'}) eq 'ARRAY'
	    ? ( '-exempt' => $args->{'-initialize'} ) ## NB ++ COUNTERINTUITIVE 
	    : () 
	    );
    }

    #########################################
    # 
    #########################################

    if ( $args->{'-hmmer'} ) {
	my ($hit) = 
	    $self->hmmer3(-application => 'hm3.hmmscan', -db => $args->{'-hmmer'});
	$self->accept('YGOB', $hit) if $hit;
    }
    
    #########################################
    # 
    #########################################

    if ( $args->{'-blast'} ) {	
	# my @bb = # phmmer still slower. extremeley so for biased composition. 	    
	#    $self->hmmer3(-application => 'hm3.phmmer', -db => $args->{'-blastdb'});	
	my $blast = 
	    $self->blast(-db => $args->{'-blast'}, -qseq => [$self]);
	my @bb = sort {$a->{EVALUE} <=> $b->{EVALUE}} @{$blast->{$self->_internal_id}};

	#
	
	if ( exists $blast->{$self->_internal_id} && $#bb > -1 ) { 
	    $self->homology( -blast => \@bb ); # store homologs etc 
	    # synteny moderated homologs. set GENE. set DESCRIPTION.   
	    $self->impute( -blast => \@bb, -reference => $args->{'-reference'} ); 
	    $self->data( $key => \@bb ) if $args->{'-store'};
	}
    }

    #########################################
    # 
    #########################################
    
    if ( $args->{'-ncbi'} ) {	
	my $blast = 
	    $self->blast(-db => 'nr', -qseq => [$self]);
	my @bb = sort {$a->{EVALUE} <=> $b->{EVALUE}} @{$blast->{$self->_internal_id}};
	
	#
	
	if ( exists $blast->{$self->_internal_id} && $#bb > -1 ) { 
	    $self->impute( -blast => \@bb, -homology => 'NCBI' ); # synteny moderated homologs 
	}
    }


    #########################################
    # 
    #########################################

    if ( $args->{'-repeats'} ) {	
	foreach my $rep ( qw(TY LTR) ) {
	    my $blast = 
		$self->blast(-db => $ENV{'ANNA_'.$rep.'_DB'}, -qmol => 'dna', -qseq => [$self]);
	    my @bb = sort {$a->{EVALUE} <=> $b->{EVALUE}} @{$blast->{$self->_internal_id}};
	    if ( exists $blast->{$self->_internal_id} && $#bb > -1 ) { 
		$self->impute( -blast => \@bb, -homology => $rep ); # synteny moderated homologs 
	    }
	}
    }

    #########################################
    # 
    #########################################

    if ( $args->{'-exonerate'} ) {
	$self->pillar( -force => 1 ); # guess the YGOB pillar
	if ( my $hom = $self->homolog( -fast => 1 ) ) { # choose a best homolog (from pillar if available)
	    $self->exonerate2( -return => 'score', -model => 'global' );  # this sets score('global')
	    $self->exonerate2( -return => 'score', -model => 'local' );   # required for HSP 
	}
	# include DNA exonerate and WISE to get a true gene score ?
    }

    #########################################
    # 
    #########################################

    if ( $args->{'-hsp'} ) {
	$self->fragment(); # tests difference between global and local scores. 
    }

    #########################################
    # 
    #########################################
    
    if ( $args->{'-synteny'} ) {
	$self->synteny( -distance => 2, -difference => 5, -spanning => 1, -set => 1 );
    }

    #########################################
    # 
    #########################################
 
    if (  $args->{'-kaks'} && $self->orthogroup ) {
	$self->throw("not implemented");
    }

    #########################################
    # 
    #########################################

    $self->evaluate(-structure => 1, -validate => 1) if $args->{'-evaluate'};
    $self->output if $args->{'-verbose'};

    return $self;
}

=head2 ibidologs(-species => [], -max => 5, -clean => 1, -trna => 1)

    Return genes at same genomic locus in other species
    that are NOT in orthogroups. Not orthologs, ibidiologs. 

    In the case where the region in the other species, is
    not bounded by OGs, we use consecutive numbering of genes
    in OGs to detemine which direction to look.

  NB: we do not return genes in calling species, even if neighbours
    of calling ORF are not in OGs. For this use $orf->gather(). 

=cut

sub ibidologs {
    my $self = shift;
    my $args = {@_};
    
    $self->throw unless my @sp = $self->up->up->bound;

    $args->{'-species'} = \@sp unless exists $args->{'-species'};
    $args->{'-max'} = 5 unless exists $args->{'-max'};
    $args->{'-og'} = 0 unless exists $args->{'-og'};
    $self->throw("Not implemented") unless $args->{'-og'} == 0;

    $args->{'-clean'} = 1 unless exists $args->{'-clean'};
    $args->{'-remove'} = [($args->{'-clean'} ? qw(GAP FEATURE REPEAT) : 'KEEPALLFEATS')] 
	unless exists $args->{'-remove'};    
    $args->{'-trna'} = 1 unless exists $args->{'-trna'};
    push @{$args->{'-remove'}}, 'TRNA' unless exists  $args->{'-trna'};
    my $excl = join('|', @{$args->{'-remove'}});    

    ######

    my ($left, $right) = $self->neighbours(-orthogroup => 1);
    return undef unless ($left || $right);
    
    my %hash;
    foreach my $sp_meth ( map { lc($_) } @{ $args->{'-species'} } ) { # names of bound species 	      
	my ($o_left,$farleft) = 
	    ( $left ? ($left->$sp_meth,$left->neighbour(-direction => 'left', -orthogroup => 1)) : (undef,undef));
	my ($o_right,$farright) = 
	    ( $right ? ($right->$sp_meth,$right->neighbour(-direction => 'right', -orthogroup => 1)) : (undef,undef));
	my $o_farleft = $farleft->$sp_meth if $farleft;
	my $o_farright = $farright->$sp_meth if $farright;
	
	# gather candidates 

	if ( ($o_left && $o_right) && $o_left->up->id == $o_right->up->id  ) {
	    push @{$hash{uc($sp_meth)}}, 
	    $o_left->intervening($o_right); # grep { $_->assign !~ /GAP|FEATURE|REPEAT/ }
	} else {
	    if ( $o_farleft && $o_farleft->up->id == $o_left->up->id ) {
		my $dir = ($o_farleft->id < $o_left->id ? 'right' : 'left' ); 
		push @{$hash{uc($sp_meth)}}, 
		$o_left->traverse(-distance =>  $args->{'-max'}, -direction => $dir); 
	    }
	    if ( $o_farright && $o_farright->up->id == $o_right->up->id ) {    
		my $dir = ($o_farright->id > $o_right->id ? 'left' : 'right' );
		push @{$hash{uc($sp_meth)}}, 
		$o_right->traverse(-distance =>  $args->{'-max'}, -direction => $dir); 
	    }
	}
    }

    # filter candidates 

    foreach my $sp ( keys %hash ) {
	$hash{$sp} = [ grep { $_->assign !~ /$excl/ } @{$hash{$sp}} ];

	if ( 1 || defined $args->{'-og'} ) {  ### FORCE TRUE 
	    $hash{$sp} = ( 
		$args->{'-og'} 
		? [ grep { $_->ogid } @{$hash{$sp}} ] 
		: [ grep { ! $_->ogid } @{$hash{$sp}} ] ## ONLY CURRENT OPTION 
		);
	}
    }

    # 

    return ( $#{$args->{'-species'}}==0 ? ( map {@{$_}} values %hash ) : \%hash);
}

=head2 bracket( -object => $other, -og => -1, -clean => 1, -trna => 1)

    Return all genes bracketed by the orthogroups self/other
    as a hash indexed by species. Returns undef if hash empty.
    
    Calls intervening().

    -clean : exclude gaps / repeats / features
    -trna : [0/-1]/1 exclude/include. 1 is default 
    -og : -1/0/1 exclude/ambivalent/require Orthogroups. -1 is default.  

=cut 

sub bracket {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-og'} = -1 unless exists $args->{'-og'};
    $args->{'-max'} = $INFINITY unless exists $args->{'-max'};
    $args->{'-min'} = 0 unless exists $args->{'-min'};    

    $args->{'-clean'} = 1 unless exists $args->{'-clean'};
    $args->{'-remove'} = [($args->{'-clean'} ? qw(GAP FEATURE REPEAT) : 'KEEPALLFEATS')] 
	unless exists $args->{'-remove'};
    $args->{'-trna'}=1 unless exists $args->{'-trna'};
    push @{$args->{'-remove'}}, 'TRNA' unless $args->{'-trna'} == 1;
    my $excl = join('|', @{$args->{'-remove'}});
    
    my $other = $args->{'-object'};
    $self->throw unless $self->isa( ref( $other ) );
    $self->throw unless $self->orthogroup && $other->orthogroup;
    $self->throw unless $self->up eq $other->up;

    # populate data structure 

    my %genes = (
	uc($self->organism) => [ $self->intervening($other) ]
	);    
    foreach my $sp ( $self->up->up->bound ) {
	my $so = $self->$sp;
	my $oo = $other->$sp;
	next unless $so->up eq $oo->up;
	my $d = $oo->distance(-object => $so, -bases => 1);
	next unless $d <= $args->{'-max'} && $d >  $args->{'-min'} ; 
	$genes{ uc($sp) } = [ grep {defined} $so->intervening($oo) ];
    }

    # 
    
    foreach my $sp ( keys %genes ) {
	$genes{$sp} = [ grep { $_->assign !~ /$excl/ } grep {defined} @{$genes{$sp}} ];
	
	if ( $args->{'-og'} ) {
	    $genes{$sp} = ( 
		$args->{'-og'} == 1
		? [ grep { $_->ogid } @{$genes{$sp}} ] 
		: [ grep { ! $_->ogid } @{$genes{$sp}} ] 
		);
	}
    }
    
    # 

    map { delete $genes{$_} if $#{$genes{$_}}==-1 } keys %genes; 

    return ( scalar(keys %genes)==0 ? undef : \%genes );
}

=head2 gather( -clean => 1, -trna => 1 )
    
    Related to the bracket method but does not require 
    OGs on either side. Gathers all genes not in OGs at 
    this locus in calling species and in all other species. 

    Calls ibidologs(). 

=cut

sub gather {
    my $self=shift;
    my $args = {@_};
    
    $args->{'-og'} = 0 unless exists $args->{'-og'};
    $self->throw("Not implemented") unless $args->{'-og'} == 0;

    $args->{'-clean'} = 1 unless exists $args->{'-clean'};
    $args->{'-remove'} = [($args->{'-clean'} ? qw(GAP FEATURE REPEAT) : 'KEEPALLFEATS')] 
	unless exists $args->{'-remove'};    
    $args->{'-trna'} = 1 unless exists $args->{'-trna'};
    push @{$args->{'-remove'}}, 'TRNA' unless exists  $args->{'-trna'};
    my $excl = join('|', @{$args->{'-remove'}});    

    ###########

    return undef if $self->orthogroup;
    
    my ($left,$right) = $self->neighbours(); 
    my (@left,@right);
    until ( ! $left || $left->ogid ) {
	push @left, $left;
	$left = $left->left;
    }
    until ( ! $right || $right->ogid ) {
	push @right, $right;
	$right = $right->right;
    }

    # 

    my %hash = %{$self->ibidologs(%{$args})};
    $hash{uc($self->organism)} = [ sort { $a->id <=> $b->id } (@left,$self,@right) ];
    
    #

    foreach my $sp ( keys %hash ) {
	$hash{$sp} = [ grep { $_->assign !~ /$excl/ } @{$hash{$sp}} ];

	if ( defined $args->{'-og'} ) {
	    $hash{$sp} = ( 
		$args->{'-og'} 
		? [ grep { $_->ogid } @{$hash{$sp}} ] 
		: [ grep { ! $_->ogid } @{$hash{$sp}} ] 
		);
	}
    }
		  
    return \%hash;
}

=head2 ohnolog(-object => newobj|undef, -score => i, -warn => )

    Get/set the genes ohnolog. If one already exists
    it may not be altered directly -- must first unset 
    by calling with -object => 'undef'. If no values are passed 
    we just return the existing ohnolog. To access the score 
    use $self->score('ohnolog');

 NB: This method utilizes _linked_pair_generic_get_set to 
    automatically handle reciprocals relationships.

=cut 

sub ohnolog {
    my $self = shift;
    my $args = {@_};
    
    my $attr = 'OHNOLOG';

    $args->{'-object'} = $_[0] if ( ! exists $args->{'-object'} && @_ );
    #$args->{'-score'} = 0 unless $args->{'-score'};
    
    if ( $args->{'-score'} ) {
	$self->throw('Depracated: Use $orf->accept(ohno) to store evidence.');
	
	if ( @_ ) { # we are setting not getting 
	    my $key =  '__'.$attr; # it is a score so gets '__'
	    $self->data( $key => ( ! $args->{'-object'} ? undef : $args->{'-score'}) );
	    if ( my $ohno = $self->ohnolog ) {
		$ohno->data( $key => ( ! $args->{'-object'} ? undef : $args->{'-score'}) );
	    }
	}
    }
    
    return $self->_linked_pair_generic_get_set(
	'-attribute' => $attr,
	(exists $args->{'-object'} ? ('-object' => $args->{'-object'}) : ()), 
	(exists $args->{'-warn'} ? ('-warn' => $args->{'-warn'}) : ())
	);
}

=head2 critique

    Return 0-1 value representing the quality of the prediction
    compared to the assocaited (SGD) gene model. 

    Returns undef if no SGD (or other gold standard) model has 
    been associated. The association process uses $genome->compare(file.gff)
    and is nontrivial. Hard even in non-fragmented genomes. 

    The score reflects base pair level overlap with penalty for overcalling. 

=cut

sub critique {
    my $self = shift;
    return undef unless $self->translatable(-fast => 1);
    return undef unless my $sgd = $self->_debug;
    $self->throw unless $self->isa(ref($sgd));
    
    my @b1 = map { $_->start..$_->stop } sort {$a->start <=> $b->start} $self->stream;
    my @b2 = map { $_->start..$_->stop } sort {$a->start <=> $b->start} $sgd->stream;
    $self->throw unless ( scalar(@b1)%3==0 && scalar(@b2)%3==0 );
    
    my %f1;
    for (my $i=0; $i <= $#b1; $i+=3) {
	$f1{$b1[$i]}++;
    }
    
    my $count;
    for (my $i=0; $i <= $#b2; $i+=3) {
	$count++ if exists $f1{$b2[$i]};
    }
    
    my $missed = @b2/3 - $count;
    my $extra =  @b1/3 - $count;
    my $frac = ($count - $missed - $extra)/(@b2/3);
    #print $count, $missed, $extra, $frac;
    return ( $frac < 0 ? 0 : $frac );
}

# _debug is used to hold another gene model that the
# present one should be compared to during debugging
# eg. if we read in SGD gff and make objects, we can 
# comepare real to .. imagined.

sub _debug {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-object'} = $_[0] if ( ! exists $args->{'-object'} && $#_ == 0 );
    
    return $self->_linked_pair_generic_get_set(
	-attribute => '_DEBUG',
	($args->{'-object'} ? ('-object' => $args->{'-object'}) : ()), 
	);    
}


# this method provides generic get/set/unset 
# methods for pairs of attributes that point to
# to other ORFs and that contain a reciprocal link. 
# ohnologs are the type example. 
# to unset pass the undefined value as the method object.
# $orf->_linked_pair_generic_get_set( -object => undef )

sub _linked_pair_generic_get_set {
    my $self = shift;
    my $args = {@_};

    $args->{'-warn'} = 1 unless exists $args->{'-warn'};

    my $ohno = $args->{'-object'};
    my $attr = $args->{'-attribute'};
    my $method = lc($attr);

    ################################################
    # test to see if we have passed 'undef' value. 
    ################################################

    # check old values if we are trying to alter
    # should not be passing undef unless ohno exists
    
    if ( exists $args->{'-object'} ) {
	if ( my $old = $self->$method ) {
	    $self->throw unless $old->$method eq $self;
	}
    }

    ################################################
    # check new values if we are tryign to update 
    ################################################
    
    # must not already be assigned etc 
    
    if ( $ohno ) {
	$self->throw unless $self->isa( ref( $ohno ) );	
	if ( my $old = $ohno->$method ) {
	    $self->throw unless $old->$method eq $ohno;	    
	    $self->throw;
	}
    }

    ################################################
    # do the work 
    ################################################

    if ( exists $args->{'-object'} && $ohno ) { # set 
	$self->throw if $self->$method;	
	$self->{$attr} = $ohno;
	$ohno->{$attr} = $self;
    } elsif ( exists $args->{'-object'} ) { # unset 
	$self->throw unless my $old = $self->$method;
	$old->{$attr} = undef;
	$self->{$attr} = undef;
    } else { # get 
	return $self->{$attr};
    }

    return $self;
}

=head2 _filter_tandems

    Takes an array of ORF objects and returns
    a set of genes that excludes possible tandem
    pairs (ie one from each tandem arrasy is chosen).

=cut 

sub _filter_tandems {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-distance'} = 20 unless exists  $args->{'-distance'};
    
    my $attr = '_filter_temp_var'.int(rand(100));
    map { $self->throw unless $self->isa(ref($_)) } @{$args->{'-object'}};
    
    my %hash;
    foreach my $o ( $self, @{$args->{'-object'}} ) {
	push @{$hash{$o->up->id}},$o;
    }
    
    my @uniq;
    foreach my $chr ( keys %hash ) {
	push @uniq, $hash{$chr}->[0] and next if $#{$hash{$chr}}==0;
	
	my %cl = map { $_ => [ $hash{$chr}->[$_] ] } 0..$#{$hash{$chr}};
	map { $cl{$_}->[0]->data($attr => $_) } keys %cl;

	for my $i ( 0..($#{$hash{$chr}}-1) ) {
	    my $g1 = $hash{$chr}->[$i];
	    $self->throw unless $g1->data($attr)==$i;
	    for my $j ( ($i+1)..$#{$hash{$chr}} ) {
		next if $cl{$i} eq $cl{$j};
		my $g2 = $hash{$chr}->[$j];
		
		if ( $g1->distance(-object => $g2) <= $args->{'-distance'} ) {
		    my @move = @{$cl{$j}};
		    foreach my $x ( @move ) { 
			$cl{ $x->data($attr) } = $cl{$i}; # update %cl so x{i} -> $i
			push @{ $cl{$i} }, $x;
		    }
		}
	    }
	}
	
	my %uniq;
	foreach my $array ( grep {!$uniq{$_}++} values %cl) {
	    my ($best) = sort { $b->score('ygob') <=> $a->score('ygob') } @{$array};
	    push @uniq, $best;
	}
    }

    map { $_->data($attr => 'delete') } ($self,  @{$args->{'-object'}});
    
    return @uniq;
}

=head2 diagnose()

    This method will look at an orthogroup and attempt to 
    identity the gene(s) responsible for at a low quality 
    score

=cut

sub diagnose {
  my $self = shift;
  my $args = {@_};
  
  my $anc = $args->{'-ancestor'};
  $args->{'-window'} = 7 unless exists $args->{'-window'};
  $args->{'-homology'} = 'YGOB' unless exists $args->{'-homology'};

  return () unless $self->_orthogroup;
  my @orfs = $self->_orthogroup;

  #####################################################
  # Baseline 
  #####################################################

  if ( 1==2 ) {
      my %hash = map { uc($_->organism) 
			   => [ $_->context(-self =>1,-distance=>$args->{'-window'}) ] } @orfs;
      my @keys = map { uc($_->organism) } @orfs;
      #print map { $_.":$#{$hash{$_}}"} @keys;
      
      my ($align,$og,$scr) = $self->dpalign(
	  # alignment 
	  -hash => \%hash, 
	  #-order => [ @keys ], #$key0 must be first 
	  -reference => $keys[0],
	  -global => 1,
	  # scoring 
	  -match => $args->{'-homology'},
	  -mismatch => -1e6,
	  -gap => 0,
	  -inversion => -100,
	  -trna => 100,
	  # 
	  -verbose => undef
	  );            
      print ">$scr";
      &_print_align( [keys %hash], $align );
  }

  #####################################################
  # Baseline 
  #####################################################

  for (my $i = 0; $i <= $#orfs; $i++ ) {
      my @set = (
	  ( $i == 0 ? () : @orfs[0..$i-1] ), 
	  ( $i == $#orfs ? () : @orfs[($i+1)..$#orfs])
	  );
      
      my %hash = map { uc($_->organism) 
			   => [ $_->context(-self =>1,-distance=>$args->{'-window'}) ] } @set;
      my @keys = map { uc($_->organism) } @set;
      #print map { $_.":$#{$hash{$_}}"} @keys;
      
      my ($align,$og,$scr) = $self->dpalign(
	  # alignment 
	  -hash => \%hash, 
	  #-order => [ @keys ], #$key0 must be first 
	  -reference => $keys[0],
	  -global => 1,
	  # scoring 
	  -match => $args->{'-homology'},
	  -mismatch => -1e6,
	  -gap => 0,
	  -inversion => -100,
	  -trna => 100,
	  # 
	  -verbose => undef
	  );         

      my ($simple_scr,$denom) = $self->_score_multiple_alignment_simple($align);

      print ">$align, $scr,$simple_scr,$denom";
      &_print_align( [keys %hash], $align, 50 );
  }

  return @issues;
}

sub _score_multiple_alignment_simple {
    my $self = shift;
    my $align = shift;
    
    my @keys = keys %{$align->[0]};
    my $unit = 1/scalar(@keys);

    my %counts;
    foreach my $k (@keys) {
	$counts{ $k } = scalar( grep { $self->isa( ref($_->{$k}) ) } @{$align} ); 
	#print $k, $counts{$k};
    }
    my ($min) = sort { $a <=> $b } values %counts;

    my $scr=0;
    foreach my $h ( @{$align} ) {
	my $lscore;
	map { $lscore+=$unit if $self->isa( ref($h->{$_}) ) } @keys;
	$scr++ if $lscore==1;
    }

    return($scr, $min);
}

=head2 ancestral_alignment()

    Similar to syntenic alignment but forces 
    alignment to the ancestral genome. 

=cut 

sub ancetral_alignment {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-ancestor'}=1 unless $args->{'-ancestor'};
    
    return $self->syntenic_alignment( %{$args} );
}

=head2 syntenic_alignment(-object => Orf, -ancestor => Anc_x.y|[Anc_x.y,Anc_w.z], 
    -homology => 'YGOB', -window => 7, -clean => 1, -score => 1)
    
    Perform a gene-level alignment of two genomeic regions centered on self/object.
    If -ancestor is true we explicitly align on the ancetral gene order.
    If -ancestor is an ancestral gene we further attempt to anchor the alignment on this.

    We return an alignment (array of hashes), score and scoring hash.
    
    dpalign() does the heavy lifting to create the alignment once the 
    datastructures have been built to represent the two regions 
    and the ancestral chromosome. 
    
=cut 

sub syntenic_alignment {
    my $self = shift;
    my $args = {@_};

    my $sister =  $args->{'-object'} || $args->{'-sister'};

    ######################################	
    # check arguments 
    ######################################	
    
    $args->{'-window'} = 7 unless exists $args->{'-window'};
    $args->{'-homology'} = 'YGOB' unless exists $args->{'-homology'};
    $args->{'-clean'} = 1 unless exists $args->{'-clean'};
    $args->{'-score'} = 'paralog' unless exists $args->{'-score'};
    
    $self->throw if $args->{'-clean'} && ! $args->{'-ancestor'};

    ######################################	
    # prepare variables 
    ######################################	

    my $key0 = 'ANC';  # used multiple places. 
    my $anc = $args->{'-ancestor'};

    ######################################	
    # Anc determines run mode --
    # Anc_x.y => Anchor alignment on this 
    # spoof => Anchor alignment on an Anc but first figure out what that is
    # undef => Do not anchor alignment on Anc chr 
    ######################################	

    # variables needed by dpalign 
    my %hash; # genes in ordered arrays by species 
    my @keys; # names of tracks(speices) to align in order of addition 

    if ( $anc ) { # explcicitly align against ancestral gene order 
	
	# we want to align to the ancestral gene order but do not know 
	# the gene to anchor on. figure it out and then proceed. 

	if ( $anc !~ /^Anc/i) {
	    my $left = $self->left;
	    $left = $left->left until ! $left || $left->ygob =~ /Anc_(\d+)\.(\d+)/o;
	    my ($l_chr, $l_index) = ($1,$2);
	    my $right = $self->right;
	    $right = $right->right until ! $right || $right->ygob =~ /Anc_(\d+)\.(\d+)/o;
	    my ($r_chr, $r_index) = ($1,$2);
	    return () unless $l_chr eq $r_chr;
	    my $index = sprintf("%d",($r_index + $l_index)/2);
	    $anc = 'Anc_'.$r_chr.'.'.$index;
	    return () unless $anc =~ /$HOMOLOGY{'YGOB'}/o;
	}   
    
	######################################	
	# create ancestral gene order using Anc and -window 
	######################################	
	
	my ($chr, $max_anc, $min_anc);
	if ( $anc =~ /^Anc_(\d+)\.(\d+)/ ) {
	    $chr = $1;
	    ($max_anc) = $2 + ( 2*$args->{'-window'} ); 
	    ($min_anc) = $2 - ( 2*$args->{'-window'} ); 
	} elsif ( ref($anc) eq 'ARRAY' ) {
	    $anc->[0] =~ /Anc_(\d+)\.(\d+)/ || $self->throw;
	    my ($chr1,$min) = ($1,$2);
	    $anc->[1] =~ /Anc_(\d+)\.(\d+)/ || $self->throw;
	    my ($chr2,$max) = ($1,$2);
	    $self->throw unless $chr1 eq $chr2;
	    ($chr, $max_anc, $min_anc) = ($chr, $min, $max);
	} else { $self->throw( $anc, ref($anc) ); }
	
	$min_anc = 1 unless $min_anc >= 1;
	$self->throw unless $min_anc <= $max_anc;	      
	
	######################################	
	# Create objects to represent ancestral gene order 
	######################################	
	
	my $genome = $self->up->up->clone;
	$genome->organism( $key0 );
	my $contig = ref($self->up)->new(SEQUENCE => 1e5 x 'ATGC', ID => 1);
	$genome->add(-object => $contig);
	#print $genome->organism, $min, $anc, $max;
	
	my @anc_gene_order;
	for my $pos ( $min_anc .. $max_anc ) {
	    my $homol = 'Anc_'.$chr.'.'.$pos;
	    my $orf = ref($self)
		->new(
		START => $pos*10,
		STOP => ($pos*10)+3,
		STRAND => 1,
		UP => undef
		);
	    $orf->data(  $args->{'-homology'} => $homol );
	    $orf->data( 'ANC_POS' => $pos );
	    $contig->add(-object => $orf);
	    push @anc_gene_order, $orf;
	}
	
	# place the ancestral data on the dpalign datastructures
	
	$hash{ $key0 } = \@anc_gene_order;
	push @keys, $key0;
    }
    
    ######################################	
    # define order and valid range. only Ancs in reference are 
    # admitted for scoring so we constrain ultimate range without 
    # having to edit the actual gene order array i_array. 
    ######################################	
    
    my @genes=($self,$sister);
    for my $i (0..1) {
	my @i_array = $genes[$i]->context(-distance => $args->{'-window'}, -self => 1);
	my @i_sort = ( ! $anc ? @i_array : grep {defined} # sort {$a <=> $b}
		       map { $_->ygob =~ /\.(\d+)/; $1 } grep { $_->ygob =~ /_$chr\./ } @i_array);    
	my @i_prune = &_prune_from_ends(\@i_sort, $args->{'-penalty'});
	
	my $i_rt; # i_rt is used to determine orientation 	
	map { $i_rt += ($i_prune[$_]>$i_prune[$_-1] ? 1 : -1) } 1..$#i_prune; # determine order on chr
	
	my $key = uc($genes[$i]->organism).($i+1);
	$hash{ $key } = [ $i_rt > 0 ? @i_array : reverse @i_array ];
	push @keys, $key;
    }

    ######################################	
    # align two regions (to the ancestor) 
    ###################################### 

    my $align = $self->dpalign(
	# alignment 
	-hash => \%hash, 
	-order => [ @keys ], #$key0 must be first 
	-reference => $keys[0],
	-global => 1,
	# scoring 
	-match =>  $args->{'-homology'} ,
	-mismatch => -1e6,
	-gap => 0,
	#-inversion => -100,
	#-trna => 100,
	-verbose => undef
	);

    ######################################	
    # strip the anc track and destroy objects. 
    ######################################	

    # all this destruction may not be necessary since fixing the 
    # problems in the DESTROY method. Mayne enough to go out of scope.. 

    if ( ! $args->{'-clean'} && ! $args->{'-score'} ) {
	delete $hash{ $key0 };
	undef( @anc_gene_order );
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @{ $align };
	$genome->DESTROY() if $genome;
	return ($align,undef);
    }
    
    ######################################	
    # trim back to only regions in YGOB ancestor 
    ######################################	
    
    my (@clean);
    if ( $args->{'-clean'} ) {
	foreach my $q ( grep { $_->{$key0} } @{$align} ) {	
	    foreach my $k ( grep {/\-/} keys %{$q} ) {
		my ($jnk,$k2) = split/\-/,$k;
		$q->{$k2}=$q->{$k};
		delete $q->{$k};
	    }		  
	    push @clean, $q;
	}
	shift(@clean) until ( ! @clean || $clean[0]->{$keys[0]} || $clean[0]->{$keys[1]} );
	pop(@clean) until ( ! @clean || $clean[-1]->{$keys[0]} || $clean[-1]->{$keys[1]} );
    } else { @clean = @{$align}; } 

    unless ( @clean ) {
	delete $hash{ $key0 };
	undef( @anc_gene_order );
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @{ $align };
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @clean;
	$genome->DESTROY() if $genome;
	return undef;
    }

    ######################################	
    # print output 
    ######################################	
    
    &_print_align( [keys %hash], \@clean, 2 ) if $args->{'-verbose'};
    
    unless ( $args->{'-score'} ) {
	delete $hash{ $key0 };
	undef( @anc_gene_order );
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @{ $align };
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @clean;
	$genome->DESTROY() if $genome;
	return (\@clean,undef);
    }

    ######################################	
    # compute scores 
    ######################################	
    
    my ($score,$hash) = 
	$self-> _score_syntenic_alignment(-align => \@clean, %{ $args } );

    my ($span,$ohno) = 
	$self->_testOhnoSpan(-align => \@clean, %{ $args } );

    ######################################	
    # final scrubbing ...
    ######################################	

    SCRUB : {
	delete $hash{ $key0 };
	undef( @anc_gene_order );
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @{ $align };
	map { delete $_->{ $key0 } } grep { exists $_->{ $key0 } } @clean;
	$genome->DESTROY() if $genome;
    }

    return (\@clean,$score,$hash,$span,$ohno);
}

sub _testOhnoSpan {
    my $self = shift @_;
    my $args = {@_};

    # 

    $self->throw unless ref( $args->{'-align'} ) eq 'ARRAY';
    $self->throw unless ref( $args->{'-align'}->[0] ) eq 'HASH';
    
    my $key0 = 'ANC';
    my ($span, $ohno)=(0,0);

    # 

    my @clean = @{$args->{'-align'}};
    my ($key1,$key2) = grep {!/$key0/} keys %{ $clean[0] };
    
    # 

    my $special_K;
    foreach my $k ( 0..$#clean ) {
	if ( $clean[$k]->{$key0} && 
	     $clean[$k]->{$key0}->data( $args->{'-homology'} ) eq $args->{'-ancestor'} ) {
	    $special_K = $k;
	    if ( $clean[$k]->{$key1} && $clean[$k]->{$key2} ) {
		$ohno=1;
	    }
	}
    }

    my @before = @clean[ 0 .. ($special_K-1) ];
    my @after =  @clean[ ($special_K+1) .. $#clean ] if $#clean >= $special_K+1;
    #&_print_align( [ keys %{ $clean[0] } ], \@before, 2 );
    #&_print_align( [ keys %{ $clean[0] } ], \@after, 2 );

    #

    my $total=0;
    foreach my $seg ( \@before, \@after ) {
	foreach my $quay ( $key1, $key2 ) {
	    $total++ if grep { $_->{$quay} } @{$seg}; 
	}
    }
    $span=1 if $total == 4;
    # print $total, $span;

    #

    return ($span, $ohno)
}

sub _score_syntenic_alignment {
    my $self = shift @_;
    my $args = {@_};
    my $align = $args->{'-align'};

    # check arguments 
    
    $self->throw unless ref($align) eq 'ARRAY';
    $self->throw unless ref($align->[0]) eq 'HASH';
    $self->throw unless $args->{'-score'} && $args->{'-score'} =~ /^[op]/i;

    $args->{'-penalty'} = 20  unless exists $args->{'-penalty'};

    my $key0 = 'ANC';
    
    # scoring matrix for alignment phase.  
    # >>>> should be possible to pass in a matrix. <<<<<
    
    my %scoring_matrix = (
	'OHNO' => 1,
	'CROSS' => ($args->{'-score'} =~ /p/i ? 1 : 0),
	'SAME' => 0,
	'GAP' => 0,
	'DELTA' => 0,
	'KC' => 0
	);

    # 
    
    my %count = (
	'OHNO' => -1, # I guess this OK ... 
	'KC' => 0,
	'GAP' => 0,
	'CROSS' => 0,
	'SAME' => 0,
	'DELTA' => 0,
	'ARRAY' => []
	);

    #######################################
    # count relationships in alignment that can later be scored.
    #######################################

    my @clean = @{$align};
    my ($key1,$key2) = grep {!/$key0/} keys %{ $clean[0] };
    
    my $oldk;
    foreach my $k ( 0..$#clean ) {
	my ($row,$old) = ($clean[$k], $clean[$oldk]);
	
	# first 2 are locus specific 
	# second 2 involve interleaving and intra locus scoring 
	
	if (! $row->{$key1} && ! $row->{$key2}) { 
	    $count{'GAP'}++;
	    next;
	} elsif ( $row->{$key1} && $row->{$key2} ) { 
	    $count{'OHNO'}++;
	    next;
	} elsif ( ($row->{$key1} && $old->{$key1}) || ($row->{$key2} && $old->{$key2}) ) {
	    $count{'SAME'}++;
	} elsif ( ($row->{$key1} && $old->{$key2}) || ($row->{$key2} && $old->{$key1}) ) {
	    $count{'CROSS'}++;
	} else {} # this is possible at the start of an alignment. 
	
	$oldk=$k; # only assigned of single-copy 
    }
    
    ################################################
    # distance approach.
    # This is similar to the old ohnologs() method. 
    ################################################	      

    my $oldk=0;
    foreach my $k ( 1..$#clean ) {
	my ($row,$old) = ($clean[$k], $clean[$oldk]);
	next unless ($row->{$key1} || $row->{$key2});
	$oldk=$k and next unless ($old->{$key1} || $old->{$key2} );
	# 
	$count{'KC'}++;

	if ( exists $old->{$key0} ) {
	    my $delta = abs( $old->{$key0}->data('ANC_POS') - $row->{$key0}->data('ANC_POS'))-1;  
	    my $pen = ($delta >= $args->{'-penalty'} ? 0 : $args->{'-penalty'} - $delta );
	    $count{'DELTA'} += $pen;
	    push @{$count{'ARRAY'}}, $pen;
	}

	$oldk=$k;
    }
    
    ################################################    
    # finalize scoring and create datastructure. 
    # we ignore delta (depracated), array (not relevant) and KC (by setting score to 0). 
    # In current implementation, only ohno and crossing contribute. 

    my $score=0;
    map { $score+=$count{$_}*$scoring_matrix{$_} } grep {!/delta|array/i} keys %scoring_matrix; 
    
    return (wantarray ? ($score, \%count) : $score);
}

sub dpalign {
    my $self = shift;
    $self->up->up->dpalign( @_ );
}

sub _print_align {
    my @keys = @{shift @_};
    my @clean = @{shift @_};
    my $lim = shift @_ || 0;
    
    foreach my $key (@keys) {
	my @print = ( 
	    scalar(@clean) >= $lim
	    ? (map { ($_->{$key} ? "X" : ".") } @clean)
	    : (map {s/Anc_//; $_} map { ($_->{$key} ? $_->{$key}->ygob : " . ") } @clean)
	    );
	print '>'.$key, join( (scalar(@clean) >= $lim ? "" : "\t"), @print);
    }
    
    return 1;
}

sub _prune_from_ends {
    my @x = @{shift @_};
    my $z = shift;

    my $i=0;
    $i++ until ( $i>=($#x-1) || $x[$i+1]-$x[$i] < $z);
    my $j=0;
    $j++ until ( $j>=($#x-1) || $x[$#x-$j]-$x[$#x-$j-1] < $z);
    
    #print $#x, $i,  $#x-$i-$j+1, "($j)";
    return splice(@x, $i, $#x-$i-$j+1)
}

=head2 distance_matrix()
=cut 

sub distance_matrix {
    my $self = shift;
    my $args = {@_};

    $args->{'-orfs'} = [ $self->_orthogroup ] unless exists $args->{'-orfs'};
    $args->{'-metric'} = 'dS' unless exists $args->{'-metric'};
    # 
    $args->{'-score'} = 0 unless exists $args->{'-score'};
    $args->{'-symmetrical'} = 0 unless exists $args->{'-symmetrical'};
    # 
    $args->{'-sort'} = 'ogid' unless exists $args->{'-sort'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};

    $self->throw unless defined $args->{'-orfs'} && ref($args->{'-orfs'}) =~ /ARRAY/;

    ########################################
    # prep vars 
    ########################################

    my $fh = STDOUT;
    my $meth = $args->{'-sort'};
    my @orfs= sort { $b->$meth <=> $a->$meth } 
    map {$self->throw($_) unless /::/; $_ } @{$args->{'-orfs'}};    
    
    ########################################
    # print pretty .. 
    ########################################

    if ( $args->{'-verbose'} >=1 ) {
	print {$fh} "\n>".$args->{'-metric'}, scalar(@orfs);
	print {$fh} (qw(Gene Anc Evalue HyperG LOSS Density OGID),
		     undef, ( map {$_->sn} @orfs ));
    }

    ########################################
    # loop + work 
    ########################################

    my %matrix;
    my @scores;
    for my $i ( 0..($#orfs-1) ) {
	my @row = ('-') x ($i+1);
	for my $j ( $i+1 .. $#orfs ) {
	    my $metric;
	    if ( $args->{'-metric'} =~ /^[pi]/i ) {
		my ($pid) = 100 - # convert to a distance 
		    $orfs[$i]->bl2seq( 
			-object => $orfs[$j],
			-identity => 1
		    );
		$metric = $pid;
	    } elsif ( $args->{'-debug'} ) {
		$metric=rand();
	    } else {
		my ($dN,$dS) = $orfs[$i]->paml( 
		    -object => [ $orfs[$j] ] ,
		    -method => 'yn00',
		    -length_ratio => 0.25,
		    -min_codons => 30
		    );
		return() unless defined $dN && defined $dS; # right behaviour?
		$metric = ($args->{'-metric'} =~ /^[dk]s/i ? $dS : $dN);
	    }
	    $matrix{ $orfs[$i]->name }{ $orfs[$j]->name } = sprintf("%.5f",$metric);
	    push @row, $matrix{ $orfs[$i]->name }{ $orfs[$j]->name };
	}
	
	########################################
	# print pretty 
	########################################
	
	if ( $args->{'-verbose'} >=1 ) {
	    print {$fh} 
	    ($orfs[$i]->sn, 
	     (map {/(\d+\.\d+)/; $1} $orfs[$i]->ygob), $orfs[$i]->logscore('ygob'), 
	     $orfs[$i]->hypergob,  $orfs[$i]->loss, $orfs[$i]->density,
	     $orfs[$i]->ogid.($orfs[$i]->ohnolog ? '*' : ''), '|', 
	     @row); 
	    print {$fh} 
	    ($orfs[$#orfs]->sn, 
	     (map {/(\d+\.\d+)/; $1} $orfs[$#orfs]->ygob), $orfs[$#orfs]->logscore('ygob'), 
	     $orfs[$#orfs]->hypergob, $orfs[$#orfs]->loss,  $orfs[$#orfs]->density,
	     $orfs[$#orfs]->ogid.($orfs[$#orfs]->ohnolog ? '*' : ''), 
	     '|', 
	     (('-') x scalar(@orfs))) if $i == ($#orfs-1);		    
	}
    }

    ######################################### 
    # Produce an overall synteny score from PW scores
    ######################################### 
    
    if ( $args->{'-score'} ) {
	my ($m,$sd) = &_calcMeanSD( @scores );
	return ( sprintf("%.1f", $m/( $sd>1 ? $sd : 1 )), $m, $sd );
    }
    
    ######################################### 
    # 
    ######################################### 

    if ( $args->{'-symmetrical'} || $args->{'-symmetric'} ) {
	foreach my $i ( keys %matrix ) {
	    foreach my $j ( keys %{$matrix{$i} } ) {
		$matrix{$j}{$i}=$matrix{$i}{$j};
	    }
	}
    }

    return \%matrix;
}

=head2 neighbour_joining( -matrix => \%sqr_dist_mat, -topology => 0|1,
    -verbose => 0|1 )

    Reconstruct a neighbour joining tree from a distance matrix. 
    Matrix must be a symmetric (square) distance matrix and can have any units. 
    
    The return value is a node object that is the center of an 
    unrooted trifurcating tree. 

    If -topology is true we force fit to the Saccharomyces sensu stricto topolgy. 

=cut

sub neighbour_joining {
    my $self = shift;
    my $args = {@_};
    
    $self->throw unless $args->{'-matrix'} && ref($args->{'-matrix'}) eq 'HASH'; #<<
    $args->{'-verbose'}=0 unless exists $args->{'-verbose'};
    $args->{'-topology'}=0 unless exists $args->{'-topology'};

    ###################################
    # Ensure matrix is well formed 
    ###################################

    my $D = $args->{'-matrix'};
    foreach my $key ( keys %{$D} ) {
	#map { print $key,$_,$D->{$key}->{$_} } keys %{$D->{$key}};
	foreach ( keys %{$D->{$key}} ) {
	    if ( $D->{$key}->{$_} <1e-10 || $D->{$_}->{$key} <1e-10 ) {
		$self->throw if $D->{$key}->{$_} > 1e-10 || $D->{$_}->{$key} > 1e-10;
		$D->{$key}->{$_} = 1e-10;
		$D->{$_}->{$key} = 1e-10;
	    } else {
		$self->throw( $D->{$key}->{$_}.':'.$D->{$_}->{$key} ) 
		    unless $D->{$key}->{$_} == $D->{$_}->{$key};
	    }
	}
    }
    
    ###################################
    # Basic calculations + early termination 
    ###################################

    my $number_taxa = scalar(keys %{$D});
    return() if $number_taxa <= 2; #<<

    ###################################
    # Set up basic vars
    ###################################

    my $fh = STDERR;
    
    # this is a hack... DEVIN 

    my @force_topology = (
	[ (grep {/Scer/} keys %{$D}), (grep {/Spar/} keys %{$D}) ],
	[ (grep {/Smik/} keys %{$D}), 'Node_1.1' ]
	);

    ###################################
    # Prepare a return data structure 
    ###################################

    my $root = Node->new(
	ROOT => 3, # trifurcating 
	NAME => 'Root'
	);
    foreach my $otu ( keys %{$D} ) {
	my $node = Node->new(
	    NAME => $otu,
	    ROOT => 0,
	    BRANCH_LENGTH => -$INFINITY,
	    BOOTSTRAP => 0
	    );
	$root->add( -object => $node );
    }

    ###################################
    # Iterate until all we attain a trifurcating root
    ###################################

    my $node_count=0;
    until ( $number_taxa <= 3 ) {
	
	#

	if ( $args->{'-verbose'} >= 2 ) {
	    print {$fh} "\nMATRIX\t$number_taxa\n", keys %{$D};
	    foreach my $i (keys %{$D} ) {
		print {$fh} $i, map { $D->{$i}->{$_} } (keys %{$D} );
	    }
	}
	
	###################################
	# calculate row and column sums 
	###################################
	
	my (%row_totals, %col_totals);
	foreach my $i ( keys %{$D} ) {
	    foreach my $j ( grep {$_ ne $i} keys %{$D} ) {		
		$row_totals{$i} += $D->{$i}->{$j};
		$col_totals{$j} += $D->{$i}->{$j};
	    }
	}
	
	# STEP 1 
	# Identify the OTUs to join by finding the minimum
	# of the Q matrix. Derive Q from D using Equation 1 
	# from  Gascuel and Steel 2006.
	
	my ($f,$g) = (
	    $args->{'-topology'} == 1 # force topology true  
	    ? @{ $force_topology[ $node_count ]}
	    : $self->_Q_matrix_minimum( $D, \%row_totals, \%col_totals, $number_taxa ) 
	    );
	
	# STEP 2
	# Join f and g through u that connects to the root. 
	# Define u by computing the distances to f,g. 
	# Gascuel and Steele, 2006. Equation 2.
	# Mol. Biol. Evol. 23(11):1997–2000. 2006
	
	my $dist_fg = $D->{ $f }->{ $g };
	my $dist_fu = sprintf("%.3f", 
	    ( 0.5 * $dist_fg ) + 
	    (( 0.5 * 1/($number_taxa-2) ) * 	     
	     ( ($row_totals{$f} -  $dist_fg ) - ($col_totals{$g} -  $dist_fg ) ))
	    );
	my $dist_gu = $dist_fg - $dist_fu;
	
	# pretty print .. 

	my $new_node = 'Node_'.(++$node_count).'.'.$node_count; # Sbay_contig.index	
	print {$fh} ">".$number_taxa, $f,$g,$dist_fg, $new_node,$dist_fu,$dist_gu 
	    if $args->{'-verbose'};
	
	# update datastructure 

	my ($f_node) = grep { $_->name eq $f } $root->stream;
	my ($g_node) = grep { $_->name eq $g } $root->stream;
	my $u_node = Node->new(
	    ROOT => 0,
	    NAME => $new_node,
	    BRANCH_LENGTH => ($dist_fg - ($dist_gu + $dist_fu))/2,
	    BOOTSTRAP => 0
	    );
	# 
	my $ancestor = $f_node->ancestor;
	$self->throw unless $f_node->ancestor eq $g_node->ancestor;
	$ancestor->add( -object => $u_node );
	#
	$ancestor->remove( -object => $f_node );
	$u_node->add( -object => $f_node );
	$f_node->branch_length( $dist_fu );
	# 
	$ancestor->remove( -object => $g_node );
	$u_node->add( -object => $g_node );
	$g_node->branch_length(  sprintf("%.5f", $dist_gu) );
	
	# STEP 3 
	# Update D matrix and iterate the process : 
	# 1. Remove f and g from the matrix
	# 2. Add u distances to the matrix 
	# Gascuel and Steele, 2006. Equation 3.
	# Mol. Biol. Evol. 23(11):1997–2000. 2006

	foreach my $k ( grep {!/^[$f|$g]$/} keys %{$D} ) {
	    $D->{$k}->{$new_node} = # calculate new distance 
		0.5 * ( $D->{$k}->{$f} - $dist_fu ) +
		0.5 * ( $D->{$k}->{$g} - $dist_gu );
		#( 0.5 * ( $D->{$k}->{$f} + $D->{$k}->{$g} - $dist_fg) );
	    $D->{$new_node}->{$k} = $D->{$k}->{$new_node}; # make symmetrical 
	}
	foreach my $remove ( $f, $g ) {
	    delete $D->{$remove};
	}
	foreach my $retain ( grep {$_ ne $new_node} keys %{$D} ) {
	    map { delete $D->{$retain}->{$_} } ($f,$g);
	}
	$number_taxa = scalar(keys %{$D});
    }

    ###################################
    # Final distances 
    ###################################
    # 1. A+B = X
    # 2. A+C = Y
    # 3. B+C = Z
    # 4. A + (Z-C) = X     [1 & 3]
    # 5. A - C = X - Z     [4]
    # 6. A = (X - Z + Y)/2 [5 & 2]
    # 7. B = X - A         [6 & 1]

    my @label = keys %{$D};
    my @dists = (3 x undef);
    
    $dist[ 0 ] = ( $D->{ $label[0] }->{ $label[1] } + 
		   $D->{ $label[0] }->{ $label[2] } - 
		   $D->{ $label[2] }->{ $label[1] } ) * .5;
    $dist[ 1 ] = $D->{ $label[0] }->{ $label[1] } - $dist[0];
    $dist[ 2 ] = $D->{ $label[0] }->{ $label[2] } - $dist[0];
    
    for my $i ( 0..$#label ) {
	my ($node) = grep { $_->name eq $label[$i] } $root->stream;
	$node->branch_length( sprintf("%.5f", $dist[ $i ]) );
    }

    ###################################
    # 
    ###################################

    $root->newick( -print => $args->{'-verbose'} );
    return $root;
}

sub _Q_matrix_minimum {
    my ($self,$D,$row,$col,$taxa) = @_;
    
    map { $self->throw unless ref($_) eq 'HASH' } ($D,$row,$col);
    $self->throw unless $taxa && $taxa >0 && $taxa <1e2;
    $self->throw if $taxa <= 3;

    ###################################
    # Initialize matrix
    ###################################

    my $Q;
    foreach my $key ( keys %{$D} ) {
	map { $Q->{$key}->{$_}=$INFINITY } keys %{$D};
    }

    ###################################
    # Equation 1
    # Mol. Biol. Evol. 23(11):1997–2000. 2006
    # http://mbe.oxfordjournals.org/content/23/11/1997.full.pdf+html
    # Gascuel and Steel
    ###################################

    foreach my $i ( keys %{$D} ) {
	foreach my $j ( grep {$_ ne $i} keys %{$D->{$i}} ) {
	    $Q->{$i}->{$j} = 
		(($taxa - 2)*$D->{$i}->{$j}) - 
		($row->{$i} - $D->{$i}->{$j}) - 
		($col->{$j} - $D->{$i}->{$j}) ; 
	}

	if ( 1==2 ) {
	    print 'Q',$i, map { 
		($Q->{$i}->{$_}==$INFINITY 
		 ? $INFINITY 
		 : sprintf("%.2f",$Q->{$i}->{$_}) ) } keys %{$D};
	}
    }

    # get the minimum of the Q matrix
    # --> nodes to join 

    my @labels=keys %{$Q};
    my ($min_i,$min_j)=@labels[0..1];
    for my $i (0..($#labels-1)) {
	for my $j ( ($i+1)..$#labels ) {
	    if ( $Q->{ $labels[$i] }->{ $labels[$j] } <
		 $Q->{ $labels[$min_i] }->{ $labels[$min_j] } ) {
		($min_i,$min_j)=($i,$j);
	    }
	}
    }
    
    return ($labels[$min_i],$labels[$min_j]);
}

=head2 sisters(-window => 10)

    Return candidate sister regions by looking at neighbouring
    ohnologs. Each candidate returned is centered on a 
    specific gene which is the return value.

=cut

sub sisters {
    my $self = shift;
    my $args = {@_};

    $args->{'-window'} = 10 unless exists $args->{'-window'};
    $args->{'-distance_test'} = 3 unless exists $args->{'-distance_test'};    
    
    # get all the genes with ohnologs left and right 
	  
    my %ohnos;
    foreach my $dir ( qw(left right) ) {
	my @context = $self->context(
	    -direction => $dir,
	    -distance => $args->{'-window'},
	    -self => -1
	    );
	push @{ $ohnos{$self->up->id}{$dir} }, grep {$_->ohnolog} @context;
    }

    # look at each CHR to see if there are ohnologs spanning the 
    # focal genes. If yes, and meet proximity criteria, then 
    # infer the location of the focal genes and store this. 
    # we will center windows on this for later analysis. 

    my @center_of_gravity;
    foreach my $chr (keys %ohnos) {
	foreach my $left ( @{$ohnos{$chr}{'left'}} ) {
	    my $d_left = $self->distance( -object => $left,-nogap => 1 )+1;		    
	    foreach my $right ( @{$ohnos{$chr}{'right'}} ) {

		# basic sanity checks -- ohnologs must be on same chr 
		
		next unless $left->ohnolog->up == $right->ohnolog->up;
		next unless $left->ygob && $right->ygob;
		next unless my @int = grep { $_->assign ne 'GAP' }  # always ordered low .. high id 
		$left->ohnolog->intervening( -object => $right->ohnolog );
		
		# require the linkage between the two ohno pairs 
		# goes through same anc genome 

		my ($sp1,$chr1,$index1) = &Annotation::Orf::_decompose_gene_name($left->ygob);
		my ($sp2,$chr2,$index2) = &Annotation::Orf::_decompose_gene_name($right->ygob);
		# print "SPAN", $left->sn, $left->ygob, $self->sn."*", $o->ygob."*", $right->sn, $right->ygob, 
		$left->up->id.'/'.$left->ohnolog->up->id;
		next unless $chr1==$chr2;
		my $d_anc = abs( $index1 - $index2 );
		next unless $d_anc > 0;
		
		# deploy basic criteria to ensure the distances makes sense 
		
		my $d_query = $left->distance( -object => $right,-nogap => 1 )+1;
		my $d_ohno = $left->ohnolog->distance( -object => $right->ohnolog,-nogap => 1 )+1;
		#next unless abs(log($d_query/$d_ohno)/log(2)) <= 2;
		map { next unless ( ($_/$d_query <= $args->{'-distance_test'}) || 
				    ($_ <= $args->{'-window'}*$args->{'-distance_test'}) ) } ($d_anc, $d_ohno);
		
		# Look for implied location of the sister gene and anc 
		
		my $cog_anc;
		if ( $index1 < $index2 ) {
		    $cog_anc = $index1 + int( ($d_left/$d_query) * $d_anc );
		    $self->throw unless $cog_anc >= $index1 && $cog_anc <= $index2; 
		} else {
		    $cog_anc =  $index1 - int( ($d_left/$d_query) * $d_anc );
		    $self->throw unless $cog_anc >= $index2 && $cog_anc <= $index1; 
		}
		
		# now the sister gene 
			
		my $cog_index = int( ($d_left/$d_query) * scalar(@int) );
		$self->throw( join(':',$#int, $cog_index) ) if $cog_index < 0 || $cog_index > $#int;
		my $cog_ohno = $int[$cog_index*( $left->ohnolog->index < $right->ohnolog->index ? 1 : -1)];
		
		# store options for this chromosome 
		
		push @center_of_gravity, [$cog_ohno, 'Anc_'.$chr1.'.'.$cog_anc]; 
	    }
	}
    }	
    
    return @center_of_gravity;
}

=head2 intervening(-object => $orf) 

    Return *all* the features that lie between the query pair as an array.
    
    NB: Calling ORFs are not included. 

=cut 

sub intervening {
    my ($self,$other) = grep { ! /\-object/ } @_;

    $self->throw() unless $self->isa(ref($other));
    #$self->throw() unless $self->up->up eq $other->up->up;

    return undef unless $self->up eq $other->up;   # DEVIN May 10
    $self->throw if $self->id eq $other->id; # should never happen 
    return() if $self eq $other;
    
    # sort by index as genes with same start can be arranged either 
    # of two ways by sort. we also do no tneed to second guess
    # the existing contig order -- should be dealt with by caller.
    
    my ($left, $right) = sort { $a->id <=> $b->id } ($self, $other);

    my @inter;
    my $next = $left->right;
    until ( $next eq $right ) {
	push @inter, $next;
	$next = $next->right;	
    }

    return @inter;
}

=head2 commoncodon

    Returns true if genes share one or more codons in same frame.
    
    Uses shortcuts but if fails, does exhaustive search. 

=cut 

sub commoncodon {
    my ($self,$other) = grep { ! /\-object/ } @_;

    ###########################################
    # basic QC
    ###########################################
    
    $self->throw unless $self->isa(ref($other));
    return undef unless $self->up eq $other->up; # desired behaviour?
    return undef unless $self->strand==$other->strand;
    map { $_->oliver(-creator => 1) and $_->throw($_->aa) 
	      unless $_->translatable(-fast => 1) } ($self, $other);

    ###########################################
    # fast answers? 
    ###########################################

    return 0 if $self->stop < $other->start || $self->start > $other->stop;
    return 1 if $self->sameframe(-object => $other, -last => 1, -strict => 1);
    return 1 if $self->exons == 1 && $other->exons == 1 && 
	$self->overlap(-object => $other) && $self->sameframe(-object => $other);

    ###########################################
    # the hard way -- check every codon 
    ###########################################
    
    my @b1 = map { $_->start..$_->stop } sort {$a->start <=> $b->start} $self->stream; # all nt in gene 
    my @b2 = map { $_->start..$_->stop } sort {$a->start <=> $b->start} $other->stream;
    $self->throw unless ( scalar(@b1)%3==0 && scalar(@b2)%3==0 );

    ###########################################
    # hash and test..
    ###########################################

    my %f1;
    for (my $i=0; $i <= $#b1; $i+=3) {
	$f1{$b1[$i]}++;
    }

    for (my $i=0; $i <= $#b2; $i+=3) {
	return 1 if exists $f1{$b2[$i]};
    }
    
    return 0;
}

=head2 sameframe( -strict => 0|1, -last => 0|1 ) 

    Returns true if genes (based on first codon) are in same 
    frame relative to contig end. We do not currently test that
    intervening codons are translatable so should be used with caution.

    With -last => 1 we also check the last codon. 

=cut 

# we use the contig terminii to test if 2 ORFs in same readng frame 

sub sameframe {
    my $self = shift; 
    my $args = {@_};

    $args->{'-last'} = undef unless exists $args->{'-last'};
    $args->{'-strict'} = undef unless exists $args->{'-strict'};

    my $other = $args->{'-object'};
    $self->throw unless $self->isa(ref($other));
    $self->throw unless $self->strand==$other->strand;
    
    return 1 if $self->start == $other->start;
    return 1 if $self->stop == $other->stop && $args->{'-last'};
    return 0 if $args->{'-strict'};

    my $f1 = $self->exons(-query => 'first')->_terminal_dist(-query => 'proximal', -method => 'start')%3;
    my $f2 = $other->exons(-query => 'first')->_terminal_dist(-query => 'proximal', -method => 'start')%3;
    return 1 if $f1==$f2;
    return 0 unless $args->{'-last'};
    
    my $f3 = $self->exons(-query => 'last')->_terminal_dist(-query => 'proximal', -method => 'stop')%3;
    my $f4 = $other->exons(-query => 'last')->_terminal_dist(-query => 'proximal', -method => 'stop')%3;
    return 1 if $f4==$f3;
    return 0;
}

=head2 _ohnolog_consistency

    Iterate through all the genes in the orthogroup, 
    lookup ohnologs and partition on the basis of their
    consistency in secondary orthogroup.
    
    -1 => no ohnolog 
    0 => ohnolog not in any OG
    >=1 => ohnolog in specified OG (key is OGID)

=cut 

sub _ohnolog_consistency {
    my $self = shift; 
    my $args = {@_};

    my %hash;
    foreach my $x ( $self->_orthogroup ) {
	my $key =  ($x->ohnolog ? $x->ohnolog->ogid || 0 : -1);
	push @{$hash{$key}},$x;
    }

    return \%hash;
}

=head2 orthogroup

    Return all orthologs in clade. Does NOT return self.

=cut 

sub orthogroup {
    my $self = shift;
#    print join("\t",%{$self})."\n" and $trigger->fail unless 
#	$self->up or (caller(1))[3] =~ /DESTROY/;
    return grep {/\w/} map { $self->$_ } map {lc($_)} grep {/\w/}  $self->up->up->bound;
}

=head2 summarize(@TEXT_TO_ADD)
    
    Print a summary of the OG.
    
=cut 

sub summarize {
    my $self = shift;
    return undef unless $self->ogid;

    print {STDOUT}
    $self->name, $self->identify, $self->data('KAKS'), $self->assign, 
    (map { sprintf("%.1f", $_) } (
	 _calcMeanSD( map {$_->length} $self->_orthogroup ),
	 _calcMeanSD( map {$_->binomialsynteny} $self->_orthogroup ),
	 _calcMeanSD( map {$_->logscore('ygob')} $self->_orthogroup )
     )),
	 @_;
    
    return $self;

    print {STDERR} 
    $self->name, $self->identify, $self->data('KAKS'), $self->assign, 
    (map { sprintf("%.1f", $_) } (
	 _calcMeanSD( map {$_->length} $self->_orthogroup ),
	 _calcMeanSD( map {$_->binomialsynteny} $self->_orthogroup ),
	 _calcMeanSD( map {$_->logscore('ygob')} $self->_orthogroup )
     )),	 
	 @_;
    
    return $self;
}

=head2 organism()

    Access genome->organism() value. 
    Baked into libraries but poor style.. 

=cut

sub organism {
    my $self = shift;
    return $self->up->up->organism;
}



=head2 intergenic(-object => $orf, -strand => 1, -seq => 0)
    
    Return the intergenic sequence between 2 genes (actually 
    any seqeunce between 2 features). By default 
    we examine the order of the genes on the contig and 
    return the top strand from the left-most gene to the rightmost. 

    If -strand -1 is supplied, we return the opposite strand. 

    The method is NOT responsible for checking the genes are 
    neighbours. It will accept any 2 features on the same contig. 

=cut

sub intergenic {
    my $self = shift;
    my $args = {@_};

    $args->{'-object'} = undef unless exists $args->{'-object'};
    $args->{'-strand'} = 1 unless exists $args->{'-strand'};
    $args->{'-seq'} = 0 unless exists $args->{'-seq'};

    $self->throw unless abs($args->{'-strand'})==1;
    $self->throw() unless exists $args->{'-object'} && $self->isa(ref($args->{'-object'}));

    # 

    $self->throw() unless $self->up == $args->{'-object'}->up;
    my ($left,$right) = sort { $a->start <=> $b->start } ( $self, $args->{'-object'} );
    return undef if $right->start < $left->stop;

    # 
    
    unless ( $args->{'-seq'} ) {
	$index = $left->index+.5;
	my $inter = Annotation::Orf
	    ->new(
	    START => $left->stop+1,
	    STOP => $right->start-1,
	    STRAND => $args->{'-strand'},
	    INDEX => ($left->index+.5),
	    UP => $left->up
	    ); 
	$inter->evidence('NONE');
	return $inter;
    }
    
    # 
    
    my $seq = substr($left->up->sequence, $left->stop, ($right->start-$left->stop)-1);
    $seq = join('',  reverse(map {tr/ATGC/TACG/; $_} split//, $seq) ) if $args->{'-strand'} == -1; 
    
    $self->warn( "Long intergenic: ".length($seq) ) 
	if ( $args->{'-verbose'} && length( $seq ) >= $args->{'-verbose'} );
    
    return $seq;
}

=head2 promoter(-length => bp) 
    
    Return upstream intergenic region. 

    This method is gap-aware so if there is a GAP feature
    (and corresponding runs of Ns) between 2 genes, they will
    be ignored. The 100Ns that is a signature of scaffold merges
    will be replaced by 10Ns. 

    If the length param is set, the method will take a fixed 
    amount and ignore all annotations. Without takes up to 
    the next annotation.     

=cut 

sub promoter {
    my $self = shift;
    my $args = {@_};

    return undef if $self->strand == 0; # DEVIN: unannotated silent refusal?? grr. 
    
    my ($start, $stop);
    if (exists $args->{'-length'}) {
	if ($self->strand == 1) {
	    if ($self->start-$args->{'-length'} > 0) {
		$start = $self->start-$args->{'-length'};
	    } else {$start=1;}
	    $stop = $self->start-1;
	} else {
	    $start = $self->stop+1;
	    if ($self->stop+$args->{'-length'} <= $self->up->length) {
		$stop = $self->stop+$args->{'-length'};
	    } else {$stop=$self->up->length;}
	}

    } else {
	if ($self->strand == -1) {
	    $start = $self->stop+1;
	    if (my $partner = $self->neighbour(-dir => 'right')) {
		$stop = $partner->start-1;
	    } else {$stop=$self->up->length;}
	} else {
	    $stop = $self->start-1;
	    if (my $partner = $self->neighbour(-dir => 'left')) {
		$start = $partner->stop+1;
	    } else {$start=1;}
	}
    }    
    return undef unless $stop >= $start;

    my $seq = $self->up->sequence;
    my $inter =  substr($seq, $start-1, ($stop - $start + 1));
    $inter =~ s/N{100}/NNNNNNNNNN/g;

    if ($self->strand == 1) {
	return($inter);
    } elsif ($self->strand == -1) {
	return(join('', reverse(map {tr/ATGCN/TACGN/; $_;} split(//, $inter))));
    }
}

=head2 neighbours(argeuments as below)

    Shortcut for: 

    ($left,$right) = $orf->neighbours(-direction => 'both', -other => args);

    Left or right (or both) may be undef. Test appropriately. 

=cut

sub neighbours { 
    my $self = shift;             
    my $args = {@_};
    $args->{'-direction'} = 'both';   
    return $self->neighbour(%{$args}); 
}

=head2 neighbour(-direction => 'left|right|both', -orthogroup => 0|-1|1, 
    -restrict => [ALL], -gap => -1, -variant => 0|-1|1 )
    
    Returns first feature of the required type in the requested direction.
    By default this is simply the first non-GAP feature. This is a modified version 
    of left/right (and should perhaps have been implemented as an override method).

    -orthogroup : require/exclude gene to be in an OG
    -variant :  require/exclude gene models that are variants of the caller
    -trna (or equiv) : require/exclude any given feature. Can only require 1.
    -restrict => [TRNA REAL]. compose complete set of allowed featurs. 

=cut

sub neighbour {
    my $self = shift;             
    my $args = {@_};
    
    $args->{'-direction'} = $args->{'-dir'} if exists $args->{'-dir'};
    $args->{'-direction'} = 'left' unless exists $args->{'-direction'};
    $args->{'-direction'} = 'both' unless $args->{'-direction'}; # undef -> both 
    $args->{'-direction'} = ['left', 'right'] if $args->{'-direction'} =~ /both/i;
    $args->{'-direction'} = [ $args->{'-direction'} ] unless ref($args->{'-direction'}) eq 'ARRAY';
    # these are orthogonal to -gap/-trna etc
    $args->{'-orthogroup'} = 0 unless exists $args->{'-orthogroup'};
    $args->{'-variant'} = 0 unless exists $args->{'-variant'};
    # 
    $args->{'-gap'} = -1 unless exists $args->{'-gap'};
    $args->{'-feature'} = -1 unless exists $args->{'-feature'};    

    ###############################################
    # define feature requirements: -1 exclude, 0 ambivalent, +1 require 
    ###############################################
    
    my %restrict = map { $_ => 0 } (
	exists $args->{'-restrict'} ? 
	@{ $args->{'-restrict'} } :
	(map { $_->{INFER} } values %EVIDENCE)
    );
    
    # count required genes

    my $test;
    foreach my $key ( map {s/^\-//; $_} keys %{$args} ) {
	$restrict{ uc($key) } = $args->{ '-'.$key } if exists $restrict{ uc($key) };
	$test++ if exists $restrict{ uc($key) } && $args->{ '-'.$key }==1;
    }

    # populate %restrict with allowed categories. 

    if ( $test > 1 ) { 
	$self->throw("May only require one feature type.");
    } elsif ( $test == 1) {
	map { delete $restrict{ $_ } unless $restrict{ $_ } == 1 } keys %restrict;
    } else {
	map { delete $restrict{ $_ } if $restrict{ $_ } == -1 } keys %restrict;
    }

    # one final wrinkle.. using -variant +/- 1 forces the use of commoncodon
    # it will blow up for tRNAs and other noncoding 
    # $self->throw if ( exists $restrict{'TRNA'} && abs( $args->{'-variant'} ) == 1);

    ###############################################
    # iterate until we meet requirements 
    ###############################################

    my @partner;
    foreach my $dir ( @{$args->{'-direction'}} ) {
	$partner = $self->$dir;
	$partner = $partner->$dir while (
	    $partner && ( 
		# not of correct type 
		! exists $restrict{$partner->assign}
		# incompatible orthogroup status 
		|| ( $args->{'-orthogroup'} == 1 && ! $partner->ogid )
		|| ( $args->{'-orthogroup'} == -1 && $partner->ogid )
		# incompatible variant status 
		|| ( $args->{'-variant'} == 1 &&            # must be variants 
		     ($self->coding && $partner->coding) && # both coding so commoncodon is valid 
		     ! $self->commoncodon($partner)         # not allowed 
		)
		|| ( $args->{'-variant'} == -1 && 
		     ($self->coding && $partner->coding) && 
		     $self->commoncodon($partner)
		) 
	    )
	    );
	push @partner, $partner;
    }

    ###############################################
    
    return ( $#{$args->{'-direction'}}==1 ? @partner : $partner[0] );
}

=head2 variants
    
    Return an array of genes that are variants of 
    the calling Orf. We use the principle that any 
    genes that share a common codon must be variants. 

    With -object argument, we perform a test and return TRUE/FALSE.
    The test differs from callign commoncodon directly in that 
    we do not throw an error for TRNAs etc.

=cut 

sub variants {
    my $self = shift;             
    my $args = {@_};

    $args->{'-variant'} = 0; # ambivalent 
    $args->{'-self'} = 0; # exclude. 
    $args->{'-distance'} = 20; # unless exists $args->{'-distance'};
    
    if ( $args->{'-object'} ) {
	$self->throw unless $self->isa( ref($args->{'-object'}) );
	$self->throw if $self eq $args->{'-object'};
	return 0 if $self->rank < -1 || $args->{'-object'}->rank < -1;
	return ( $self->commoncodon( $args->{'-object'} ) ? 1 : 0 );
    } 

    return grep { $self->commoncodon($_) } $self->context( %{$args} ); # 
}

=head2 context(-distance => 5, -direction => undef, 
    -self => 0|1, -neighbour => args)

    Returns an ordered array of up/down-stream genes.
    All hard work done by neighbour() and by default 
    we exclude variants, tRNAs, GAPs. 

=cut 

sub context {
    my $self = shift;
    my $args = {@_};

    $args->{'-direction'} = undef unless exists $args->{'-direction'};
    $args->{'-distance'} = int($args->{'-window'}/2) if $args->{'-window'};
    $args->{'-distance'} = 5 unless exists $args->{'-distance'};
    # 
    $args->{'-self'} = undef unless exists $args->{'-self'}; # include self on return array. 
    $args->{'-all'} = undef unless exists $args->{'-all'}; # accept all objects. unrestricted. 
    # 
    $args->{'-orthogroup'} = 0 unless exists $args->{'-orthogroup'};
    $args->{'-variant'} = -1 unless exists  $args->{'-variant'};
    #
    $args->{'-trna'} = -1 unless exists  $args->{'-trna'};
    $args->{'-gap'} = -1 unless exists  $args->{'-gap'};

    # 
    
    my %restrict = map { $_ => $args->{$_} } grep {!/distance|window|self|direction|all/} keys %{$args};
    map { $restrict{$_}=0 } keys %restrict if $args->{'-all'};

    #
    
    my @left;
    unless ( $args->{'-direction'} =~ /^r/i ) {
	my $left = $self->neighbour( -direction => 'left', %restrict );
	until ( ! $left || scalar(@left) == $args->{'-distance'} ) {
	    unshift @left, $left;
	    $left = $left->neighbour( -direction => 'left', %restrict );
	}
    }
    
    #

    my @right;
    unless ( $args->{'-direction'} =~ /^l/i ) {
	my $right = $self->neighbour( -direction => 'right', %restrict ); 
	until ( ! $right || scalar(@right) == $args->{'-distance'} ) {
	    push @right, $right;
	    $right = $right->neighbour( -direction => 'right', %restrict );
	}
    }
    
    unshift @right, $self if $args->{'-self'} == 1;
    return (@left, @right);
}

=head2 show( -window => 15 ) 

    Show a -window sized genomic region centered on locus 
    and print out the BLAST information for each: Hit and Evalue. 
    We also print the synteny delta for each gene from the caller.

=cut

sub show {
    my $self = shift;             
    my $args = {@_};
    
    $args->{'-window'} = 15 unless exists $args->{'-window'};
    $args->{'-species'} = [sort grep {!/SGD|YGOB|MITO|SPOM|KWAL|YLIP/} keys %HOMOLOGY] unless exists $args->{'-species'};

    my @genes = $self->context( -self => 1, -variants => -1, -trna => -1, %{$args} );
    my $fh = *STDOUT;

    print {$fh} '#id', undef, map { $_->id } @genes;
    print {$fh} '#length', undef, map { $_->length } @genes;
    print {$fh} '#assign', undef, map { $_->assign } @genes;
    print {$fh} undef, undef, map { ($_ eq $self ? '*****' : undef ) } @genes;

    foreach my $sp ( grep {defined} (qw(YGOB), @{$args->{'-species'}}) ) {
	# name
	print {$fh} $sp, undef, map { ($_->data($sp) ? join('.', (_decompose_gene_name($_->data($sp)))[1..2] ) : undef) } @genes;
	# score 
	print {$fh} undef, undef, map { $_->score($sp) } @genes;
	# syn
	print {$fh} undef, undef, map { ($_ eq $INFINITY ? '.' : $_ ) } 
	map { ( $self->data($sp) && $_->data($sp) ? _compute_gene_distance($self->data($sp), $_->data($sp)) : undef) } @genes;
	# len 
	#print {$fh} undef, undef, map {  $_->length } @genes;
    }
    
    return $self;
}

=head2 glyph(-print => KEY, -tag => undef)

    Draft version. Prints basic info and the chosen data('KEY').
    Copies to allele/Anna/current.svg for viewing. Use -tag to 
    append addtional identifier info to file name.

=cut

sub glyph {
    my $self = shift;             
    my $args = {@_};

    $args->{'-print'} = 'global' unless exists $args->{'-print'};
    $args->{'-distance'} = 5 unless exists $args->{'-distance'};
    $args->{'-tag'} = undef unless exists $args->{'-tag'};

    my $rowheight = 18;
    my $offset = 20;
    my $mod = 1;

    # 

    my @orfs = $self->context( -self => 1, -all => 1, -distance => $args->{'-distance'} );
    
    my ($left) = sort { $a <=> $b } map { $_->start} @orfs;
    my ($right) = sort { $b <=> $a } map { $_->stop} @orfs;

    my $width = (( $right - $left) + 1)/$mod + (2*$offset);
    my $height = ( (2+@orfs) * $rowheight) + (2*$offset);
    #$self->throw if $width > 50000;
    #$self->throw if $height > 1000;

    my $im = GD::SVG::Image->new($width, $height);
    my $gray = $im->colorAllocate(127,127,127);
    my $black = $im->colorAllocate(0,0,0);
    $im->rectangle($offset/2 ,$offset/2, $width-($offset/2), $height-($offset/2), $gray);
    my $font = gdLargeFont;
    
    my %col = (
	'REAL' => $im->colorAllocate(255,0,0),
	'SELF' => $im->colorAllocate(0,0,255),
	'TRNA' => $im->colorAllocate(0,255,0),
	'GAP' => $im->colorAllocate(0,0,0)
	);

    foreach my $i ( 0..$#orfs ) {	
	my $obj = $orfs[$i];
	my $h_off = ($i*$rowheight) + $offset;
	my $w_off = $offset - $left;       

	my $col = ( $obj eq $self ? $col{SELF} : $col{ $obj->assign} );
	$col = $col{ 'REAL' } unless $col;	

	my $oldex;
	foreach my $ex ( $obj->stream ) {
	    $im->filledRectangle(
		($ex->start + $w_off)/$mod, ($h_off)+2, 
		($ex->stop + $w_off)/$mod, ($h_off+$rowheight)-2, 
		$col 
		);
	    $im->line(
		($oldex->stop + $w_off)/$mod, ($h_off+($rowheight/2)), 
		($ex->start + $w_off)/$mod, ($h_off+($rowheight/2)), 
		$col 
		) if $oldex;
	    $oldex=$ex;
	}

	# prettify 

	my $poly = new GD::SVG::Polygon;
	if ( $obj->strand == 1 ) {
	    $poly->addPt( ($obj->stop + $w_off)/$mod, ($h_off-($rowheight/2)) );
	    $poly->addPt( ($obj->stop + $w_off)/$mod, ($h_off+($rowheight*(3/2))) );
	    $poly->addPt( ($obj->stop + $w_off + 2*$offset)/$mod, ($h_off+($rowheight/2)) );	    
	} elsif ($obj->strand == -1) {
	    $poly->addPt( ($obj->start + $w_off)/$mod, ($h_off-($rowheight/2)) );
	    $poly->addPt( ($obj->start + $w_off)/$mod, ($h_off+($rowheight*(3/2))) );
	    $poly->addPt( ($obj->start + $w_off - 2*$offset)/$mod, ($h_off+($rowheight/2)) );	    	    
	}
	$im->filledPolygon($poly,$col);

	# add info 
	my $text = join(" / ", ($obj->name, $obj->gene, $obj->ygob, $obj->score( $args->{'-print'} ), $obj->_creator) );
	$im->string($font, ($obj->start + $w_off)/$mod, $h_off, $text, $black);
    }

    # make a file in the working dir 

    my $fileid = $self->name.($args->{'-tag'} ? '_'.$args->{'-tag'} : undef).".svg";
    open(SVG, ">$fileid");
    print SVG $im->svg;
    close(SVG);

    # if on a non-local server, put in web dir so I can look at it 

    chomp(my $host = `hostname`);
    if ( $host !~ /local/ ) {
	my $svg = "/Library/WebServer/Documents/Anna/".$fileid;
	my $current = "/Library/WebServer/Documents/Anna/current.svg";
	system("cp $fileid $svg");
	system("cp $fileid $current");
    }

    return ( $host !~ /local/ ? $host."/Anna/$fileid" : undef);
}

=head2 combinations( -object => [orfs], -overlap => undef ) 

    Return all combinations of the supplied orfset
    subject to an overlap constraint. Overlap (redundancy) 
    is assessed either using commoncodon() or parametrically
    using the -overlap parameter. 

=cut 

sub combinations {
    my $self = shift;
    my $args = {@_};

    $self->throw unless $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->throw unless $self->isa(ref($_)) } @{$args->{'-object'}};

    $args->{'-overlap'} = undef unless exists $args->{'-overlap'};

    my @comb = ([$self]);
    foreach my $add ( @{$args->{'-object'}} ) {
	foreach my $existing ( @comb ) {
	    if ( $args->{'-overlap'} ) {
		map { next if $_->overlap( -object => $add ) > $args->{'-overlap'} } @{$existing};
	    } else {
		map { next if $_->commoncodon($add) } @{$existing};
	    }
	    push @comb, [ @{$existing}, $add ];
	}
    }

    return @comb;
}


=head2 kaks(-object => [$orf, $objects], -method => 'codeml|yn00', -range => [from..to])

    Run either codeml or yn00 on the calling ORF and 
    populate KA, KS and KAKS evidence parameters. 

    NB: Currently only codeml implemented, objects automatically drawn from
    orthogroup (-object ignored) and tree (required by codeml) taken
    from genome object (reference genome). 

=cut 

sub kaks {
    my $self = shift;             
    my $args = {@_};

    $args->{'-object'} = [ $self->orthogroup ] unless exists $args->{'-object'};
    $self->throw unless $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->throw unless $self->isa(ref($_)) } @{$args->{'-object'}};

    my %sp = map { $_->organism => 1} @{  $args->{'-object'} };     
    $self->throw unless scalar(keys %sp) == scalar($self->up->up->bound);
    
    my ($ka, $ks) =  $self
	->paml(
	-method => 'codeml', 
	-object => [  @{$args->{'-object'}} ], 
	-range => $args->{'-range'} # use a subset of the alignment 
	);    
    
    return $self unless $ka >= 0 && $ka < 100;
    return $self unless $ks >= 0 && $ks < 100;

    foreach my $x ( $self, @{$args->{'-object'}} ) {
	$x->data('KA' => $ka);
	$x->data('KS' => $ks);
	if ( defined $ks && $ks > 0 ) {
	    $x->data('KAKS' => sprintf("%.3f", $ka/$ks) );	
	    $x->data('KAKS2' => $x->data('KAKS') ); # both of these are used for testing 
	}
	$x->evaluate;
    }
    $self->output if $args->{'-verbose'};

    return $self;
}

#########################################
#########################################

# manipulation methods 

=head2 optimise(-feature => 'terminii|junctions|both', 
    -M_min => 90, -M_max => 300, -M_frac => .5,
    -gap_min => 15, -gap_max => 120)

    This method should be replaced by one that uses alignments. 
    We currerntly try to maximise the length of genes (with some 
    allowances for structure) and this can result in loss of info
    that was previously obtained by BLASTing or other sequence 
    methods. 

    3 run modes. terminii and/or junctions: 	
    - terminii will "optimise" the terminii of each exon
    - junctions will ensure that all exon-exon junctions are healthy
    - defualt is both

=cut 

sub optimise {
    my $self = shift;             
    my $args = {@_};
    # $self->throw;

    # uses labels + (regretably) goto.  
    # 3 main labels: TERM, JUNC, FINISH 
    # TERM uses 3 inner labels: EXON, START, STOP 
    
    ################################################
    # do not inherit somebody elses problems.

    return $self if $self->rank < 0; # RNA (-Inf) or pseudo (-1)	

    $self->index; # we may NOT call ->index again until method is done. see note below. 
    $self->throw("ORF does not validate\n".$self->debug)
	unless $self->translatable;

    ################################################
    # 2 run modes. terminii and/or junctions.	

    $args->{'-feature'} = 'both' unless exists $args->{'-feature'};    
    # terminii will optimise the terminii of each exon
    $args->{'-M_min'} = 120 unless exists $args->{'-M_min'};  # if exon <, do not seek M       
    $args->{'-M_frac'} = .5 unless exists $args->{'-M_frac'}; # if > exon lost, do not accept M 
    $args->{'-M_max'} = 300 unless exists $args->{'-M_max'};  # if > exon lost, do not accept M 
    # junctions will ensure that all exon-exon junctions are healthy
    $args->{'-gap_min'} = 15 unless exists $args->{'-gap_min'};
    $args->{'-gap_max'} = 120 unless exists $args->{'-gap_max'};
    

    ################################################
    # TERMINII: expand and contract exons to maximise sequence coverage
    ################################################       

  TERM: # there are 3 inner labels EXON|START|STOP 
    goto JUNC if $args->{'-feature'} =~ /junc/;

    my $first = $self->exons(-query => 'first');
    my $last = $self->exons(-query => 'last');

  EXON: foreach my $exon ($self->stream) { # every exon in turn.... 
      next if $exon->_creator =~ /reconstitute|exonify/ ||
	  $exon->_depracate; # this is depracated ... 

      ################################################ 	    
      # expand exon start unless protected or at contig end
    START:
      {
	  goto STOP if $exon->intron(-direction => 'left') > 0 ||
	      $exon->_terminal_dist(-method => 'start') < 3 ||
	      $args->{'-feature'} =~ /stop/i; # ignore 5' 

	  $exon->morph(
	      -method => 'expand',
	      -terminus => '5 prime',
	      -stop => '\*'
	      );

	  $exon->start(-adjust => +3, -R => 1) #  remove leading * 
	      unless ($self->translatable && $exon->start > 0); 

	  # $exon->start == 0 should not be possible. 
	  # there seems to be an off-by-one bug in morph.
	  
	  # Special treatment for first exons. 
	  # In general, we roll back to find an M unless
	  # 1. exon looks like terminated by contigs end 
	  # 2. exon very short-- assume truncated by error. 
	  # 3. we lose too much gene body. 
	  
	  goto STOP unless ($exon eq $first);
	  goto STOP if $exon->_terminal_dist(-method => 'start') < 3
	      || $exon->length < $args->{'-M_min'};
	  
	  my $M = $exon->start(-R => 1);
	  my $L = $exon->length;
	  
	  $exon->morph(
	      -method => 'contract',
	      -terminus => '5 prime',
	      -stop => 'M'
	      );

	  # if fractional or absolute loss of length too great, revert to old STOP 

	  if ( (($L - $exon->length)/$L) > $args->{'-M_frac'} || # what price M?
	      ($L - $exon->length) > $args->{'-M_max'} ) {
	      $exon->start(-R => 1, -new => $M);
	  }
      }
     
      ################################################ 	    
      # expand exon stop unless protected or at contig end	          
    STOP:
      {
	  next EXON if $exon->intron(-direction => 'right') > 0 ||
	      $exon->_terminal_dist(-method => 'stop') < 3;
	  
	  if ( $args->{'-feature'} =~ /stop/i ) {
	      next EXON unless $exon eq $last;
	  }

	  $exon->morph(
	      -method => 'expand',
	      -terminus => '3 prime',
	      -stop => '\*'
	      );		

	  $exon->stop(-adjust => -3, -R => 1) # remove trailing * 
	      unless ($self->translatable || $exon eq $last);	      
      }
  }

    ################################################
    # NOTE
    ################################################

    # Exon order can change due to terminus optimization and calling of ->index 
    # This can cause: 
    # 1. (previously) internal exons to end with * 
    # 2. can cause previoulsy OK intron score pairs to disordered + illegal
    # =>      1       85622   85647   26      1       1e+100  6.91
    # =>      3       85713   86394   682     1       6.91    0
    # =>      0       86563   87093   531     1       0       0
    # becomes: 
    # =>      1       85622   85644   23      1       1e+100  6.91
    # =>      0       85696   87090   1395    1       0       0
    # =>      3       85713   87093   1381    1       6.91    0

    #
    # ->index is not allowed until the end of the method. just priot to exonerate. 
    # 

    ################################################
    # JUNCTIONS 
    ################################################	
    
  JUNC: 
    $self->evaluate(-structure => 1,-validate => 1);
    goto FINISH if $args->{'-feature'} =~ /ter|stop/ || $self->exons == 1;

    $self->throw("-gap/-overlap params not set ".%{$args}) 
	unless defined $args->{'-gap_min'} && defined $args->{'-gap_max'};	

    my @add; # these are novel exons recovered by exonerate. 
             # they are added at the end becuase of the ban on calling 
             # index within this method. 
    
  PAIR: foreach my $exon ($self->stream) {
      next unless my $lexon = $exon->left;     

      ################################################ 	    
      # 1. check exon boundary scores 
      # $INF -> feature. boundaries are locked.  
      # 1e-5 < intron < $INF -> true intron. boundaries are locked.   
      # 1e-5 -> artificially paired intron. hacked from ->merge. boundaries are locked.
      # 0 -> moveable boundary 
      
      my $ir = $exon->intron(-direction => 'left');
      my $il = $lexon->intron(-direction => 'right');
      
      if (($il > 0 && $il < $INFINITY) && 
	  ($ir > 0 && $ir < $INFINITY)) {
	  next PAIR; # true introns 
      } elsif ($il == $INFINITY && $ir == $INFINITY) {
	  next PAIR; # locked. intervening feature. probably gap.
      } elsif ($il == 0 && $ir == 0) {
	  # proceed ....
      } else {
	  $self->output();
	  $self->throw("Hybrid introns not allowed: $il, $ir");
      }
      
      ################################################ 	    
      # 2) check exon boundary sequences 
      # gap=0 : merge 
      # gap>0 : either fill in region between HSPs or ignore 
      # gap<0 : complicated. frame preservation main criterion. 
      
      my $gap = (
	  $self->strand == 1 
	  ? ( $exon->start(-R => 1) - $lexon->stop(-R => 1) )
	  : ( $lexon->stop(-R => 1) - $exon->start(-R => 1) )
	  ) -1;
      
      #$lexon->output;
      #$exon->output;
      #print $gap;

      if ($gap == 0) { # PERFECT: just merge 
	  #$exon->_merge_data_structures(-object => $lexon);
	  $exon->start(-new => $lexon->start(-R=>1), -R => 1);
	  $self->remove(-object => $lexon, -force => 1);
	  
      } elsif ($gap > 0) { # GAP: ignore or recover  	  
	  # This is very questionable. we are randomly adding seqeunce. ditch it? 
	  # 
	  # ... if so small it does not matter or so big it cannot be right 
	  # just ignore it. otherwise try recover the internal sequence. 
	  # should we add a method to check for tandems and split? 
	  
	  if ($gap > $args->{'-gap_min'} && $gap < $args->{'-gap_max'}) {
	      my $intron = $exon->intron(-direction => 'left', -object=>1);
	      @add = $intron->exonify(-minimum => $args->{'-gap_min'}, -overlap => 0);
	      $intron->DESTROY;
	      # $self->index; # we use @add, to add these in later. 
	  }
	  
	  } elsif ($gap < 0) { # OVERLAP 

	      # the internal junctions have score 0 (checked above) and can be moved at will. 
	      # we have the option to make a single exon if possible. 
	      # the external junctions we may not have freedom with however. 
	      # we therefore begin by examining these.
	      # 1. we choose either the coord required by frame + score constraints
	      # or the coords which maximize the exon span. 
	      # 2. we then go back and decide what to do with the internal junction. 
	      # This usually amounts to trimming back any overlapping regions. 

	      # It is also possible an exon overlaps >= 1 other.   
	      # In this case we need to iterate. 

	      until ( ! $lexon || $exon->overlap(-object => $lexon) <= 0 ) { 

		  my ($start_exon, $start_other);
		  if ( $lexon->frame(-first => 1) != $exon->frame(-first => 1) ||
		       $lexon->intron(-direction => 'left') > 0 ) {
		      # we are obliged to use the start position of the upstream exon 
		      ($start_exon, $start_other) = ($lexon, $exon);
		  } else {
		      # we can use either start but if using one from downsteam exon,
		      ($start_exon, $start_other) = sort { $a->_terminal_dist(-method => 'start') 
							       <=> $b->_terminal_dist(-method => 'start') } ($exon, $lexon);
		  }
		  
		  my ($stop_exon, $stop_other);
		  if ( $lexon->frame(-last => 1) != $exon->frame(-last => 1) ||
		       $exon->intron(-direction => 'right') > 0 ) {
		      # we are obliged to use the stop position of the downstream exon 
		      ($stop_exon, $stop_other) = ($exon, $lexon);
		  } else {
		      ($stop_exon, $stop_other) = sort { $a->_terminal_dist(-method => 'stop') 
							     <=> $b->_terminal_dist(-method => 'stop') } ($exon, $lexon);
		  }

		  # this can only happen due to length not structure concerns. 
		  # it is probably a sign that somethign has gone wrong. 
		  # (when index is called results in swapping of exon order).
		  # simplest solution is just to restore exon order and 
		  # trim back sections. 

		  if ( $start_exon eq $exon && $stop_exon eq $lexon ) { 
		      ($start_exon, $start_other) = ($lexon, $exon);
		      ($stop_exon, $stop_other) = ($exon, $lexon);
		  } 

		  # 

		  if ( $start_exon eq $stop_exon ) {
		      $self->remove(-object => $start_other);
		      
		  } elsif ( # we know they overlap so iff they are in frame (relative to contig) can merge 
			    ($start_exon->_terminal_dist(-method =>'start')%3 == 
			     $stop_exon->_terminal_dist(-method =>'start')%3) && 
			    ($start_exon->_terminal_dist(-method =>'stop')%3 == 
			     $stop_exon->_terminal_dist(-method =>'stop')%3) ) {	

		      $start_exon->stop(-new => $stop_exon->stop(-R => 1), -R => 1);
		      $start_exon->intron(-direction => 'right', -new => $stop_exon->intron(-direction => 'right') );
		      $self->remove(-object => $start_other);
		      
		  } elsif ( # not in frame or upstream not permissive. must maintain 'intron'. 
		      $stop_exon->_terminal_dist(-method => 'stop') < # stop_exon downstream 
			   $stop_other->_terminal_dist(-method => 'stop') ) {
		      if ( $start_exon->_terminal_dist(-method => 'start') < # start_exon upstream 
			   $start_other->_terminal_dist(-method => 'start') ) {
			  # A--------------
			  #         --------------B 
			  #         >>>>>>>		      
		      } else {
			  #       A--------------
			  #  --------------------------B 
			  #  >>>>>>>>>>>>>>>>>>>>>         
		      }		  
		      until ( $stop_exon->_terminal_dist(-query => 'distal', -method => 'start') < 
			      $start_exon->_terminal_dist(-query => 'proximal', -method => 'stop') ) {
			  $stop_exon->start(-R => 1, -adjust => +3); 
		      }
		      
		  } else {
		      if ($stop_exon->_terminal_dist(-method => 'stop') < # stop_exon downstream 
			  $stop_other->_terminal_dist(-method => 'stop') ) {
			  #         <<<<<<<<<<<<<<<<<<<<<<<<
			  # A-------------------------------
			  #         --------------B 
		      } else {
			  #                   <<<<<<<<<<<<<<<<<<<<<<<<
			  #           A-------------------------------
			  #    ---------------------B 
			  #    >>>>>>>>>>>>>>>>
			  
			  my $mid = int( ($start_exon->start(-R => 1) + $stop_exon->stop(-R => 1) )/2 );
			  if ( $self->strand == 1 ) {
			      $stop_exon->start(-R => 1, -adjust => 3) until ( $stop_exon->start(-R => 1) > $mid );
			  } else {
			      $stop_exon->start(-R => 1, -adjust => 3) until ( $stop_exon->start(-R => 1) < $mid );
			  }
		      }

		      until ( $start_exon->_terminal_dist(-query => 'distal', -method => 'stop') < 
			      $stop_exon->_terminal_dist(-query => 'proximal', -method => 'start') ) {
			  $start_exon->stop(-R => 1, -adjust => -3); 
		      }
		  }  

      		  $lexon = $exon->left; # exons can overlap >= 1 other. iterate. 
	      } # Iterate multiple exons 

	  } # OVERLAP 
  } # each exon pair
    
    # unfinished business from the GAP case above 

    map { $self->add(-object => $_) } @add;

    ################################################ 	        
    # finsih up 
    ################################################ 	    

  FINISH:

    $self->index;
    $self->evaluate(
	-structure => 1,
	-validate => 1
	); 

    return $self;
}

=head2 merge(-object => $orf, -reference => homolog, -extend => 500)

    Merge two ORF objects into one. Return caller on success 
    and destroy other. Else return undef.

    Rather than deal with how to integrate the existing models, 
    we generate a range of new models and choose the best one. 
    Since homology is invariably the basis for merging two
    genes, we use this (-reference) to pick best new model. It can 
    be either an object with a sequence(-molecule => 'aa') method
    or it can be a attribute/gene/file name.

    The models to choose from are generated by exonerate and 
    search_contig in an area around the merged gene (-extend). 

    We use the update method to regenerate what evidence we can. 

=cut

sub merge {
    my $self = shift;             
    my $args = {@_};
    
    #$args->{'-reference'} = 'SGD' unless exists $args->{'-reference'};   

    my ($other,$ref) = ($args->{'-object'},$args->{'-reference'});

    # same type, same contig, same strand, neighbours (checked below)
    
    $self->throw unless $self->isa(ref($other));	
    $self->throw unless $other->up eq $self->up;
    $self->throw unless $self->strand == $other->strand;
    map { $_->throw unless $_->translatable } ($self, $other);

    # basic QC 

    map { return $self unless $_->assign eq 'GAP' } $self->intervening(-object => $other);
    my ($ortho, $two) = grep {defined} map {$_->orthogroup || $_->ogid} ($self, $other);
    $self->throw and return undef if $two;
    my ($ohno, $ohtwo) = grep {defined} map {$_->ohnolog} ($self, $other);
    $self->warn("Two Ohnologs!", $other) and return undef if $ohtwo;

    # 

    return undef unless my $hom = $self->homolog( -object => [$other] );

    # we want to compare a few possibilities 
    # 1. just sticking two halves together and not thinking 
    # 2. take 2 parts and "optimise" 
    # 3. repredict the whole region using exonerate 
    # 4. reannotate region using a lower intron threshold 

    # make a proxy merge: this can be really messed up. vorsicht! 
    # this does not have to be exact since we only look in region around model
    # it is not used to guide choice. 

    my $merge = $self->clone;
    my $temp = $other->clone;
    map { $_->transfer(-from => $temp, -to => $merge, -warn => 0) } $temp->stream;
    $temp->DESTROY;
    map { $_->stop(-R => 1, -adjust => -3) if $_->length > 3  } $merge->stream; # remove STOP codons
    
    # optimize around the merged model 
    
    $merge->update( -evaluate => 0 ); # call update on $merge 
    my $success =     # calls update on $success 
	$merge->reoptimise(-reference => $hom, -verbose => 0); # can return  merge proposal 
    return undef unless $success && $success->translatable;
    
    # compare to self? pointless? 
    # my ($best,$delta) = $self->choose(-object => [$merge], -reference => $ref);
    # not implemented. or desirable? 
    
    map { $self->remove(-object => $_, -warn => 0); $_->DESTROY; } $self->stream;
    map { $_->transfer(-from => $merge, -to => $self, -warn => 0) } $merge->stream;
    print $self->aa and $self->throw unless $self->translatable;
    $merge->DESTROY;

    # update does not handle the following: 
    # EXTRA is depracated and no longer set. 
    # $self->data('EXTRA' => join(';', (grep {defined} $self->data('EXTRA'),$other->data('EXTRA')))  );

    if ( $ohno && ($ohno ne $self) ) { # ohno == other 
	$other->ohnolog($self) unless $other->ohnolog eq $self;
	#$self->ohnolog($other) unless $self->ohnolog eq $other;
    }

    if ( $ortho && ($ortho ne $self) && $ortho->orthogroup ) { # if non ref memeber cannot do transfers. 
	map { $self->$_( $other->$_ ); } map { lc($_) } $self->up->up->bound;
	my $test = ($self->up->up->bound)[0];
	$self->throw unless $self->$test->id == $other->$test->id;
    }

    $other->DESTROY;
    $self->_creator( (caller(0))[3] );
    $self->update();
    #$self->output(-creator => 1);

    return $self;
}


=head2 overlap(-object => $orf, -compare => 'coords|seq|homolog', 
    -evalue => 1e-10, -contig => 1)

    Compare caller and object in terms of either their genomics coords, alignability
    (to each other; bl2seq) or their alignment to a shared homoog. In the last case
    we return a 0-1 scale indicating what fraction of the aligned positions in the 
    homolog covered twice vs once. Consequently, 0 is usually the desired outcome 
    and return values must be carefully chekced -- undef indicates failure. 
    
    -evalue : specify alignment cutoff for direct sequence comparison 
    -object : may be a fasta file if compare => 'seq' 
    -contig : set to zero to skip 'same contig' test

=cut

sub overlap {
    my $self = shift;             
    my $args = {@_};	
    
    $args->{'-contig'} = 1 unless exists $args->{'-contig'};
    $args->{'-compare'} = 'coords' unless exists $args->{'-compare'};
    $args->{'-evalue'} = 1e-10 unless exists $args->{'-evalue'};
    my $other = $args->{'-object'};

    $self->throw unless $self->isa(ref($other)) || -e $other;
    return undef unless ! $args->{'-contig'} || $self->up eq $other->up; # correct behaviour?

    ######################

    goto HOMOLOG unless $args->{'-compare'} =~ /coords|seq/;

    # get shorty
    
    my $objlen;
    if ( -e $other ) {
	chomp($objlen = `grep -v '>' $other | sed s/[[:space:]]//g | sed s/[[:space:]]+//g | wc -c `);
	$objlen =~ s/\s+//; # required for score normalization 
    } else {
	$objlen = $other->length;
    }        
    my $bar = ( $self->length > $objlen ? $objlen : $self->length );
    
    # full sequence comparison or just coords

    my $overlap;
    if ($args->{'-compare'} =~ /coord/) {
	foreach my $x ($self->stream) {
	    foreach my $y ($other->stream) {
		$overlap += $x->overlap(-object => $y);
	    }
	}
    } elsif ($args->{'-compare'} =~ /seq/) {

	# change bl2seq to exonerate? 
	$overlap = $self->bl2seq(
	    -object => $other, 
	    -evalue => $args->{'-evalue'}
	    );
    }
    
    return 0 unless $bar;
    $overlap /= $bar;
    return $overlap;

    ######################

  HOMOLOG: 

    # get a homolog unless one supplied 

    my ($file);
    if ( $args->{'-compare'} =~ /^hom/i ) {
	return unless my @r = $self->homolog(-object => $other);
	$file = $r[-1];
    } else {
	($file) = $self->fetch(-id => $args->{'-compare'}, -molecule => 'aa', -file => 1);
    }
    $self->throw unless -e $file;
    
    # ensure queries both overlap the homolog 

    if ( defined $args->{'-evalue'} ) {
	my $o1 = $self->overlap(-object => $file, -compare => 'seq', -evalue => $args->{'-evalue'} );
	my $o2 = $other->overlap(-object => $file, -compare => 'seq', -evalue => $args->{'-evalue'});
	return undef unless $o1 >= .5 && $o2 >= .5; # weird hard coding .. 
    }

    # 
    
    my @res = $self->exonerate2(
	-object => $other, 
	-homolog => $file, 
	-model => 'local'
	#,-verbose => 2
	);

    return ( @res && exists $res[2] ? $res[2]->[0] : undef);
}     

=head2 adjust(adjustment)

    Change coordinates by sliding gene model along sequence.

=cut 

sub adjust {
    my $self = shift;
    map { $_->_adjust(@_) } $self->stream;
    return $self;
}

=head2 locus( -extend => 5000 ) 

    Return the gene plus -extend up and downstream 
    as a _contig_ object. We also return the offset
    relative to the original coord scheme. It is a 
    wrapper on _pseudocontig().

=cut 

sub locus {
    my $self = shift;
    my $args = {@_}; 
   
    $args->{'-extend'} = 5000 unless exists $args->{'-extend'};

    $self->throw unless wantarray;

    return $self->up->_pseudocontig( -object => $self , -extend => $args->{'-extend'} );
}

=head2 reoptimise(-orfmin => undef, -hmm => ygob(), -protein => homology(),
    -reference => obj/gene/file, -extend => 500, -safe => 0, -update => 1, 
    -atg => 0|i)

    Repredict gene using hmology and/or more lenient intron params.
    Use homology to select among the pool of new models and the original. 
    ORF-based prediction is _off_ by default. 

    By default returns the reoptimised version of self only. With -safe 
    returns the new version and does not alter original. Both are attached to contig array.
    Setting -update => 0 disables default behaviour of updating evidence hash.

    NB: can return self unaltered if this is best match to reference. Returns undef
    on failure (no translatable models found). 

=cut 

sub reoptimise {
    my $self = shift;
    my $args = {@_};
    
    # basic prediction 
    $args->{'-extend'} = 5000 unless exists $args->{'-extend'};
    $args->{'-orfmin'} = undef unless exists $args->{'-orfmin'};
    # homology 
    $args->{'-hmm'} = ( $self->pillar('YGOB') )[0] unless exists $args->{'-hmm'};  # 
    $args->{'-hmm'} = $self->ygob unless defined $args->{'-hmm'};
    $args->{'-protein'} = $self->homolog(-fast => 1) unless exists $args->{'-protein'}; # also uses pillar()
    #$args->{'-reference'} = undef unless exists $args->{'-reference'}; 

    # post processing 
    $args->{'-update'} = 1 unless exists $args->{'-update'};
    $args->{'-atg'} = 300 unless exists $args->{'-atg'};
    # generic 
    $args->{'-safe'} = 0 unless exists $args->{'-safe'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    
    $self->throw unless $args->{'-reference'};

    ##################################################################
    # Repredict and generate a panel of alt mdoels 
    ##################################################################

    # 0. get relevant region 

    my ($locus,$adjust) = 
	$self->locus( -extend => $args->{'-extend'} );

    # 1. reannotate based on open reading frames 

    my @anno = 
	$locus->search( -min => $args->{'-orfmin'}, -hmm => undef ) 
	unless ! $args->{'-orfmin'};

    # 2. use Exonerate. 

    my @exonerate = #  -intron => -10,
	$locus->exonerate( -protein => $args->{'-protein'}, -verbose => ($args->{'-verbose'}>=2 ? 2 :0 ) ) 
	unless ! $args->{'-protein'};

    # 3. run GeneWise  
    
    my @wise = 
	$locus->wise( -hmm => $args->{'-hmm'}, -verbose => ($args->{'-verbose'}>=2 ? 1 :0 ) ) 
	unless ! $args->{'-hmm'};

    # post process pseudocontig effects:: all orfs to real contig  
    # my @a2 = grep { $_->homology(-object => $self) } grep {defined} @anno;
    
    my @newmodels = grep { $_->commoncodon($self) } #grep { $_->strand == $self->strand }
    map { $_->transfer(-from => $locus, -to => $self->up, -warn => 0) }
    map { $_->adjust($adjust) } grep {defined} (@anno,@exonerate,@wise);
    $locus->DESTROY;    

    ##################################################################
    # choose best model from panel and destroy all others 
    ##################################################################
    
    my ($call,@valid) = grep {$_->translatable} grep {defined} ($self, @newmodels);
    my ($best,$delta) = $call->choose( # choose uses exonerate protein homology 
	-object => \@valid,            # to the supplied reference to choose "best"
	-verbose => $args->{'-verbose'}, 
	-reference => $args->{'-reference'}
	) if $call;

    # transfer scores from WISE etc. not 100% accurate but still useful. 

    my ($wisest) = sort { $b->score('global') <=> $a->score('global') } grep {$_->_creator =~ /wise/i} grep {defined} @wise;
    if ( $wisest && $wisest ne $best ) {
	$best->data('WISE' => $wisest->data('WISE') );
	$best->data('_WISE' => $wisest->data('_WISE') );
    }
    my ($exo) = sort { $b->score('global') <=> $a->score('global') } grep {$_->_creator =~ /exo/i} grep {defined} @exonerate;
    if ( $exo && $exo ne $best ) {
	$best->data('EXONERATE.DNA.LOCAL' => $exo->data('EXONERATE.DNA.LOCAL') );
    }
    
    # DESTROY 
    
    map { $_->DESTROY } grep { $_ ne $best } grep {defined} (@anno,@exonerate,@wise); # clean up 

    ##################################################################
    # simple return cases 
    ##################################################################

    return undef unless $best; # no valid models? caller must deal with this return code 
    $self->throw unless $best->translatable;
    return $self if $self eq $best;

    ##################################################################
    # new gene structure: double check boundaries 
    ##################################################################

    # some genes have legitimate extensions in the stricto but are 
    # penalized by Exonerate. we try to extend 
    
    my ($tp1,$tl1) = $self->_top_tail;
    my ($tp2,$tl2) = $best->_top_tail;
    if ( $tp2 ne 'ATG' && $tp1 eq 'ATG' ) {
	my $fes = $self->firstexon;
	my $feb = $best->firstexon;	
	my $val = -1*$fes->strand*($fes->start(-R => 1) - $feb->start(-R => 1));
	print 'ATG=>', $val, $args->{'-atg'} if $args->{'-verbose'};
	if ( $fes->sameframe(-object => $feb) && $val > 0 && $val <= $args->{'-atg'} ) {
	    $best->atg( -safe => 0, -verbose => $args->{'-verbose'} );
	}
    }
    
    ##################################################################
    # modify coordinates 
    ##################################################################

    my $retorf;
    if ( $args->{'-safe'} ) {

	$retorf = $best;

    } else {

	# transfer _exons_ between orfs so we do not have to deal with ohnologs/orthogroups etc 

	map { $self->remove(-object => $_, -warn => 0); $_->DESTROY; } $self->stream;
	map { $_->transfer(-from => $best, -to => $self, -warn => 0) } $best->stream;
	$self->_creator( $best->_creator );

	# destroy unrequired shell 

	$self->up->remove(-object => $best, -warn => 0, -log => 0);
	$best->DESTROY;
	
	# ugh ####################################
	my %str; map { $str{$_->strand}++ } $self->stream;
	$self->throw if scalar(keys %str)>1;
	$self->strand( (keys %str)[0] );
	$self->index;
	##########################################
	
	$retorf = $self;
    }

    ##################################################################
    # one last check ... then wipe info and recalculate 
    ##################################################################

    $self->throw unless $retorf->translatable;
    
    if ( $args->{'-update'} ) { 

	# exempt some things things that we have just calculated 
	# and some that are not easily recalculated 
	# NB: listed items are NOT reinited. all others are. 	
	
	$retorf->update(
	    -initialize => [qw(LTR TY LDA HMM RNA WISE EXONERATE.DNA.LOCAL)]
	    );
    }

    return( $retorf ); # repopulates %HOMOLOG and some of %EVIDENCE 
}

sub _eliminate_frameshifts {
    my $self = shift;
    my $args = {@_};
    return $self unless $self->exons > 1;

    foreach my $r ( $self->stream ) {
	next unless my $l = $r->left;
	if ( ($r->start(-R => 1) - $l->stop(-R => 1))*$l->strand == 1 ) {
	    $r->start(-R => 1, -new => $l->start(-R => 1));
	    $r->intron(-direction => 'left', -new => $l->intron(-direction => 'left') );
	    $self->remove(-object => $l);
	    $l->DESTROY;
	}
    }
    $self->throw unless $self->translatable;

    return $self;
}

=head2 choose(-object => [@models], -reference => 'homolog()|YAL017W')

    Choose amongst alternative gene models by comparing all to a homolog
    with exonerate and selecting best score. We always return a choice 
    (assuming protein coding genes with homology criterion met) but the 
    second return value (delta) may be undef, signifying a random choice. 
    
    Return choice and delta to next model. All scores are on ->data('_TMP').
    Alternatives are neither removed from contig nor destroyed. 

=cut 

sub choose {
    my $self = shift;
    my $args = {@_};

    # basic param checking 

    $self->throw unless exists $args->{'-object'};
    $args->{'-object'} = [$args->{'-object'}] unless ref($args->{'-object'}) eq 'ARRAY';
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    $args->{'-reference'} = 0 unless exists $args->{'-reference'};

    ############################################
    # set up vars 
    ############################################

    # my $key = 'EXONERATE.PROTEIN.GLOBAL'; # **** linked to var in exonerate2 ****    
    # 'global' is a shorcut for 'EXONERATE.PROTEIN.GLOBAL'    

    my $key = 'global';
    my @all = ($self, @{$args->{'-object'}}); # useful... 
    my @ref=();

    ############################################
    # basic QC 
    ############################################

    foreach my $obj ( @{$args->{'-object'}} ) {
	$self->throw unless $self->isa(ref($obj));
	$self->throw if $self eq $obj;
	$self->throw unless $self->commoncodon(-object => $obj);
    }
    map { $_->evaluate; $self->throw if $_->rank < -1 } @all;
    
    ############################################
    # fast return trivial case 
    ############################################

    # return inf not 0 as we are highly confident
    # that choice does not matter. 0 expresses no conf in choice. 

    return ($self, $INFINITY) if $#{$args->{'-object'}}==0 && 
	$self->sequence eq $args->{'-object'}->[0]->sequence;

    ############################################
    # sort out homolog
    # 1. everyone uses there own -- 0|undef 
    # 2. defined by caller -- provided 
    # 3. choose consensus -- 1 
    ############################################
    
    if ( $args->{'-reference'} ) {
	@ref = ( '-homolog', 
		 ( $args->{'-reference'} == 1 
		   ? $self->homolog(-object => $args->{'-object'}).'' # make consensus 
		   : $args->{'-reference'}) # use provided 
	    );
	@ref = () unless $ref[1] =~ /\w+/; # in case there is a failure 
    }
    
    ############################################    
    # finally, do the work.. 
    ############################################
    
    map { $_->exonerate2( -model => 'global', -return => 'score', @ref ) } @all;
    my ($best,@order) = sort { $b->score('global') <=> $a->score('global') } @all;
    my $delta =  $best->score('global') - ( $order[0] ? $order[0]->score('global') : 0 );    

    ############################################    
    # output 
    ############################################
    
    if ( $args->{'-verbose'} ) {
	print ">";
        $self->oliver(-prepend => ['1.SELF'], -append => [(caller(1))[3]]);
	map { $_->oliver(-prepend => ['2.MODEL'], -creator =>1) } 
	grep {$_ ne $self} (reverse(@order),$best);
	$best->oliver(-prepend => ['3.BEST'], -creator => 1);
    }

    return ( wantarray ? ($best, $delta) : $best );
}

=head2 purge2(-object => gene, -safe => 0)

    Choose between 2 genes. They do not need to be 
    alternative gene models as is required by choose(). 

=cut 

sub purge2 {
    my $self = shift;             
    my $args = {@_};
    
    my $other = $args->{'-object'};
    $self->throw unless $self->isa(ref($other));

    $self->evaluate(-structure => 1);
    $other->evaluate(-structure => 1);
    return undef if $self->rank < -1 && $other->rank < -1;

    my ($win,$lose,$delta);
    if ( $self->rank != $other->rank ) {

	($win,$lose) = sort {$a->rank <=> $b->rank} ($self,$other); 

	##################################################
	# below here ranks == and >= -1 (Pseudo or worse)
	##################################################

    } elsif ( $self->rank >= 2 && $other->rank >= 2 ) {
	
	if ( $self->evidence eq $other->evidence ) {
	    ($win,$lose) = sort {$a->data($self->evidence) <=> 
				     $b->data($self->evidence)} ($self,$other);
	    ($win,$lose) = ($lose, $win) if 
		$EVIDENCE{$self->evidence}->{'TEST'} eq '>';
	} else {
	    ($win,$lose) = sort {$b->length <=> $a->length} ($self,$other);				     
	}

	##################################################
	# below here ranks == and >= -1 and <= -1.5 (have homology)
	##################################################

    } elsif ( # if YGOB evidence, try shortcut 
	      ($self->evidence eq $other->evidence) && 
	      ($self->evidence eq 'YGOBHMM') && 
	      (abs(log($self->logscore('ygob')/$other->logscore('ygob'))/log(2)) > 1 )
	) {

	($win,$lose) = sort { $a->logscore('ygob') <=> $b->logscore('ygob') } ($self, $other);

    } else { # otherwise do full comparison to best homolog 
	($win,$delta) = $self->choose(-object => $other, -reference => 'SGD'); # should be GENE in general 
	$lose = ( $win eq $self ? $other : $self );
    }

    $win->up->remove(-object => $lose,-force => 1) and $lose->DESTROY 
	unless $args->{'-safe'};

    return $win;
}

=head2 distance(-object => , -bases => 1, -strand => 1, 
    -rank_min => x, -rank_max => y)

    Get the distance between two orfs in terms of the number of genes.
    Can choose to count only genes on same strand or with certain 
    types of evidence (-rank_min / -rank_max). With -bases => true
    returns the number of bases instead. 

=cut

sub distance {
    my $self = shift;
    my $args = {@_};
    my $dist = 0;
    
    # fails to return when clustering HSPs becasue may have 
    # many HSPs with same coords. just return -- dont need anyway. 
    return undef if defined $args->{'-pseudo'}; 
    
    $args->{'-rank_max'} = $INFINITY unless exists $args->{'-rank_max'};
    $args->{'-rank_min'} = -$INFINITY unless exists $args->{'-rank_min'};
    $args->{'-strand'} = undef unless exists $args->{'-strand'};
    $args->{'-nogap'} = undef unless exists $args->{'-nogap'};
    
    $self->throw("Class") unless ref($args->{'-object'}) eq ref($self);
    return $INFINITY unless $args->{'-object'}->up eq $self->up;
    return -1 if ($args->{'-object'} eq $self)
	|| ($args->{'-object'}->start == $self->start);
    
    # order genes
    
    my ($x,$y) = sort {$a->start <=> $b->start}
    ($self, $args->{'-object'});
    
    # distance in bases (minimum)
    
    if ($args->{'-bases'} || $args->{'-nucleotides'}) {
	return $y->start - $x->stop; # can be negative 	
    }
    
    # distance in number of genes 
    
    $self->throw unless my $right = $x->right;
    until ($right eq $y) {
	$self->throw unless $right = $right->right;
	$dist++ unless 
	    ($y->assign eq 'GAP' && $args->{'-nogap'}) ||
	    $y->rank > $args->{'-rank_max'} ||
	    $y->rank < $args->{'-rank_min'} ||
	    ($args->{'-strand'} && $right->strand != $self->strand);
    }
    
    return $dist;
}

#########################################
#########################################

# evidence accumulation methods  

=head2 homology(-blast => 'array_of_blast_hits', -object => orf_obj, 
    -species => [keys %HOMOLOGY], -direct => 1, -evalue => 10)

    -blast : An array of BLAST hits as returned by ->blast().
    Populate DATA hash with top hits from each species using
    regexes from %HOMOLOGY in GlobalVars.pm
    -object : Orf object. Initiate a comparison between caller and
     arg for each species in DATA hash. Return number of shared hits. 
    -species : restrict comparison to named species. ['KLAC', 'KPOL', ect]    
    -direct : do direct (using ->compare) comparison of orfs. 
    -evalue : only consider hits where both species are -evalue. 

=cut

sub homology {    
    my $self = shift;
    my $args = {@_};

    # 
    $args->{'-direct'}=0 unless exists $args->{'-direct'};
    $args->{'-align'}=0 unless exists $args->{'-align'};
    $args->{'-species'}=[ keys %HOMOLOGY ] unless exists $args->{'-species'};
    $args->{'-evalue'}=10 unless exists $args->{'-evalue'};
    # process BLAST 
    $args->{'-sort'} = 'score' unless exists $args->{'-sort'};

    $self->throw("Require Orf obj or BLAST report")
        unless exists $args->{'-blast'} || exists $args->{'-object'};
    $self->throw("Supply Orf obj OR BLAST report")
        if exists $args->{'-blast'} && exists $args->{'-object'};

    goto COMPARE if $args->{'-object'};

    #########################################        
  BLAST: # extract homology info from BLAST report
    
    $self->throw("Expected BLAST report object: $arg")
	unless ref($args->{'-blast'}) eq 'ARRAY';

    # we sort by evalue but break ties by score 
    # because we get extensive evalue saturation
    # -- many hits recording evalue of 0 

    $self->_init_homology_data( -exempt => ['YGOB'] );
    foreach my $hsp ( sort {$a->{ 'EVALUE' } <=> $b->{ 'EVALUE' } } @{ $args->{'-blast'} } ) {
	next unless my $key = _species_key( $hsp->{'HIT'} );	
	$self->accept($key, $hsp) and next unless $self->data($key);	
	next unless $hsp->{EVALUE} <= $self->evalue($key); # both 0, saturated 
	$self->accept($key, $hsp) if $hsp->{SCORE} > $self->score($key);	
    }
    return $self;
    
    #########################################        
  COMPARE: # comapre 2 objects

    my $obj = $args->{'-object'};
    $self->throw unless $self->isa( ref( $obj ) );
    
    my %scr;
    foreach my $k ( @{$args->{'-species'}} ) {
	$k .= '*' if $self->data($k.'*'); # use synteny moderated homology if available
	next unless my $sh = $self->data($k);
	next unless my $ah = $obj->data($k);
	next unless $ah =~ /\w/ && $sh =~ /\w/;
	$scr{$ah}++ if ($ah eq $sh) &&
	    ($self->evalue($k) <= $args->{'-evalue'}) &&
	    ($obj->evalue($k) <= $args->{'-evalue'});
    }

    if ($self->data('GENE') && $self->data('GENE') eq $obj->data('GENE')) {	
	$scr{ $self->data('GENE') }++ unless exists $scr{ $self->data('GENE') };
    }
    return scalar(keys %scr) unless $args->{'-direct'};

    #########################################        

    my $bl2 = $self->overlap(-object => $obj, -compare => 'seq');
    my $scr = (scalar(keys %scr) || 1) if $bl2;
    return $scr; # $bl2 is fractional overlap 
} 

=head2 coding

    Return true if the object is a protein coding gene.
    We return false for transposable elements (eg Ty ORFs).

=cut 

sub coding {
    my $self=shift;
    return 1 if $self->evidence =~ /HCNF/;
    return 1 if $self->ygob eq 'Anc_1.380' && $self->hypergob >=5; # HAP1 exception. this is crazy...
    return ( $self->assign =~ /FEATURE|GAP|TELOMERE|CENTROMERE|RNA|PSEUDO|REPEAT/ ? 0 : 1); # INTER?
}

=head2 noncoding

    Returns true/false.

=cut 

sub noncoding {
    my $self=shift;
    return ( $self->coding ? 0 : 1);
}

=head2 syntenic(-distance => 5, -hyper => 5, -loss => 25)
    
    Return 1/0 value for whether gene is at the syntenic locus by 
    runnning both synteny(loss()) and synteny(hypergeometric()).
    Examines -distance genes either side of caller. 

    Defaults to *liberal* hardcoded score cutoffs from distribution
    plots obtained using complete genome assemblies. 
    
=cut 

sub syntenic {
    my $self = shift;
    my $args = {@_};

    return undef if $self->rank < 0;

    $args->{'-distance'} = 5 unless exists $args->{'-distance'};
    $args->{'-hyper'} = 5 unless exists $args->{'-hyper'};
    $args->{'-loss'} = 25 unless exists $args->{'-loss'};

    my @scores = (
	$self->synteny(-loss => 1, -distance => $args->{'-distance'}),  
	$self->synteny(-hyper => 1, -distance => $args->{'-distance'}, -restrict => ['YGOB'])
	);
    my $res = ( ($scores[0] >= $args->{'-loss'} || $scores[1] >= $args->{'-hyper'}) ? 1 : 0);

    return (wantarray ? ($res, (map { sprintf("%.1f", $_) } @scores) ) : $res );
}

=head2 outgroup(-object => [], -species => [qw(PRE WGD SPECIES)],
    -evalue => 1e-10)

    Return a credible outgroup by using _homolog() to look up
    pillar and other info and returning best hit from a default 
    panel of Pre-WGD yeasts. 

    This method is entirely dependent on -species array -- 
    it does not perform any genuine tests for outgroup status. 

=cut

sub outgroup {
    my $self = shift;
    my $args = {@_};

    $args->{'-object'} =  [$self->orthogroup]  unless exists $args->{'-object'};
    $args->{'-species'} = [qw(SKLU KLAC ZROU EGOS KTHE KWAL)] unless exists $args->{'-species'};
    #$args->{'-synteny'} = 0 unless exists $args->{'-synteny'};
    $args->{'-evalue'} = 1e-10 unless exists $args->{'-evalue'};
    my $postwgd = [qw(KPOL SCAS CGLA)];

    # 
    
    my $outG=undef;
    unless ( $outG = $self->_homolog( %{$args} ) ) {
	$args->{'-species'} = $postwgd;
	$outG = $self->_homolog( %{$args} );
    } 
    return () unless $outG;

    # do some light QC. outgroups with low %ID are not worthwhile. 
    # NB: even if the returned homolog is not the one stored by $self->data
    # the fact that it was preferred, implies similar qulaity. 
    
    my ($sp) = _species_key($outG);
    return () unless $self->evalue($sp) <= $args->{'-evalue'};
    
    #
    
    my ($aafile) = $self->fetch(-id => $outG, -file => 1, -molecule => 'aa' );
    return () unless -e $aafile && -s $aafile > 20;
    chomp(my $seqlen = `grep -v '>' $aafile | sed s/[[:space:]]//g | sed s/[[:space:]]+//g | wc -c `);
    $seqlen =~ s/\s+//;
    
    return ( wantarray ? ($sp, $outG, $seqlen, $aafile) : $outG );
}

=head2 homolog(-object => $orf, -species => [ keys %HOMOLOGY ], 
    -delta => 20, -model=> 'local', -fast => 0|1) 

    Accept either a single object or a set of objects and return 
    the highest confidence (shared) homolog among -species. 

    If pillar() returns a value and its score exceeds -delta the
    homolog will be chosen from this pillar. Else we weight all 
    available proposals [pillar(), ygob(), data()] by the exonerate
    score (using -model) and choose the best. 

    Once a pillar has been chosen, we sort the genes in the pillar
    by exonerate score to obtain a final homolog. 
    
  NB: This method utilizes synteny information indirectly
    through pillar() and synteny moderated BLAST. Both use LOSS
    data so if not computed, homolog should be regarded carefully.
    
    Returns (species, id, length, aaFsaFile) OR id in scalar context

    -fast : return precomputed homolog if set. if not defualts to full
            and _slow_ homolog inference. 
    
=cut 

sub homolog {
    my $self = shift;
    my $args = {@_};
    
    # [ grep {!/YGOB|SPOM|MITO/} keys %HOMOLOGY ]
    $args->{'-object'} = [ $args->{'-object'} ] if $args->{'-object'} && $self->isa(ref( $args->{'-object'} ));
    $args->{'-species'} = [qw(SGD SKLU KLAC EGOS CGLA KPOL SCAS KWAL)] unless exists $args->{'-species'};
    
    return undef if $self->rank < 0;
    $self->throw if $args->{'-fast'} && $args->{'-object'}; 
    $self->throw if $args->{'-fast'} && wantarray; 

    my $key = 'HOMOLOG';

    if ( $args->{'-fast'} ) { 
	if ( my $hom = $self->data($key) ) {
	    return $hom;
	}
    } 
    
    return $self->_homolog( %{$args} );
}

sub _homolog {
    my $self = shift;
    my $args = {@_};
    
    ###################################################
    # preamble
    ###################################################

    $args->{'-model'} = 'local' unless exists $args->{'-model'};    
    $args->{'-delta'} = 20 unless exists $args->{'-delta'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    $args->{'-pillar'} = undef unless exists $args->{'-pillar'};

    my $fh = *STDERR;
    my $key = 'HOMOLOG';
    my $choice;

    # basic checking 
    
    my @orfs = grep {defined} ($self, @{$args->{'-object'}});
    map { $self->throw unless $self->isa(ref($_)) } @orfs;
    map { $self->throw($_) unless exists $HOMOLOGY{$_} } @{$args->{'-species'}};

    goto PILLAR if $args->{'-pillar'};

    ###################################################
    # speed up for consensus case 
    ###################################################

    my $homolog;
    if ( $#orfs > 0 ) {
	my %hom;
	# the grep function removes homologs that are of form 6319628
	# these are hits from Genbank
	map { $hom{ $_ }++ } grep {!/^\d+$/} grep {/\w/} grep {defined} 
	map { $_->homolog(-fast => 1).'' } @orfs;
	my ($first) = sort { $hom{$b} <=> $hom{$a} } keys %hom;

	if ( $first ) {
	    my ($sp) = _species_key($first);
	    if (grep {/$sp/} @{$args->{'-species'}} ) {
		$homolog = $first and goto FINISH if 
		    $first && $hom{$first} > scalar(@orfs)/2;
	    }
	}
    }
    
    ###################################################
    # grab data from trusted sources 
    ###################################################
    
    my (%pillar, %ygob, %blast);
    foreach my $o ( @orfs ) {
	# pillar() -- best YGOB pillar from BLAST/impute/GENE and decided using LOSS
        my ($p,$s) = $o->pillar;
	#print $o->name,$p,$s;
	if ( ! $p && $o->gene ) { # eg no pillar in YGOB but gene exists in SGD 
	    ($p,$s) = ($o->gene,$o->synscore('gene')) # swamped by _PILLAR values 
		unless $o->gene =~ /[a-zA-Z]{3,}\:[a-zA-Z]{3,}/; # exclude tRNAs
	}
	$pillar{ $p } += $s if $p;
	print {$fh} 'PILLAR', $p, $s if $args->{'-verbose'} >= 2; #######

	# ygob() -- direct look up of best YGOB-HMM hit 
	my $py = $o->lookup( -query => $o->ygob ) if $o->ygob;
	$ygob{ $py->{ID} } += $o->synteny(-hyper => 1, -restrict => ['YGOB']) if $py;	
	print {$fh} 'YGOB', $py->{ID}, $o->synteny(-hyper => 1, -restrict => ['YGOB']) 
	    if $args->{'-verbose'} >= 2; #######
	
	# tedious blast() -- reexamine all synteny moderated BLAST best hits 
	foreach my $sp ( grep {$o->data($_.'*')} @{$args->{'-species'}} ) {
	    
	    my ($chr,$index,$sd,$bigN) = $o->chromosome( $sp );
	    next unless my $syn_mod_hit = $o->data( $sp.'*' );
	    next unless $bigN >= 3 && $sd < 10;
	    
	    my ($sp2, $chr2, $index2) = _decompose_gene_name( $syn_mod_hit );
	    next unless $chr && $chr2;
	    next unless "$chr" eq "$chr2";       
	    next unless abs($index - $index2) < 10;
	    
	    my $pill = $o->lookup( -query => $syn_mod_hit );
	    $blast{ ($pill ? $pill->{ID} : $syn_mod_hit ) } += $o->synscore($sp.'*'); # if no pillar, just use ID 

	    print {$fh} ($sp, ($pill ? $pill->{ID} : $syn_mod_hit ), $o->synscore($sp.'*')) 
		if $args->{'-verbose'} >= 2; #######
	}
    }
    
    ###################################################
    # choose best by each metric
    ###################################################

    my %uniq;
    my ($pil) = sort { $pillar{$b} <=> $pillar{$a} } keys %pillar;
    my ($gob) = sort { $ygob{$b} <=> $ygob{$a} } keys %ygob;
    my ($blast) = sort { $blast{$b} <=> $blast{$a} } keys %blast;  
    my $nn = scalar( grep { !$uniq{$_}++ } grep {defined} ($pil,$gob,$blast) );

    # 
    
    #print {$fh} $self->_method, $self->identify, 
    #($pil || 'NA'), $pillar{$pil}, ($gob || 'NA'), $ygob{$gob}, ($blast || 'NA'),  $blast{$blast}
    #if $args->{'-verbose'} >= 2; 

    ###################################################
    # we have candidates. choose a pillar based on local aligment score. 
    # but use a voting based shrotcut where possible. 
    # NB: we cannot use synteny since 
    ###################################################

    # 0. if pil is good, just use 
    # 1. if pil weaker do comparison 
    # 2. if no choice, just use 
    # 3. if nothing, just return 

    if ( $nn == 0 ) {
        ($choice) = undef;
    } elsif ( $nn == 1 ) {
	($choice) = grep {defined} ($pil,$gob,$blast);
    } elsif ( $nn > 1 ) {
	($choice) = grep {defined} ($pil,$gob,$blast);
	
	unless ( $pil && ( $pillar{$pil}/scalar(@orfs) >= $args->{'-delta'} ) ) {
	    foreach my $sp ( @{$args->{'-species'}} ) {
		my %scr;
		foreach my $p ( grep {defined} ($pil,$gob,$blast) ) {
		    next unless my $px = $self->lookup( -query => $p );
		    $scr{ $px->{ID} }{GENE} = $px->{ $sp }->[0] if $px->{ $sp }->[0]; 
		}
		next unless scalar( keys %scr ) == $nn; # same sp present in all
		
		# compare pillars using this species 
		
		foreach my $p ( keys %scr ) {
		    map { $scr{ $p }{SCR} += 
			      $_->exonerate2( 
				  -homolog => $scr{$p}{GENE}, 
				  -model => $args->{'-model'}, # local adequate? certianly faster... 
				  -return => 'score'
			      ) || 0 } @orfs;
		    #print $self->name, $sp, $p, $scr{$p}{GENE}, $file, $scr{$p}{SCR};
		}
		($choice) = sort { $scr{$b}{SCR} <=> $scr{$a}{SCR} } keys %scr;
		last;
	    }
	}
    }
    return undef unless $choice;

    ###################################################
    # if we have a pillar, choose a species using local align scores. 
    ###################################################

  PILLAR:

    # we have supplied a pillar but do not know the correct 
    # homolog to use. it may be either an Anc or a real gene. 
    if ( $args->{'-pillar'} && ! $choice ) {
	$self->throw($args->{'-pillar'}) unless 
	    my $pil = $self->lookup( -query => $args->{'-pillar'} );
	$choice = $pil->{ID};
    } elsif ( $args->{'-pillar'} || ! $choice ) { 
	$self->throw(" $args->{'-pillar'}/$choice ");
    }

    # 

    if ( $choice =~ /YGOB_\d+_DS/ ) {

	$self->throw unless 
	    my $pill = $self->lookup( -query => $choice );
    
	# we have a pillar. choose the best species to use as an outgroup 
	# from the pre-specificed list by using exonerate scores.  
	
	my %mat;
	foreach my $i ( 0..$#{$args->{'-species'}} ) {
	    my $sp = $args->{'-species'}->[$i];
	    foreach my $gene ( grep {defined} @{ $pill->{$sp} } ) {
		next unless my $file = $self->fetch( -id => $gene, -file => 1, -molecule => 'aa' );
		map { $mat{ $gene } += $_->exonerate2( 
			  -homolog => $file, -model => $args->{'-model'}, -return => 'score') || 0 } @orfs;
		$self->_cleanupfile( $file );
	    }
	    # speed up hack. do not look at all. assumes @species is ordered 
	    last if $i >= 2 && scalar( grep { $_>0 } values %mat ) > 1;
	}
	return undef unless %mat;
	
	# 
	
	($homolog) = sort { $mat{$b} <=> $mat{$a} } keys %mat;
    } else { 
	$homolog = $choice;	    
    }

    ###################################################    
  FINISH: # standard and shortcut share remaining code (return values)
    ###################################################    
    
    print {$fh} $self->name, $homolog if $args->{'-verbose'} >=2;    
    $self->data( $key => $homolog );
    return $homolog unless wantarray;

    ###################################################    
    # epilogue: print results. long return form 
    ###################################################
    
    my ($species) = _species_key($homolog);
    my ($aafile) = $self->fetch(-id => $homolog, -file => 1, -molecule => 'aa' );
    return () unless -e $aafile && -s $aafile > 20;
    chomp(my $seqlen = `grep -v '>' $aafile | sed s/[[:space:]]//g | sed s/[[:space:]]+//g | wc -c `);
    $seqlen =~ s/\s+//;
    
    # 

    print {$fh} $self->name, $self->sgd, $self->gene, $self->ygob, $species, $homolog, $seqlen, $aafile, $pill->{ID}
	if $args->{'-verbose'};
    return ($species, $homolog, $seqlen, $aafile);
}

=head2 impute(-evalue => 1e-5, -blast => [], -homology => 'AA|NCBI', 
    -reference => 'SGD', -safe => 1)

    Impute the gene name based on the supplied BLAST results and synteny. 
    We also set the relevant homology vars. 
    
    Returns new gene name if GENE updated otherwise returns undef.

    -reference : GENE must be from this species. 
    -safe : return gene name. do not set GENE.

=cut 

sub impute {
    my $self  = shift;
    my $args = {@_};

    $args->{'-homology'} = 'AA' unless exists $args->{'-homology'}; # nonoptional ?
    $args->{'-evalue'} = 1e-5 unless exists $args->{'-evalue'};
    $args->{'-reference'} = 'SGD' unless exists $args->{'-reference'};
    #
    $args->{'-drop'} = 1e-5 unless exists $args->{'-drop'};
    $args->{'-percdrop'} = 10 unless exists $args->{'-percdrop'};
    # 
    $args->{'-safe'} = 0 unless exists $args->{'-safe'};
    $args->{'-store'} = 0 unless exists $args->{'-store'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    
    #my $blast = $self->blast(-db => $args->{'-blastdb'}, -qseq => [$self]);
    #goto FINISH unless exists $blast->{$self->_internal_id};
    #my @bb = sort {$a->{EVALUE} <=> $b->{EVALUE}} @{$blast->{$self->_internal_id}};

    ################################################################
    # do some standard processing. 
    ################################################################

    $self->throw unless ref($args->{'-blast'}) eq 'ARRAY';
    return undef unless $#{$args->{'-blast'}} >= 0; # an empty aray is OK 
    $self->throw unless exists $args->{'-blast'}->[0]->{'HIT'};

    # simplify and check safety 

    my $key = $args->{'-homology'};
    my @bb = sort {$b->{SCORE} <=> $a->{SCORE}} @{ $args->{'-blast'} }; # resort array 

    return undef unless ! $args->{'-safe'} || # -safe=0 is default 
	$bb[0]->{'EVALUE'} < $self->data($key);
    
    # update homology evidence
    # evidence is retained irrespective of significance 

    $self->accept($key, $bb[0]);
    
    ################################################################
    # Based on homolgy type, do we do any extra processing?
    # mainly relevant for AA where we use synteny to impute a most likely 
    # syntenic ortholog where similarity is ambiguous. 
    ################################################################
    
    my $best;
    if ( $key eq 'NCBI' ) { # best set by NBCI only if no local DB hit 
	
	$self->accept($key => $bb[0]); 	
	unless ( $self->data('GENE') ) {
	    $bb[0]->{'HIT'} =~ /gi\|(\d+)\|/ || $self->warn("No NCBI name!");
	    $best = { HIT => $1, EVALUE => $bb[0]->{EVALUE}, SCORE => $bb[0]->{SCORE} } if $1;
	}
	
    } elsif ( $key eq 'AA' ) {	# best is ALWAYS set by AA 

	( $best ) = map { shift(@{$_}) } grep { $_->[1] eq $args->{'-reference'} } map { [$_, _species_key($_->{HIT})] } @bb;
	$best = $bb[0] unless $best; # we set best to any other sp if ref not present. 
	
	# 2 scenarios. best hit > others or ~ equal. in latter case we use synt. 
	# is synteny (ie LOSS) is not available we do no imputation. 

	if ( defined $self->synteny( -lod => 1 ) ) { # 0 OK, undef not 
	    
	    foreach my $sp ( keys %HOMOLOGY ) {		
		my @rest = map { shift(@{$_}) } grep { $_->[1] eq $sp } map { [$_, _species_key($_->{HIT})] } @bb;
 		next unless @rest;

		my @top = shift( @rest );
		foreach my $bx ( @rest ) {
		    last if _logdropbig(
			$top[-1]->{EVALUE}, $bx->{EVALUE}, 
			$args->{'-percdrop'}, # these are minimum 
			$args->{'-drop'}      # required values ****
			);
		    push @top, $bx;
		}
		
		# is this family just too complicated? 
		# max top->bottom drop 1e-20 and 20% from height 
		
		my $complex =
		    # if either scores 1, then too complex. 
		    # return value of 1 means non-zero value has been exceeded.  
		    _logdropbig($top[0]->{EVALUE}, $top[-1]->{EVALUE}, 20, 1) +   # these are max 
		    _logdropbig($top[0]->{EVALUE}, $top[-1]->{EVALUE}, 0, 1e-20); # allowed values **** 
		next unless ($#top==0 || ! $complex);

		# compute synteny for each and choose the best one 

		my %val;
		foreach my $i ( 0..$#top ) {
		    my $restore = $self->data($sp);
		    $self->data( $sp => $top[$i]->{HIT} );
		    $val{ $i } = $self->synteny( -lod => 1, -distance => 15, -restrict => [$sp]);
		    $self->data( $sp => $restore);
		}
		my ($max,$sec) = sort { $val{$b} <=> $val{$a} } keys %val;
		my $spbest = $top[$max];
		
		# add a synteny score based on LOSS delta to second hit 
		# and accept as usual. 

		$spbest->{'SYNTENY'} = $val{$max} - ($val{$sec} || 0);
		$self->accept($sp.'*', $spbest); 
		$best = $spbest if $sp eq $args->{'-reference'}; # we reset best if we have a better guess 
		
	    } # %homology      
	}
	
    } elsif ($key =~ /LTR|TY/) { # best is NEVER set by LTR/TY 
	$best=undef;
    } else { $self->throw; }

    #####################################################
    # can we set GENE using imputed best hits? 
    #####################################################

    if ( $best && $best->{EVALUE} <= $args->{'-evalue'} ) {
	$self->accept('GENE', $best); 

	# set the DESCRIPTION atrribute. this is an ORF property 
	# NOT stored on DATA hash. do not rely. experimental.
	
	my ($org,$def);
	if ( $homology eq 'NCBI' ) {
	    if ( my $content = get($eutils.$self->data('GENE')) ) {
		#$content =~ /DEFINITION\s+([^\n]+)\n/;
		$content =~ /GBSeq_definition\>([^\<]+)/;#   <GBSeq_definition>Spp381p [Saccharomyces cerevisiae S288c]</GBSeq_definition>
		$def = $1;
		#$content =~ /ORGANISM\s+([^\n]+)\n/;
		$content =~ /GBSeq_organism\>([^\<]+)/;# Saccharomyces cerevisiae S288c</GBSeq_organism>
		$org = $1;
	    }
	}

	##########################################################################
	# DESCRIPTION is only set here and by contig->fuse() to label GAP.
	# redundat with GENE. think about this. 
	$self->_description( join('/', grep {defined} ($self->data('GENE'), $org, $def)) );
	##########################################################################
    }

    return $self->description;
}

sub description {
    my $self = shift;
    $self->throw if shift;
    return $self->{DESCRIPTION};
}

sub _description {
    my $self = shift;
    $self->throw unless exists $_[0];
    $self->{DESCRIPTION} = shift;
    return $self->{DESCRIPTION};
}

sub _logdropbig {
    my ($e1, $e2, $pd, $ed) = @_;
    #print @_;

    die unless defined $e1 && defined $e2;
    die unless $e1 >= 0 && $e1 <= 10;
    die unless $e2 >= 0 && $e2 <= 10;    
    die unless $e1 <= $e2; # can be same 
    
    return 0 if $e1 >= 1;  # OK with this? 
    
    $pd = 20 unless defined $pd;
    die unless $pd >= 0 && $pd <= 100; # min 1% drop 
    $ed = 1e-5 unless defined $ed;
    die unless $ed >= 0 && $ed <= 1;

    # 

    my $le1 = - ($e1 == 0 ? -200 : log($e1)/log(10) ); # 1e-150  --> 150 
    my $le2 = - ($e2 == 0 ? -200 : log($e2)/log(10) ); # 1e-20   --> 20 
    my $led = - ($ed == 0 ? -200 : log($ed)/log(10) ); # 1e-5   --> 5

    #

    my $delta = $le1 - $le2; # --? 130 
    my $perc = $le1 * ( $pd/100); # --> 15 
    #print $delta, $perc, $led;
    
    return( (($delta >= $perc) && ($delta >= $led) ) ? 1 : 0); # >= required fro complexity test in blast2gene 
}

##########################################
##########################################

=head2 synteny(-object -> orf, -neighbour => 1, 
    -logodds => 1, -hypergeometric => 1,
    -difference => 5, -distance => 2, -spanning => 1, 
    -restrict => ['Zbai'], -set => 1)
    
    This method computes a variety of synteny scores for either the
    query gene alone or between the query and a specified object in the 
    same genome. To compares orfs from diff genomes see: conserved_synteny().

  Return: undef if the desired metric cannpt be computed, otherwise 0 .. X.

    There are 6 run modes though only 1,2,3,4 are truly distinct 
    (5,6 are iterations of 3): 
    
    1. -loggodds => 1 (-lod also accepted)
    Compute a log-odds synteny score (in bits) for the caller given the
    relationships of its homologs to the homologs of neighbouring genes. 
    Given a 10 gene window and 10 species, this tends to be a continuously 
    distributed value between -100 .. +400. It is unique in that it does 
    not require (or accept) a -difference argument.

    2. -hypergeometric => 1 (-hyper also accepted)
    Compute hypergeometric based score 
    For 10 gene window, score is quasi-continuous in the range 0..50.

    3. -object => $orf 
    Core method for remaning call forms (non-statistical counting based). 
    Compares both orfs best hits in all other species (based on previoulsy 
    run BLASTs) and asks whether the hits are on the same contig and 
    are <= -difference genes apart. The returned score is the number of 
    speices in which this condition is met. If -restrict is set only 
    species with keys matching the supplied regex are used to determine 
    synteny (e.g. Z will match ZBAI/ZKOM/ZBIS but not EGOS/SCER/KWAL). 
    
    4. -neighbour => 'all|left|right', -restrict => 'Sbay|Zlen|YGOB'
    Calls -object form on $left and/or $right. If return value is positive
    (ie determined to be syntenic), score is incremented and returned. 
    Score is either 0|1 or 0|1|2 depending on call param. 
    This method would be equivalent to -distance => 1 (and still can be) 
    except that with appropriate values for -restrict it knows how to 
    compare synteny for genomes in bound genome objects -- not just using 
    homology from BLAST etc.
    
    5. -direction => 'left|right', -distance => i
    Compares self to -distance genes in the specified -direction
    by calling $self->synteny(-object => $obj) on each gene in set.
    Returns number of valid genes considered and the number of genes 
    that have non-zero score in array form. Else latter only. 
    
    6. none of the above. This is default. 
    The genes both sides of the gene of interest are queried up to 
    -distance. Equivalent to summing -direction => 'left' and -direction => 'right'.
    If -spanning is set (default) then the specified level of support must
    be met on both sides for a non-zero synteny score to be returned. 
    
=cut

sub querysynteny { my $self = shift; return $self->synteny(@_)+0 } # depracated :: check all uses and wantarray returns 
sub hypergob { my $self = shift; return sprintf("%.2f",$self->synteny( -hypergeometric => 1, -restrict => ['YGOB'], -distance => 7 )); }
sub loss { my $self = shift; return sprintf("%.2f",$self->synteny( -logodds => 1, -distance => 3 )); } 

sub synteny {
    my $self = shift;
    my $args = {@_};

    $self->throw if exists $args->{'-query'}; # depracated 
    
    $self->throw if 
	scalar(grep {$args->{$_}} qw(-neighbour -object -hyper -lod -logodds -hyerpgeometric)) > 1;
    
    if ( $args->{'-hyper'} || $args->{'-hypergeometric'} ) {
	return $self->_hypergeometric_synteny_score(@_)+0;
    } elsif ( $args->{'-lod'} || $args->{'-logodds'} || $args->{'-loss'} ) {
	# test whether LOSS available. now tested in _logodds_synteny_score;
	# return undef unless defined $self->up->up->_access_synteny_null_dist(); 
	return $self->_logodds_synteny_score(@_); # +0
    } else {
	return $self->_synteny(@_);
    }
}

sub _synteny {
    my $self = shift;
    my $args = {@_};

    # special synteny params 
    #$args->{'-neighbour'}; # we test existence 
    #$args->{'-object'};    # we test existence 
    $args->{'-set'} = 0 unless exists $args->{'-set'};
    # generic synteny params 
    $args->{'-restrict'} = [ keys %HOMOLOGY ] unless exists $args->{'-restrict'};
    $args->{'-distance'} = 2 unless exists $args->{'-distance'};
    $args->{'-difference'} = 5 unless exists $args->{'-difference'};
    $args->{'-spanning'} = 1 unless exists $args->{'-spanning'};
    # generic params 
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    
    # 4 possibilities 
    
    my ($num,$denom) = (0,0);  
    if (exists $args->{'-neighbour'}) {
	$args->{'-neighbour'} = 
	    [ ($args->{'-neighbour'} =~ /left|right/ ? $args->{'-neighbour'} : qw(left right)) ];

	my %bound = map {uc($_) => 1} $self->up->up->bound;
	foreach my $dir (@{$args->{'-neighbour'}}) {
	    next unless my $other = $self->neighbour(-dir => $dir);
	    foreach my $key ( map {uc($_)} @{$args->{'-restrict'}}) {
		if ( exists $HOMOLOGY{$key} ) {
		    $num++ and last if 
			$self->synteny(-object => $other, -restrict => [$key], -difference => 1);
		} elsif ( exists $bound{$key} ) {
		    my ($sp,$op) = map { $_->$meth } ($self,$other);
		    next unless $sp && $op;
		    $num++ and last if 
			($sp->neighbour(-dir => 'left') eq $op || $sp->neighbour(-dir => 'right') eq $op);
		} else { $self->throw("$dir,$key"); }
	    }
	}
	return $num;
	
    } elsif (exists $args->{'-object'}) { # here we do the bulk of the work. return a score only. 
	$self->throw unless $self->isa(ref($args->{'-object'}));
    	
        foreach my $key ( grep {!/SPOM/} @{$args->{'-restrict'}}) {
            next unless my $sh = $self->data($key);
            next unless my $ah = $args->{'-object'}->data($key);
	    unless ( $args->{'-difference'} == 0 ) { next if $sh eq $ah; } # used in contig::merge 
	    next unless my $delta = _compute_gene_distance($ah,$sh);
            $num++ if $delta < $args->{'-difference'};
	    #print $key, $delta, $args->{'-difference'}, $num; 
        }
	return $num;
	
    } elsif ($args->{'-direction'}) {
	
	# We recruse all the real work to -object variant.
	# return both numerator and denominator. 

	foreach my $obj ( grep {$_ ne $self} $self->traverse(
			    -direction => $args->{'-direction'},
			    -distance => $args->{'-distance'})) {
	    $denom++;
	    $num++ if $self->synteny(
		-object => $obj,
		-difference => $args->{'-difference'},
		-restrict =>  $args->{'-restrict'}
		) > 0;
	}

    } else { # recurse work to -direction, then to -object .. 
	
	my %hash;
	foreach my $dir ( qw(left right) ) {
	    my ($n2,$d2) = $self->synteny(
		-direction => $dir,
		-distance => $args->{'-distance'},
		-difference => $args->{'-difference'},
		-restrict =>  $args->{'-restrict'}
		);
	    $hash{$dir}{'NUM'} = ($n2 || 0);
	    $hash{$dir}{'DENOM'} = ($d2 || 0);
	}
	
	$num = $hash{'left'}{'NUM'} + $hash{'right'}{'NUM'};
	$denom = $hash{'left'}{'DENOM'} + $hash{'right'}{'DENOM'};
	$num = 0 if ( $hash{'left'}{'NUM'} < $args->{'-spanning'} || $hash{'right'}{'NUM'} < $args->{'-spanning'});
	
	if ( $args->{'-set'} ) {
	    $self->data( 'SYNT' => $num );
	    $self->evaluate;
	}
    } 
    
    return (wantarray ? ($num,$denom) : $num );
}

=head2 _hypergeometric_synteny_score(-mode => 'position|species', 
    -difference => 7, -distance => 10, -restrict => [ALL]) 
    
    Compute a synteny score based on the hypergeometric
    probability of observing X homologous genes in a Y 
    gene window around the gene of interest given Z genes
    in the genome. 
    
    eg we look at YGOB hits for 10 genes up- and down-stream
    and 14 of these are within W genes of the query YGOB hit:

    -log10( Phyper(14|20,4700) )
    
    We automatically adjust the window size for genes at contig
    ends (ie if only 15 genes examined, we write P(14|15,4700)).

    The major weakness of this method is that it requires the
    fudge param W and does not distinguish proximity within the 
    window. 

    Two modes: 
    1. By default we consider positions in turn, iterate over 
    species and accept any position with a syntenic gene in any
    species as positive. This tends to mildly inflate synteny.
    2. We iterate over species, compute hypergeomettric independently 
    for each and average the results. This tends to deflate. 

=cut

sub _hypergeometric_synteny_score {
    my $self = shift;
    my $args = {@_};
    
    delete $args->{'-hyper'};
    delete $args->{'-hypergeometric'};
    
    # 

    $args->{'-mode'} = 'position' unless exists $args->{'-mode'};
    # 
    $args->{'-distance'} = 7 unless exists $args->{'-distance'};
    $args->{'-difference'} = $args->{'-distance'}*(1.5) unless exists $args->{'-difference'};
    $args->{'-restrict'} = [ keys %HOMOLOGY ] unless exists $args->{'-restrict'};
    # 
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};

    # 
    
    if ( $args->{'-mode'} =~ /^p/i ) {

	my $genomesize = ( ($#{$args->{'-restrict'}} == 0 && 
			    $args->{'-restrict'}->[0] eq 'YGOB') ? 4704 : 5500 );
	my $whiteballs = $args->{'-difference'}*2;
	my $blackballs =  $genomesize - $whiteballs;
	my ($whitedraws,$totaldraws) = $self->synteny(%{$args}); # total based on -distance	
	
	my $dhyper = _hypergeom( $whiteballs, $blackballs, $totaldraws, $whitedraws );
	print $whiteballs, $blackballs, $totaldraws, $whitedraws, $dhyper if $args->{'-verbose'};
	return -1*(log($dhyper)/log(10)); # convert to positive score 
	
    } elsif ( $args->{'-mode'} =~ /^s/i ) {
	
	my (%raw, $average);
	foreach my $sp (grep {defined $self->data($_) } @{$args->{'-restrict'}} ) {
	    my $syn = $self->_hypergeometric_synteny_score(-restrict => [$sp]);	    
	    next unless $syn > 1;
	    $raw{$sp} = $syn;
	    $average += $syn;
	}
	return 0 unless %raw;

	$average/=scalar(keys %raw);
	return (wantarray ? ($average,\%raw) : $average );

    } 
    $self->throw("Must specify mode.");
}

=head2 _logodds_synteny_score(-distance => 3, -restrict => [! qw(MITO SPOM YGOB)]) 
    
    Synteny method that composes a neighbours x species
    matrix around the gene of interest and computes the 
    log-odds score of the observed data under a model of 
    synteny and a model of no synteny (based on randomizations).
    Following BLAST we consider every position in the matrix to 
    be independent and sum the log-odds scores over all 
    positions. It does not account for the phylogenetic relationships 
    between species (this can be moderated with -restrict).

    In future we will distinguish and automatically choose between 
    subtelomeric and chromosome internal models since these
    regions exhibit radically different synteny regimes. 

=cut 

sub _logodds_synteny_score {
    my $self = shift;
    my $args = {@_};   
    
    $args->{'-distance'} = 3 unless exists $args->{'-distance'};
    $args->{'-restrict'} = [ grep {!/MITO|SPOM|YGOB/} keys %HOMOLOGY ] 
	unless exists $args->{'-restrict'};
    #
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};

    $self->throw if exists $args->{'-difference'};
    
    $self->throw unless my $genome = $self->up->up;
    
    # syn can be 0 
    # but undef means that scoring matrices not available 
    return undef unless defined $genome->_access_synteny_null_dist(); 
    
    return $genome->_synteny( -object => $self, %{$args} );
}

sub _compose_synteny_delta_matrix { 
    my $self = shift;
    my $args = {@_};

    #$args->{'-restrict'} = [qw(KLAC KPOL SKLU EGOS SGD)] unless exists $args->{'-restrict'};
    #$args->{'-distance'} = 10 unless exists $args->{'-distance'};
    $args->{'-normalize'} = 1 unless exists $args->{'-normalize'};

    # 
    
    my %restrict = map { $_ => 1 } (
	ref($args->{'-restrict'}) eq 'ARRAY' ? 
	@{ $args->{'-restrict'} } :
	keys %HOMOLOGY
    );
    
    if ( $args->{'-weights'} ) {
	$self->throw;
	$self->throw unless ref($args->{'-weights'}) eq 'HASH';
	map { $self->throw unless exists $restrict{$_} } keys %{$args->{'-weights'}};
	map { $restrict{$_} = $args->{'-weights'}->{$_} } keys %{$args->{'-weights'}};
    }
    
    # scoring policy is to penalize ($INF) for homologs that are undefined in the caller
    # but to ignore (undef) if no homolog in comparator

   my @mat;
    foreach my $orf ( $self->context( -distance => $args->{'-distance'}, -self => 1, -variant => 0) ) {
	my $hash = { OBJ => $orf };
	foreach my $key ( keys %restrict) { # grep {$self->data($_)} 	
	    if ( my $other = $orf->data($key) ) {
		if ( my $ref = $self->data($key) ) {
		    $hash->{$key} = _compute_gene_distance($ref,$other);
		    #($other ne $ref ? _compute_gene_distance($ref,$other) : undef);
		} else { $hash->{$key} = $INFINITY; }
	    } else { $hash->{$key} = undef; } 
	    print $key, $self->name, $orf->name, $ref, $other, $hash->{$key} if $args->{'-verbose'} >= 3;
	}
	push @mat, $hash;
    }
    
    # pad matrix if at telomere 
    
    foreach my $i ( 0..$#mat ) {
	if ($mat[$i]->{OBJ} eq $self ) {
	    my $left = $args->{'-distance'} - $i;
	    my $right = $args->{'-distance'} - ( $#mat - $i );
	    push @mat, (undef) x $left;
	    unshift @mat, (undef) x $right;
	    last;
	}
    }
    $self->throw() unless $#mat == $args->{'-distance'}*2;
    
    return \@mat unless $args->{'-normalize'};

    # subtract expectation (assumes coding only genome but contrls assume same)
    
    foreach my $i ( 0..$#mat ) {
	my $exp_delta = abs( $i - $args->{'-distance'});
	foreach my $sp ( grep {!/OBJ/} keys %{$mat[$i]} ) {
	    next if $mat[$i]->{$sp} == $INFINITY || ! defined $mat[$i]->{$sp}; 
	    $mat[$i]->{$sp} -= $exp_delta;
	}
    }

   return \@mat;
}

=head2 synteny_conserved(-object => $orf, -homology =>'YGOB', 
    -genome => 4703, -distance => 10) 
    
    Returns -log10(Phyper) for the number of conserved neighbours
    between caller and object. Uses specified # of genes up- and 
    down-stream and assumes a genome size of -genome. 
    
    Hypergeometric probability density is a pure perl 
    implementation borrowed from perlmonks. 

=cut 

sub synteny_conserved {
    my $self  = shift;
    my $args = {@_};

    $args->{'-genome'} = 4703 unless exists $args->{'-genome'};
    $args->{'-distance'} = 10 unless exists $args->{'-distance'};
    $args->{'-homology'} = 'YGOB' unless exists $args->{'-homology'};

    # 

    my ($obj, $homolog) = ($args->{'-object'}, $args->{'-homology'});
    $self->throw unless $self->isa(ref($obj));
    
    # make data arrays 

    my %hash;
    foreach my $set ( [$self, 'Q'], [$obj, 'S'] ) {
	foreach my $dir ('left', 'right') {
	    my ($x, $y) = @{$set}; 
	    my $key = join('_', uc($y), uc($dir));
	    next unless $x = $x->$dir;
	    until ($#{$hash{$key}} == $args->{'-distance'}-1) {
		push @{$hash{$key}}, $x if ( $x->rank > -1 && $x->assign ne 'REPEAT' );
		last unless $x = $x->$dir;
	    }
	}
    }

    my %count;
    foreach my $q ('Q_LEFT', 'Q_RIGHT') {
	foreach my $s ('S_LEFT', 'S_RIGHT') {
	    my $ori = ( (($q =~ /LEFT/ && $s =~ /LEFT/) || ($q =~ /RIGHT/ && $s =~ /RIGHT/)) ? 1 : -1);
	    my ($max) = sort {$a <=> $b} (
		scalar(grep {$_->data($homolog) =~ /\w/} @{$hash{$s}}), 
		scalar(grep {$_->data($homolog) =~ /\w/} @{$hash{$q}})
	    );
	    $count{$ori}{DENOM}+=$max;

	    # do comparisons 

	    my %uniq; # screen out tandems 
	    foreach my $i (grep {!$uniq{$_->data($homolog)}++} 
			   grep {$_->data($homolog) =~ /\w/} @{$hash{$q}}) {
		my %uniq2;	    
		foreach my $j (grep {!$uniq2{$_->data($homolog)}++}  
			       grep {$_->data($homolog) =~ /\w/} @{$hash{$s}}) {
		    $count{$ori}{NUMER}+=( $i->data($homolog) eq $j->data($homolog) ? 1 : 0 );
		}
	    }	    
	}
    }

    $self->throw if $count{'-1'}{NUMER} > $count{'-1'}{DENOM};

    # where are the stats mo-fo?

    $count{'1'}{DENOM}=1 unless $count{'1'}{DENOM};
    $count{'-1'}{DENOM}=1 unless $count{'-1'}{DENOM};
    
    my $ori;
    if ($count{'1'}{NUMER} / $count{'1'}{DENOM} > $count{'-1'}{NUMER} / $count{'-1'}{DENOM} ) {
	$ori = 1;
    } elsif ($count{'1'}{NUMER} / $count{'1'}{DENOM} < $count{'-1'}{NUMER} / $count{'-1'}{DENOM} ) {
	$ori = -1;
    } else {return (0);}
    
    my $dhyper = _hypergeom( 
	$count{ $ori }{DENOM}, # white balls 
	($args->{'-genome'} - $count{ $ori }{DENOM}),  # black balls 
	$count{ $ori }{DENOM}, # draws 
	$count{ $ori }{NUMER}  # drawn white balls 
	);
    
    return sprintf("%.1f", -1*(log($dhyper)/log(10)) );
}

=head2 density( -window => 7 ) 
    
    Wrapper for ancestralSyntenyDensity that returns the scaled (0-1)
    value only i.e. the density of syntenic ancestral homologs in the 
    region. 

=cut 

sub density { my $self = shift; return sprintf("%.2f", $self->ancestralSyntenyDensity( @_ )+0); }

=head2 ancestralSyntenyDensity( -window => 10 )

    Another synteny scoring method. Requires only a single
    parameter and returns nice smooth values from 0 - 1. 

    Unlike other synteny methods it is fairly graceful in 
    how it handles missing data (e.g. at contig ends etc).  
    Should probably be merged with synteny() method and used
    to improve synteny(). 
    
=cut 

sub ancestralSyntenyDensity {
    my $self = shift;
    my $args = {@_};

    $args->{'-window'} = 10 unless exists $args->{'-window'};   
    
    $self->throw unless my $win = $args->{'-window'};
    
    #######################################
    # 
    #######################################

    return 0 unless $self->ygob;
    $self->ygob =~ /Anc_(\d+)\.(\d+)/ || $self->throw; 
    my $ref_chr = $1;
    my $ref_q = $2;

    #######################################
    # 
    #######################################

    my $factor = ( $self->up->up->wgd ? 2 : 1 );
    my $limit = $win*$factor;
    my $min = $limit*($win + 1);
    my $max = $limit*$win*2;
    my $scaleF = $limit*($win - 1); # max - min
    my @best = map { $_*$factor } (1..$win);
    
    #######################################
    # 
    #######################################

    my $raw;
    my @raw;
    foreach my $dir ( qw(left right) ) {

	my @context = grep {$_->ygob} $self->context(
	    -distance => $args->{'-window'},
	    -direction => $dir,
	    -self => -1
	    );
	splice(@context, 0, $#context, reverse(@context)) if $dir eq 'left';

	# compute distances 

	my @dist;
	foreach my $pos ( 1..$win ) {
	    my $index = $pos-1;	    
	    my $dist=undef;
	    if ( exists $context[$index] ) {
		$context[$index]->ygob =~ /Anc_(\d+)\.(\d+)/ || $self->throw;
		my ($chr,$q)=($1,$2);
		$dist = ( ($ref_chr eq $chr) && (abs($q-$ref_q) <= $limit) ? abs($q-$ref_q) : $limit);
		$dist = $best[$index] if $dist < $best[$index]; # && $factor==1; # non-WGD cannot be better than best
	    } 
	    push @dist, $dist;
	}

	# compute expected ratio across all actual data 

	my @exp;
	foreach my $pos ( 1..$win ) {
	    my $index = $pos-1;	
	    next unless defined $dist[$index];
	    push @exp, ( $dist[$index] > $best[$index] ? $dist[$index]/$best[$index] : 1);
	} 
	my ($exp) = ( $#exp>1 ? &_calcMeanSD( @exp ) : undef);

	# estimate the missing distances based on observed data 

	foreach my $pos ( 1..$win ) {
	    my $index = $pos-1;	
	    next if defined $dist[$index];
	    $dist[$index] = $limit and next unless $exp;
	    $dist[$index] = ( $exp*$best[$index] < $limit ? $exp*$best[$index] : $limit);	    
	}
	map { $raw += $_ } @dist;
	push @raw,@dist;
	#print $dir, sprintf("%.2f",$exp), @dist;
    }

    # normalize to a linear score 

    my $norm = $max - $raw; # higher is better 
    my $scaled = $norm / $scaleF;
    
    #print $self->name, map { sprintf("%.1f", $_) } ($min, $raw, $max, $norm, $scaled);
    
    return ( wantarray ? ($scaled,$norm,$raw) : $scaled );
}

##########################################
##########################################

=head2 quality

    Return orthogroup quality metrics (Z-scores) 
    by calling genome method and comaring caller
    to genomewide statistics. 

    Returns:
    Percentile, overall quality (synteny+homology), synteny, homology.
    
=cut 

sub quality {
    my $self = shift;
    my $args = {@_};

    return undef unless $self->orthogroup;
    
    unless ( $args->{'-calculate'} || $args->{'-force'} ) {
	return (wantarray ? @{$self->{QUALITY}} : $self->{QUALITY}->[0]) if 
	    defined  $self->{QUALITY} && ref($self->{QUALITY}) =~ /ARRAY/;   
    }
    
    my $G = $self->up->up;
    $self->{QUALITY} = [$G->quality(-object => $self)];
    return (wantarray ? @{$self->{QUALITY}} : $self->{QUALITY}->[0]);
}

=head2 direction( -object => orf, -numeric => 0|1)
    
    Return direction of object relative to caller.
    
    By default return 'left' or 'right' but with 
    -numeric switch will return -1 and 1 respectively. 

=cut

sub direction {
    my $self = shift;
    my $args = {@_};

    my %map = ( 'right' => 1, 'left' => -1 );

    $self->throw unless $self->isa( ref( $args->{'-object'} ) );
    $self->throw unless $self->up eq $args->{'-object'}->up;
    $self->throw if $self eq $args->{'-object'};

    my $dir = ( $self->index < $args->{'-object'}->index ? 'right' : 'left' );
    
    return ( $args->{'-numeric'} ? $map{$dir} : $dir );
}

=head2 codons(-query => 'coding|stop')
    
    Return codon usage frequencies as an array.
    Coding returns coding freqs only (ie no STOPs).

=cut 

sub codons {
	my $self = shift;
    my $args = {@_};

    $args->{'-query'} = 'coding' unless exists $args->{'-query'};
    my $seq = $self->sequence;

	my (%hash,$count);
    for (my $x = $TRIPLET; $x <= length($seq)-$TRIPLET; $x += $TRIPLET) {
        my $codon = substr($seq,$x,$TRIPLET);                           
        $self->throw("Codon not 3 nt long!") unless 
            length($codon) == $TRIPLET;
            
        if ($args->{'-query'} =~ /stop/i) {
        	next unless $codon =~ /TAG|TAA|TGA/;
        } else {next if $codon =~ /N|TAG|TAA|TGA/;}

        $hash{$codon}++;                   
        $count++;
    }

    return $count if $args->{'-query'} =~ /stop/i;
    
    my @r;
    foreach my $c (keys %CODONS) {
        next if $c =~ /N|TAG|TAA|TGA/;        
        if (exists $hash{$c}) {   
            push @r, $hash{$c}/$count;
        } else {push @r, 0;}
    }
    
    $self->throw("wrong number of codons") unless $#r == 60;
    return \@r;     
}


#########################################
#########################################

# evidence evaluation methods  

=head2 evaluate()

    Evaluate accumulated evidence according to the hierarchy 
    defined in %EVIDENCE and update the inferred gene type if
    required. 

=cut

sub evaluate {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-structure'} = 1 unless exists $args->{'-structure'};
    $args->{'-validate'} = 1 unless exists $args->{'-validate'};
    $args->{'-force'} = 0 unless exists $args->{'-force'};

    ############################################################
    # test structure and translatibility for protein coding 
    ############################################################

    unless ( $self->rank < -1 ) {
	$self->structure unless $args->{'-structure'} == 0;
	($self->output && $self->throw("Bad structure!\n".$self->debug))
	    unless $args->{'-validate'} == 0 || $self->translatable;
    }

    ############################################################
    # go thourgh the evidence in defined order. 
    ############################################################

    $self->evidence('NONE') if $args->{'-force'}; # hard reset 

    foreach my $evidence ( @EVIDENCE ) {

	##############################
	# these should all be removed.
	# rank based tests rely on existing structure however. 	
	if ($evidence eq 'INTRONS') { 
	    next unless $self->data('STRUCT') == # only relevant if
		$EVIDENCE{'STRUCT'}->{'EVAL'};   # structure valid 
	} elsif ( $evidence eq 'LENGTH' && ! $self->data('MANUAL') ) {
	    next unless $self->translatable;
	}
	##############################

	# this is scary far from a properly integrated model.... 
	
	##############################
	# evalue is basis for test and is stored on _KEY for homology based evidence.
	# all others stored on KEY itself. 
	# this is really awful... at the very least should modify
	# %EVIDENCE so that evidence knows if it is of type
	# homology or not...
	my $test_key = ($evidence =~ /AA|NCBI|TY|LTR|YGOB/ ? '_'.$evidence : $evidence);
	# NB: HSP does not contain an E-value but a ratio (from ->fragment() ?)
	##############################

	# core comparisons 
	
	if ($EVIDENCE{$evidence}->{'TEST'} eq '=') {
	    next unless $self->data($test_key) == $EVIDENCE{$evidence}->{'EVAL'};
	} elsif ($EVIDENCE{$evidence}->{'TEST'} eq '>') {		
	    next unless $self->data($test_key) >= $EVIDENCE{$evidence}->{'EVAL'};
	} elsif ($EVIDENCE{$evidence}->{'TEST'} eq '<') {
	    next unless $self->data($test_key) <= $EVIDENCE{$evidence}->{'EVAL'};
	} else {$self->throw("Evidence unasociated with test: $evidence");}

	# this is the big step.. 
	$self->evidence($evidence) if $EVIDENCE{$evidence}->{'RANK'} <= $self->rank;
    }
    
    ############################################################
    # update evidence/assignment/rank 
    ############################################################

    # evaluate is the only method that can set evidence so. 
    # we can just set these here. we call rank/assign A LOT so it is taking up cycles. 

    $self->{ASSIGN}=$self->_evidence( -evidence => $self->evidence, -query => 'infer' );
    $self->{RANK}=$self->_evidence( -evidence => $self->evidence, -query => 'rank' );

    return($self->assign);     
}

=head2 structure(-exons => 1|0)

    Returns a numeric value based on whther the gene has good structure: 
    2: M,* 
    1: * 
    -1: M
    0 

    Also recalculates the number of introns, frameshits etc by querying exons

=cut

sub structure {
	my $self = shift;
	my $args = {@_};
	
	$args->{'-exons'} = 1 unless exists $args->{'-exons'};	
	
	if ($self->rank < -1) {
		$self->data('STRUCT', 0); 		
		$self->data('STOP', 0);
		$self->data('INTRONS', $self->exons -1);
		return 0;		
	}
	
	# get START and STOP 
	
	my ($start, $stop) = $self->_top_tail;
	$start = $CODONS{$start};
	$stop = $CODONS{$stop};
	
	# this scale is not very good but probably workable. 	
	# tests should call ->structure but the actual value can be 
	# accessed as usual via ->data

	if ($start eq 'M' && $stop eq '*') {
		$self->data('STRUCT', 2); 
	} elsif ($stop eq '*') {
		$self->data('STRUCT', 1); 	
	} elsif ($start eq 'M') {
		$self->data('STRUCT', -1); 
	} else {
		$self->data('STRUCT', 0); 
	}

	# 
	
	$self->data('LENGTH' => $self->length);

	# go through exons and ensure intron boundaries are consistent
	# recalculate the number of frameshifts and the number of
	# real introns and update 
	
	unless ($args->{'-exons'} == 0) {
		my ($in, $fr, $gap) = (0,0,0);
		foreach my $ex ($self->stream) {
			next unless my $left = $ex->left;
			my $ir = $left->intron(-direction => 'right');
			my $il = $ex->intron(-direction => 'left');					
			
			if ($il == $INFINITY || $ir == $INFINITY) {
				$gap++;
			} elsif ($il == 0 && $ir == 0) {
				$fr++;
			} elsif (($il > 0 && $il < $INFINITY) && 
  					 ($ir > 0 && $ir < $INFINITY)) {	
				$in++;			
			} else {$self->throw("Hybrid introns not allowed: $il, $ir");}
		}

		$self->throw("Exons do not validate: $in, $fr".$self->exons)
			unless $self->exons-1 == $in+$fr+$gap;

		$self->data('STOP', $fr);
		$self->data('INTRONS', $in);
		$self->data('_GAP', $gap);		
	}
	
	return abs($self->data('STRUCT'));
}

=head2 bounded

    Retrun TRUE (self) if orthologs of left and right neighbours 
    are neighbours of orthologs in all species. 

=cut 

sub bounded {
    my $self = shift;
    my $args = {@_};
    
    return undef unless $self->orthogroup;

    my ($lx,$rx) = $self->neighbours;
    return undef unless $lx->orthogroup && $rx->orthogroup;
    
    foreach my $og ($self->orthogroup) {
	my ($l,$r) = $og->neighbours;
	return undef unless $l && $r;
	my $meth = $og->organism;
	return undef unless ($l eq $lx->$meth || $r eq $lx->$meth) &&
	    ($l eq $rx->$meth || $r eq  $rx->$meth);
    }
    
    return $self;
}

=head2 lookup(-query => , -file => , -species => )

    Use DB (YGOB pillar file) to look up orthologs for query. 
    Returns gene name(s) if a species specified else hash of Pillar data. 

    Use -file to over ride the default database. 

=cut 

sub lookup {
    my $self = shift;
    my $args = {@_};

    $args->{'-species'} = undef unless exists  $args->{'-species'};
    $args->{'-file'} = undef unless exists $args->{'-file'};
    $args->{'-query'} = $self->pillar unless exists $args->{'-query'};
    
    $self->throw unless defined $args->{'-query'};
    $self->throw unless ! $args->{'-file'} || -e $args->{'-file'};
    $self->throw unless ! $args->{'-species'} || exists $HOMOLOGY{ $args->{'-species'} };

    $YGOB = YGOB->new( $args->{'-file'} ) unless $YGOB;

    return $YGOB->access($args->{'-query'}, $args->{'-species'});
}

=head2 pillar(-query => undef|species, -distance => 10, -evalue => 1e-10,
    -force => 0|1)

    Identify the YGOB Pillar - and hence homologs - for use in other 
    methods. We use the log-ods synteny score (LOSS) to choose among
    a panel of pillars implicated by pre-computed BLAST results.

    NB: retuns undef if LOSS not available.

    -distance : compute LOSS based on +/- distance genes
    -force : force re-compute. else we return existing value if it exists. 
    
=cut 

sub pillar {
    my $self = shift;
    my $args = {@_};

    # pillar uses LOSS matrices to choose among candidates. 
    # test if available 
    return undef unless defined $self->up->up->_access_synteny_null_dist();

    # mode 
    $args->{'-query'} = undef unless exists $args->{'-query'};
    $args->{'-query'} = shift if $#_ == 0;
    # params 
    $args->{'-distance'} = 10 unless exists $args->{'-distance'};
    $args->{'-evalue'} = 1e-5 unless exists  $args->{'-evalue'};
    # generic 
    $args->{'-force'} = 0 unless exists $args->{'-force'};  
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    
    #

    my $fh = *STDOUT;
    my $key = 'PILLAR'; # assumed by score('pillar')

    # just want to access data? great! 
    
    if ( $args->{'-query'} ) {
	$self->throw unless exists $HOMOLOGY{ $args->{'-query'} };	
	return undef unless my $id = $self->pillar;
	return undef unless my $pil = $self->lookup( -query => $id );
	return undef unless exists $pil->{ $args->{'-query'} };
	return (grep { defined } @{$pil->{ $args->{'-query'} }});
    } elsif ( ! $args->{'-force'} && $self->data($key) ) {
	return ( wantarray ? ($self->data($key) ,$self->data( '_'.$key)) : $self->data($key) );
    }
    
    # initialize... 

    $self->data( $key => undef );
    $self->data( '_'.$key => undef);

    ###################################################    
    # specifically include GENE. chosen by impute. 
    # single best candidate. 
    ###################################################

    my %count;
    if ( my $gene = $self->gene ) {
	my $hash = $self->lookup( -query => $gene );
	#print $gene, $hash, $hash->{'ID'};
	if ( $hash && $hash->{'ID'} ) {	
	    $self->throw unless my $sp = _species_key( $gene );
	    my $store = $self->data($sp);
	    $self->data($sp => $gene);
	    my $loss = $self->synteny( -lod => 1, -distance => $args->{'-distance'}, -restrict => [$sp] );	    
	    $self->data($sp => $store);    
	    $hash->{LOSS} = $loss;
	    $count{ $hash->{ID} } = $hash;
	}
    }
    
    ###################################################
    # use BLAST homology and sort by synteny. 
    # we generate as many candidates as possible 
    ###################################################
    
    foreach my $sp ( grep {!/MITO|SPOM/} keys %HOMOLOGY ) {
	foreach my $mod ( '*', undef ) {
	    next if $self->data( $sp ) eq $self->data( $sp.'*' ) && $mod eq '*';
	    
 	    next unless my $hit = $self->data( $sp.$mod );
	    next unless $self->data( '_'.$sp.$mod ) <= $args->{'-evalue'};
	    next unless my $loss = 
		$self->synteny( -lod => 1, -distance => $args->{'-distance'}, -restrict => [$sp] );		

	    my $hash = $self->lookup( -query => $hit );
	    #print $hit, $hash, $hash->{'ID'};
	    if ( $hash && $hash->{'ID'} ) {
		$count{ $hash->{'ID'} } = $hash unless exists $count{ $hash->{'ID'} };	
		$count{ $hash->{'ID'} }->{'LOSS'} += $loss;
	    }
	}
    }
    return undef unless %count;

    # 

    my ($first,@others) = sort { $b->{LOSS} <=> $a->{LOSS} } values %count;
    my $delta = sprintf("%.2f", $first->{LOSS} - (@others ? $others[0]->{LOSS} : 0));
    my $zed = ($first->{SD}>0 ? ($self->length - $first->{MEAN}) / $first->{SD} : 'NA');
    #print $first, $delta, $zed;

    # 
    
    #if ( $delta > $args->{'-delta'} || ($first->{SD}>0 && abs($zed) < $args->{'-sd'} ) ) {
    $self->data( $key => $first->{'ID'} );
    $self->data( '_'.$key => $delta );
    #}
    
    # 
    
    my @sort = sort @{$first->{SGD}};
    print {$fh} $self->_method, $self->name, $self->gene, @sort,
    $self->length, (map { sprintf("%.2f",$_) } ($first->{MEAN}, $first->{SD}, $zed) ),
    $self->data( '_'.$key ), $self->data( $key ) if $args->{'-verbose'} >= 1;
    
    # 

    return ( wantarray ? ( $self->data( $key ), $self->data( '_'.$key)) :  $self->data( $key ) );
}

=head2 pillarscore

    Return delta between the best pillar choice and the second best. 
    Value is in LOSS units and significance therefore depends on the number 
    of species considered. 50 is a relatively safe cutoff. 

  NBB: This is not conventionally consistetn with other 
    evidence / score systems. Currently accessed through 
    $o->score('evidence') and calling __EVIDENCE.
    By convention we are storing the P-/E-values on _EVIDENCE 
    and the score on __EVIDENCE. Need to unify all evidence. 

=cut 

sub pillarscore { my $self = shift; return $self->data( '_PILLAR' ); } 

=head2 chromosome( -species => '' ) 

    Guess the homologous chromosome in another species
    by looking at (+/-5) linked genes.
    
    We do not attempt to assess uncertainty.

    If the caller requests an array we also return the number
    of genes on the chr, the mean index on the chr and the sd.

    Values are stored as ORG_CHR .

=cut 

sub chromosome {
    my $self = shift;
    my $args = {@_};

    $args->{'-species'} = $_[0] if $#_ == 0;
    $args->{'-force'} = undef unless exists $args->{'-force'};
    $args->{'-distance'} = 10 unless exists  $args->{'-distance'};
    $args->{'-self'} = -1 unless exists $args->{'-self'};
    $self->throw() unless exists $HOMOLOGY{ $args->{'-species'} };
    
    ###############################################
    # are these variables sanctioned? 
    ###############################################

    my $attr = $args->{'-species'}.'_CHR';

    ###############################################
    # look for a fast return 
    ###############################################

    return $self->data($attr) unless 
	(wantarray || $args->{'-force'} || ! $self->data( $attr ));

    ###############################################
    # look up data 
    ###############################################
    
    # by default ->context() excludes variants/tRNAs/GAPs/features
    # we also exclude REPEATs since they will lead to errors in 
    # chromosome guessing. 

    my %chr;
    map { push @{$chr{$_->[1]}}, $_->[2] } grep {defined $_->[0]} 
    map { [ _decompose_gene_name($_) ]  }  grep {defined}
    map { $_->data($args->{'-species'}) } 
    $self->context(
	-distance => $args->{'-distance'}, 
	-repeat => -1, 
	-self => ($args->{'-self'}==1 ? 1 : 0)
	);

    my ($max) = sort { $#{$chr{$b}} <=> $#{$chr{$a}} } keys %chr;
    return undef unless $#{$chr{$max}} >= 0; # 2 is the min to call 

    ###############################################
    # compute 
    ###############################################
    
    my @sort = sort {$a <=> $b} @{$chr{$max}};    
    my @meansd = ( $#sort >= 2 ? @sort[ 1..($#sort - 1) ] : @sort );
    my ($mean,$sd) = _calcMeanSD( @meansd ) ; # trimmed mean 

    ###############################################
    # 
    ###############################################
    
    $self->data( $attr => $max );
    
    return (wantarray ? ($max, $mean, $sd, scalar(@{$chr{$max}}) ) : $max);
}


=head2 telomere( -bases => 1 ) 

    Returns number of genes (or array if requested) between
    the caller and the closest contig end. With -bases => 1
    returns bp instead. 

=cut 

sub telomere {
    my $self=shift;
    my $args = {@_};

    if ( $args->{'-bp'} || $args->{'-bases'} ) {
	my ($min) = sort {$a <=> $b} ($self->start, $self->up->length - $self->stop);
	return $min;
    }
    
    my (@l,@r);
    if ( my $next = $self->left ) {
	unshift @l, $next while ( $next = $next->left );
    }
    if ( my $next = $self->right ) {
	push @r, $next while ( $next = $next->right );
    }
    my @min = ( $#l > $#r ? @r : @l);

    return ( wantarray ? @min : scalar(@min) );
}

=head2 subtelomeric
    
    Does no computation, just reports results of the 
    $genome->subtelomeres() method. Returns undef for 
    non-subtelomeric genes and left|right for others. 
    
=cut 

sub subtelomeric {
    my $self=shift;    
    return ( $self->data('_SUBTEL') =~ /SUBTEL/ ?  $self->data('_SUBTEL') : undef );
}

=head2 introns

    Return introns as a stream.

=cut 

sub introns {
    my $self = shift;
    my $args = {@_};
    
    return unless $self->exons>1;

    my @int;
    foreach my $ex ( grep {$_ ne $self->fex} $self->stream ) {
	my $ls = $ex->intron( -direction => 'left');
	my $rs = $ex->left->intron( -direction => 'right');
	next unless $ls >=1  && $rs >= 1;
	next unless my $int = $ex->intron( -object => $ex->left );
	push @int, $int if $int->length => 30;
    }

    return @int;
}

sub _species_key {
    print caller() and die($gene) 
	unless my $gene = shift;

    if ($gene =~ /([YA])[A-P][LR]\d{3}.+/) { # SGD 
    } elsif ($gene =~ /([A-Z]{4})\d+[A-Z]\d+[grst]?/) { # Genolevures 
    } elsif ($gene =~ /(\w{3,4})_\d+\.\d+[a-z]?/) { # Standard
    } else { return undef; } 
    #} else { die($gene); } 
    
    my $sp = uc($1);

    return( exists $HOMOLOGY{$sp} ? $sp : $SPECIES{$sp} );
}

sub _decompose_gene_name {
    die unless my $gene = shift;
    
    my ($sp,$chr,$index,$other,$arm,$strain)=( 6 x undef);
    if ($gene =~ /([YA])([A-P])([LR])(\d{3})(.+)/) { # SGD 
	($sp,$chr,$index,$other,$arm,$strain)=( $1, $2, $4, $5, $3, undef);
    } elsif ($gene =~ /([A-Z]{4})(\d+)([A-Z])(\d+)([grst])?/) { # Genolevures 
	($sp,$chr,$index,$other,$arm,$strain)=( $1, $3, int($4/11), $5, undef,$2);
    } elsif ($gene =~ /(\w{3,4})_(\d+)\.(\d+)([a-z])?/) { # Standard
	($sp,$chr,$index,$other,$arm,$strain)=( $1, $2, $3, $4, undef, undef);
    #} else { return undef; } 
    } else { print (caller(1)); die($gene); } 

    $sp = uc($sp);
    $sp = $SPECIES{$sp} unless ( exists $HOMOLOGY{$sp} );

    return($sp,$chr,$index,$other,$arm,$strain);
}

sub _compute_gene_distance {
    my ($x,$y) = @_;
    return undef if $y =~ /_YGOB_/ || $x =~ /_YGOB_/;

    die("$x,$y") unless $x && $y;
    my @x = _decompose_gene_name($x);
    my @y = _decompose_gene_name($y);

    die("$x,$y,$x[0],$y[0]") unless $x[0] && $x[0] eq $y[0];
    return $INFINITY unless $x[1] eq $y[1];
    
    my $delta;
    if ( defined $x[4] && $x[4] ne $y[4] ) { # arm 
	$delta = $x[2] + $y[2] -1; # counts from 0 not 1 
    } else {
	$delta = abs($x[2] - $y[2]);
    }
    
    # we ignore suffixes -- no way to order. 

    return $delta;
}

sub _compute_average_synteny_distance {
    my $self = shift;
    my $other = shift;
    $self->throw unless $self->isa( ref($other) );

    my @delta;
    foreach my $sp ( keys %HOMOLOGY ) {
	next unless $self->data($sp) && $other->data($sp);
	next unless $self->data('_'.$sp) <= 1e-5 && $other->data('_'.$sp) <= 1e-5;
	push @delta,  grep {defined} _compute_gene_distance( $self->data($sp), $other->data($sp) );
    }
    return undef unless @delta;

    my @sort = sort {$a <=> $b} @delta;
    my $median = $sort[ ($#sort%2==0 ? $#sort/2 : ($#sort-1)/2 ) ];
    
    #print @sort;

    return $median;
}

=head2 _loss( $gene1,$gene2, $offset )
    
    Return the log-odds synteny score (LOSS) for 
    a pair of gene *names*. Offset (number of separating)
    genes in the alignemnt is required to calculate.     

    Do not call directly. 

=cut 

sub _loss {
    my $self = shift;
    my ($g1,$g2,$offset)=@_;
    
    $self->throw unless my $genome = $self->up->up;
    
    # $g1 ~~ @REGEX is broken 
    $self->throw("$g1,$g2") unless (grep {$g1 ~~ $_} @REGEX) && (grep {$g2 ~~ $_} @REGEX);
    $self->throw unless $offset > 0;
    
    # 

    my $sp = _species_key($g1);
    my $sp2 = _species_key($g2);
    $self->throw unless $sp && $sp eq $sp2;
    $sp = uc($sp);

    # 

    return undef unless my $d = _compute_gene_distance($g1,$g2);
    $self->throw unless my $limit = $genome->_access_synteny_null_dist('LIMIT');
    $d = $limit if $d > $limit;

    # 
    
    $self->throw unless my $null = 
	$genome->_access_synteny_null_dist('NULL_DIST', $offset, $sp, $d );
    $self->throw unless my $syn = 
	$genome->_access_synteny_null_dist('SYN_DIST', $offset, $sp, $d );    
    return log( $syn/$null );   
}

=head2 identify 
    
    Return the most common SGD (based on gene not SGD attribute) 
    and YGOB hits for an *orthogroup*. Former may be tRNA. 
    
    (SGD,YGOB,GENE) = $self->identify

=cut 

sub identify {
    my $self = shift;
    my $args = {@_};

    my @og = grep {defined} ($self, $self->orthogroup);

    my %gene;
    my %ygob;
    my %sgd;
    foreach my $o ( @og ) {
	$sgd{$o->data('GENE')}++ if $o->data('GENE') =~ /^$HOMOLOGY{'SGD'}|\w+\:\w+/;
	$ygob{$o->ygob}++ if $o->ygob =~ /^$HOMOLOGY{'YGOB'}/;
	$gene{$o->data('GENE')}++;
    }
    
    my ($sgd) = sort { $sgd{$b} <=> $sgd{$a} } keys %sgd;
    my ($ygob) = sort { $ygob{$b} <=> $ygob{$a} } keys %ygob;
    my ($gene) = sort { $gene{$b} <=> $gene{$a} } keys %gene;
    
    return (($sgd || undef) , ( $ygob|| undef), ( $gene || undef ) );
}


sub _genericlinks {
    my $self = shift;
    my ($gene) = grep { /^Y/ } $self->identify;
    return(
	'http://www.yeastgenome.org/cgi-bin/locus.fpl?locus='.$gene,
	'http://wolfe.gen.tcd.ie/cgi/browser/ygob.pl?gene='.$gene,
	'http://www.yeastgenome.org/cgi-bin/FUNGI/FungiMap?locus='.$gene
	);
}

=head2 fex/firstexon 

    Return first exon of gene.
    
=cut 

sub firstexon { return $_[0]->exons(-query => 'first'); }
sub fex { return $_[0]->firstexon; } 

=head2 dump

    Override method for parent dump().
    Print all values on on %{DATA}.

=cut 
sub dump { return $_[0]->output(-dump => 1); }

=head2 atg( -safe => 1, -mode => 'FIRST|last' )

    Return the *first* upstream in-frame ATG that does not 
    conflict with other models. 

    Returns new start coord in safe mode, else return modified self.

=cut 

sub atg {
    my $self = shift;
    my $args = {@_};

    $self->throw if $self->rank < -1;
    return $self if $self->rank == -1;

    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    $args->{'-safe'} = 1 unless exists $args->{'-safe'};
    $args->{'-mode'} = 'first' unless exists $args->{'-mode'};
    $self->throw("not implemented") unless $args->{'-mode'} eq 'first';

    # find a sensible upstream boundary at which to limit expansion
    # we want to stop at any upstream feature (ORF/TRNA etc) that is not a 
    # variant of the caller. stop for GAP?

    my $init = $self->exons(-query => 'first')->start(-R => 1); 
    my $dir = ( $self->strand == 1 ? 'left' : 'right');
    my $obs = $self->$dir; # obs = obstacle 
    $obs = $obs->$dir while $obs && $self->variants(-object => $obs) && $self ne $obs;

    # 

    my $clone = $self->clone() || die;
    my $fex = $clone->exons(-query => 'first');
    
    my ($top,$tail) = $clone->_top_tail;
    until ( $top eq 'ATG' || $top =~ /TAA|TAG|TGA/i ||
	    $fex->start(-R => 1) < $TRIPLET || 
	    ($fex->start(-R => 1) > ($self->up->length - $TRIPLET)) ||
	    ($obs && $self->overlap(-object => $obs))
	)  {
	
	$fex->start( -R => 1, -new => ($fex->start(-R=>1) + -1*$self->strand*$TRIPLET) );
	($top,$tail) = $clone->_top_tail;
	print $fex->start(-R=>1), $top if $args->{'-verbose'};
    }
    
    ($top,$tail) = $clone->_top_tail;
    my $new = ( $top eq 'ATG' ? $fex->start(-R => 1) : $init);

    print $self->_method, $self->name, $init, $new, $top,$tail 
	if $args->{'-verbose'};
    
    return $new if $args->{'-safe'};

    my $fex = $self->exons(-query => 'first');
    $fex->start(-R => 1, -new => $new );
    $self->throw unless $self->translatable;

    return $self;
}

=head2 atg_og(-reference => undef, -verbose => undef, -_comb_limit => 10)

    Optimise the location of the start Codon based on comparative data. 

=cut

sub atg_og {
    my $self = shift;
    my $args = {@_};
    
    return undef unless $self->orthogroup;

    $args->{'-_comb_limit'} = 10 unless exists $args->{'-_comb_limit'};
    $args->{'-reference'} = undef unless exists $args->{'-reference'};
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};

    my ($m,$s) = _calcMeanSD( map { $_->length } ($self, $self->orthogroup) );
    print "\n>", (map { $_->length } ($self, $self->orthogroup)), $m, $s if $args->{'-verbose'};
    
    ################################################################
    # clone genes and expand start position up to first STOP codon 

    my @clones = map { $_->clone } ($self,$self->orthogroup);

    foreach my $x ( @clones ) {
	my ($fex) = $x->exons(-query => 'first');
	my $init = $fex->start(-R => 1);
	#map { $x->remove(-object => $_, -warn => 0) } grep { $_ ne $fex } $x->stream;
	#$fex->stop(-R => 1, -adjust => -1 ) until $fex->length%3==0;
	$fex->morph(
	    -method => 'expand',
	    -step => 3,
	    -stop => '\*',
	    -terminus => '5'
	    );
	$fex->start(-R => 1, -adjust => 3) unless $fex->start(-R => 1) == $init;
    }
    map { return undef unless $_->aa =~ /M/ } @clones; # we require M 
    
    ################################################################
    # align expanded clones and parse 
    
    my $aln = 'data/'.$self->_internal_id.'_'.$self->gene;
    unless ( -e $aln && -s $aln > 100 ) {
        my $newaln = $clones[0]->align(-object => [@clones[1..$#clones]], -molecule => 'aa', -format => 'fasta');
	copy($newaln, $aln) || $self->throw("$newaln, $aln");
	#print `cat $aln`;
    }

    local $/ = ">";
    open(my $fh, $aln);
    <$fh>;
    my %hash;
    while ( my $s = <$fh> ) {
	chomp($s);
	my ($id, @r)=split/\n/, $s;
	$id =~ /^(\w+)/ || die($s);
	$hash{ lc($1) } = [split//, join('', @r)];
    }
    close $fh;
    local $/ = "\n";
    my $alnLen = $#{$hash{lc($self->organism)}};

    print $aln.'/'.$self->name.'/'.$self->gene.`cat $aln` and $self->throw
	unless ($alnLen+6) >= $clones[0]->length/3;
    map { $self->throw unless $#{$hash{lc($_->organism)}} == $alnLen } @clones;

    ################################################################
    # we define an acceptable window based on the initial predicted lengths.
    # we want to exclude perfectly aligned Ms at 3' of gene 

    my @start;
    foreach my $sp ( keys %hash ) {
	my ($or) = grep { $_->organism =~ /$sp/i } ($self,$self->orthogroup);
	my ($cl) = grep { $_->organism =~ /$sp/i } @clones;
	$self->throw("$sp $or $cl ") unless $or && $cl;
	$self->throw("$sp $or $cl ") if $or->down eq $cl->down;
	my $ext = ($cl->length - $or->length)/3;

	# 

	my $q=0;
	foreach my $i ( 0..$#{$hash{$sp}} ) {
	    push @start, $i and last if $q==($ext+1);
	    $q++ unless $hash{$sp}->[$i] eq '-';
	    #print $i, $q, $ext;
	}
    }
    $self->throw unless scalar(keys %hash) == scalar($self->orthogroup)+1;

    # define region around gene starts 

    my @sort =  sort { $a <=> $b } @start;
    my $median = $sort[ ($#sort%2==0 ? $#sort/2 : ($#sort+1)/2 ) ];
    my ($m,$s) = _calcMeanSD( @sort );
    my ($lo,$hi) = ( int($median-$s)-10 , int($median+$s)+10 );
    $lo = 0 if $lo < 0;
    $hi = $alnLen if $hi > $alnLen;
    print ">>$alnLen", @sort, $m, $s, $median, $lo, $hi if $args->{'-verbose'};

    # use min variance approach in case NO homologous M.
    # gather up to 10 M in each species -- 100K paths for stricto. 
    # combinatorics are a bitch otherwise.

    my %M;
    for my $i ( $lo .. $hi ) {
	foreach my $sp ( grep { $#{$M{$_}} < $args->{'-_comb_limit'} } keys %hash ) {
	    #print $i, $sp,  $hash{$sp}->[$i];
	    push @{$M{$sp}}, {M => $i, SP => $sp} if $hash{$sp}->[$i] eq 'M';		
	}
    }
    map { print `cat $aln` and return undef unless $#{$M{$_}} >= 0 } keys %hash;

    # make all possible combinations and compute summary stats 

    my @matrix;
    foreach my $ref_m ( @{$M{ lc($self->organism) }} ) {
	my @array=[$ref_m];
	foreach my $og ( $self->orthogroup ) {	      	    
	    foreach (my $i=$#array; $i >= 0; $i--) {
		my $stub = splice(@array, $i, 1);
		#print $ref_m->{M}, $og->organism, $i,$stub->[0]->{M};
		foreach my $og_m ( @{$M{ lc($og->organism) }} ) {
		    push @array, [ @{$stub}, $og_m ]
		}
	    }
	}

	foreach my $path ( @array ) {
	    $self->throw unless ref($path) eq 'ARRAY' && $#{$path}==scalar($self->orthogroup);
	    $self->throw("@{$path}") unless exists $path->[0]->{'M'};
	    my ($m1,$s1) = _calcMeanSD( map { $_->{M} } @{$path} );
	    my %data = map { $_->{SP} => $_->{M} } @{$path};
	    $data{MEAN} = $m1;
	    $data{SD} = $s1;
 	    push @matrix, \%data;
	}
    }
    $self->throw unless $#matrix >= 0;

    # choose favourite: lowest variance, longest align 

    $self->throw unless my ($best) = sort { $a->{SD} <=> $b->{SD} } @matrix;
    my @top = grep { $_->{SD} <= $best->{SD} } @matrix; 
    ($best) = sort { $a->{MEAN} <=> $b->{MEAN} } @top;
    
    # 

    foreach my $sp ( keys %hash ) {
	my @cut = grep {!/\-/} splice(@{$hash{$sp}}, 0, $best->{$sp});
	my $delta = 3*scalar(@cut);

	my ($cl) = grep { $_->organism =~ /$sp/i } @clones;
	my ($or) = grep { $_->organism =~ /$sp/i } ($self,$self->orthogroup);

	my $fex = $cl->exons(-query => 'first');
	my $fx2 = $or->exons(-query => 'first');
	
	if ( $delta >= $fex->length ) {
	    $cl->oliver;
	    $or->oliver;
	    $or->warn("Mulit-exon adjustments not implemented! $delta");
	    next;
	}

	$fex->start(-R => 1, -adjust => $delta);
	$or->remove(-object => $fx2, -warn => 0);
	$fex->transfer(-from => $cl, -to => $or, -warn => 0);
	$fx2->DESTROY;
	$or->index;
	$or->evaluate;

	if ($args->{'-verbose'}) {
	    print ">>>".$or->name,$best->{$sp}, $delta, $or->length, substr($or->aa, 0, 30);
	}
    }
    print if $args->{'-verbose'};
    return $self;

    # print out summary 
    # my ($m,$s) = _calcMeanSD( map { $_->length } ($self, $self->orthogroup) );
    # print '>>', (map { $_->length } ($self, $self->orthogroup)), $m, $s if $args->{'-verbose'};

    # return ? or do additional search against outliers ?
    # yileds very poor results in test. 

    my ($outlier, @others) = sort { abs($b->length - $m) <=> abs($a->length - $m) } ($self, $self->orthogroup);
    my ($m1,$s1) = _calcMeanSD( map { $_->length } @others );

    if ( $args->{'-reference'} && ! $self->alignedStartCodon && abs($outlier->length - $m1 ) > 2*$s1 ) {
	print '>>>>',$outlier->length,(map { $_->length } (@others)), $m1, $s1 if $args->{'-verbose'};
	$outlier->oliver(-append => ['OUTL']);
	$outlier->reoptimise(-reference => $args->{'-reference'} );
	$outlier->oliver(-append => ['FIX']);
	my ($m,$s) = _calcMeanSD( map { $_->length } ($self, $self->orthogroup) );
	print '>>>>', (map { $_->length } ($self, $self->orthogroup)), $m, "$s\n" if $args->{'-verbose'};
    }

    return $self;
}

sub _calcMeanSD {
    # return undef if $#_ <= 0;
    my ($mean, $sd);
    map { $mean += $_ } @_;
    $mean /= ($#_+1);
    map { $sd += ($mean-$_)**2 } @_;
    return (  sprintf("%.3f", $mean),  sprintf("%.3f", sqrt($sd/($#_+1))  ));
}

=head2 alignedStartCodon

    Returns true (self) if both conditions hold :
    1. all genes in OG begin with M 
    2. first -align AAs align without gaps. 
    
=cut 

sub alignedStartCodon {
    my $self = shift;
    my $args = {@_};

    return undef unless $self->orthogroup;
    my @og = ($self,$self->orthogroup);

    $args->{'-align'} = 10 unless exists $args->{'-align'};
    my $adj =  $TRIPLET*$args->{'-align'} - 1;

    map {return undef unless $_->exons(-query => 'first')->length >= $TRIPLET*$args->{'-align'} } @og;

    foreach my $seq ( @og ) {
	my ($atg) = $seq->_top_tail;
	return undef unless $atg eq 'ATG';
    }

    # prepare 5' clones 
    
    my ($head, @clones) = map { $_->clone } @og; 
    foreach my $cl ( $head, @clones ) {
	my $fex = $cl->down;
	map { $cl->remove(-object => $_, -warn => 0); $_->DESTROY; } grep { $_ ne $fex }  $cl->stream;
	my $start = $fex->start(-R => 1);
	my $stop = ($fex->strand == 1 ? $start + $adj : $start - $adj );
	$fex->stop(-R => 1 , -new => $stop);
    }

    # make an alignmente and count gaps 
    
    my $aln = $head->align(-object => \@clones, -molecule => 'aa', -format => 'hash') || die;
    map { $_->DESTROY} ($head,@clones);

    map {return undef if $_ eq '-'} map {@{$_}} values %{$aln};

    return $self;
}

=head2 pseudo()

    Call pseudogene using : 
    1. Fshift/STOP density (2 SDs from mean..)
    2. Ka/Ks ~ 1
    3. Kaessmann simulation based method (reEVOLVER)
    http://www.pnas.org/content/103/9/3220
    4. this wold be better but no implementation:
    http://mbe.oxfordjournals.org/content/19/1/110.long

    None of these account for assembly problems which seems to 
    be the greatest source of disruptions. 

=cut 

sub pseudo {
    my $self = shift;
    $self->throw;
}

=head2 fragment() 
    
    Test whether model meets the criteria for a full
    model or a fragment using the local - global score test.
    
    Set the FRAG evidence attribute. 

=cut 

sub fragment {
    my $self = shift;
    my $args = {@_};

    my $gs = $self->score('global');
    my $ls = $self->score('local');
    return undef unless defined $gs && defined $ls;
    return undef unless $ls;    
    
    return $self->data( 'HSP' => sprintf("%.2f", ($ls-$gs)/($ls+$gs)) );
}

=head2 partial()

    Compare length to lengths of YGOB homologs 
    and compute Z-score. 

=cut

sub partial {
    my $self = shift;
    my $args = {@_};

    my $key = 'FRAG';

    # 

    $self->data($key => undef); # unless $self->data($key);
    return undef if $self->rank < -1; # /FEATURE|RNA|GAP/;

    # too tolerant here?

    my $id;
    if ( $id = $self->pillar() ) {
    } elsif ( $id = $self->gene ) {
    } elsif ($id = $self->ygob ) {
    } elsif ( $self->evalue('aa') <= 1e-5 ) {
	$id = $self->data('AA');
    }

    return undef unless $id;

    # 
    
    return undef unless my $pil = $self->lookup( -query => $id );
    return undef unless $pil && $pil->{ID};
    return undef unless $pil->{MEAN};
    return undef unless $pil->{SD}>0;

    # 
    
    my $zscr = ($self->length - $pil->{MEAN} )/ $pil->{SD};
    $self->data($key => $zscr);
    
    print $self->name,$self->gene,$self->length, $pil->{ID}, $pil->{YGOB}->[0], 
    $pil->{MEAN},$pil->{SD}, $self->data($key) if $args->{'-verbose'};

    return $self->data($key);
}

=head2 fragments( -object => , -overlap => , -score => )

    We require evidence from direct alignment to a homolog to 
    determine that two models represent fragments of the same gene. 
    However, we have a variety of methods to rapidly reject this idea:
    1. If the genes are confidently assigned to different pillars. 
    2. If assigned to different pillars and have matching length profiles.
    3. If both preferentially align to the same region of the homolog.

=cut 

sub fragments {
    my $self = shift;
    my $args = {@_};

    my $obj = $args->{'-object'};
    $self->throw unless $self->isa( ref($obj) );
    $self->throw unless $self->strand == $obj->strand;
    $self->throw unless $self->up eq $obj->up;
    #$self->throw if $self->commoncodon($obj);
    
    $args->{'-verbose'} = 0 unless exists  $args->{'-verbose'};
    $args->{'-length'} = undef unless exists $args->{'-length'};
    $args->{'-overlap'} = -0.5 unless exists $args->{'-overlap'};
    $args->{'score'} = 100 unless exists $args->{'score'};

    # cannot be confirmed w/o homology 

    return undef unless my $hom = $self->homolog( -object => [$obj] );

    # but can be excluded 
    # does not require computatoin 
    
    if ( my $i1 = $self->pillar() ) {
	if ( my $i2 = $obj->pillar() ) {
	    #print 'PIL',$i1, $i2, $self->pillarscore, $obj->pillarscore, $args->{'score'} 
	    # if $args->{'-verbose'};
	    if ( $i1 ne $i2 ) {
		return 0 if $self->pillarscore > $args->{'-score'} && 
		    $obj->pillarscore > $args->{'-score'};	
	    } elsif ( $i1 eq $i2 ) {
		my $h1 = $self->homolog( -fast => 1 );
		# +ve == unique, -ve double tiling  
		my $olap = $self->overlap( -object => $obj, -compare => $h1, -evalue => undef );	
		print 'OLAP',$h1, $h2, $len, $olap if $args->{'-verbose'};
		return 0 if $olap <= $args->{'-overlap'};
	    }
	}
    }

    ######################################
    #
    ######################################
    
    if ( defined $args->{'-length'} ) {
	$self->throw;
	if ( my $i1 = $self->pillar() ) {
	    if ( my $i2 = $obj->pillar() ) { 
		my $p1 = $self->lookup( -query => $i1 );
		my $p2 = $obj->lookup( -query => $i2 );	    
		if ( $p1->{ID} && $p2->{ID} ) {
		    if ( $p1->{SD}>0 && $p2->{SD}>0 ) {
			my $z1 = ($self->length - $p1->{MEAN})/$p1->{SD};
			my $z2 = ($obj->length - $p2->{MEAN})/$p2->{SD};
			print 'LEN',@{$p1->{SGD}}, $z1, @{$p2->{SGD}}, $z2 if $args->{'-verbose'};  
			return 0 if abs($z1) <= $args->{'-length'} && abs($z2) <= $args->{'-length'};
		    }
		}
	    }
	}
    }    

    ######################################
    # positively confirm a candidate 
    ######################################
    
    # we use a weaker gap extension penalty to permit gaps to be 
    # crossed in merged seqeunces. 
    # --joinrangeext / --spanrangeext ??
 
    my @params=( -homolog => $hom, -return => 'score', -model => 'local', -extend => -1 );
    return 0 unless my $sc1 = $self->exonerate2( @params, -verbose =>0 );
    return 0 unless my $sc2 = $obj->exonerate2( @params, -verbose =>0 );
    
    my ($hi,$lo) = grep {defined} sort {$b <=> $a} ($sc1,$sc2);
    return 0 unless $lo > 50;

    # stop codons not dealt with properly here. 
    
    my @dummy = map { $_->stop(-R => 1, -adjust => -3) unless $_->length <= 3; $_ } 
    map { $_->clone() } map { $_->stream } ($obj,$self); 
    my $fuse = ref($self)->new( 
	START => 1, 
	STOP => 2, 
	STRAND => $self->strand, 
	EXONS => \@dummy, 
	UP => $self->up 
	);
    $self->throw unless $self->stream + $obj->stream == $fuse->stream; # exon count    
    $fuse->index;
    
    # get exonerate score for the fused gene 
    
    my $sc3 = $fuse->exonerate2( @params, -verbose =>0 );
    $fuse->oliver(-prepend => [$hom], -append => [$hi, $lo, $sc3] ) if $args->{'-verbose'};
    $fuse->DESTROY;

    # 
    #$orf->merge(-object => $lorf, -reference => $file);
    return ( $sc3 > ($hi + $lo/2) ? 1 : 0 );
}

=head2 validate 

    Verify object integrity by checking all reciprocal linkages
    and verifiying data hash is appropriately populated. 

=cut 

sub validate {
    my $self = shift;
    $self->throw("Exon recursion not implemented. Method not tested.");

    # 

    $self->throw unless $self->up && $self->up->contains($self);
    $self->throw unless ! $self->left || $self->left->right eq $self;
    $self->throw unless ! $self->right || $self->right->left eq $self;

    # descendents 
    
    $self->throw unless $self->down;
    map {  $self->throw unless $_->up eq $self } $self->stream;
    map { $_->validate } $self->stream; 

    # orthogroups
    
    # ohnologs
    
    # coordinate system 

    $self->throw unless $self->stop > $self->start;
    map { $self->throw unless $_->stop > $_->start } $self->stream; 
    map { $self->throw unless $self->strand == $_->strand } $self->stream;

    # data hash. check to ensure no dangling object refs 

    map { $self->throw unless defined $self->data($_) } @EVIDENCE;
    map { $self->throw if ref($self->data($_)) =~ /\:\:/ } 
    grep {!/EXCLUSIONS/} keys %{ $self->_data() };
    
    return $self;
}

=head2 translatable(-fast => 0|1)

    Validate that reading frame is open up to last codon.
    Returns true (self) on success, false (undef) on fail.

    Fails on:
    1. rank < -1 (TRNA/GAP/TELOMERE/MANUAL)
    2. length%3 != 0 
    3. premature STOP 

    -fast returns TRUE after 2. 
    
=cut

sub translatable {
    my $self = shift;
    my $args = {@_};

    return undef if $self->rank < -1; # not codon-based : TRNA / TELOMERE / GAP / MANUAL  
    return undef unless $self->length%3 == 0; # codon-based  

    return $self if $args->{'-fast'};

    my @aa = split //, $self->sequence(-molecule => 'aa');
    pop(@aa);
    map { return undef if $_ eq '*' } @aa;

    return $self;
}

=head2 rank

    Return rank (ie evidence code) associated with object. 

=cut

sub rank {
    my $self = shift;
    $self->{RANK}=$self->_evidence( -evidence => $self->evidence, -query => 'rank' ) unless $self->{RANK};
    return $self->{RANK};

    $self->throw;
    return $self->_evidence(
	-evidence => $self->evidence,
	-query => 'rank'
	);
}

=head2 assign

    Return inferred (assigned) sequence type based on evidence. 

=cut

sub assign {
    my $self = shift;
    $self->{ASSIGN}=$self->_evidence( -evidence => $self->evidence, -query => 'infer' ) unless $self->{ASSIGN};
    return $self->{ASSIGN};

    $self->throw;
    return $self->_evidence(
	-evidence => $self->evidence,
	-query => 'infer'
	);
}

#########################################
#########################################

# getters and setters 

=head2 exons(-query => 'number|initial|term|real')

    Return exons of specified type or number of exons. 

=cut

sub exons {
	my $self = shift;
    my $args = {@_};
              
	$args->{'-query'} = 'number' unless exists $args->{'-query'};

	######################################################
	# NNB exons : default behaviour reversed (DEVIN 4/2/9) 
	# needs proper checking. (very) cursory look indicates
	# no reprecussions but careful checking required. 
	$self->index if exists $args->{'-index'}; 
	######################################################

	# first is 'first exon' as opposed to first on contig

	my $return=1;
	if ($args->{'-query'} eq 'number') {
	    $return = scalar($self->stream);
	} elsif ($args->{'-query'} =~ /fir|init/ || $args->{'-query'} == 1) {	
	    $return = $self->down(-direction => 1);
	} elsif ($args->{'-query'} =~ /las|term/ || $args->{'-query'} == -1) {	
	    $return = $self->down(-direction => -1);
	} elsif ($args->{'-query'} =~ /internal/) {	
	    $self->warn("Not IMplemented");
	} elsif ($args->{'-query'} =~ /\d+/) {			
	    $self->warn("Not IMplemented");
	} elsif ($args->{'-query'} eq 'real') {		
	    foreach my $exon ($self->stream) {
		$return++ if $exon->intron(-direction => 'left') > 0;
	    }
	} else {$self->die;}
	
	$self->throw("Exons must return a value: $return (@_)") 
		unless $return || $return =~ /\d/;

	return($return);    
}

=head2 length

    ...

=cut 

sub length {
    my $self = shift;
    my $args = {@_};

    my $length;
    foreach my $ex ($self->stream) {
    	$length += $ex->length;
    }
    return $length;
}

=head2 sequence(-molecule => 'aa|dna', -format => 'string|fasta', 
    -decorate => '1|0', -id => 'internal', -nostop => '1|0')
    
    Return sequence as string. 

=cut 

sub sequence {
    my $self = shift;
    my $args = {@_};

    # exons
    
    my $seq;
    foreach my $ex ($self->stream) {	
    	if ($self->strand == -1) {
	    $seq = $ex->sequence.$seq;
    	} else {$seq .= $ex->sequence;}
    }
    
    # revomp
    
    $seq = join('', reverse(map {tr/ATGCN/TACGN/; $_;} split(//, $seq)))    
   	if ($self->strand == -1);
    
    # prot 
    
    if ($args->{'-molecule'} =~ /aa|prot|amino/i) { 
	my $prot;
	for (my $x = 0; $x < length($seq); $x += $TRIPLET) {
	    my $nirnberg = substr($seq, $x, $TRIPLET);
	    $prot .= (exists $CODONS{$nirnberg} ? $CODONS{$nirnberg} : 'X');
    	}
    	$seq = $prot;
	$seq =~ s/\*$// if ($args->{'-nostop'} && $seq =~ /\*$/);
    }

    # fasta output ?
    
    if ($args->{'-format'} eq 'fasta') {
	my $id;
	if ($args->{'-id'} eq 'internal') {
	    $id = $self->_internal_id;
	} else {$id = $self->name}		
	$id .= " ".join('; ', '['.$self->assign.'/'.$self->evidence, 			
			'YGOB:'.$self->ygob,
			'HMM:'.$self->evalue('ygob'),
			'SCER:'.$self->data('SGD'), # .'/'.$self->data('SCER'),

			'SYNT:'.($self->left ? $self->left->data('SGD') : 'NoLeft').
			'/'.$self->data('SYNT').'/'.
			($self->right ? $self->right->data('SGD') : 'NoRight'),

			'LEN:'.$self->data('LENGTH'),
			'LDA:'.$self->data('LDA'),
			'Ka/Ks:'.$self->data('KAKS'),
			'OHNO:'.($self->ohnolog ? $self->ohnolog->name : undef)
			.']') if $args->{'-decorate'};	
	$seq = ">$id\n".$seq;
    }
    
    return($seq);
}


=head2 fasta
    
    Return sequence fasta as formatted string. 

=cut 

sub fasta {
	my $self = shift;
	return ">".$self->name."\n".$self->sequence;
}

=head2 name
    
    Return ORF name. 

=cut 

sub name {
    my $self = shift;
    my $args = { @_ };
    $self->throw($self->debug) unless $self->up;
    my $contig = $self->up->id;
    my $organism = $self->up->up->organism;
    my $name = $organism."_".$contig.".".$self->id;

    if ( $args->{'-genbank'} ) {
	my $replace;
	if ( $self->evidence =~ /YGOB|KAKS|NCBI|AA|HSP|HCNF|SYNT|LDA|LENGTH|INTRONS|STRUCT|STOP/ ) {
	    $replace = 'p';
	} elsif ( $self->evidence =~ /LTR|TY/ ) {
	    $replace = 'x';
	} elsif ( $self->evidence =~ /RNA/ ) {
	    $replace = 'r';
	} elsif ( $self->evidence =~ /HMM/ ) {
	    $replace = 'f';
	} elsif ( $self->evidence =~ /MANUAL/ ) {
	    $replace = 'm';
	} else { return $name; }

	$name =~ s/\./$replace/ || die;
    }

    return( $name );
}

=head2 shortname
    
    Orf name that fits within one tab : ChrSpOrf
    Alias 'sn' also works. 
    
=cut 

sub shortname {
    my $self = shift;
    $self->throw($self->debug) unless $self->up;
    my $contig = $self->up->id;
    my $organism = substr($self->organism,1,1);
    return($contig.$organism.$self->id);
}

sub sn { return $_[0]->shortname; }

=head2 unique_id

    Return a unique stable id.

=cut

sub unique_id {
    my $self = shift;
    return $self->organism.'+'.$self->_internal_id;
}

=head2 family(-set => 'new|undef')

    Get/set family param. 
    
=cut 

sub family {
    my $self = shift;
    my $args = {@_};
    
    $self->throw unless (! @_) || exists $args->{'-set'};
    
    if ( exists $args->{'-set'} ) {
	$self->throw($args->{'-set'}) unless 
	    (! defined $args->{'-set'}) || ($args->{'-set'} =~ /^Anc_/);
	$self->throw("Must unset first.") if defined $args->{'-set'} && $self->family;
	$self->{'FAMILY'} = $args->{'-set'};
    }
    
    return $self->{'FAMILY'};
}

sub oliver { my $self = shift; $self->output(@_, -oliver=>1); }

=head2 output(-recurse => 1|0, -substr => 50, -string => 1|0,
    -ontology => 'P|F|C', -separator => "", -append => [], -prepend => [])
    
    Name, start, stop, length, strand, exons, structure...

=cut 

sub output {
    my $self = shift;
    my $args = {@_};

    # output format 
    my $fh = (exists $args->{'-fh'} ? $args->{'-fh'} : STDOUT);
    $args->{'-all'} = undef unless exists $args->{'-all'};
    $args->{'-oliver'} = undef unless exists $args->{'-oliver'};
    $args->{'-simple'} = undef unless exists $args->{'-simple'};
    $args->{'-quality'} = undef unless exists $args->{'-quality'};
    # string formatting 
    $args->{'-substr'} = 50 unless exists $args->{'-substr'};
    $args->{'-string'} = undef unless exists $args->{'-string'};
    $args->{'-separator'} = "" unless exists $args->{'-separator'}; # "\t*"
    # append extra data
    $args->{'-creator'} = undef unless exists $args->{'-creator'};
    $args->{'-internal'} = undef unless exists $args->{'-internal'};
    $args->{'-ontology'} = undef unless exists  $args->{'-ontology'};
    $args->{'-ohnolog'} = undef unless exists $args->{'-ohnolog'};
    $args->{'-og'} = undef unless exists $args->{'-og'};
    # recurse to exons etc
    $args->{'-recurse'} = 1 unless exists $args->{'-recurse'};
    $args->{'-recurse'} = undef if $args->{'-oliver'};

    if ( $args->{'-quality'} || $args->{'-qc'} ) {
	$args->{'-quality'}=1;
	$args->{'-recurse'}=undef;
    }

    #########################################
    # preparation 
    #########################################

    my $neighbour = 
	join('/', map {$_->gene if defined }		
	     ( ($self->neighbour(-direction => 'left', -orthogroup=>0, -variant => 0) || undef), 
	       ($self->neighbour(-direction => 'right', -orthogroup=>0, -variant => 0) || undef) )
	);

    #########################################
    # create output value array 
    #########################################
    
    my @r;
    if ( $args->{'-simple'} ) {
	@r=( $self->name,$self->start,$self->stop,$self->strand,$self->assign,$self->gene);
	
    } elsif ( $args->{'-quality'} ) {

	@r=(
	    $self->sn, $self->ogid, $self->data('KAKS'),(map {/Anc_(\d+\.\d+)/; 'A'.$1} $self->family),
	    (map {/Anc_(\d+\.\d+)/; 'A'.$1} $self->ygob), $self->logscore('ygob'),
	    $self->gene, $self->logscore('gene'),
	    $self->sgd, $self->logscore('sgd'),
	    $self->loss, $self->hypergob, $self->pillarscore,
	    ($self->ohnolog ? ('*'.$self->ohnolog->sn, $self->score('ohno'),$self->ohnolog->ogid) : (('NoOhno')x3) ),
	    (defined $self->quality ? ($self->quality) : (('NoOG')x4) )
	    );
	
    } elsif ( $args->{'-dump'} ) {
	$self->oliver;
	map { print "$_\t".$self->data($_) } sort keys %{ $self->_data };

    } elsif ( $args->{'-oliver'} ) {	
	@r = ( 
	    $self->name, $self->gene, ($self->score('global') || 0), #$self->data('__GENE'),
	    $self->start,$self->stop,$self->strand,# Coordinates
	    $self->data('STRUCT'),$self->data('STOP'),$self->data('_GAP'),
	    substr($neighbour, 0, 18),
	    $self->ygob, $self->score('ygob'), 
	    $self->evidence.'/'.$self->assign,
	    join('/', map {$_->length} $self->stream),
	    $self->data('_SUBTEL')
	    );		
	map { $r[$_] = '.' unless defined $r[$_] } 0..$#r;
	map { $r[$_].=' ' until length($r[$_])>=12; $_ } (1,10);
	map { $r[$_].=' ' until length($r[$_])>=10; $_ } (11..13);
	map { $r[$_].=' ' until length($r[$_])>=18; $_ } (9,12);
    } else {
	my ($name) = grep {defined} ($self->description, $self->gene, $self->sgd, $self->data('AA'), '.');
	$name .=' ' until length($name)>=18;
	$name = substr($name, 0, 18) if $name >18;
	my $scername = ($self->gene =~ /^Y[A-P][LR]/ ? $self->gene : $self->sgd) || '.';
	$scername .=' ' until length($scername)>=9;
	my $ygobname = $self->ygob || '.';
	$ygobname .=' ' until length($ygobname)>=9;

	@r = ( 
	    #$self->name, ($self->score('global') || 0), $name.$args->{'-separator'},	    
	    $self->name, ($self->data('EXONERATE.PROTEIN.GLOBAL') || 0), $name.$args->{'-separator'},	    
	    $self->start,$self->stop,$self->strand,$self->length.$args->{'-separator'}, # Coordinates 
	    $self->exons,$self->data('STRUCT'),$self->data('STOP').$args->{'-separator'}, # ORF Structure
	    $scername,$ygobname.$args->{'-separator'}, # Similarity 
	    #$self->data('SYNT'), $args->{'-separator'}, # Synteny
	    $self->logscore('ygob'),$self->logscore('aa'),$self->logscore('ncbi'),
	    $self->logscore('ty'),$self->logscore('ltr').$args->{'-separator'}, # Homology 
	    #$self->data('HCNF'),$self->data('RNA').$args->{'-separator'}, # Repeat Potential 
	    $self->data('LDA'),#$self->data('KAKS').$args->{'-separator'}, # Coding Potential
	    $self->evidence,$self->assign.$args->{'-separator'}, # Conclusion  
	    $self->data('HSP'), $self->data('HOMOLOG')
	    );
    }
  
    #########################################
    # fries with that? 
    #########################################
    
    push @r, ($self->ohnolog ? $self->ohnolog->name : undef)
	if $args->{'-ohnolog'};
    #push @r, substr($self->description, 0, $args->{'-substr'}) 
    #if $self->description;
    push @r, $self->_creator if $args->{'-creator'};
    push @r, substr(
	join("; ", grep {defined} $self->ontology($args->{'-ontology'})), 0 , $args->{'-substr'}) 
	if $args->{'-ontology'};

    #########################################    
    # 
    #########################################

    push @r, @{$args->{'-append'}} if $args->{'-append'};
    unshift @r, @{$args->{'-prepend'}} if $args->{'-prepend'};
    unshift @r, ($self->ogid || 'NA') if $args->{'-og'};
    unshift @r, $self->_internal_id if $args->{'-internal'};

    #########################################
    # pritn array / return as string 
    #########################################

    # temp: force columsn to line up 
    return join("\t", @r) if $args->{'-string'};
    print $fh @r;

    #########################################
    # recursion 
    #########################################

    delete $args->{'-append'};
    delete $args->{'-prepend'};
    map($_->output(%{$args}), $self->stream) if ($self->exons > 1 && $args->{'-recurse'});
    
    return $self;
}

=head2 _dissolve_ohnolog
=cut 

sub _dissolve_ohnolog {
    my $self = shift;
    my $args = {@_};
   
    $args->{'-verbose'}=1 unless exists  $args->{'-verbose'};
    
    $self->throw unless my $oh = $self->ohnolog;
    
    # electing to leave the data on the data hash for now. 
    
    $self->ohnolog( -object => undef );
    $self->throw if $self->ohnolog || $oh->ohnolog;
    return $self;
}

=head2 _dissolve_orthogroup
=cut 

sub _dissolve_orthogroup {
    my $self = shift;
    my $args = {@_};
   
    $args->{'-verbose'}=1 unless exists  $args->{'-verbose'};

    map { $self->throw unless $self->isa(ref($_)) && $_->up } @_;
    
    my @meth = map { uc($_) } $self->up->up->bound;

    if ( $args->{'-verbose'} ) {
	map { $_->oliver(-prepend => ['DISSOLVE']), 
	      -append => [$_->ogid] } ($self, $self->orthogroup);
	print;
    }
    
    foreach my $m ( @meth ) {
	$self->throw unless my $sp = $self->$m;
	$sp->_set_ogid(undef);	
	foreach my $k (qw(KA KS KAKS)) {
	    $sp->data($k => undef);
	}
	delete $self->{$m};
    }

    foreach my $k (qw(KA KS KAKS)) {
	$self->data($k => undef);
    }
    $self->_set_ogid(undef);
    return $self;
}

=head2 _define_orthogroup(-kaks => 0, -sowh => 0, -matrices => 0)

=cut 

sub _define_orthogroup {
    my $self = shift;
    my $args = {@_};

    ################################
    # Set defautl argnuments ... 
    ################################
    
    my @og = ( exists $args->{'-object'} ? @{$args->{'-object'}} : (grep {/\:\:/} @_) );

    $args->{'-verbose'}=1 unless exists  $args->{'-verbose'};
    $args->{'-kaks'}=0 unless exists $args->{'-kaks'};
    $args->{'-matrices'}=0 unless exists $args->{'-matrices'};
    $args->{'-sowh'}=0 unless exists $args->{'-sowh'};
    $args->{'-debug'}=0 unless exists $args->{'-debug'};

    ################################
    # Basic param checking 
    ################################

    if ( $args->{'-debug'} ) {
	my $fh = STDERR;
	print $fh caller(1);
	map { $_->oliver(-og => 1, -fh => $fh) } ($self,@og);
    }

    $self->throw( $#og ) unless scalar(@og) == scalar($self->up->up->bound);
    map { $self->throw unless $self->isa(ref($_)) && $_->up } @og;
    map { $self->throw if ($_->ogid || $_->orthogroup) } @og;
    my %sp = map { uc($_->organism) => $_ } @og;
    map { $self->throw($_) unless exists $sp{uc($_)} } $self->up->up->bound;
    $self->throw if exists $sp{uc($self->organism)};

    ################################
    # bad asignments -- this needs to be changed ... ropey. 
    ################################
    
    my %assign;
    map { $assign{$_->assign}++ } ($self, values %sp);
    print %assign and $self->throw if exists $assign{'GAP'} || 
	(scalar(keys %assign) > 1  && (exists $assign{'TRNA'} || exists $assign{'FEATURE'}) );  
    # exists $assign{'REPEAT'} || -- good decision? required to get HAP1.. Anc_1.380

    ################################
    # last chance to reject 
    ################################

    if ( $args->{'-sowh'} && $#og >= 1 && $self->coding ) {
	if ( my $sowh = $self->phyml( -object => [values %sp], -sowh => $args->{'-sowh'}, -verbose => 1 ) ) {
	    print ">>>".++$success, #$o->identify, $o->outgroup(), $o->outgroup(-ygob => 1),
	    $sowh->{DELTA}, $sowh->{MEAN}, $sowh->{SD}, $sowh->{PVAL}, $sowh->{RESULT} if $args->{'-verbose'};
	    return undef unless $sowh->{PVAL} < 0.05;
	}
    }

    ################################
    # make OG 
    ################################

    # we should be blessing into a new class here...
    # to be added later. we can then move some methods 
    # e.g. the paml() call below to that class and store
    # the data on an appropriate data strucutre. 

    my ($newogid) = (map { $_->ogid } sort { $b->ogid <=> $a->ogid } $self->up->up->orthogroups)[0] + 1;
    $self->_set_ogid($newogid);
    
    foreach my $sp ( keys %sp ) {
	my $meth = uc($sp);
	$self->$meth($sp{$sp});
	$sp{$sp}->_set_ogid($newogid);
    }
    
    ################################
    # optional compute Ka/Ks 
    ################################
    
    if ( $args->{'-kaks'} && $self->assign !~ /RNA/ ) {
	$self->kaks;
    }

    if ( $args->{'-matrices'} ) {
	$self->paml(
	    -object => [$self->_orthogroup],
	    -method => 'yn00'
	    );
    }

    ################################
    # verbose 
    ################################
    
    if ( $args->{'-verbose'} ) {
	print ">$newogid",$self->identify, 
	( map { $self->data($_) } qw(KA KS KAKS) ),
	( $self->data('KAKS') > .5 ? '*' : undef ).
	    ( $self->data('KAKS') > 1 ? '*' : undef );
	map { $_->oliver(-prepend => [$_->ogid] ) } $self->_orthogroup;
    }
    
    return $self;
}

=head2 _orthogroup 

    Return ($self, $self->orthogroup)

=cut 

sub _orthogroup { my $self = shift; return ( $self->orthogroup ? ($self, $self->orthogroup) : () ); }

=head2 ogid 

    Return orthogroup id. 

=cut 

sub ogid {
    my $self = shift;
    $self->throw if @_;
    return $self->{'OGID'};
}

sub _set_ogid {
    my $self = shift;    
    $self->throw unless $#_ == 0;
    $self->{OGID} = shift;
    return $self->{OGID};
}


=head2 structure_conserved

    Analog of synteny_conserved but returns self or undef. 

=cut 

sub structure_conserved {
    my $self = shift;        
    my $args = {@_};
    my $other = ( $args->{'-object'} ? [ $args->{'-object'} ] : [ $self->orthogroup ] );

    # 

    map { $self->throw unless $self->isa(ref($_)) } @{$other};      
    map { $_->structure } ( $self, @{$other} );

    # 

    foreach my $cmp ( @{$other} ) {
	return undef unless 
	    $self->data('INTRONS') == $cmp->data('INTRONS') && 
	    $self->data('STOP') == $cmp->data('STOP') && 
	    $self->exons == $cmp->exons &&
	    $self->assign eq $cmp->assign;
    }

    my ($m,$sd) = _calcMeanSD( map { $_->length } ($self, @{$other}) );
    return ( $sd/$m < .05 ? $self : undef);
}

=head2 ontology(-attribute => 'process') 
=cut 

sub ontology {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-attribute'} = 'P' if ! $args->{'-attribute'};
    foreach my $attr ('P', 'F', 'C') {
	$args->{'-attribute'} = $attr if $args->{'-attribute'} =~ /$attr/i;
    }

    my $OntologyBeast = $self->up->up->ontology;
    return undef unless $OntologyBeast && $self->ygob && 
	( $self->evalue('ygob') <= 1e-5 );

    return @{ $OntologyBeast->{$self->ygob}->{$args->{'-attribute'}} };
}

=head2 genbank_evidence

    Return appropriate Genbank feature for ASN1 format.

=cut

sub genbank_evidence {
    my $self = shift;
    my $args = {@_};

    my $string;
    if ( $self->evidence eq 'YGOB' ) {
	$string = "protein motif:HMMER:3.0b3:YeastGeneOrderBrowser:".$self->ygob;	    
    } elsif ( $self->evidence =~ /KAKS|NCBI|AA|HSP|HCNF/ ) {
	$string = "similar to AA sequence:blastall:2.2.17:GenBank:".$self->gene;	    

    } elsif (  $self->evidence eq 'LTR' ) {
	$string = "similar to DNA sequence:blastall:2.2.17:SGD:".$self->hit('ltr');
    } elsif (  $self->evidence eq 'TY' ) {
	$string = "similar to AA sequence:blastall:2.2.17:SGD:".$self->hit('ty');

    } elsif (  $self->evidence eq 'RNA' ) {
	$string = "nucleotide motif:tRNAscan-SE:1.23:".$self->gene;	    
    } elsif (  $self->evidence eq 'HMM' ) {
	$string = "nucleotide motif:HMMER:1.8.4:".$self->gene;

    } elsif (  $self->evidence eq 'NNNN' ) {
	$string = "alignment"; # slight abuse of the ontology ... 

    } elsif ( $self->evidence =~ /SYNT|LDA|LENGTH|INTRONS|STRUCT/ ) {
	$string = 'non-experimental evidence, no additional details recorded';
    } elsif (  $self->evidence eq 'STOP' ) {
	$string = 'non-experimental evidence, no additional details recorded';
    } elsif (  $self->evidence eq 'MANUAL' ) {
	$string = 'non-experimental evidence, no additional details recorded';
    }

    if ( $args->{'-existence'} == 1 ) {
	unless ( $self->evidence =~ /NNNN|MANUAL|STOP|SYNT|LDA|LENGTH|INTRONS|STRUCT/ ) {
	    $string = 'EXISTENCE:'.$string;
	}
    }

    return $string;
}

=head2 genbank_feature_key

    Return appropriate Genbank feature for ASN1 format.

=cut

sub genbank_feature_key {
    my $self = shift;       
    my $sofa = $self->_evidence(
	-evidence => $self->evidence,
	-query => 'genbank_asn1'
        );
    my $subsofa = $self->_evidence(
	-evidence => $self->evidence,
	-query => 'genbank_asn1_partof'
        );

    return($sofa,$subsofa);
}

=head2 sofa

    Return sofa ontology term.

=cut

sub sofa {
    my $self = shift;       
    my $sofa = $self->_evidence(
	-evidence => $self->evidence,
	-query => 'sofa'
        );
    my $subsofa = $self->_evidence(
	-evidence => $self->evidence,
	-query => 'sofa_partof'
        );

    return($sofa,$subsofa);
}

=head2 gff(-fh => $fh, -string => 1|0) 

    From http://www.sequenceontology.org/gff3.shtml
    
    Column 1: "seqid"
    Column 2: "source"
    Column 3: "type" http://www.sequenceontology.org/wiki/index.php/Category:SO:SOFA
    Columns 4 & 5: "start" and "end"
    Column 6: "score"
    Column 7: "strand"
    Column 8: "phase" {0,1,2} for CDS 
    Column 9: "attributes" tag=value pairs are separated by semicolons

=cut 

sub gff {
    my $self = shift;
    my $args = {@_};
    
    my $fh = (exists $args->{'-fh'} ? $args->{'-fh'} : STDOUT);
    
    $self->throw unless $args->{'-source'};

    #

    my %data = (
        ID => $self->_internal_id,
	#Parent => $self->name,
        Name => $self->name,
        gene => $self->gene,
        ygob => $self->ygob,
	sgd => $self->data('SGD') || $self->data('SCER')
        );

    #
    
    foreach my $key (grep {!/kaks2|none|nnnn/i} keys %GlobalVars::EVIDENCE) {
        $data{lc($key)} = $self->data($key);
	$data{lc($key).'_evalue'} = $self->evalue($key) if defined $self->data('_'.$key);
    }

    $data{'synteny_loss'} = sprintf("%.2f",$self->loss);
    $data{'synteny_hyper'} = $self->hyper;
    if ( $self->coding ) {
	map { $data{ $_.'_score' } = $self->score( $_ ) } qw(wise pillar local global ygob);
    }
    $data{'orthogroup'} = ($self->ogid ? 'OG'.$self->ogid : undef);
    $data{'status'} = $self->_evidence(
	-evidence => $self->evidence,
	-query => 'status'
        ) || '.';

    #
    
    my ($sofa,$subsofa) = $self->sofa;
    $sofa = 'centromere' if $self->gene =~ /^CEN/;   

    push my @gff, 
    $self->up->id, 
    $args->{'-source'},#( map {s/\s+/_/; $_} $self->up->up->source),
    $sofa,
    $self->start,
    $self->stop,
    (defined $self->score ? $self->score : '.'), # SCORE ? 
    ($self->strand > 0 ? '+' : '-' ),
    '.',
    join(';', map { $_."=".(defined $data{$_} ? $data{$_} : '.') } keys %data);

    return join("\t", @gff) if $args->{'-string'};

    print $fh @gff;
    if ( $sofa eq 'gene' ) {
	$gff[2] = 'mRNA';
	print $fh @gff;
    }

    next unless $subsofa;
    
    map { $_->gff(@_, -relationship => $subsofa) } $self->stream;
    map { $_->gff(@_, -relationship => 'intron') } grep {$_->length >= 30} $self->introns;
    
    return $self;
}

=head2 genbank_tbl()
=cut 

sub genbank_tbl {
    my $self = shift;
    my $args = {@_};

    #################################
    # validate args and organize variables 
    #################################

    $self->throw unless exists $args->{'-fh'};
    my $fh = $args->{'-fh'};
    
    my @bump = (3 x undef);

    #################################
    # exclude certain features 
    #################################

    my ($asn,$asn_sub) = $self->genbank_feature_key;
    return undef unless $asn;
 
    $self->evaluate( -force => 1);

    #################################
    # standard components 
    #################################
    
    my $start = $self->strand == -1 ? $self->stop : $self->start;
    my $stop = $self->strand == -1 ? $self->start : $self->stop;
    
    my ($top,$tail) = $self->_top_tail;
    if ( $asn eq 'CDS' ) {
	unless ( $top eq 'ATG' ) {
	    $start = $self->strand == -1 ? '>'.$start : '<'.$start;
	}
	unless ( $tail =~ /TAA|TAG|TGA/ ) {
	    $stop = $self->strand == -1 ? '<'.$stop : '>'.$stop;
	}
    }

    print {$fh} $start, $stop, $asn, undef, undef;
    print {$fh} @bump, 'citation', $args->{'-reference'} if $args->{'-reference'};
    unless (  $asn eq 'assembly_gap' ) {
	print {$fh} @bump, 'gene', $self->name( -genbank => 1 );
	print {$fh} @bump, 'inference', $self->genbank_evidence( -existence => 1 );
	print {$fh} @bump, 'locus_tag', (map { s/\+//; $_ } $self->unique_id);
    }

    #print {$fh} @bump, 'standard_name', $self->name; # diff from gene ?

    #################################
    # standard components 
    #################################

    if ( $asn eq 'CDS' ) {

	# stable database Id

	#print {$fh} @bump, 'protein_id', (map {uc($_)} map { /^(\w{3})/ } $self->organism).'1';

	# homology 

	print {$fh} @bump, 'product', 'hypothetical protein';
	print {$fh} @bump, 'product', 'homology to '.$self->gene if 
	    $self->gene && $self->evalue('gene') <= 1e-5;
	print {$fh} @bump, 'product', 'homology to '.$self->ygob if
	    $self->ygob && $self->evalue('ygob') <= 1e-5;

	# function 
	
	my @func = join(';', grep {!/molecular_function/} $self->ontology( -attribute => 'F'));
	print {$fh} @bump, 'function', @func, 
	'[general function prediction based on homology to S. cerevisiae]'  if $func[0];	
	
	# translation 

	print {$fh} @bump, 'codon_start', $self->fex->frame+1;
	if ( $self->pseudogene ) {
	    print {$fh} @bump, 'pseudogene', 'unitary';	 # unknown
	} elsif ( $self->data('STOP') ) {
	    # use pseudo if not an actual pseudogene but gene is disrupted by sequencing error etc 
	    print {$fh} @bump, 'pseudo', undef;
	} else {
	    print {$fh} @bump, 'translation', $self->sequence( -molecule => 'aa' );
	    #print {$fh} @bump, 'transl_table', 1;
	}

    } elsif ( $asn eq 'assembly_gap' ) {

	print {$fh} @bump, 'estimated_length','unknown';
	print {$fh} @bump, 'gap_type','within_scaffold';
	print {$fh} @bump, 'linkage_evidence','align genus'; # unspecified 

    } elsif ( $asn eq 'tRNA' ) {
	
	print {$fh} @bump, 'pseudogene', 'unitary' if $self->gene =~ /pseudo/i;
	
    } elsif ( $asn eq 'misc_feature' ) {

	print {$fh} @bump, 'function', $self->description;

    } elsif ( $asn eq 'mobile_element' ) {

	print {$fh} @bump, 'mobile_element_type', 'retrotransposon';
	print {$fh} @bump, 'rpt_family', 'TY';
	print {$fh} @bump, 'rpt_type', 'dispersed';

    } elsif ( $asn eq 'LTR' ) {

	print {$fh} @bump, 'rpt_family', 'TY';
	print {$fh} @bump, 'rpt_type', 'dispersed';

    } else { $self->throw( $asn ); }

    #################################
    # include exon and 'intron' features 
    #################################

    return 1;
}


=head2 _upgrade_orf_structure() 

    Method to change the modify the Orf structure used for 
    the 2011 paper to the one expected by the current code base. 
    Almost all changes are to the data attribute ($Orf->data())
    but a few key elements are now stored on the main object
    reference : 
    $orf->_debug() 
    $orf->ohnolog()
    $orf->ogid()

=cut

sub _upgrade_orf_structure {
    my $self = shift;

    ##############################
    # Catalog things to change 
    ##############################

    my %repair_map = (
	'YGOBHMM' => '_YGOB',
	'NCBI' => '_NCBI', # there is no '_NCBI' attribute -- stored directly on 'GENE'. bad. 
	'LTR' => 'swap',
	'TY' => 'swap',
	'AA' => 'swap',
	#'_SGD_GENE' => undef, # this is now stored on the Orf object as _DEBUG
	#'_OG' => undef,  # this is now stored on the Orf object 
	#'_OHNOLOG' => undef, # this is now stored on the Orf object 
	#'_GAP' => undef, # Honor. I think we just track this for no reason 
	#'_TMP'=> undef, # Do not honor 
	#'_TEMP'=> undef, # Do not honor 
	);
    
    ##############################
    # Test all elements of the array to make sure we understand them 
    ##############################

    my $old = $self->_data;
    
    foreach my $key (keys %{$old} ) {
	my $key_bkp = $key;
	my ($prefix,$suffix)=();

	# see _init_homology_data

	if ( $key =~ /(_CHR|\*)$/ ) {
	    $suffix = $1;
	    $key =~ s/$suffix//;
	} else { $suffix=undef; }

	#  see _init_homology_data

	if ( $key =~ /^(_{1,3})/ ) {
	    $prefix=$1;
	    $key =~ s/$prefix//;
	    #$self->throw( $key_bkp, $old->{$key_bkp} ) if $old->{$key_bkp} =~ /[a-z]/i;
	} else { $prefix=undef; }
	
	#print $key_bkp, $prefix, $key, $suffix, $old->{ $key };

	# 

	if ( $key eq 'OG' || $key eq 'GAP' || 
	     $key eq 'TMP' ||  $key eq 'TEMP' || 
	     $key eq 'SGD_GENE' || $key eq 'OHNOLOG' ) {
	    next;
	} elsif ( $suffix ) {
	    $self->throw( $key_bkp ) unless exists $HOMOLOGY{ $key }; 
	} elsif ( $prefix ) {
	    $self->throw( $key_bkp ) unless exists 
		$HOMOLOGY{ $key } || $EVIDENCE{ $key } || $key eq 'GENE'; 
	}
    }

    ##############################
    # Edit the data -- special ones first 
    ##############################

    if ( exists  $old->{'_OG'} ) {
	$self->throw if $self->{'OGID'};
	$self->{'OGID'} = $old->{'_OG'};
	delete $old->{ '_OG' };
    }

    if ( exists  $old->{'_OHNOLOG'} ) {
	$self->throw if $self->ohnolog;
	$self->ohnolog( -object => $old->{'_OHNOLOG'} );
	delete $old->{ '_OHNOLOG' };
        $self->ohnolog->data( '_OHNOLOG' => 'delete' )
	    if $self->ohnolog->data( '_OHNOLOG' ) ;
    }
    
    if ( exists  $old->{'_SGD_GENE'} ) {
	$self->throw if $self->_debug;
	$self->_debug( $old->{'_SGD_GENE'} );
	delete $old->{ '_SGD_GENE' };
    }

    delete $old->{ '_TMP' };
    delete $old->{ '_TEMP' };

    ##############################
    # Edit the data -- 
    ##############################
    
    foreach my $k ( keys %repair_map ) {
	if ( $repair_map{ $k } eq 'swap' ) {
	    my ($x,$y) = ($old->{ $k }, $old->{ '_'.$k });
	    $old->{ $k } = $y;
	    $old->{ '_'.$k } = $x;
	} else {
	    $old->{ $repair_map{ $k } } = $old->{ $k };
	    delete $old->{ $k };
	}
    }
    
    ##############################
    # Edit the data -- 
    ##############################
    
    #print %{ $old };
    $self->_data( $old ); 
    #print %{ $self->_data };
    $self->_init_data_structure( $old );

    # 

    $self->evaluate;

    return $self;
}

=head2 data(KEY => val)
    
    Get/Set values from DATA hash. 

    KEY => undef     -- set value to undef
    KEY => 'delete'  -- remove attribute from DATA hash

=cut

# retrieve/alter the entire data hash 

sub _data { 
    my $self = shift;
    $self->throw if @_ && ref($_[0]) ne 'HASH';
    $self->{'DATA'} = shift if @_;
    return $self->{'DATA'}; 
}

# access elements 

sub data {
    my $self = shift;             
    my ($k, $v) = @_;
    
    $k =~ s/^\-//;	
    $self->throw("No key value - $k, $v") unless $k;		

    # this is a protective mechanism 
    # LENGTH is never actually set and just stores the default 
    return $self->length if ($k eq 'LENGTH'); 
    # 

    my $data = $self->_data();
    if ( $v eq 'delete' ) {
	return delete $data->{$k};
    } elsif ( ref($v) eq 'CODE' ) { # can handle either { delete } or { undef } 
	$v->( $data->{$k} );
    } elsif (scalar(@_) == 2) { # can handle undef or defined value 
	$data->{$k} = $v;
    }

    return $data->{$k};
}

=head2 accept('YGOB|AA|..', $hit)

    Add hit, evalue and score to the DATA hash
    in appropriate fashion.

=cut 

sub accept {
    my $self = shift;
    my $cat = shift;
    my $hit = shift;
    
    my $catstar = $cat;
    $catstar =~ s/\*$//;
    
    $self->throw($cat) unless $cat && 
	($catstar eq 'GENE' || $catstar eq 'OHNO' || # will OHNO break something?
	 exists $EVIDENCE{$catstar} || exists $HOMOLOGY{$catstar});
    $self->throw unless $hit && ref($hit) eq 'HASH';
    $self->throw unless (exists  $hit->{HIT} && exists $hit->{EVALUE} && exists $hit->{SCORE});
    
    $self->data($cat => $hit->{HIT});
    $self->data('_'.$cat => $hit->{EVALUE}*1);
    $self->data('__'.$cat => $hit->{SCORE}*1);
    $self->data('___'.$cat => $hit->{SYNTENY}*1) if defined $hit->{SYNTENY};    
    
    return $self;
}

=head2 _init_homology_data
    
=cut 

sub _init_homology_data {
    my $self = shift;
    my $args = {@_};

    my $exempt = join( '|', ( exists $args->{'-exempt'} ? @{$args->{'-exempt'}} : ('_acceptall_') ) );

    my $data = $self->_data();
    foreach my $sp ( grep {!/$exempt/} (keys %HOMOLOGY, 'GENE') ) {
	foreach my $prefix (undef, '_', '__', '___') {
	    foreach my $suffix (undef, '*', '_CHR') {
		delete $data->{ $prefix.$sp.$suffix } if 
		    exists $data->{ $prefix.$sp.$suffix };
	    }
	}	
    }

    return $self;
}

=head2 hit('aa')

    Return best hit for a given DB search. 

=cut 

sub hit {
    my $self = shift;
    my $key = uc(shift);
    $key = 'AA' if $key eq 'BLAST';
    $key = 'YGOB' if $key =~ /^HMMER/;
    return $self->data($key);
}

=head2 evalue('aa')

    Return E-value for a given DB search. 

=cut 

sub evalue {
    my $self = shift;
    my $key = uc(shift);
    $key = 'AA' if $key eq 'BLAST';
    $key = 'YGOB' if $key =~ /^HMMER/;
    return $self->data('_'.$key);
}

=head2 score('ygob')
    
    Retrieve score (not evalue) of specified type: 
    blast/ygob/global/local/gene/wise

    Default is to gene (ie imputed gene and BLAST score). 
    
    NB: functions as a getter - not setter - only. 

    See logscore() to return the exponent of the evalue. 

=cut 

sub score {
    my $self = shift;
    my $score = lc(shift || 'dna');
    
    if ( $score eq 'gene' ) { # imputed gene BLAST score.  
	return $self->data('__GENE');
    } elsif ( $score eq 'blast' || $score eq 'aa' ) { # best BLAST hit 
	return $self->data('__AA'); 
    } elsif ( $score eq 'ygob' || $score =~ /^hmmer/ ) { # best YGOB hit 
	return $self->data('__YGOB');
    } elsif ( $score eq 'pillar' ) {   # confidence for YGOB pillar assignment 
	return $self->data('_PILLAR'); # also pillarscore(). for now. 

    } elsif ( $score =~ /ohno$/i ) {   # ohnolog scoring algorithm ohnolog2()
	return $self->data('__OHNO');  
    } elsif ( $score =~ /ohno/i ) {   # changed in ohnologs2()
	$self->throw("Depracated: Use 'ohno' not 'ohnolog'.");
	
    } elsif ( $score eq 'wise') {    # 
	return $self->data('WISE');
    } elsif ( $score eq 'local' ) {  # exonerate local protein align
	$self->exonerate2( -model => 'local', -return => 'score' ) if 
	    ! defined $self->data('EXONERATE.PROTEIN.LOCAL') && $self->rank > 0 && $self->homolog(-fast => 1);
	return $self->data('EXONERATE.PROTEIN.LOCAL');
	
    } elsif ( $score eq 'global' ) { # exonerate global protein align
	$self->exonerate2( -model => 'global', -return => 'score' ) if 
	    ! defined $self->data('EXONERATE.PROTEIN.GLOBAL') && $self->rank > 0 && $self->homolog(-fast => 1);
	return $self->data('EXONERATE.PROTEIN.GLOBAL');	
	
    } elsif ($score eq 'dna') {
	return $self->data('EXONERATE.DNA.LOCAL');
    } elsif ( exists $HOMOLOGY{uc($score)} ) {
	return $self->data('__'.uc($score)); # species specific scores 
    } else { $self->throw($score); }

} 

=head2 logscore('ygob')

    Return exponent of the evalue. See 'score' for the raw score. 
    
    eg 2.1e-300 => -330 

    The original motivation for this was to not return 0 
    - which is hard to work with - but to return the lower 
    limit supported by HMMER3, ~ 1e-450.

=cut

sub logscore { 
    my $self = shift;
    
    my $attr = lc(shift) || 'ygob';
    die unless $attr;
    my $key = uc($attr);

    my %limits = (
	'perl' => {'min' => 1e-300, 'use' => -300},
	'ygob' => {'min' => 1e-450, 'use' => -450},
	'aa' => {'min' => 1e-180, 'use' => -180},
	'ncbi' => {'min' => 1e-180, 'use' => -180},
	'ty' => {'min' => 1e-180, 'use' => -180}, # ty 
	'ltr' => {'min' => 1e-180, 'use' => -180}, # ltr  
	);
    
    if ( $self->evalue($key) eq '0' || $self->evalue($key) == 0 || $self->evalue($key) eq '0.0' ) { #sting compare reqd.
	return $limits{$attr}{'use'};
    } elsif ( $self->evalue($key) <= $limits{'perl'}{$min} ) { # limit of what perl can handle 
	$self->evalue($key) =~ /e\-(\d+)/;
	return -$1;
    } else {
	return sprintf("%.2f", log( $self->evalue($key) )/log(10));
    }
}

sub synscore { 
    my $self = shift;
    my $key = uc(shift);
    $self->throw unless $key;
    return $self->data('___'.$key);
}

###############################################
# still used by _process_BLAST. ugh.  
sub blastscore { return $_[0]->score('aa'); }
###############################################

=head2 ygob
    
    Return YGOB hit. 

=cut

sub ygob {return $_[0]->data('YGOB');}

=head2 gene
    
    Return gene name imputed from BLAST and synteny 
    in the preferred species. 

=cut

sub gene { return $_[0]->data('GENE'); }

=head2 sgd

    Return best SGD hit. 

=cut

sub sgd { return $_[0]->data('SGD'); }

=head2 aa

    Return amino acid sequence as string. 

=cut

sub aa { return $_[0]->sequence(-molecule =>'aa'); }

#########################################
#########################################

# override methods 
# this is correct! 
# NB the 'relative' notation is NOT suported for ORFs.
 
sub start {
    my $self = shift;
    if (my $ex = $self->down(-direction => $self->strand)) {
	return $ex->start(@_);
    } else {$self->throw;}
}

sub stop {
    my $self = shift;
    if (my $ex = $self->down(-direction => $self->strand*-1)) {
	return $ex->stop(@_);
    } else {$self->throw;}
}

#########################################
#########################################

=head2 _init_data_structure(-exempt => []) 

    A DATA hash contains all information about an ORF 
    that can be used to establish its vearicty, nature
    or quality. These include scores for RNAfolding/coding 
    potential and repeat similarity that argue for a specific 
    function. These are typically assessed by the evidence/
    evaluate method set and are listed in the %EVIDENCE hash.

    In addition, protein coding sequences have available scores
    for gene quality and similairty as well as a full set of 
    homologs in other species.

    Information on ORF coordinates etc are not stored on the 
    DATA hash. They are devolved to exon objects. In addition, 
    information such as neighbours on the contig, ohnologs and
    orthologs in other annotated genomes are the preserve of 
    the Orf object itself. They are not stored in the DATA hash.

    The DATA hash is initialized when the ORF is first created.
    All evidence values are created as defaults and all other 
    properties are created as undef except homologs in other 
    speices (these are inferred from BLAST hits as opposed to 
    orthologs which are computed) which are created as required. 
    
    The DATA hash is reinitialized when update() is called though
    certain pieces of info are permitted to persist. 
    
=cut

sub _init_data_structure {
    my $self = shift;

    #########################################
    # preamble 
    #########################################

    # gather values to exempt from overwriting. 

    my $exempt;
    if ( @_ && $#_ == 0 && ref($_[0]) eq 'HASH' ) { # pass an existing data hash ref
	$exempt = $_[0];
    } elsif ( @_ && ref($_[0]) ne 'HASH' ) { # pass an array 
	map { $exempt->{$_} = $self->data($_) } @_;
    } elsif (@_) {$self->throw("@_");}

    # create a new DATA hash. 

    my $data = $self->_data( {} ); 

    #########################################
    # set all evidence to defaults 
    #########################################
    
    # evaluate can be called at any time so defaults are required. 
    # since update does not recall LTR methods and some others,
    # we allow these to persist -- though they become inaccurate. 

    foreach my $db ( keys %EVIDENCE) {
	# if this section is changed, corresponding changes must be made in 
	# evaluate() which also treats homology methods differently 
	if ( $db =~ /AA|YGOB|TY|LTR|NCBI/ ) { # these are multivalued (effectively "HITs") and handled differently 
	    $self->accept($db, {'HIT' => undef, 'SCORE' => 0, 'EVALUE' => $EVIDENCE{$db}->{'DEF'}} );
	} else { $self->data( $db => $EVIDENCE{$db}->{'DEF'} ); }
    }

    #########################################    
    # delete all homology data
    #########################################

    # these keys are not assumed. delete unless needed.     
    # $self->_init_homology_data();
    
    #########################################
    # undef all remaining vars 
    #########################################

    #this is a mess. why undef? not delete?     
    #map { $data->{ $_ } = undef if exists $data->{ $_ } } 
    #grep {!/$exempt/} (
    #qw(EXTRA TANDEM _GAP _TMP _CLUSTER_VAR) , # junk vars 
    #qw(HOMOLOG _HOMOLOG _SUBTEL SUBTEL PILLAR _PILLAR FRAG _FRAG) , # computed vars 
    #qw(EXONERATE.PROTEIN.LOCAL EXONERATE.PROTEIN.GLOBAL EXONERATE.DNA.LOCAL WISE _WISE) # external scores 
    #);

    #########################################    
    # restore exempted values 
    #########################################

    map { $data->{$_} = $exempt->{$_} } keys %{$exempt};

    #########################################
    # evaluate 
    #########################################

    $self->{'EVIDENCE'} = 'NONE'; # no key piece of evidence yet   
    $self->evaluate(-validate => 0, -structure => 0);    
    
    return $self; 
}

sub _init_exon_structure {
    my $self = shift;

    $self->throw if $self->exons;

    if ( @_ ) {
	
	map { $self->throw unless $self->_down->isa(ref($_)) } @_;
	map { $self->add(-object => $_) } @_;
	
    } else {

	my $exon = ref($self->_down)	    
	    ->new(
	    START => $self->SUPER::start, # autoload in 'Annotation' used  
	    STOP => $self->SUPER::stop,   # because overridden in 'Orf'
	    STRAND => $self->strand,
	    INTRON => [0, 0],
	    UP => $self
	    );	
	$self->add(-object => $exon);
	
    }
 
    $self->index;
    return $self;
}

sub _invert_coords {
    return map { $_->_invert_coords } $_[0]->stream;
}

sub _top_tail {
    my $self = shift;
    my $first = $self->exons(-query => 'first');
    my $last = $self->exons(-query => 'last');
    my ($ftop, $ftail) = $first->_top_tail;
    my ($ltop, $ltail) = $last->_top_tail;
    return ($ftop, $ltail);
}

sub _terminal_dist2 {
    my $self = shift;
    my $args = {@_};

    if ( $#_==0 ) {
	my $val = shift;
	$self->throw unless abs($val)==1;
	$args->{'-direction'} = ( $val == 1 ? 'right' : 'left');
    }

    if ( $args->{'-direction'} ) {
	return( 
	    ($args->{'-direction'} eq 'left' || $args->{'-direction'} == -1) ?  
	    $self->start : ($self->up->length - $self->stop +1) 
	    );
    }

    return( $self->start, ($self->up->length - $self->stop +1) );
}

# this method is incomprehensible. 
# gut it and exon->morph() ASAP. 
# rework in terms of contig_term => , exon => start|stop 

sub _terminal_dist {
    my $self = shift;

    $args->{'-query'} = 'proximal' unless exists $args->{'-query'};
    $args->{'-method'} = 'start' unless exists $args->{'-method'};

    my @val; 
    foreach my $ex ($self->stream) {
	push @val, $ex->_terminal_dist(@_);
    }
    return $val[0] unless $#val > 0;

    # values are in exon order 
    # ex1, ex2 .. exn 

    my $val;

    ################################
    # Mistake? 
    # proximal does not necessarily refer to first exon 
    # DS. This is the same problem as calling Orf start/stop with -R. 
    if ($args->{'-query'} =~ /^p/) {
	$val = shift(@val);
    } elsif ($args->{'-query'} =~ /^d/) {
	$val = pop(@val);
    ################################

    } elsif ($args->{'-query'} =~ /^m/) {
	my ($min, $max) = sort {$a <=> $b} ($val[0], $val[$#val]);
	if ($args->{'-query'} =~ /^min/) {
	    $val = $min;
	} else {$val = $max;}
    } else {$self->throw;}

    return $val;
}

##########################################
##########################################
# Borrowed Code. Lightly Tested. 
# http://www.perlmonks.org/?node_id=466599 

sub _logfact {
   return _gammln(shift(@_) + 1.0);
}

sub _hypergeom {
   # There are m "bad" and n "good" balls in an urn.
   # Pick N of them. The probability of i or more successful selections:
   # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
   my ($n, $m, $N, $i) = @_;
   my $loghyp1 = _logfact($m)+_logfact($n)+_logfact($N)+_logfact($m+$n-$N);
   my $loghyp2 = _logfact($i)+_logfact($n-$i)+_logfact($m+$i-$N)+_logfact($N-$i)+_logfact($m+$n);
   return exp($loghyp1 - $loghyp2);
}

sub _gammln {
  my $xx = shift;
  my @cof = (76.18009172947146, -86.50532032941677,
             24.01409824083091, -1.231739572450155,
             0.12086509738661e-2, -0.5395239384953e-5);
  my $y = my $x = $xx;
  my $tmp = $x + 5.5;
  $tmp -= ($x + .5) * log($tmp);
  my $ser = 1.000000000190015;
  for my $j (0..5) {
     $ser += $cof[$j]/++$y;
  }
  -$tmp + log(2.5066282746310005*$ser/$x);
}

##########################################
##########################################

sub time {
    my $self = shift;
    my $args = {@_};

    my $delta = 0;
    my $stamped_time = time;
    if ( $args->{'-init'} || ! $global_time_tracker ) {
	$delta = 0;
    } else {
	$delta = $stamped_time - $global_time_tracker;
    }
    $global_time_tracker=$stamped_time;

    return $delta;
}

##########################################
##########################################

# Depracated 
sub _splice_exons {
    my $self = shift;
    my $args = {@_};
    
    # 'splice' overlapping exons together to minimise overlap.
    # if this routine createas a STOP at an exon splice junction
    # will move juntion until STOP is removed
    
    # do not inherit somebody elses problems
    
    $self->throw("ORF does not validate\n".$self->debug)
	unless $self->translatable;
    
    # argument checking 
    
    $self->throw("Must supply valid object: $args->{'-object'}") 
	unless ref($args->{'-exon1'}) eq ref($args->{'-exon2'});
    $self->throw("Exons do not belong to this Orf") 
	unless $args->{'-exon1'}->up eq $self && 
	$args->{'-exon2'}->up eq $self;
    #
    ($self->output && $self->throw("Objects must overlap")) unless  
	$args->{'-exon1'}->overlap(-object => $args->{'-exon2'});
    
    # order exons -- do not take anything on trust 
    
    my ($x, $y) = sort {$a->_terminal_dist <=> # sort in exon order 
			    $b->_terminal_dist} ($args->{'-exon1'}, $args->{'-exon2'});	
    
    # y too small. probably won't happen. if y is very small
    # but accounts for a large % of the exon (as tested in optimise
    # and must be case to call _splice_exons) exon will be destroyed. 
    
    return $self if $TRIPLET >= 
	abs($y->stop(-R => 1) - $x->stop(-R => 1));
    
    # make the cut. since there is no principaled way to do this
    # we will at least be consistent -- always favour the leading exon 
    # and cut just before the STOP codon 	
    
    my ($xs,$ys); # dist to end of gene 
    if ($self->strand == -1) {
	$ys = $y->start(-R => 1) - 1;
	$xs = $x->stop(-R => 1) - 1;
    } else {
	$ys = $self->up->length - $y->start(-R => 1);
	$xs = $self->up->length - $x->stop(-R => 1);				
    }
    
    until ($ys < $xs) {
	$y->start(-R => 1, -adjust => +3);
	if ($self->strand == -1) {
	    $ys = $y->start(-R => 1) - 1;
	    $xs = $x->stop(-R => 1) - 1;
	} else {
	    $ys = $self->up->length - $y->start(-R => 1);
	    $xs = $self->up->length - $x->stop(-R => 1);				
	}
    }
    
    return $self;		
}

########################################
# I AM A MODULE 
########################################

1;
