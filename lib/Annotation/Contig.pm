#!/usr/bin/perl

package Annotation::Contig;

use Annotation;
@ISA = qw(Annotation);

use Annotation::Orf;
use GlobalVars;

#########################################
#########################################

# autoloaded methods 

our $AUTOLOAD;

my %AUTOMETH = (
    sequence => undef,
#    scaffold => undef,
#    kmercoverage => undef
);

sub AUTOLOAD {
    my $self = shift;    
    my $name = $AUTOLOAD;
    my $name = (split(/\:/, $AUTOLOAD))[-1];
    #$name =~ s/.*://;   # strip fully-qualified portion
    $name = uc($name);

    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

sub DESTROY {
    my $self=shift;

    $self->_initialise_features;
    map { $_->DESTROY } grep {defined} $self->_get_features; 

    return $self->SUPER::DESTROY;
}

#########################################
#########################################

=head2 new()

    SEQUENCE => seq, 
    SCAFFOLD => undef

=cut

sub new {
    my $class = shift;        
    my $args = {@_};

    $args->{'KMERCOVERAGE'} = -1; # DEVIN 20110303 

    my $base = Annotation::Orf->new(ID => 0);
    my $self  = $class->SUPER::new(@_, DOWN => $base);
    $base->up($self);

    foreach my $k (map {uc($_)} keys %AUTOMETH) {
	if ($self->id == 0) {
	    $self->{$k} = $AUTOMETH{$k};
	} else {
	    $self->throw("$k : $args->{$k} ") 
		unless exists $args->{$k} && defined $args->{$k};
	    $self->{$k} = $args->{$k};
	}
    }    
    
    $self->{'_FEATS'} = [];    
    bless $self, $class;
}

#########################################
#########################################

=head2 predict(-orf_min => 300, -intron_min => x, -term_min => y, -intron_max => z)

    Enumerate all ORFs that satisfy constraints. 
    
=cut 

sub predict {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-stop'} = '\*' unless exists $args->{'-stop'};
    $args->{'-optimise'} = 1 unless exists $args->{'-optimise'};
    $self->throw("Cannot search for Orfs this small: @_")
	unless $args->{'-orf_min'} >= 3;
    
    ###########
    # INITIALISE
    ###########
    
    # prep. set _orf_min in  %EVIDENCE (accessed by evaluate).
    # remove all extants orfs. this permits annotate to be called
    # multiple times if necessary. intron models etc are still 
    # on _FEATS. tRNAs are recreated if desired and annotation
    # begins from scratch. 
    
    $self->_set_default(	
	-evidence => 'LENGTH',
	-set => $args->{'-orf_min'},
	);
    
    map {$_->DESTROY} $self->stream; # blank slate 
    $self->_initialise_features;  # destroy all refs except 'link'
    
    my %strand = (
	1 => $self->sequence,
	-1 => $self->revcomp
	);
    
    ###########
    # NONCODING
    ###########

    # link intron ends etc 

    foreach my $cw (keys %strand) {
	$self->_link_noncoding(
	    -strand => $cw,
	    -feat_min => $args->{'-feat_min'},
	    );
    }

    # run tRNAscan + other tools 
			
    foreach my $method (keys %{$args}) {
	next if $method =~ /min|max|frame|strand|op/; # st[op]timise
	my $param = $args->{$method};
	$method =~ s/^-//;
	$self->throw($method) unless $self->can($method);
	$self->$method(-score => $param); # I guess this is trnascan etc but no idea whst going on.... 
    }	

    map {$_->evaluate(-validate => 0, -structure => 0)} $self->stream;
    
    #########
    # CODING
    #########
    
    # the basic outline: 
    # 1) all relevant features (introns/STOPs/etc) are placed
    # 2) neighbouring features are linked to produce exons
    # 3) exons are linked to produce full gene structures
    # we use a triple loop to do one reading frame at a time
    
    foreach my $cw (keys %strand) {      
	next if defined $args->{'-strand'} && 
	    $args->{'-strand'} != $cw;
	
	foreach my $frame (0, 1, 2) {		
	    next if defined $args->{'-frame'} && 
		$args->{'-frame'} != $frame;
	    
	    # use a pseudo contig to hold all features
	    
	    my $pseudo = Annotation::Contig
		->new(
		SEQUENCE => $strand{$cw},
		SCAFFOLD => -$INFINITY,
		UP => $self->up		
		);		
	    
	    ###########
	    # FEATURES
	    ###########
	    
	    # ok. add all features that can affect gene prediction
	    # (centromers/telomers etc) -- in correct reading frame. 
	    # strand zero objects must be flipped when on - strand
	    # so they are available to both strands (eg gap objects -- NNNNNN)
	    
	    foreach my $feat ($self->_get_features) {
		next if $feat->strand == -1*$cw;

		# correct orientation of strand 0 objects 

		my $clone = $feat->clone;				
		if ($clone->strand == 0 && $cw == -1) {
		    $clone->_invert_coords;
		    $clone->_invert_feature;
		}
		
		# micro-adjust coords for correct reading frame
		
		if ($clone->feature =~ /INTRON/) {
		    # do nothing. already adjusted to coding frame 
		} elsif ($clone->feature =~ /_1/)  {
		    # used as right (x) feature only
		    until ($clone->coord%$TRIPLET == $frame) {
			$clone->coord($clone->coord-1);
		    }
		} elsif ($clone->feature =~ /_2/)  {
		    until ($clone->coord%$TRIPLET == $frame) {
			$clone->coord($clone->coord+1);
		    }
		}
		
		# add unless micoradjust has moved us off contig
		# does not affect feature annotation (done above)  
		# eg if a TELOMERE_2 gets adjusted off the end of the contig
		# cannot affect 3' orfs -- there are none. 

		$pseudo->add(-object => $clone) unless 
		    $clone->coord > $self->length ||
		    $clone->coord < 0;		
	    }
	    
	    ###########
	    # PLACE STOPS
	    ###########
	    
	    # begin annotation proper 
	    # mark terminii 			
	    
	    my $init = Annotation::Feature # 5' end 
		->new(
		STRAND => $cw,
		COORD => $frame, # not in gene. corrected in _link_feats
		SCORE => 0,
		FEATURE => 'TERM'
		);
	    $pseudo->add(-object => $init);
	    
	    # place all STOP codons 
	    
	    my $coord;
	    for (my $x = 0; ($x+$frame+$TRIPLET) <= $self->length; 
		 $x += $TRIPLET) {                
		
		$coord = $frame + $x;
		my $codon = substr($strand{$cw}, $coord, $TRIPLET);
		next unless $CODONS{$codon} =~ /$args->{'-stop'}/;
		
		my $stop = Annotation::Feature
		    ->new(
		    STRAND => $cw,
		    COORD => $coord + $TRIPLET, # final letter 
		    SCORE => 0,
		    FEATURE => 'STOP'
		    );						
		$pseudo->add(-object => $stop);
	    }               
	    
	    my $term = Annotation::Feature # 3' end 
		->new(
		STRAND => $cw,
		COORD => $coord + $TRIPLET, # not in gene (+1)
		SCORE => 0,
		FEATURE => 'TERM'
		);
	    $pseudo->add(-object => $term);		
	    
	    #print 'STOP', $cw, $frame;
	    #$time = gettimeofday;    
	    
	    ###########
	    # MAKE EXONS
	    ###########
	    
	    # link feature pairs to make exons			

	    my @r = $pseudo->_link_features(
		-strand => $cw,
		-orf_min => $args->{'-orf_min'},
		-exon_min => $args->{'-exon_min'},
		-term_min => $args->{'-term_min'},
		-optimise => $args->{'-optimise'}				
		);

	    # transfer all the created exons to self
	    map { $self->add(-object => $_) } @r;
	    $pseudo->DESTROY; # DESTRYO is calling feature->output on all features. not sure how. 

	} # FRAME LOOP 

	###########
	# MAKE ORFS
	###########

	# all frames have been covered. start sticking things together.
	# coords are converted to absolute coords here.
	
	$self->_link_exons(
	    -strand => $cw,		
	    -intron_max => $args->{'-intron_max'},
	    -intron_min => $args->{'-intron_min'}
	    );

    } # STRAND LOOP

    # remove all unnecessary feaures 
    
    $self->_cleanup_features;
    map { $_->evaluate(-structure => 1, -validate => 1)  } $self->stream;
    $self->index;
    return $self;	
}

# an over-ride method 

sub contains {
    my $self = shift;
    my $obj = shift;
    $self->throw unless $self->_down->isa( ref($obj) ) || ref($obj) =~ /Feature/;
    return( grep {$_ eq $obj} $self->stream );
}

#########################################
#########################################

=head2 merge(-gap => 500, -tandem => .5)

    Walk along contig and merge neighbouring objetcs which look like 
    parts of same gene. Designed for coding genes only. 
    
    Do not merge genes if distnace >= gap.
    Treat as tandems if >= tandem of seqs align.
    Many more options... 

=cut

sub merge {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};

    # find candidates 
    $args->{'-cluster'} = 'homology' unless exists $args->{'-cluster'};
    $args->{'-param'} = 1 unless exists $args->{'-param'};
    $args->{'-distance'} = 0 unless exists $args->{'-distance'};    
    # eliminate candidates 
    $args->{'-align'} = .5 unless exists $args->{'-align'};
    $args->{'-overlap'} = .5 unless exists $args->{'-overlap'};
    $args->{'-gap'} = 500 unless exists $args->{'-gap'};
    # 
    $args->{'-force'} = undef unless exists $args->{'-force'};
    $args->{'-tandem'} = .5 unless exists $args->{'-tandem'};
    $args->{'-overwrite'} = 2.5 unless exists $args->{'-overwrite'};

    my $fh = \*STDERR;
    
    ######################################### 
    # initial clusters by homology. consider only neighbouring genes on same strand. 
    #########################################
    
    my $clusters = ( 
	ref($args->{'-cluster'}) eq 'ARRAY' ? { 1 => $args->{'-cluster'} } :
	$self->cluster(
	    -cluster => $args->{'-cluster'},   # by homology ?
	    -param => $args->{'-param'},       # how much homology. be lax. we recehck later.   
	    -distance => $args->{'-distance'}, # max distance (features except gaps)
	    -strand => 1,	    
	    -frame => 0 # do not enforce frame criterion 
	)
	);
    
    #########################################
    # consider each cluster in turn 
    ######################################### 
    
    foreach my $cl (grep {$#{$_} > 0} values %{$clusters}) {
	
	# order candidates and process pairs left -> right 
	
	my @order = sort {$a->start <=> $b->start} @{$cl};
	map { $self->throw unless $_->up eq $order[0]->up } @order;
	next if scalar(@order) == scalar(grep { $_->assign eq 'REPEAT' } @order); # exlude repeats 

	# ouptut ...

	if ($args->{'-verbose'}) {
	    print {$fh} '>';
	    map { $_->output( -fh => $fh, -recurse => 0 ) } @order;
	}

	my %lookup; # handles merge order issues 
	for my $i (1..$#order) {
	    my $lorf = (exists $lookup{$order[$i-1]->_internal_id} ? 
			$lookup{$order[$i-1]->_internal_id} : $order[$i-1] );
	    my $orf = (exists $lookup{$order[$i]->_internal_id} ? 
		       $lookup{$order[$i]->_internal_id} : $order[$i] );
	    $self->throw unless $orf && $lorf;
	    
	    #print  {$fh} 0, $lorf->name, $orf->name, $orf->distance(-object => $lorf, -bases => 1); 
	    
	    #########################################
	    # do some basic proximity tests
	    #########################################

	    map { $self->throw unless $_->assign eq 'GAP' } $orf->intervening($lorf);
	    next unless $orf->distance(-object => $lorf, -bases => 1) <= $args->{'-gap'};
	    
	    #########################################
	    # find the specific homolog and ensure it is also supported. 
	    #########################################

	    my %count;
	    foreach my $key ( grep {!/YGOB/}  ( keys %HOMOLOGY, 'GENE') ) {
		foreach my $append (undef, '*') {
		    $count{ $key.$append }++ if $orf->data($key.$append) && 
			$orf->data($key.$append) eq $lorf->data($key.$append);
		}
	    }
	    my ($homolog,@other) = sort { $count{$b} <=> $count{$a} } keys %count;	    
	    next unless $homolog && $count{$homolog} > 0;

	    # use SGD if all other things equal 
	    
	    $homolog = ( $count{$homolog} == $count{'SGD'} && 
			 $homolog !~ /Y[A-P][LR]/ ? 'SGD' : $homolog);
	    
	    # gathter contextual information 

	    my $ref = (
		$homolog eq 'YGOB' ?
		$orf->homolog(-object => $lorf, -pillar => $orf->data($homolog)) : 
		$orf->data($homolog)
		);
	    my $pil = $orf->lookup( -query => $orf->data($homolog) );
	    my $exp_len = $pil->{'MEAN'} || 100; # we test onlt for shortness realitive to exp 
	    
	    my $syn = $lorf->synteny(-direction => 'left', -restrict => [$homolog]) +
		$orf->synteny(-direction => 'right', -restrict => [$homolog]) || 0;	    
	    
	    my @names = ($lorf->name, $orf->name);
	    #print {$fh} 2, @names, $homolog, $ref, $syn, $orf->data('SGD'); 
	    
	    #########################################
	    # Generate required test values 
	    #########################################
	    
	    my $frame = $orf->commoncodon(-object => $lorf);

	    my $frag1 = $lorf->fragment; # >0 means fragment (local align > global) 
	    my $frag2 = $orf->fragment;  # >0 means fragment (local align > global) 
	    my $plen1 = 100*$lorf->length/$exp_len;
	    my $plen2 = 100*$orf->length/$exp_len;

	    my $percent_overlap = $orf->overlap(-object => $lorf, -compare => 'seq');
	    my $additive_score = $orf->fragments( -object => $lorf ); # yes or no 	    
	    my $additive_tiling = # -ve for double tiling; +ve for single tiling 
		$orf->exonerate2( # normalized to total tiled region of homolog (untiled excluded)
		    -object => $lorf, 
		    -homolog => $ref,
		    -model => 'local',
		    -return => 'overlap'
		);

	    #########################################
	    # OK. lets start making some decisions. 
	    #########################################

	    # what are the possibilities ? 
	    # 2. fragments of same gene (merge) -> min align, min overlap
	    # -
	    # 1. tandems (no merge, no choose) -> sig align, min overlap 
	    # 6. partial tandem duplication 
	    # -
	    # 3. alt version of same gene (choose not merge) -> sig align, sig overlap 
	    # - 
	    # 4. chance homology 
	    # 5. repeats 
	    
	    my ($call,$reason, $master,$slave);
	    if ( $additive_score || ( $overlap < 0.2 && $additive_tiling > 0.5 ) ) {
		($call,$reason) = qw(merge additive);
		
	        ($master, $slave) = sort { $b->ogid <=> $a->ogid } ($orf, $lorf);
		$reason='fail' unless $master->merge(
		    -object => $slave, 
		    -reference => $ref, 
		    -index => $args->{'-index'} 
		    );

	    } elsif ( $frame ) {
		($call,$reason) = qw(reject altmodel);

		($master,$delta) = $orf->choose(-object => $lorf, -reference => $ref);
		($slave) = ( $master eq $orf ? $lorf : $orf );
		$master->integrate( -object => $slave, -index => $args->{'-index'} );
		
	    } elsif ( $percent_overlap > $args->{'-tandem'} && $syn ) {
		($call,$reason) = qw(reject tandem);
		map { $_->data('TANDEM' => ($_ eq $lorf ? 'right' : 'left' ) ) } ($lorf, $orf);
		
	    } else {
		($call,$reason) = qw(reject passthrough);
		# do not fit together nicely 
		# do not have overlapping codons 
		# do not align to each other (or only a little) 
	    }

	    #########################################
	    # Cleanup 
	    #########################################

	    if ($master && $reason ne 'fail') {
		$lookup{$slave->_internal_id}=$master;
		$slave->DESTROY if $slave;
	    }

	    #########################################
	    # Output 
	    #########################################
	    
	    if ( $args->{'-verbose'} ) {
		$master->output(-fh => $fh, -prepend => [ $call ], -recurse => 0) if $master;
		print {$fh} "$call ($reason):", @names, $homolog, $syn, $exp_len,
		(map { sprintf("%.2f", $_) } 
		 ($plen1, $plen2, $percent_overlap, $frag1, $frag2, 
		  $additive_score,$additive_tiling));
	    }
	}
    }
    
    $self->index;
    return $self;	
}	

=head2 simplify( -object => [], -overlap => undef ) 

    Accept a set of ORFs and return highest-scoring
    non-overlapping subset. Defaults to commoncodon() 
    but will use overlap() to determine overlaps if 
    -overlap defiend. 

    Related to choose() but differs in 2 ways:
    1. the models do not have to be gene variants -- any set OK 
    2. we automatically DESTROY unfavoured models to simplify region 

=cut 

sub simplify {
    my $self = shift;
    my $args = {@_};

    $args->{'-overlap'} = undef unless $args->{'-overlap'};

    $self->throw unless $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->throw unless $_->translatable(-fast=>1) }  @{$args->{'-object'}};

    my $attr = 'global'; # global protein level exonerate score 
    
    ############################################
    # generate all compatible model combinations and compute
    # a _combination_ score as linear sum of constitutent gene scores

    my ($orf,@orfs) = sort { $b->score($attr) <=> $a->score($attr) } @{$args->{'-object'}};
    my @comb = $orf->combinations( 
	-object => \@orfs, 
	( $args->{'-overlap'} ? (-overlap => $args->{'-overlap'} ) : () ) 
	); 

    my @scr;
    foreach my $i ( 0..$#comb ) {
	map { $scr[$i] += $_->score($attr) } @{ $comb[$i] };
    }
    
    # 
    
    my ($hi) = sort { $scr[$b] <=> $scr[$a] } 0..$#scr;
    my %keep = map { $_->_internal_id => 1 } @{$comb[$hi]};

    map { $self->remove(-object => $_); $_->DESTROY; } 
    grep { !$keep{$_->_internal_id} }  @{$args->{'-object'}};
    
    return $self;
}

=head2 purge2( -overlap => 0.8 )

    Replacement for older purge method. 

=cut 

sub purge2 {
    my $self = shift;
    my $args = {@_};

    $args->{'-overlap'} = 0.8 unless exists $args->{'-overlap'};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};
    $args->{'-merge'} = 1 unless exists $args->{'-merge'};
    $args->{'-debug'} = 0 unless exists $args->{'-debug'};
  
    my $fh = \*STDOUT;

    # 0. group genes that share codons -- alternative gene versions 
    # 1. deal with hybrid/compound model clusters -- make discrete 
    # 2. use all members of cluster to identify the best homolog (prefer YGOB)
    # 3. do homology based prediction -- exonerate + genewise 
    # 4. do model selection using homology. Retain gene scores. 

    $self->index;

    ############################################ 
    # group variant models and choose best one 
    ############################################ 
    
    my $clusters = # --> hinges on commoncodon method 
	$self->cluster( -frame => 1, -param => $NONZERO, -cluster => 'overlap' ); 	
    
    # generate all compatible model combinations and compute
    # a _combination_ score as linear sum of constitutent gene scores
    # map { $_->oliver } map {@{$_}}  ( grep { $#{$_}>0 } values %{$clusters} );
    
    map { $self->simplify( -object => $_ ) }  ( grep { $#{$_}>0 } values %{$clusters} );
    
    $self->index;
    
    ############################################ 
    # spurious overlaps 
    ############################################ 
    
    my $clusters = $self->cluster(
	-frame => 0, 
	-stranded => 0, # either strand OK
	-param => $args->{'-overlap'}, 
	-cluster => 'overlap'
	); 

    # this is logically identical to above but incoporates a 
    # parameter for overlap rather than using commoncodon 

    map { $self->simplify( -object => $_, -overlap => $args->{'-overlap'} ) }  
    ( grep { $#{$_}>0 } values %{$clusters} );

    $self->index;
    
    return $self unless $args->{'-merge'} == 1;

    ############################################     
    # merge fragments 	
    ############################################     

    # this is a relatively conservative method as
    # it requires direct neighbours.
    # the contig->merge method is an alternative. 
    
    foreach my $o ( $self->orfs ) {
	next unless my $l = $o->neighbour(-direction => 'left');
	next unless $l->strand == $o->strand;
	next unless $o->homology( -object => $l );	    
	next if $o->assign eq 'REPEAT' || $l->assign eq 'REPEAT';
	# $o->gene =~ /YNL054W-B|YPL060C-A/; # these are Ty elements. occur in multiple broken fragments.
	
	# this is the key step.
	# fragments compares fusion alignment score to a homolog 
	# to the sum of the individual scores 
	
	next unless $o->fragments( -object => $l, -verbose => 1 );  # 1 | 0 | undef   

	if ( $args->{'-debug'} ) {
	    print {$fh} ">>".(++$ds9);
	    $o->ouptut(-fh => $fh, -creator =>1, -prepend => [$ds9] );
	    $o->glyph( -tag => 'current'.$ds9 );
	}
	
	$o->merge( -object => $l );
	$o->oliver if $args->{'-verbose'};
    }
    
    $self->index;
    return $self;
}

=head2 purge(-overlap => .9, -simple => 1, -method => 'length|evalue')
    
    Force choice between overlapping gene models. Consider purge2() instead. 
    
    Choose between two deletion models:
    1) Clustering (default) -- assemble clusters of overlapping orfs 
    by single linkage clustering, then select best gene in each cluster.
    2) Directional (simple=>1) -- sorts genes and make 'safe' one at a time.
    genes that conflict with a 'safe' gene are deleted. the 
    sorting method is crucial in this case. May be any method. 
    
=cut

sub purge {
    my $self = shift;
    my $args = {@_};

    # ordered deletions 
    $args->{'-overlap'} = .9 unless exists $args->{'-overlap'};
    $args->{'-exempt'} = -1 unless exists $args->{'-exempt'};	
    $args->{'-ygob'} = undef unless exists $args->{'-ygob'};
    # -method, -order -> orfs 
    
    # greedy clustering 
    $args->{'-cluster'} = undef unless exists $args->{'-cluster'}; # defaults to ordered if undef 
    $args->{'-param'} = $args->{'-overlap'} unless exists $args->{'-param'}; # to capture accidents 
    $args->{'-reference'} = 'GENE' unless exists $args->{'-reference'}; 
    # -cluster, -param, -frame, -strand, -distance, -bases -> cluster

    # shared 
    $args->{'-force'} = undef unless exists $args->{'-force'};
 
    #######

    # choose between greedy and stingy purging.
    # 1) greedy -- assemble clusters of overlapping orfs 
    # by single linkage clustering, then select one per cluster.
    # This is used for alternative model selection (with -frame argument) 
    
    # 2) stingy -- sort genes and make 'safe' one at a time.
    # genes that conflict with a 'safe' gene are deleted. the 
    # sorting method is crucial in this case.
    
    if ( defined $args->{'-cluster'} ) {
	my %cluster = %{$self->cluster(@_)}; # arrays of overlapping features 
	
	#print '>'.$self->id, scalar($self->orfs), scalar(keys %cluster);
	foreach my $clst (grep {$#{$_}>0} values %cluster) {	    
	    my $pop = pop(@{$clst});
	    my ($best, $delta) = $pop->choose( # make choice using exonerate 
		-object => $clst, 
		-reference => $args->{'-reference'}
		);
	    $self->throw unless $best; # always returns a choice 

	    # select best orf from cluster, re-init _TMP. delete others.  	   

	    map { $best->up->remove(-object => $_); $_->DESTROY;} grep { $_ ne $best } ($pop,@{$clst});	    
	    $best->data('_TMP' => 'delete');
	}
	
    } else {
	
	# uses sorting to establish a sensible order of 
	# eliminations. recommend 'evalue' and a moderate 
	# value for -overlap: .3 - .5
	
	# we collect genes that have been 'made safe' on %safe.
	# init with RNAs etc. and use to eliminate rev strand orfs
	
	my %safe = map { $_->_internal_id => $_ } grep { $_->rank < $args->{'-exempt'} } $self->stream;	
	if ( defined $args->{'-ygob'} ) {
	    map { $safe{$_->_internal_id} = $_ } grep { $_->evalue('ygob') <= $args->{'-ygob'} } $self->stream;
	}
	if ( defined $args->{'-og'} ) {
	    map { $safe{$_->_internal_id} = $_ } grep { $_->ogid } $self->stream;
	}

        #map { $_->oliver(-prepend => ['PRE'])  } $self->orfs;
	
      ORF:foreach my $cand ( grep {! exists $safe{$_->_internal_id} } $self->orfs(@_)) {	  
	DIR: foreach my $dir ( qw(left right) ) {
	    next unless my $neighbour = $cand->$dir;
	    
	    while ( $cand->distance(-object => $neighbour, -bases => 1) <= 20000 ) { # do not want to check whole chr .. 
		my $olap = $cand->overlap( -object => $neighbour );
		if ( exists $safe{ $neighbour->_internal_id } && $olap > $args->{'-overlap'} ) { 
		    # delete and go to next candidate 		      
		    if ( $args->{'-verbose'} ) {
			$neighbour->oliver(-prepend => ['SAFE'], -append => [$olap]);
			$cand->oliver(-prepend => ['PURGE'], -append => [$olap]);
		    }
		    $self->remove(-object => $cand,-warn => 0);
		    $cand->DESTROY;
		    next ORF;
		}
		last unless $neighbour = $neighbour->$dir;
	    } # 
	} # DIR 
	  $safe{$cand->_internal_id} = $cand;
      }	# ORF 
	
        #map { $_->oliver(-prepend => ['POST'])  } $self->orfs;	
    }
    
    $self->index;
    return $self;
}

=head2 cluster(-cluster => 'overlap|homology|commoncodon', -param => .8, 
    -distance => $INF, -bases => $INF, -frame => 0, -stranded => 1,
    -accelerate => 20)

    Cluster genes on a contig by single-linkage clustering. Can cluster using any
    method ( or data('KEY')) supplied as -cluster and using the value supplied on -param. 
    Link formation between features/clusters can be subject to several constraints:
    -distance/-strand/-frame/-bases. 

    Homology/sequence overlap are most likely clustering criteria.     
    tRNAs/GAPS are ignored by default.

=cut

sub cluster {
    my $self = shift;
    my $args = {@_};

    # 1. orf selection (all passed down to orfs()) 
    # start / stop / strand ... but not method. we sort by START.
    $self->throw if $args->{'-method'} || $args->{'-order'};
    # 2. additional orfs to include 
    $args->{'-object'} = undef unless exists $args->{'-object'};
    # 3. clustering conditions 
    $args->{'-cluster'} = 'overlap' unless exists $args->{'-cluster'};
    $args->{'-param'} = .8 unless exists $args->{'-param'};
    # 4. additional constraints :: off by default 
    $args->{'-distance'} = $INFINITY unless exists $args->{'-distance'};
    $args->{'-bases'} = $INFINITY unless exists $args->{'-bases'};
    $args->{'-frame'} = ($args->{'-cluster'} eq 'commoncodon' ? 1 : 0) unless exists $args->{'-frame'};

    # 5. additional constrains :: on by default 
    $args->{'-stranded'} = 1 unless exists $args->{'-stranded'};
    # 6. technical parameters  
    $args->{'-accelerate'} = 20 unless exists $args->{'-accelerate'};

    # 

    my $method = $args->{'-cluster'};
    my $key = '_CLUSTER_VAR';
    my $attr = '_INDEX_VAR';
    my $flag = '_FLAG_STORE_CTG'; # could be N diff contigs 

    ######################################
    # get orfs and index for accelearated search 
    ######################################

    # are there extra genes ? 
    # label so we can temporarily place them on 
    # the contig array and later remove 

    my @current = $self->orfs(@_);
    #print {STDERR} map { $_->assign } @current;

    my @extra;   
    my $base = $self->down;
    foreach my $o ( @{ $args->{'-object'} } ) {
	$self->throw unless $base->isa(ref($o));
	$o->data($flag => $o->up);
	$o->transfer( -from => $o->up, -to => $self, -fast => 1, -warn => 0 );
	push @extra, $o;
    }

    # 1. we pass down -start, -stop, -strand, -method, -order etc to orfs().     
    # 2. then sort all the ORFs. 
    # 3. index the orfs and do some QC 

    $self->index; # index by start order. must then sort by index below as  
    my @orfs =    # genes with same start can be returned out of order. 
	sort {$a->id <=> $b->id} (@current, @extra);
    foreach my $i ( 0..$#orfs ) {
	$orfs[$i]->data($key => undef); 
	$orfs[$i]->data($attr => $i); 
    }

    ######################################
    # compare all pairs of genes and make groups 	
    ######################################

    my $clusterid;
    my %cluster;
    foreach my $orf ( @orfs ) {
	$self->throw if $orf->data($key); # already on a cluster 
	$self->throw if exists $cluster{ $orf->_internal_id }; # ditto 
	
	$orf->data( $key => 'c'.(++$clusterid) ); # name new cluster 
	$cluster{ $orf->data($key) } = [$orf];    # make new cluster
	$cluster{ $orf->_internal_id } = $orf;    # object tracker. does not change.  	  
	
	######################################
	# Acceleration: determine local search space 
	######################################

	my ($start,$stop)=(0,$#orfs);
	if ( $args->{'-accelerate'} && 
	     ($args->{'-cluster'} eq 'overlap' || $args->{'-cluster'} eq 'commoncodon') ) {
	    ($start,$stop) = map { $orf->data($attr)+($_*$args->{'-accelerate'}) } (-1,1);
	    $start = 0 unless $start >= 0;
	    $stop = $#orfs unless $stop <= $#orfs;
	}

	my @searchspace = 
	    grep { $orf ne $_ } grep { exists $cluster{$_->_internal_id} } @orfs[$start..$stop];
	
	# can we link to existing cluster?	
	
	foreach my $x ( @searchspace ) { # non-self ORFs in clusters	    
	    
	    ######################################
	    # CLUSTER 0: trivial case ..
	    ######################################
	    
	    next if $orf->data($key) eq $x->data($key); # already on same cluster 

	    ######################################
	    # CONSTRAINT 1 : strand 
	    ######################################

	    next unless ! $args->{'-stranded'} || $x->strand == $orf->strand;

	    ######################################
	    # CLUSTER 1: veryify cluster data structure 
	    ######################################

	    $self->throw unless $cluster{$x->_internal_id} eq $x;
	    $self->throw unless my $new_cl = $x->data($key);
	    $self->throw unless ref($cluster{$new_cl}) eq 'ARRAY';
	    $self->throw unless ref($cluster{$new_cl}->[0]) eq ref($orf);
	    $self->throw unless $cluster{$new_cl}->[0]->data($key) eq $new_cl;

	    ######################################
	    # CLUSTER 2: test clustering condition (-cluster/-param)
	    ######################################

	    my $v = 0;
	    if ($orf->can($method)) { # homology|overlap
		$v = $orf->$method(-object => $x);
	    } elsif ($self->data($method)) {
		$v = 1 if $self->data($method) == $x->data($method);
	    } else {$self->throw($orf->debug);}
	    next unless $v >= $args->{'-param'};	    
	    #print $method, $v, $args->{'-param'};
	    
	    ######################################
	    # CONSTRAINT 2 : frame
	    ######################################

	    # commoncodon is slow so we try avoid.. do last.

	    next unless ( (! $args->{'-frame'}) || $x->commoncodon($orf) );

	    ######################################
	    # CONSTRAINT 3 : distance 	    
	    ######################################

	    # distance constraint ? allow to merge across certain objects?

	    next unless 
		scalar(grep { $_->assign ne 'GAP' } $orf->intervening(-object => $x))
		<= $args->{'-distance'};
	    next unless # mainly used for HSP clustering 
		$orf->distance(-object => $x, -bases => 1) <= $args->{'-bases'};	        

	    ######################################
	    # CLUSTER 3: just $orf to deal with or does it have friends? 
	    ######################################

	    my $old_cl =  $orf->data($key);
	    map { $_->data($key => $new_cl) } @{$cluster{ $old_cl }};
	    push @{ $cluster{$new_cl} }, @{$cluster{ $old_cl }};
	    delete $cluster{$old_cl};
	}
    }

    ######################################
    # Unwind everything ....
    ######################################

    map { delete $cluster{$_} } grep {! /c/} keys %cluster;
    map { $_->transfer( -from => $self, -to => $_->data($flag) ) } @extra; # restore to contig...
    foreach my $o (@orfs) {
	map { $o->data( $_ => 'delete' ) } ($flag,$key,$attr);
    }

    $self->index;
    return \%cluster;
}

#grep {defined} map { ( (! $args->{'-rank_max'} || $_->rank <= $args->{'-rank_max'}) ? $_ : undef ) }
#grep {defined} map { ( (! $args->{'-strand'} || $_->strand == $args->{'-strand'}) ? $_ : undef ) }

sub _common_codon_cluster {
    my @collect = @{$_[0]};

    my %clst;
    for my $i (0..$#collect-1) {
	$clst{ $collect[$i]->_internal_id } = [ $collect[$i] ] unless 
	    exists $clst{ $collect[$i]->_internal_id };
	for my $j ($i+1..$#collect) {
	    $clst{ $collect[$j]->_internal_id } = [ $collect[$j] ] unless 
		exists $clst{ $collect[$j]->_internal_id };
	    next unless  $collect[$i]->commoncodon(-object =>  $collect[$j]);
	    foreach my $mv ( @{$clst{ $collect[$j]->_internal_id }} ) {
		next if $clst{ $mv->_internal_id } eq $clst{ $collect[$i]->_internal_id };
		push @{$clst{ $collect[$i]->_internal_id }}, $mv;  
		$clst{ $mv->_internal_id } = $clst{ $collect[$i]->_internal_id };		      
	    }
	}
    }

    return \%clst;
}

=head2 _validate_assembly_gaps()

    Examine all GAP objects and ensure they are properly documented
    and their borders are true. We do ab initio GAP detection to 
    guarantee complete detection. 

    To ensure GAPs have correct properties, we correct any issues with 
    strand (0 for GAPs) and GAP type. 

=cut 

sub _validate_assembly_gaps {
    my $self = shift;
    my $args = {@_};

    $args->{'-gap_len'} = 10 unless $args->{'-gap_len'}; # what is NCBI number? 

    ######################################
    # 0. vars 
    ######################################

    my $fh =  \*STDERR;

    # important for control over start/stop 

    map { $_->strand( 0 ) } grep { $_->assign =~ /GAP|FEATURE/ } $self->stream;

    ######################################
    # 1. look for mislabelled gaps 
    ######################################
    # Is the order right here? 

    foreach my $gap ( grep { $_->assign ne 'GAP' }  $self->stream ) { 
	my $gap_type;
	if ( $gap->sequence =~ /^N{100}$/ ) {
	    $gap_type='Homology - Unknown length';
	} elsif ( $gap->sequence =~ /^N+$/ ) {
	    $gap_type='Mate pairs - Estimated length';
	}
	
	# 
	
	if ( $gap_type ) {
	    $gap->output(-prepend => [ $self->_method, __LINE__, $gap_type ], 
			 -fh => $fh );
	    $gap->throw;
	}
    }

    ######################################
    # 2. remove erroneus gaps + fix borders 
    ######################################

    my @remove;
    foreach my $gap ( grep { $_->assign eq 'GAP' }  $self->stream ) { 

	my $path;
	if ($gap->sequence !~ /N/) {
	    $path='excess';
	    push @remove, $gap;

	} else {
	    $path='border';
	    my ($start,$stop)=($gap->start, $gap->stop);

	    # greedy expand -- end with one base excess 
	    
	    $gap->down->start( -adjust => -1 )
		until $gap->first_base ne 'N' || $gap->start == 1;
	    $gap->down->stop( -adjust => +1 )
		until $gap->last_base ne 'N' || $gap->stop == $self->length;
	    
	    # greedy retract -- should end in the right place 
	    
	    $gap->down->start( -adjust => +1 )
		until $gap->first_base eq 'N';
	    $gap->down->stop( -adjust => -1 )
		until $gap->last_base eq 'N';
	    $self->throw if $gap->sequence =~ /[^N]/;

	    # logging 

	    if ( $gap->start != $start || $gap->stop != $stop ) {
		$gap->history( $self->_method.':'.$path );
	    }
	}
	
	# pretty pretty 
	
	$gap->output(
	    -prepend => [ $self->_method, __LINE__, 'CURRENT' ], 
	    -fh => $fh, -append => [$path, $gap->_top_tail] 	    
	    ) if $args->{'-verbose'};
    }
    map { $self->remove( -object => $_ ) } @remove;
    $self->index;
    
    ######################################
    # 3. look for gaps that were missed completely 
    # -- by reannotating all gaps above size x..
    ######################################

    my @fake_gaps = ();
    my @seq = split(//, $self->sequence);

    my $seq;
    my $changepoint=1;
    for my $i (0..$#seq) {
	
	if ( ! $seq ) { #initialize 
	    # set the change point to 1 and sequence to $seq[0] (below)

	} elsif ( $seq =~ /$seq[$i]/ || # more of the same 
		  ($seq[$i] ne 'N' && $seq !~ /N/) ) { # short DNA sequences / polyX tracts 
	    $seq .= $seq[$i];
	    next;
	    
	} elsif ( $seq[$i] ne 'N' ) { # transition -- ending a GAP seqeunce 
	    if ( length( $seq ) >= $args->{'-gap_len'} ) {
		my $fake_gap = ref($self->down)
		    ->new(
		    START => $changepoint,
		    STOP => $i,
		    STRAND => 0,
		    UP => $self
		    );
		
		$fake_gap->data('NNNN' => $INFINITY);
		$fake_gap->evidence('NNNN');
		$fake_gap->evaluate(-structure => 0, -validate => 0);
		$self->throw( $fake_gap->assign ) unless $fake_gap->assign eq 'GAP';
		
		$self->add( -object => $fake_gap );		
		$fake_gap->output( 
		    -prepend => [ $self->_method, __LINE__, 'FOUND' ], 
		    -fh => $fh, ) if $args->{'-verbose'};
		push @fake_gaps, $fake_gap;
	    }

	} elsif ( $seq[$i] eq 'N' ) {
	    # transition : ending a non-gap --> No action 
	    # print {$fh} ">$changepoint - $i\n$seq";

	} else { $self->throw(); }

	$changepoint = $i+1;
	$seq = $seq[$i];
    }
    $self->index;
    
    ######################################
    # 4. Make non-redundant set 
    ######################################

    my $gap_overlaps = $self->cluster(
	# clustering 
	-cluster => 'overlap',
	-param => .001,
	# orfs params 
	-start => 1,          # order by start coord
	-noncoding => 1,      # enable non-coding sequences
	#-rank => -$INFINITY, # enforce non-coding ...
	-restrict => ['GAP']  # enforce GAP // exclude FEATURE 
	);
    
  GAP: foreach my $cl ( grep { $#{ $_ } > 0 } values %{ $gap_overlaps } ) {
      
      # we take the maximum GAP span. 
      
      map { $self->throw unless $_->exons==1 } @{ $cl };
      my ($start) = map { $_->start } sort { $a->start <=> $b->start } @{ $cl };
      my ($stop) = map { $_->stop } sort { $b->stop <=> $a->stop } @{ $cl };
      my $gap = shift( @{ $cl } );

      # logging 

      if ( $gap->start != $start || $gap->stop != $stop ) {
	  $gap->history( $self->_method.':merge' );
      }

      # execution 

      $gap->start( $start );
      $gap->stop( $stop );
      map { $self->remove( -object => $_ ) }  @{ $cl };
      map { $_->DESTROY() }  @{ $cl };
  }
    $self->index;

    ######################################
    # Ugly print .... 
    ######################################
    
    if ( $args->{'-verbose'} )  {
	map {$_->output(-debug => 1, -fh => $fh, -append => ["\n".$_->sequence])} 
	grep {$_->assign eq 'GAP' || $_->sequence =~ /N/ } $self->stream;
    }

    ######################################
    # Housekeeping .... 
    ######################################
   
    map { $_->strand( 0 ) } grep {$_->assign eq 'GAP'} $self->stream;    
    return $self;
}

=head2 _validate_overlapping_features()

    Run two procedues to ensure overlaps are consistent with expectations: 
    1. Find and reconcile REDUNDANT *coding* features 
    2. Cluster all overlapping features and attempt to merge 

=cut 

sub _validate_overlapping_features {
    my $self = shift;
    my $args = {@_};
    
    my $fh = *STDERR;

    my @out = ('-fh' => $fh, '-debug' => 1, '-og' => 1, '-recurse' => 0);

    ################################
    # 
    ################################

    $self->throw unless my $index = $args->{'-index'};
    $self->throw unless ref($index) =~ /HASH/;

    ################################
    # 1. Find and remove REDUNDANT CODING features 
    # Based on commoncodon test so should be identical (+ errors) for the most part 
    ################################

    my $clref = $self->cluster( 
	-cluster => 'commoncodon',
	-stranded => 1,
	-noncoding => 0, # enforce CODING via orfs() call in cluster() 
	-accelerate => 20
	);
    
    foreach my $cl ( values %{$clref} ) {
	my ($pop, @others) = @{ $cl };       
	next unless @others; 
	    
	# choose best seq...
	
	my ($best,$scr) = $pop->choose( 
	    -object => \@others, 
	    -verbose => $args->{'-verbose'},
	    -strict => 0
	    );
	
	# ouptut 
	
	if ( $args->{'-verbose'} >= 2 ) {
	    print {$fh} "\n--------------------------\n";
	    map { $_->output(-prepend => ['OPTION'], @out, 
			     -append => [$best->overlap(-object => $_) ] ) } @{$cl}; 
	    $best->output(-prepend => ['BEST'], @out);	    
	    print {$fh} undef;
	}    
	
	# manage orthogroups 
	
	my ($other_og, @ogid) = grep { $_->ogid } grep { $_ ne $best} @{ $cl };
	$self->throw if @ogid;
	
	my $path; # ='Best?/Other-';
	if ( $other_og ) {
	    
	    if ( $best->ogid ) {
		$path='collide'; #Best+/Other+
		
		# 1. Both are in OGs
		
		my $og1 = $index->{ $best->ogid };
		my $og2 = $index->{ $other_og->ogid }; 
		map { $self->throw unless $best->isa( ref($_->[0]) ) } ( $og1, $og2 );
		
		# Iterate through species / OGs and show relationships between models
		# Shit's ugly. 

		if ($args->{'-verbose'}) {		
		    foreach my $pos ( 0..$#{$og1} ) {
			map { $_->[$pos]->output(
				  -prepend => [ $_->[$pos]->organism ], @out,
				  -append => [ ($_->[0] eq $best ? '**' : '') ] )
			} ( $og1, $og2 );
			print { $fh } (
			    $og1->[$pos] eq $og2->[$pos] 
			    ? ('Same','Same') 
			    : ($og1->[$pos]->commoncodon( $og2->[$pos] ), 
			       $og1->[$pos]->overlap( -object => $og2->[$pos] )) 
			),"\n";
		    }
		}
		
		# Use collide to build best OG 
		
		$og1->[0]->collide( -object => $og2->[0], -cleanup => 1 );
		
	    } elsif ( ! $best->ogid ) {
		
		# 2. $best is not in an OG. Either use other_og or port over OG relationships
		
		if ( ! $scr ) { # $scr=0 indicates a random choice by choose. just swap.  
		    $path='swap'; #"Best-/Other+ (other by OG) > swap";
		    $best = $other_og;
		} else {
		    $path='integrate'; # "Best-/Other+ (best by score: $scr) > integrate()";
		    $best->integrate( -object => $other_og, @_ );
		}
	    } else { $self->throw; }

	}#other_og 
	
	# delete redundant genes
	
	foreach my $del ( grep { $_->up } # may have delelted in the call to collide() above 
			  grep { $_ ne $best } @{$cl} ) {
	    $del->output(@out, -recurse => 0, -prepend => ['DEL']) 
		if $args->{'-verbose'} >= 2;		
	    $self->remove( -object => $del );
	    $del->DESTROY;
	}

	# logging 

	$best->history( $self->_method.':'.$path ) if $path;
	print {$fh} "\n$path" if $args->{'-verbose'};	

    }# cluster loop

    ################################
    # 2. Handle any other overlapping features 
    # Note: Overlaps not expressly forbidden 
    ################################

    my $clref = $self->cluster(
	-cluster => 'overlap',
	-stranded => 0,
	-frame => 0,
	-accelerate => 20,
	-param => 0.2
	);
    
    my $clx;
    foreach my $cl ( grep { $#{$_} >0 } values %{$clref} ) {
	map { $_->output( -fh => $fh, -debug => 1 ) } @{$cl}
	if $args->{'-verbose'} >= 2;
	
	$self->merge(
	    -cluster => $cl,
	    -index => $args->{'-index'}, 
	    -verbose => $args->{'-verbose'}
	    );
    }
    
    $self->index;
    return $self;
}

=head2 _genbank_gap_overlaps


=cut 

sub _genbank_gap_overlaps { 
    my $self = shift;
    my $args = { @_ };

    #################################	        
    # Basic QC 
    #################################	    

    $self->throw unless my $index = $args->{'-index'};
    $self->throw unless ref($index) =~ /HASH/;
    
    #################################	        
    # Some vars. 
    #################################	    

    my $fh = *STDERR;
    my $loop=1;
    my %path;
    my $old_path_internal = $path{ 'internal' };

    #################################	        
    # Set up comparisons: GAPS vs all else 
    #################################	    

  GAP: foreach my $gap ( grep { $_->assign eq 'GAP' } $self->stream ) {
      #next GAP if $gap->genbank_exclude;
      
      $gap->output(
	  -debug => 1, 
	  -prepend => [$self->_method, __LINE__, "GAP", $gap->_top_tail], 
	  -fh => $fh, -recurse => 0, -append => [$gap->genbank_exclude]	  
	  ) if $args->{'-verbose'}; #=1;
      
      foreach my $nei ( $gap->context( -distance => 10, -all => 1, -self => -1 ) ) {
	  next if $nei->genbank_exclude;
	  
	  ################################
	  # identify overlaps we need to deal with  
	  ################################

	  next unless my $olap = $gap->overlap( -object => $nei, -compare => 'gross' );
	  my $bad_trans_ratio = ( $nei->composition( -aa => 'X' ) >= $GENBANK_GAP_LIMIT ? 1 : 0);

	  $nei->output(
	      -debug => 1, 
	      -prepend => [$self->_method, __LINE__, "TEST", $nei->_top_tail], 
	      -fh => $fh, -recurse => 0,
	      -append => [$olap, $bad_trans_ratio]
	      ) if $args->{'-verbose'};
	  
	  # Not worth fighting for Celine's centromeres 
	  $nei->genbank_exclude( $self->_method ) and next if 
	      $nei->evidence eq 'MANUAL';
	  
	  ################################
	  # There are three options here : 
	  ################################
	  # 1. It is a partial overlap -->try tweak gene coordaintes 
	  # -- We are not careful about how much we trim 
	  # -- The only guarnatee is that you get something translateable
	  # -- We assume that _make_genbank_gene_terminii() cleans up the mess. 
	  # 2. Full overlap that fails Genbank criteria for GAP overlap 
	  # -- We adjust coords or split gene as required 
	  # 3. Full overlap that passes Genbank criteria for GAP overlap 
	  # -- We do not try to change. Just add a note.

	  my $path;
	  if ( $olap < 1 ) {
	      
	      if ( $gap->start < $nei->start && $gap->stop < $nei->stop ) {
		  $path='start';# gap overlaps start of gene 
		  
		  my $edit; 
		  foreach my $exon (  sort { $a->start <=> $b->start} $nei->stream) {
		      $edit = $exon and last if $exon->stop > $gap->stop;
		      $nei->remove( -object => $exon );
		      $exon->DESTROY;
		  }
		  $edit->start(-adjust => +1) until 
		      $edit->start > $gap->stop && $nei->length%3==0;
		  
	      } elsif ( $gap->start > $nei->start && $gap->stop > $nei->stop ) {
		  $path='back';# gap overlaps backend of gene 

		  my $edit;
		  foreach my $exon (sort { $b->start <=> $a->start} $nei->stream) {
		      $edit = $exon and last if $exon->start < $gap->start;
		      $nei->remove( -object => $exon );
		      $exon->DESTROY;
		  }
		  $edit->stop(-adjust => -1) until 
		      $edit->stop < $gap->start && $nei->length%3==0;
		  
	      } else {$self->throw;}
	      
	  } elsif ($gap->start < $nei->start || $gap->stop > $nei->stop) {
	      map { $_->output(-fh => $fh, -prepend => [$self->_method, __LINE__] ) } ($nei,$gap);
	      map { print {$fh} $_->sequence(); }  ($nei,$gap);
	      $self->throw("Gap encompasses feature");
	      
	  } elsif ( $gap->gap_type == 2 || $bad_trans_ratio ) { 
	      $path='internal'; # complete overlap of gap in coding 
	      
	      # cannot cross GAP of unknown length or have >50% in GAP of known length 
	      # http://www.ncbi.nlm.nih.gov/genbank/wgs_gapped#both
	      
	      # eliminate all Ns from coding gene 

	      $nei->excise_gaps( -override => 1 ); # ugly hardcode but OK 
	      
	      # collect all the exons on either side of the gap 
	      # and split gene into to separate fragments 

	      my @left = grep { $_->stop <= $gap->start } $nei->stream;
	      my @right = grep { $_->start >= $gap->stop } $nei->stream;
	      #print { STDERR } $self->_method, __LINE__, @left+0, @right+0, $nei->exons;
	      $nei->throw(-output => 1) unless @left+@right == scalar($nei->stream);
	      $nei->throw(-output => 1) unless 
		  scalar(keys %{{ map { $_ => 1 } (@left,@right) }}) == scalar($nei->stream);
	      
	      # the longer (presumably better supported) fragment keeps OG etc. 
	      # the other exons become the $new ORF object. 
	      
	      my ($left,$right);
	      map { $left += $_->length } @left;
	      map { $right += $_->length } @right;
	      my ($new) = $nei->fission( -exons => [ ($left < $right ? @left : @right) ] );
	      $self->index;
	  } else { 
	      $path='tolerate'; 
	      $nei->genbank_note('Coding sequence crosses assembly gap');
	  }
	  
	  $path{ $path }++;
	  $nei->history( $self->_method.':'.$path );

	  $nei->output(
	      -prepend => [$self->_method, __LINE__, $path.$loop, $nei->_top_tail], 
	      #-append => [$nei->genbank_exclude || 0, $nei->sequence],
	      -recurse => 0, -fh => $fh 
	      ) if $args->{'-verbose'};
      }
      print { $fh } "\n===========================" if $args->{'-verbose'}
  }
    $self->index;

    ################################################
    # Iternate as needed 
    ################################################

    # we need to run until steady state as new genes are created by fission 
    # and will not be in initial test array. 
    
    if ( $path{ 'internal' } > $old_path_internal ) {
	$old_path_internal =  $path{ 'internal' };
	$loop++;
	goto GAP;
    }

    ################################################
    # Wrap up 
    ################################################

    #print { $fh } '>', ( map { $_.":".($path{$_} || 0) } keys %path ); 
    map { $_->output(-fh => $fh) and  $_->throw unless $_->translatable } grep { $_->coding } $self->stream;
    
    $self->index;
    return $self;
}

=head2 _genbank_gene_terminii()

    There is an important decision to be made here about how to implement 
    Genbank compliance. Should we calc coords when this method is called
    and update all the features? or compute an alternative set of coords
    and store those separately (to be used in over-ride manner by genbank_tbl())? 
    Or should we write just-in-time methods that calc adjusted coords on the fly 
    but never alter original and are never stored? 

    Implementing #2 but have doubts about design.

=cut 

sub _genbank_gene_terminii {
    my $self = shift;
    my $args = {@_};

    ################################# 
    # Defaults 
    #################################    
    
    $args->{'-index'}=undef unless exists $args->{'-index'};
    $args->{'-verbose'}=0 unless exists $args->{'-verbose'};
    
    ################################# 
    # Set some variables 
    #################################    

    my $fh = \*STDERR;
    
    ################################# 
    # Main loop...
    #################################    
        
    # We run 2 sets of tests -- for START codon and for STOP codon.
    # M more interesting than * . True start could be :
    # 1. upstream of initial location 
    # 2. downstream of initial location 
    # 3. we might be starting the search in the wrong place completely
    # due to challenging upstream exon or addition of 'exon' to 
    # step over frameshift. 
    # By contrast it is fairly trivial to always find a STOP codon .. 

    foreach my $orf ( grep { $_->coding( -pseudo => 1 ) } $self->stream) {	
	# next unless $orf->up; # fixed below. No longer needed. 	

	#################################	    
	# Debug nonsense 
	#################################

	#next unless $orf->ygob eq 'Anc_8.270';
	#$orf->output();
	# $args->{'-verbose'} = ( $orf->start > 260000 && $orf->stop < 280000 ? 2 : 0);

	my ($start_mem, $stop_mem) = ($orf->start, $orf->stop);
	
	#################################	    
	# Useful test cases 
	#################################
	
	# 5' degraded in Sbay -- we find good ORF 
	# next unless $orf->gene eq 'YOL053W'; 	
	# Fully degraded in Sbay -- we shrink to absurdity 
	# Not really properly resolved. 
	# Points to need for a proper stopping condition when shrinking. 
	# next unless $orf->gene eq 'YEL035C'; 

	################################# 
	# First codon 
	my ($poison, $force_reoptimise)=(0,0);
      FINDSTART:
	#################################
	
	unless ($orf->translatable) {
	    $orf->output(-prepend => [$self->_method, __LINE__], -recurse => 1, -fh => $fh, -creator => 1);
	    $self->throw;
	}

	# get surrounding features to conduct necessary tests
	
	my $upstream =
	    $orf->neighbour(-direction => ($orf->strand==-1 ? 'right' : 'left'));	
	my $gap = 
	    $orf->neighbour(-direction => ($orf->strand==-1 ? 'right' : 'left'), -gap => +1); 
	$self->throw if $gap && $gap->assign ne 'GAP';

	# start tests 
	
	my $fex = $orf->firstexon;
	unless ( $orf->first_codon eq $START_CODON ) {	    
	    $orf->output(
		-prepend => [$self->_method, __LINE__, $orf->_top_tail], 
		-recurse => 0, -fh => $fh ) if $args->{'-verbose'};
	    
	    my $path=undef;

	    #################################	    
	    # 1A - look upstream for M 
	    #################################
	    
	    $fex->start( -adjust => -$TRIPLET, -R => 1 ) until 
		( $orf->first_codon eq $START_CODON ||            # OK outcome 
		  $orf->_terminal_dist2( $orf->strand*-1 ) < 3 || # Almost OK (contig end) 
		  ($gap && $orf->overlap( -object => $gap )) ||   # Almost OK 
		  $orf->first_codon =~ $STOP_CODON);              # NOT OK 		  

	    #################################	    
	    # 1B - Look for downstream M 
	    #################################
	    # It would neater to do 1B before 1A but we prefer upstream M 
	    # E.g., where *XXXMXXX+XXXXMXXXXXXX* and + is starting position. 
	    
	    my $mark = $fex->start( -R => 1 );
	    $fex->start( -adjust => $TRIPLET, -R => 1) until 
		( $orf->first_codon eq $START_CODON || $fex->length < 4 );

	    if ($orf->first_codon eq $START_CODON ) { 
		$path='move';
		$orf->history( $self->_method.':'.$path );
		goto FINDSTOP;
	    } else { $fex->start( -R => 1, -new => $mark ); } # failed so revert 
	    
	    #################################
	    # 1C - deal with failures -- contig ends, GAPs, STOP etc 
	    #################################
	    
	    my $start_gbk;
	    if ( $gap && $orf->overlap( -object => $gap ) ) {
		$self->throw if $gap->start <= $fex->start && $gap->stop >= $fex->stop;

		$fex->start( -adjust => +1, -R => 1) until 
		    ! $orf->overlap(-object => $gap) || $fex->length == 1;
	        $start_gbk = $fex->start( -R => 1 );
		$fex->start( -adjust => ( $fex->length < 3 ? -1 : +1 ), -R => 1) until 
		    $orf->length % $TRIPLET == 0;
		$path='gap';
		
	    } elsif ( $orf->_terminal_dist2( $orf->strand*-1 ) == 0 ) {
		$start_gbk = $fex->start( -R => 1 ); # do nothing 
		$path = 'contig_end';
		
	    } elsif ( $orf->_terminal_dist2( $orf->strand*-1 ) < 3 ) {
		$start_gbk = ($orf->strand == 1 ? 1 : $self->length);
		$path = 'contig_end2';

	    } elsif ( $orf->first_codon =~ $STOP_CODON )  {
		$fex->start( -adjust => $TRIPLET, -R => 1);		
		
		##############################################################
		# Merge possible here but not appropriate in _genbank_gene_trminii() 
		# Functionality has been moved. 
		##############################################################

		$path='stop';

		# Run reoptimize (really exonerate) with relaxed params 
		# where we think missing intron may be the problem. 
		# Mising intron fail mode is indicated by uninterrupted reading frame. 
		# Interrupted reading frame usually a sign of pseudgenization. 	       

		if ( $orf->first_codon ne $START_CODON ) {
		    $orf->structure();

		    # pretty pretty 
		    
		    if ($args->{'-verbose'}) {
			#$orf->glyph( -print => 'SGD', -tag => 'pre' );		    
			$orf->output(
			    -prepend => [$self->_method, __LINE__, $orf->_top_tail], 
			    -append => [$orf->data('STOP'), $force_reoptimise], 
			    -recurse => 1, -fh => $fh
			    );
		    }

		    # We run reoptimize if either of two conditions hold : 
		    # 1. we have been told to do it ($force_reopt) 
		    # becuase the rest of the gene is a mess. 
		    # 2. we have iterated -pruned- through each 5' exon and 
		    # and there is only a single exon left. 
		    # What we do not have is a nice way to say stop pruning. 
		    # i.e., the exons are more important than finding M. 
		    # IF we do ultimately run reoptimise() this should be less 
		    # of an issue because we should reover whetever the best structure is. 

		    # The fun here is on things like S. bayanus YOL053W
		    # The 5' part of the gene has degraded in Sbay (or extended in Scer?) 
		    # So we initially find 5' segments that have homology (and score!)
		    # but which are not part of the true cding structure. 
		    
		    if ( $orf->data('STOP') == 0 || $force_reoptimise ) {
			my ($best) = sort {$orf->logscore($a) <=> $orf->logscore($b)} qw(gene sgd);
			if ($best && $orf->logscore($best) < -10) {
			    $orf->reoptimise(
				-reference => $orf->data( uc($best) ),
				-atg => undef, # let's not complicate 
				-hmm => ($orf->evalue('ygob') < 1e-5 ? $orf->hit('ygob') : undef),
				-intron => 1, # make introns favorable 
				-filter => 5, # but not to many introns ... coarse 
				-verbose => 0
				) unless $poison;
			}
			$poison++;
			
		    } else {

			# Iterative test exons to find an M. 
			# We prune one (or the the minimum number of) exons at each step.
			# Then send back to FINDSTART to try find an M. 
			# If it fails we prune another exon. 
			# If we fail to remove an exon, we go directly to reoptimise()
			# We can fail to remove an exon because : 
			# 1) The gene is a mess with no hope of creating a sensible ORF. 
			# 2) The exon we want to strip is linked by a well supported intron. 

			$force_reoptimise=0;
			my @removed_exons = ( $orf->prune_fiveprime_exon() );
			$force_reoptimise=1 unless @removed_exons;
			
			# what do we do if we messed up? restore exons? 
		    }		    
		    $orf->structure();
		    
		    # pretty pretty
		    
		    if ( $args->{'-verbose'} ) {
			# $orf->glyph( -print => 'SGD', -tag => 'reopt' );
			$orf->output(
			    -prepend => [$self->_method, __LINE__, $orf->_top_tail], 
			    -recurse => 0, -fh => $fh);
		    }

		    if ($orf->first_codon ne $START_CODON && $poison<=1) {
			$orf->history( $self->_method.':'.$path );
			goto FINDSTART;
		    }
		}
	    } else { $self->throw( $orf->first_codon ); } 

	    $orf->update();
	    $fex->genbank_coords( -mode => 'start', -R => 1, -set => $start_gbk ) if $start_gbk;
	    $orf->history( $self->_method.':'.$path );
	}
	
      FINDSTOP:
	
	$orf->output(
	    -prepend => [$self->_method, __LINE__, $orf->_top_tail], 
	    -recurse => 1, -fh => $fh ) if $args->{'-verbose'};

	################################# 
	# Last codon 
	#################################
	
	# get downstream gap (if there is one) 
	# this is needed to suport one of the tests below 
	
	my $gap = $orf->neighbour( 
	    -direction => ( $orf->strand==-1 ? 'left' : 'right'),
	    -gap => +1 # require 
	    );
	
	if ( $gap ) {

	    # hygiene criteria. 

	    $self->throw if $gap->assign ne 'GAP';
	    $self->throw if $gap->spans( $orf ); 
	    #$gap->output( -fh => $fh ) and $orf->throw( -output => 1 ) if 
	    #	$orf->spans( $gap ) && # should be clean ( e.g., excise_gaps() ) 
	    #	$orf->overlap( -object => $gap, -bases => 1 ) > 6; # but allow messy borders
	    
	    $gap->output(
		-prepend => [$self->_method, __LINE__, $gap->_top_tail], 
		-recurse => 1, -fh => $fh ) if $args->{'-verbose'};
	}

	# OK. get into loop and start expansion  
	
	my $lex = $orf->exons( -query => 'last' );
	unless ( $orf->last_codon =~ $STOP_CODON ) {	    
	    
	    # adjust last codon to downstream STOP (or other stopping condition)
	    
	    $lex->stop( -adjust => +$TRIPLET, -R => 1 ) until 
		( $orf->last_codon =~ $STOP_CODON ||            # OK outcome 
		  $orf->_terminal_dist2( $orf->strand ) < 3 ||  # almost OK outcome (must be exact) 
		  ($gap && $gap->gap_type==2 && # Gap of N=100
		   $orf->overlap( -object => $gap )) || # NOT OK by Genbank 
		  ($gap && $gap->gap_type==1 && # Other gap
		   $orf->overlap( -object => $gap ) && $orf->composition( -aa => 'X' ) >= 0.5 ) # NOT OK by Gbk
		  # ^ need to fix this below as cannot have stop in gap region 
		);
	    
	    $orf->output(
		-prepend => [$self->_method, __LINE__, $orf->_top_tail], 
		-recurse => 1, -fh => $fh ) if $args->{'-verbose'};
	    
	    # Further adjust coords if no STOP to satisfy Genbank conventions 
	    
	    my ($stop_gbk, $path)=(undef, 'contig_end');	    
	    if ( $gap && $orf->overlap(-object => $gap) ) { 
		$path='gap';

		if ( $gap->spans( $lex ) ) { # should pop last exon in this case. not throw. 
		    $path .= '0';
		    $self->throw( $path );

		} elsif ( $lex->spans( $gap ) ) { # went _through_ a gap. probably found STOP. 
		    
		    # Criteria to this point means that _orf_ can span gap if 
		    # minimal sequence overlap (see hygiene criteria above). 
		    # Last exon could not span gap. Until iterator above.  
		    # If lex has gone clean through a gap, we need to spool back. 
		    
		    if ( $orf->composition( -aa => 'X' ) >= $GENBANK_GAP_LIMIT ) {
			$lex->stop( -adjust => -1, -R => 1) until ! $orf->overlap(-object => $gap);
			$path .= '1';
		    } else { $path .= '2'; }

		} elsif ($gap->spans( $lex->stop(-R => 1) ) ) { # stop in a gap 
		    $lex->stop( -adjust => -1, -R => 1) until ! $orf->overlap(-object => $gap);
		    $path .= '3';

		} elsif ($gap->spans( $lex->start(-R => 1) ) ) { # stop in a gap 
		    $lex->start( -adjust => +1, -R => 1) until 
			! $orf->overlap(-object => $gap) && $orf->length % $TRIPLET == 0;
		    $path .= '4';
		    
		} else { # diff exon overlaps 
		    $orf->excise_gaps();
		    $path .= '5';
		    #map { $_->output( -debug => 1, -fh => $fh,-recurse => 1) } ( $gap, $orf);
		    #$self->throw; 
		}
		
	        $stop_gbk = $lex->stop( -R => 1 );
		$lex->stop( -adjust => -1, -R => 1) until $orf->length % $TRIPLET == 0;
		$self->throw if $lex->length < 1; # can do better than this surely.. 
		
	    } elsif ( $orf->last_codon =~ $STOP_CODON ) {
	        $stop_gbk = $lex->stop( -R => 1 ); # do nothing 
		$path='stop';

	    } elsif ( $orf->_terminal_dist2( $orf->strand ) < 3 ) {
		$stop_gbk = ($orf->strand == 1 ? $self->length : 1);
		
	    } elsif ( $orf->_terminal_dist2( $orf->strand ) == 0 ) {
		$stop_gbk = $lex->stop( -R => 1 ); # do nothing 
		
	    } else { $self->throw(); }	    

	    # wrap up 

	    $lex->genbank_coords( -mode => 'stop', -R => 1, -set => $stop_gbk );
	    $orf->history( $self->_method.':'.$path );	    
	    $orf->output( 
		-prepend => [$self->_method, __LINE__, $orf->_top_tail], 
		-append => [ $path.':'.$lex->stop(-R => 1).' -> '.$stop_gbk ], 
		-fh => $fh,-recurse => 1) if $args->{'-verbose'}; 

	} # STOP 

	$orf->evaluate( -force => 1 );
    } # orf loop 
    
    map { $_->structure } $self->stream;

    return $self;
}

=head2 _genbank_quality_filter

=cut 

sub _genbank_quality_filter {
    my $self = shift;
    my $args = {@_};

    $args->{'-orf_min'} = 15 unless exists $args->{'-orf_min'};

    #################################	
    # make sure it looks like a pseudogene 
    #################################
    
    my $fh = \*STDERR;
    my $fhide = $self->up->_out_file_handle(
	-fh => 1, 
	-suffix => 'hide', 
	-tag => undef, 
	-append => 1
	);

    #################################	
    # make sure it looks like a pseudogene 
    #################################

    my %path;
    foreach my $orf ( grep { $_->coding } $self->stream ) {	
	
	my $path;
	if ( $orf->first_codon ne $START_CODON ) {
	    
	    #################################	
	    # make sure it looks like a pseudogene 
	    #################################
	    
	    if ( $orf->data('STOP') == 0 || $orf->hypergob > 5 ) {
		$path='hide';

		$orf->genbank_misc_feature( $self->_method );
		$orf->genbank_note('Could not find start codon');

	    } elsif ( $orf->ogid ) { # break OGs -- no synteny anyway. 
		$path='og-delete';
		
		$self->throw unless $args->{'-index'} && ref( $args->{'-index'} ) eq 'HASH';
		my @og = @{$args->{'-index'}->{ $orf->ogid }};
		my $og_master = shift(@og);
		$og_master->_dissolve_orthogroup( -verbose => 1 );
		delete $args->{'-index'}->{ $orf->ogid };
		
	    } else {$path='delete';}
	    
	} elsif ($orf->length <= $args->{'-orf_min'}) {

	    #################################	
	    # Hide corner cases ... GBK does not accept <=3 nt.     
	    #################################
	    
	    $path='hide-small';
	    $orf->genbank_exclude( $self->_method );
	}
	
	#################################
	# finish up 
	#################################

	# logging etc
	
	$path{ $path }++;
	$orf->history( $self->_method.':'.$path );
	$orf->output(-prepend => [$self->_method, $path, $orf->_top_tail], -fh => $fhide, -recurse => 1);
	$orf->output(-prepend => [$self->_method, $path, $orf->_top_tail], -fh => $fh, -recurse => 0)
	    if $args->{'-verbose'};

	# clean up 

	if ( $path =~ /delete/ ) {
	    $self->remove( -object => $orf );
	    $orf->DESTROY();
	}
    }

    #################################
    # print summary and clean up 
    #################################
    
    print { $fh } '>'.$self->organism.$self->id, 
    ( map { $_.":".($path{$_} || 0) } qw(delete og-delete hide hide-small) ) ;    
    map { close($_) } ( $fh, $fhide);

    #################################
    # fin
    #################################

    $self->index;
    return $self;
}


=head2 orfs(-method => 'start', -order => 'undef|up|down', -restrict => ['REAL']
    -strand => 1, -rank => 3, -noncoding => '1|0')

    Get specified genes on contig as a stream. Unlike generic ->stream 
    method allows considerable flexibility in how the set of genes are 
    selected and sorted. 

    -method : method to sort by -- start/stop/exons/evalue etc
    -order : determines the sort order (up|down)
    -strand : consider only genes on specified strand. default both. 
    -restrict : consider only ORFs with assignments listed in -restrict => []
    -rank : consider only ORFs with evidencen rank <= specificed number
    -noncoding : by default does not return genes with rank < -1. request. 
    -intergenic : return models for intergenics between selectetd features 

    $self->orfs (no args) returns all ORFs up to and including pseudos in 
    ascending start order. 

    In scalar context returns a count. 

=cut
    
sub orfs {
    my $self = shift;
    my $args = {@_};

    # defaults 
    
    $args->{'-method'} = 'start' unless exists $args->{'-method'};
    $args->{'-order'} = undef unless exists $args->{'-order'};
    # select orfs on the DNA 
    $args->{'-strand'} = undef unless exists $args->{'-strand'};
    $args->{'-start'} = undef unless exists $args->{'-start'};
    $args->{'-stop'} = undef unless exists $args->{'-stop'};
    # restrict by rank/inference/coding potential. 
    $args->{'-noncoding'} = undef unless exists $args->{'-noncoding'};
    $args->{'-intergenic'} = undef unless exists $args->{'-intergenic'};
    $args->{'-rank'} = $INFINITY unless exists $args->{'-rank'};
    # $args->{'-restrict'} = [] unless exists $args->{'-restrict'};

    # restrictions 

    map {push @{$args->{'-restrict'}}, $EVIDENCE{$_}->{'INFER'}} (keys %EVIDENCE)
	unless exists $args->{'-restrict'};
    my %hash = map {$_ => 1} @{$args->{'-restrict'}};    
    
    # checks 
    
    my $method = $args->{'-method'};
    $self->throw() unless $self->_down->can( $method );
    
    ####################################
    # apply restrictions 
    ####################################
    
    my @rank;
    foreach my $x ($self->SUPER::stream) {  # sort {$a->start <=> $b->start}
	my $rank = $x->rank;
	next unless $rank >= -1 || $args->{'-noncoding'}; # we ignore tRNAs and GAPs -- coding related only 
	next unless $rank <= $args->{'-rank'}; # we can choose to ignore more ! 	
	next unless exists $hash{$x->assign};  # or do it by assignment ? 
	push @rank, $x if ( ! $args->{'-strand'} || $x->strand == $args->{'-strand'} );
    }

    ####################################
    # impose coordinate restrictions ?
    ####################################
    
    my @r = ( 
	$args->{'-start'} && $args->{'-stop'} 
	? grep { $_->stop >= $args->{'-start'} && $_->start <= $args->{'-stop'} } @rank
	: @rank
	);
    
    ####################################
    # get intergenics? 
    ####################################

    if ( $args->{'-intergenic'} ) {
	my @sort = sort {$a->start <=> $b->start} @r;
	splice(@r);
	push @r, 
	grep {$_->length >= $args->{'-intergenic'} }  
	grep {defined} $sort[$_]->intergenic(-object => $sort[$_-1]) for (1..$#sort);
    }
    
    ####################################
    # sort
    ####################################
    
    if (! wantarray) {
	return(scalar(@r));
    } elsif ( $args->{'-order'} ) {
	if ($args->{'-order'} eq 'up' || $args->{'-order'} == +1) {
	    return(sort {$a->$method <=> $b->$method} @r);
	} else {
	    return(sort {$b->$method <=> $a->$method} @r);                              
	}
    } else { return @r; }
}

#########################################
#########################################

=head2 _pseudocontig(-start => , -stop => , -object => , -extend => )

    Create a pseudo-contigs object for localized search etc. 

    Can either supply absolute start/stop coords on the contig or can 
    supply an object and an optinal extension argument that is used to 
    adjust its start/stop coords. 

    NB: We return the pseudo-object AND the offset required to convert 
    coords to that of the calling contig. The only exception is if
    the lengths are equal between caller and pseudo AND if the calling
    context expects a scalar value (ie wantarray is FALSE).

    As with all pseudo-objects, the created contig is not added to the genome
    however, it is essential to destroy the object  post-use (as well as any
    ORFs that are not transferred to off other objects).    
    
=cut 

sub _pseudocontig {
    my $self = shift;
    my $args = {@_};

    $args->{'-extend'} = 0 unless exists $args->{'-extend'};
    $self->throw if $args->{'-extend'} != 0 && ! exists $args->{'-object'}; # extend requires object 
    $self->throw if ( exists $args->{'-start'} || exists $args->{'-stop'} ) 
	&& exists $args->{'-object'}; # object over-rides stat/stop 
    
    $args->{'-start'} = 1 unless $args->{'-start'};
    $args->{'-stop'} = $self->length unless $args->{'-stop'};
    
    if ( exists $args->{'-object'} ) {
	$self->throw unless $self->_down()->isa(ref($args->{'-object'}));
	$self->throw unless $args->{'-object'}->up eq $self;
	$args->{'-start'} = $args->{'-object'}->start - $args->{'-extend'};
	$args->{'-start'} = 1 if $args->{'-start'} <= 1;
	$args->{'-stop'} = $args->{'-object'}->stop + $args->{'-extend'};
	$args->{'-stop'} = $self->length if $args->{'-stop'} > $self->length;
    }

    $self->throw unless $args->{'-start'} =~ /^\d+$/;
    $self->throw unless $args->{'-stop'} =~ /^\d+$/;
    $self->throw if $args->{'-start'} < 1 || $args->{'-stop'} > $self->length;
    $self->throw if $args->{'-start'} > $args->{'-stop'};    

    my $pseudo = ref($self)
	->new(	
	UP => $self->up,
	SEQUENCE => substr(
	    $self->sequence, 
	    $args->{'-start'}-1, 
	    ($args->{'-stop'}-$args->{'-start'}+1) 
	),
	SCAFFOLD => -1,
	KMERCOVERAGE => -1
	);

    return ( (wantarray || $self->length != $pseudo->length ) 
	     ? ($pseudo, ($args->{'-start'}-1)) : $pseudo );
}

=head2 search(-start => , -stop => , -extend => , -safe => 1, -fast => 0,
    -db => 'undef|h3m', min => 120, -intron => [model.hmm, cutoff])
    
    In normal mode (-fast => 0) this is method runs predict()
    on a user defined sub region of a contig and selects from 
    the ORFs based on the HMM homology. It allows you to seach a 
    region for a model of interest. 
    
    The details: 
    Search a defined genomic region (-start, -stop, -extend) for 
    orfs (-min,-intron) with homology to an hmm of interest (-hmm). 
    By default (-safe) qualifying Orfs are not  attached to the 
    calling object and are returned as an array. If no HMM is provided
    we use update() to distinguish real from spurious Orfs. 
    
    The -fast option changes the run mode from performing de novo 
    ORF prediction (potentially with introns)  and susbsequently 
    selecting on HMM homology to running WISE directly on the DNA. 
    Returns fewer options but is *much* faster usually. 

=cut

sub search {
    my $self = shift;
    my $args = {@_};

    $args->{'-extend'} = 0 unless exists  $args->{'-extend'};
    $args->{'-start'} = 1 unless exists $args->{'-start'};
    $args->{'-stop'} = $self->length unless exists $args->{'-stop'}; 
    # 
    $args->{'-min'} = 120 unless exists $args->{'-min'};
    $args->{'-intron'} = [$ENV{'ANNA_HMMER2_LIB'}.'/INTRON.hmm', 2] 
	unless exists $args->{'-intron'};
    # 
    $args->{'-hmm'} = undef unless exists $args->{'-hmm'};
    # 
    $args->{'-safe'} = 1 unless exists $args->{'-safe'};
    $args->{'-fast'} = 0.1 unless exists $args->{'-fast'};
    $args->{'-top'} = undef unless exists $args->{'-top'};

    $self->throw if $args->{'-fast'} && ! $args->{'-hmm'};

    ###########################################
    # get HMM 
    ###########################################

    my ($hmmfile, $minlen);
    if ($args->{'-hmm'}) {
	($hmmfile, $minlen) = $self->_locate_hmm( $args->{'-hmm'} );
	return () unless $hmmfile;
    }
    
    ###########################################
    # create a localized region to search 
    # we pass down the -start, -stop arguments 
    ###########################################

    my $start = ($args->{'-start'} - $args->{'-extend'} < 1 ? 
		 1 : $args->{'-start'} - $args->{'-extend'});
    my $stop = ($args->{'-stop'} + $args->{'-extend'} > $self->length ? 
		$self->length : $args->{'-stop'} + $args->{'-extend'});
    my ($range, $adjustment) = 
	$self->_pseudocontig(-start => $start, -stop => $stop);
   
    ###########################################
    # use WISE to do focused accelerated search 
    ###########################################

    my @cand;
    if ( $args->{'-fast'} ) {

	# should I run twice? once with a light-weight model?
	# and then with heavier model if success? 4:21 only small saving. 
	# occasionally get very small models from wise. ignore... 	  
  
	@cand = grep { $_->length > $minlen } grep {defined}	
	$range->wise( -hmm => $hmmfile, -verbose => 0, -safe => 0 );

    } else {

	###########################################
	# place introns, make orfs 
	###########################################
	
	$range->features(-intron => $args->{'-intron'}) if $args->{'-intron'};
	
	$range->predict(
	    -orf_min => $args->{'-min'},
	    -term_min => int($args->{'-min'}/2), 
	    -exon_min => int($args->{'-min'}/2),
	    -intron_min => 30,
	    -intron_max => 1500
	    );
	
	###########################################
	# get best hit to model 
	###########################################
	
	foreach my $o ( sort {$b->length <=> $a->length} $range->stream ) {
	    if ($hmmfile) {
		next unless my ($hit) = sort { $a->{EVALUE} <=> $b->{EVALUE} }  
		$o->hmmer3(-application => 'hm3.hmmsearch', -db => $hmmfile);
		$hit->{'HIT'}=$1; # clean up ...
		$o->accept('YGOB' => $hit);
	    } else { $o->update( -exonerate => undef ); } # need enough info to do rank test
	    
	    # tandems complicate the insert of a break (-accelerate) here.
	    # if we search a large area that contains homologs, then we cannot be sure that
	    # we have recovered the _right_ homolog and cannot terminate the search. 
	    
	    $o->evaluate;
	    push @cand, $o if $o->rank < 2;
	}
    }

    ###########################################
    # transfer all good candidates 
    ###########################################

    if ( @cand ) {    
	my $lim = ($args->{'-top'} ? $args->{'-top'}-1 : $#cand);
	my @sort = sort {$a->logscore('ygob') <=> $b->logscore('ygob') } @cand;
	
	foreach my $o (@sort[0..$lim]) {
	    $range->remove(-object => $o, -warn => 0);
	    $o->adjust($adjustment);
	    $o->up($self); # must have a UP =>contig 
	    $self->add(-object => $o) if $args->{'-safe'} == 0; 
	    # $_->transfer(-from => $locus, -to => $contig, -warn => 0)
	    $o->evaluate(-validate => 1, -structure => 1);
	}
    }
    
    ###########################################
    # Clean up detritus on the _pseduocontig object 
    ###########################################

    map {  $_->DESTROY; } map { $range->remove(-object => $_, -warn => 0); } $range->stream;
    $range->DESTROY;    
    return @cand; 
}

=head2 _exons2orf( @exons )

    We make an orf object from the supplied exons.
    We do NOT validate it and we do NOT add to contig. 

=cut

sub _exons2orf {
    my $self = shift;
    my @exons = @_;
    
    foreach my $ex ( @exons ) {
	$self->throw unless $#{$ex} == 2 || $#{$ex} == 4;
	$self->throw unless $ex->[0] =~ /^\d+$/;
	$self->throw unless $ex->[1] =~ /^\d+$/;
	$ex->[2] .= '1' unless $ex->[2] =~ /1/;
	$self->throw unless abs($ex->[2]) == 1;
    }

    # 

    my $orf = Annotation::Orf
	->new(
	START => $exons[0]->[0],
	STOP => $exons[-1]->[1],
	STRAND => $exons[0]->[2]*1,
	UP => $self,
	_CREATOR => ( (caller(1))[3] )
	);
    my $fex = $orf->down;
    
    # 
    
    foreach  my $ex (@exons) {
	my $exon = Annotation::Exon
	    ->new(
	    START => $ex->[0],
	    STOP => $ex->[1],
	    STRAND => $ex->[2]*1,
	    INTRON => [ ( $#{$ex}==4 ? ($ex->[3],$ex->[4]) : (1e100, 1e100) ) ]
	    );
	$orf->add(-object => $exon);
    }
    
    $orf->remove(-object => $fex);
    $fex->DESTROY;
    $orf->index;

    return $orf;
}

=head2 features(-intron => ['INTRON.hmm', 4.5], -telomere => ['TEL.hmm', 10])

    Places featuers that are well described by HMMs but not by BLAST. 

    Models MUST be declared in the %FEATURES hash (GlobalVars.pm).

=cut

sub features {
    my $self = shift;
    my $args = {@_};	
    
    my %hash;
    foreach my $f (keys %{$args}) {
	my $model = uc($f);
	$model =~ s/^-//;
	
	next unless ref($args->{$f}) eq 'ARRAY';
	my ($hmm, $p) = @{$args->{$f}};	
	$args->{'-application'} = 'hmmer.'.($model =~ /intron/i ? 'hmmfs' : 'hmmls');
	
	$self->hmmer(
	    -application => $args->{'-application'},
	    -hmm => $hmm,     # actual hmm file
	    -model => $model, # just a name 
	    -score => $p      # log likelihood cutoff
	    );
    }

    return $self;
}

#########################################
#########################################

# override methods

=head2 remove(-log => 1|0)

    An override method to allow logging of all Orf
    removals. use -log => 0 to suppress. 

=cut 

sub remove {
    my $self = shift;
    my $args = {@_};

    $args->{'-log'} = 0; #  unless exists $args->{'-log'};

    # we log all removals for Orf models 

    unless ( $args->{'-log'} == 0 ) { # not a real part of annotation 
	$args->{'-creator'}=1;
	$args->{'-debug'}=1;
	#$args->{'-append'}=[ top and tail ];
	$args->{'-prepend'}=['REMOVE'];
	$self->up->log(%{$args}); # passed to orf->output
    }

    return $self->SUPER::remove(@_); # only accepts -object 
}

#########################################
#########################################

# getters and setters 

=head2 length()

    ....

=cut

sub length {
    my $self = shift;
    return length($self->sequence);
}

=head2 revcomp()

    ....

=cut

sub revcomp {
    my $self = shift;
    my $seq = $self->sequence;
    
	my $rc;
   	for ($x = length($seq) - 1; $x >= 0; $x--){
      	my $base = substr($seq, $x, 1);
      	if ($base =~ /a/i) {$rc .= 'T';}
        elsif ($base =~ /t/i) {$rc .= 'A';}
        elsif ($base =~ /g/i) {$rc .= 'C';}
        elsif ($base =~ /c/i) {$rc .= 'G';}
        else {$rc .= uc($base);}	
 	}
 	
	return $rc;
}


#########################################
#########################################

# make pairs from features describing NNNN, Telomeres, Centromeres etc.
# introns must be considered with coding segments so are dealt with elsewhere.
# tRNAs/snoRNAs are placed later using a different method. 

sub _link_noncoding {
    my $self = shift;
    my $args = {@_};
	
	$args->{'-feat_min'} = 30 unless defined $args->{'-feat_min'};
	$self->throw("Must specify strand: @_")
		unless abs($args->{'-strand'}) == 1;	

	# go through one-by-one
 	
	my %hash;
	foreach my $x (sort {$b->coord <=> $a->coord} $self->_get_features) {

	    next if $x->feature =~ /TERM|STOP|INTRON|RNA/;
	    next unless my $y = $x->link;
	    next unless $y->link eq $x;
	    
	    next unless $x->strand == $args->{'-strand'} ||
		($args->{'-strand'} == 1 && $x->strand == 0); # NNNNN 
	    
	    next if exists $hash{$y->_internal_id};           
	    $hash{$x->_internal_id} = 1;       
	    
	    my $orf = Annotation::Orf
		->new(
		START => $y->coord,
		STOP => $x->coord,
		STRAND => $args->{'-strand'},
		UP => $self # must point to self for morph method
		);			
	    my $exon = $orf->down;				
	    
	    # add (and add data) or destroy 
						
	    if ($exon->length >= $args->{'-feat_min'}) {	
		if ($x->feature =~ /NNN/) {
		    $orf->data('NNNN', $x->score);
		} else {$orf->data('HMM', $x->score);}
		
		my $feat = $x->feature;
		$feat =~ s/_\d//;
		$feat = substr($feat, 0, 5);
#		$orf->data('GENE', $feat);
		
		$exon->intron(-direction => 'left', -new => $INFINITY);
		$exon->intron(-direction => 'right', -new => $INFINITY);
		
		$orf->_invert_coords if $orf->strand == -1;
		$self->add(-object => $orf);	
	    } else {$orf->DESTROY;}
	}
    
    return $self;
}

# private methods. the prediction engine. 

# this method binds exon objects into full genes  
# it makes *all* possible combinations of exons 
# return ORF objects with correct (absolute) coords. 

sub _link_exons {
    my $self = shift;
    my $args = {@_};
    
    # checking that they are defined. we know they will exist. 
    
    $args->{'-intron_max'} = 1200 unless defined $args->{'-intron_max'};
    $args->{'-intron_min'} = 30 unless defined $args->{'-intron_min'};
    
    $self->throw("Must specify strand: @_")
	unless abs($args->{'-strand'}) == 1;
    
    # build up index of locations	
    # but skip genes that are single exon genes 

    my %hash;		
    foreach my $orf ($self->stream) {
	
	my $mark;
	if ($orf->strand != $args->{'-strand'}) {
	    next;
	} elsif ($orf->rank < -1) { # RNA genes
	    next;
	} elsif ($orf->down->type == 0) { # one exon 
	    $orf->_invert_coords if $orf->strand == -1;
	    next;
	} elsif ($orf->down->type == 1) {
	    $mark = $orf->down->stop(-R => 1);			
	} elsif ($orf->down->type == 2) {
	    $mark = $orf->down->start(-R => 1);			
	} elsif ($orf->down->type == 3) {
	    $mark = $orf->down->start(-R => 1);			
	}

	$hash{$mark}{'ASSIGN'} = $orf->down->type;	 
	push @{$hash{$mark}{'EXONS'}}, $orf->down; # prep exon 
	$self->remove(-object => $orf, -warn => 0); # remove orf 	
    }
    
    # stage 2: go along chromosome and build up all possible 
    # gene structure combinations 
    
    my @genes;
    foreach my $coord (sort {$a <=> $b} keys %hash) { # move along in order  	
	
	my $counter = 0;
	foreach my $x (@{$hash{$coord}{'EXONS'}}) { # foreach exon 
	    if ($hash{$coord}{'ASSIGN'} == 1) { # first exon -> new gene
		
		# make new gene 
		
		my $new = Annotation::Orf
		    ->new(
		    START => $coord,
		    STOP => $coord,
		    STRAND => $x->strand,
		    UP => $self						
		    );
		
		# strip exon and add clone of correct one
		
		$new->remove(-object => $new->down, -warn =>0);     	
		$new->add(-object => $x->clone);				
		push @genes, $new;
		
	    } elsif ($hash{$coord}{'ASSIGN'} == 2 || # -> add to all
		     $hash{$coord}{'ASSIGN'} == 3) {
		
		# foreach existing gene...
		
		my @temp;
	      M: foreach my $orf (@genes) {				
		  
		  # do not add if either gene is complete or
		  # exons overlap
		  
		  foreach my $ex ($orf->stream) {
		      next M if $ex->type == 3;				
		      next M if $x->overlap(
			  -object => $ex,
			  -fraction => 1
			  );
		  }
		  
		  # duplicate gene
		  
		  my $new = Annotation::Orf
		      ->new(
		      START => $coord,
		      STOP => $coord,
		      STRAND => $x->strand,
		      UP => $self
		      );		
		  
		  # remove dummy exon
		  
		  $new->remove(-object => $new->down, -warn =>0);
		  
		  # copy over all exons but do not add new exon 
		  
		  foreach my $exon ($orf->stream) {
		      $new->add(-object => $exon->clone);	
		  }
		  push @temp, $new;
		  
		  # add new exon to the original copy of the gene
		  # this mimics exon skipping 
		  
		  $orf->add(-object => $x->clone);  
	      }
		push @genes, @temp;
		
	    } else {$self->throw;}		
	    
	    # dont want these hanging around 
	    
	    $x->DESTROY;
	}
    }
    
    # all possible combinations should now be on array
    # go through each one and test whether we want it
    # 3 tests: gene length, intron lengths and can it be translated
    
  G: foreach my $gene (@genes) {
      $gene->_invert_coords if $gene->strand == -1;
      $gene->index;
      
      # single-exon genes will be recreated 
      
      $gene->DESTROY and next G unless $gene->exons > 1;		
      
      # length ok?
      
      # $gene->DESTROY and next unless $gene->length >= $args->{'-orf_min'};			
      
      # translatable ?
      
      $gene->DESTROY and next G unless $gene->translatable;

      # unrealistic intron ??
      
      foreach my $ex ($gene->stream) {
	  next unless my $lex = $ex->left;
	  my $len = abs($ex->start(-R => 1) - $lex->stop(-R => 1));
	  $gene->DESTROY and next G if 
	      $len > $args->{'-intron_max'} || $len < $args->{'-intron_min'};
      }
      
      $self->add(-object => $gene);
  }
    
    return $self;
}

# this method binds feature objects into Exon objects. 
# it makes *all* possible combinations of features. 

# receives contig loaded with feature objects (STOP|TERM)
# and INTRON_1|INTRON_1|ARS_1|ARS_1|TRNA_1 etc 

# all coords are relative at this point and all features think they
# are on the + strand. this makes programming much easier. 

# returns a contig loaded with orf objects, each of which has exactly
# one exon. 

sub _link_features {
    my $self = shift;
    my $args = {@_};
    	
	$args->{'-orf_min'} = 300 unless defined $args->{'-orf_min'};
	$args->{'-term_min'} = 60 unless defined $args->{'-term_min'};	
	$args->{'-exon_min'} = 60 unless defined $args->{'-exon_min'};
	$args->{'-optimise'} = 1 unless defined $args->{'-optimise'};

	# no need for heavy checking but need to know strand and 
	# objects MUST be in right order 
	
	$self->throw("Must specify strand: @_")
		unless abs($args->{'-strand'}) == 1;	
	$self->index(-method => 'coord');	


	#############
	# GET STOPS
	#############
	
	# go through one-by-one and make all possible valid structures
 	
 	my @q;
    R: foreach my $x (sort {$a->coord <=> $b->coord} $self->stream) {
		next R unless my $z = $x->left;

		# deal with $x 		           
		           
		my ($feat,$type) = (undef,undef);
		if ($x->feature =~ /STOP|TERM/) { 
			# simple case.
		} elsif ($x->feature =~ /_2/) {	
			# will be got by next STOP if INTRON_2
			# will be got by _link_noncoding otherwise
			next R; 			
		} elsif ($x->feature =~ /(\w+)_1(.*)/) {
			$feat = $1; # INTRON | ARS | TEL etc 
			$type = $2;	# species or nothign. eg SCER|AGOS INTRON
		} else {$self->throw($x->debug);} 
			
		#############
		# GET START(S) 
		#############

		# coding. x={TERM|STOP|*_1}. all _2 thrown out  		
		# collect y={INTRON_2|STOP|TERM|FEAT_2}
		# if find FEAT_1 am inside a noncoding feature. 
		
		my @r = ($z) unless $z->feature =~ /INTRON_1/;
		until ($z->feature !~ /INTRON/) { # introns transparent
			$self->throw($x->debug) unless $z = $z->left;			
			push @r, $z unless $z->feature =~ /_1/; # INTRON_1
		}
		next R if $z->feature =~ /_1/; # FEAT_1
		
		#############
		# MAKE PAIRS
		#############
		
		# derive start coord for each possible start feature 
		# will be sorted out later 
		
		L:foreach my $y (@r) {			
		
			my $adjust = 1;			
			if ($y->feature =~ /TERM|STOP/) {
				# this leads to a very nice speed up
				next if ($x->feature =~ /STOP/ && $y->feature =~ /STOP/ &&
					$x->coord - $y->coord < $args->{'-orf_min'});
			} elsif ($y->feature =~ /(\w+)_2(.*)/) {			
				if ($1 eq 'INTRON') {
					# this may change. allow 'mixed' species introns?  
					next L if ($feat eq 'INTRON' && $2 ne $type);
					$adjust = 0; # already in coding frame 
				}		 
			} else {$self->throw($y->debug);}

			#############
			# MAKE OBJECTS 
			#############
			
			my $pseudo = Annotation::Orf
				->new(
					START => $y->coord+$adjust,# moving from STOP->coding
					STOP => $x->coord,
					STRAND => $args->{'-strand'},
					UP => $self # must point to self for morph method
				);						

			my $exon = $pseudo->down;			
			$exon->intron(-direction => 'left', -new => $y->score);
			$exon->intron(-direction => 'right', -new => $x->score);
				
			# find M if 1st exon. Not using absolute coords
			# at this point so flip strand and then reflip.

			my ($top, $tail) = (); # used below. 
			unless ($y->feature =~ /TERM|NNN|INTRON/ || # intron|feature
					$args->{'-optimise'} == 0) {
				$exon->strand(1);
				$exon->morph(
					-method => 'contract',
					-terminus => '5 prime',
					-stop => 'M',
					-frame => [0,($exon->stop - $exon->start)%3]
				);					
				($top, $tail) = $exon->_top_tail;				
				$exon->strand($args->{'-strand'});				
			}

			#############
			# MAKE TESTS
			#############

			# does exon pass basic tests? morph is run at a 
			# later stage to get precise coords	if necesary 	
		
			if ($exon->type == 0) { # regular orf 

				# is it interrupted by a contig end ?
			    # does it run into an RNA gene/telomere? 
				# is optimise == 0? ie do not require M 
				# otherwise require successful M seek

				if ($x->feature =~ /TERM|NNNN/ || $y->feature =~ /TERM|NNNN/) {
					$pseudo->DESTROY and next L
						unless $exon->length >= $args->{'-term_min'};		
				} elsif ($x->feature =~ /_/ || $y->feature =~ /_/) {
					$pseudo->DESTROY and next L
						unless $exon->length >= $args->{'-orf_min'};
				} elsif ($args->{'-optimise'} == 0) {
					$pseudo->DESTROY and next L
						unless $exon->length >= $args->{'-orf_min'};
				} elsif ($top eq 'ATG') {
					$pseudo->DESTROY and next L
						unless $exon->length >= $args->{'-orf_min'};
				} else {$pseudo->DESTROY and next L;}

			} elsif ($exon->type == 1) { # first exon 
			
				# the hard work is done above -- unless () {}	
				# donor will be used in 3 frames. assume start 0 each
				# time and adjust the frame of the last nt. 
					
				# this is the only mechanism we have for preferrign
				# correct RP gene structure. finding an M - no matter
				# how close to the slpice junction - is better than a 
				# longer first exon with no M
				
				if ($top eq 'ATG') {
					$exon->intron(
						-direction => 'left',
						-new => $INFINITY
					);		
				} elsif ($y->feature =~ /TERM/) {
					# exempted from M requirement 
				} else {$pseudo->DESTROY and next L;}

			} elsif ($exon->type == 2) { # internal exon 

				$pseudo->DESTROY and next L
					unless $exon->length >= $args->{'-exon_min'};	

			} elsif ($exon->type == 3) { # terminal exon 		   
				if ($x->feature =~ /TERM/ || $y->feature =~ /TERM/) {
					my ($i,$j) = sort {$a <=> $b} # use more lenient 
						($args->{'-exon_min'}, $args->{'-term_min'});  
					$pseudo->DESTROY and next L 
						unless $exon->length >= $i; 
				} else {
					$pseudo->DESTROY and next L
						unless $exon->length >= $args->{'-exon_min'};		
				}
	
			} else {$self->throw;}

			# coords are still relative when they leave this routine
			# but going next to _link_exons. 
			
			push @q, $pseudo;
		}
	}

	return @q;
}

sub _invert_contig {
    my $self = shift;

    my %test = map { $_->_internal_id => $_->translatable } $self->stream;

    $self->sequence($self->revcomp);

    foreach my $o ($self->stream) {
	$o->_invert_coords;
	$o->strand($o->strand*-1);
	map { $_->strand($o->strand) } $o->stream;
	$o->oliver and $self->throw 
	    unless $o->translatable == $test{$o->_internal_id};
    }

    $self->index;
    return $self;
}

=head2 fuse(-object => $ctg, -orientation => +/-1)

    Create a single contig by making an end-to-end fusion of two
    contigs. We join -object to the right end of self according 
    to -orientation, adjust all feature coordinates and insert a 
    GAP object at the fusion site. -object is destroyed. 

=cut 

sub fuse {
    my $self = shift;
    my $args = {@_};

    $args->{'-orientation'} = 1 unless exists $args->{'-orientation'};
    $self->throw unless $self->isa(ref($args->{'-object'}));
    my $other = $args->{'-object'};

    ########################################
    # orient other and define new sequence 
    ########################################

    my $offset = $self->length + 100;
    $other->_invert_contig() unless $args->{'-orientation'}==1;
    $self->sequence( $self->sequence.('N' x 100).$other->sequence );

    ########################################
    # adjust all coordinates && parents 
    ########################################
    
    map { $_->throw unless $_->noncoding || $_->translatable } 
    map { $_->adjust($offset); $_ }
    map { $_->transfer(-from => $other, -to => $self, -warn => 0); $_ }
    $other->stream;
    $self->index;

    ########################################
    # insert gap object at junction 
    ########################################
    
    my $gap = Annotation::Orf
	->new(
	START => $offset-99,
	STOP => $offset,
	STRAND => 0
	);
    
    my $ex = $gap->exons(-query => 'first');
    $ex->introns(-direction => 'left', -new => $INFINITY);
    $ex->introns(-direction => 'right', -new => $INFINITY);
    
    $self->add(-object => $gap);
    $gap->data('NNNN' => $INFINITY);
    $gap->evaluate(-structure => 0, -validate => 0);
    $self->index;

    $gap->_description( join('_', 'FUSE', $self->id, $other->id*$args->{'-orientation'}) );
    
    ########################################
    # record manipulation history 
    ########################################
    
    my @append = (
	ref($other->scaffold) eq 'ARRAY' ?
	( $args->{'-orientation'} == -1 ? reverse( @{$other->scaffold} ) :  @{$other->scaffold} ) : 
	$other->scaffold
	);
    map { $_*= -1 } @append if $args->{'-orientation'};
    
    my @prepend = (ref($self->scaffold) eq 'ARRAY' ? @{$self->scaffold} : $other->scaffold);

    $self->scaffold( [ @prepend, @append] );

    ########################################    
    # clean up 
    ########################################
    
    # we have to be a little careful here as sometimes
    # we pass pseudo-objects which need - should - not be removed. 

    my $genome = $self->up;
    $genome->remove( -object => $other, -warn => 0 ) if 
	$genome->contains($other);
    $other->DESTROY;
    
    # 

    print join(',',@{$self->scaffold}) if $args->{'-verbose'};
    return (wantarray ? ($self, $gap) : $self);
}

=head2 break( -object => $gap ) 

    Break a contig object into 2 separate objects at the junction
    defined by the supplied GAP object. The objects to the right of
    GAP are placed on a new contig that is the return object.     
    The new object is automatically added to the genome object and
    all necessary modifications to coords etc are made. 

=cut 

sub break {
    my $self = shift;
    my $args = {@_};

    my $gap = $args->{'-object'};
    $self->throw unless $self->down->isa( ref($gap) );
    $self->throw unless $gap->sequence =~ /^N+$/;
    $self->throw unless $gap->up eq $self && $self->contains($gap);
    
    # we cannot split if an object spans the gap. 
    # we need to add an ORF split method. 
    # choose objects that will be on each contig.  

    map { return undef if $gap->overlap(-object => $_) } grep { $_ ne $gap } $self->stream;
    my @left = grep { $_->stop < $gap->start } $self->stream;
    my @right = grep { $_->start > $gap->stop } $self->stream;
    $self->throw unless scalar(@left,@right)+1 == scalar($self->stream);
    
    # divvy up sequences

    my @seq = split//, $self->sequence();
    my $seq1 = join('', @seq[0..($gap->start-2)]);
    my $seq2 = join('', @seq[($gap->stop)..$#seq]);
    $self->throw unless length($seq2) + $gap->length + length($seq2) == $self->length;

    # new contig. transfer all objects and verify coding. 

    my ($max) = sort { $b->id <=> $a->id } $self->up->stream;  # grr.. 
    my $new = ref($self)->new(SEQUENCE => $seq2, ID => $max->id+1);

    my $offset = -$gap->stop;
    map { $self->throw unless $_->translatable } grep { $self->rank >= 0 } 
    map { $_->adjust( $offset ) }  
    map { $_->transfer(-from => $self, -to => $new) } @right;

    # 
    
    $self->sequence( $seq2 );
    $self->remove($gap);
    $gap->DESTROY;
    map { $self->throw unless $_->translatable } grep { $self->rank >= 0 } $self->stream;

    # 

    $self->up->add( -object => $new );
    return $new;
}

=head2 fasta(-fh => STDOUT|*handle, -genbank => 0|1, 
    -organism => , -strain => )

    -genbank option triggers different fasta line decorations
    suitable for genbank submission. Uses -orgnaism, -strain. 
    
=cut

sub fasta {
    my $self = shift;
    my $args = {@_};

    my $fh = (exists $args->{'-fh'} ? $args->{'-fh'} : STDOUT);
    
    my $decoration;
    if ( $args->{'-genbank'} ) {
	$decoration = 
	    join(" ", 
		 $self->id, 
		 '[molecule=dna]', 
		 '[moltype=genomic]', 
		 '[gcode=1]', 		 
		 '[organism='.($args->{'-organism'} || 'unknown').']',
		 '[strain='.($args->{'-strain'} || 'unknown').']'
	    );
    } else {
	my @info;
	foreach my $type ( qw(GAP TRNA REAL) ) {
	    push @info, ($type eq 'REAL' ? 'CODING' : $type )."=".
		scalar(grep {$_->assign eq $type } $self->stream);	
	}
	
	$decoration = 
	    join(' ', 	    
		 $self->up->organism.'_'.
		 $self->id,
		 $self->length.'bp',
		 '['.join('; ', @info).']'
	    );
    }
    
    print {$fh} ">".$decoration;
    
    if ( $args->{'-strict'} ) {
	my @seq = split//, $self->sequence;
	
        my $line_length=60;
        until ($#seq==-1) {
            my $lim = (scalar(@seq) > $line_length ? $line_length : scalar(@seq));
            my $tseq = join('', splice(@seq, 0, $lim));
            print {$fh} $tseq;
            # print $#seq, $lim, $tseq;
        }
    } else { print {$fh} $self->sequence; }    
    
    return $self;
}

=head2 genbank_tbl( -fh => ) 

    Output information needed by tbl2asn1 for submission to genbank.

=cut 

sub genbank_tbl{
    my $self = shift;
    my $args = {@_};

    #################################
    # defaults and variables 
    #################################
    
    $self->throw unless exists $args->{'-organism'};
    $self->throw unless exists $args->{'-fh'};
    my $fh = $args->{'-fh'};

    my @bump = (3 x undef);

    #################################
    # scaffold output 
    #################################

    print {$fh} ">Feature", $self->id;
    print {$fh} 1, $self->length, 'REFERENCE';
    print {$fh} @bump, 'mol_type', 'genomic DNA';
 
    foreach my $key ( grep { defined $args->{$_} } qw(-organism -strain -taxid -reference) ) {
	my $printkey = $key;
	$printkey =~ s/^\-//;
	$printkey = 'db_xref' if $printkey eq 'taxid';
	$printkey = 'PubMed' if $printkey eq 'reference';
	print {$fh} @bump, $printkey, $args->{ $key };	
    }

    #################################
    # orf output 
    #################################

    foreach my $feat ( $self->stream ) {
	next if $feat->genbank_exclude;
	$feat->genbank_tbl( %{$args} );
    }

    return 1;
}

#########################################
#########################################

# _feature_management 

sub _add_feature {
	my $self = shift;
	my $args = {@_};

	$self->throw("Must supply a valid feature object - $args->{'-object'}")
		unless exists $args->{'-object'} && 
			$args->{'-object'} =~ /::/; 

	$args->{'-object'}->up($self);
	push @{$self->{'_FEATS'}}, $args->{'-object'};

	return $self;
}

sub _remove_feature {
	my $self = shift;
	my $args = {@_};

	$self->throw("Must supply a valid feature object - $args->{'-object'}")
		unless exists $args->{'-object'} && 
			$args->{'-object'} =~ /::/; 

	for my $i (0..$#{$self->{'_FEATS'}}) {
		next unless $self->{'_FEATS'}->[$i]->_internal_id
			== $args->{'-object'}->_internal_id;
		$self->throw unless $self->{'_FEATS'}->[$i]->coord
			== $args->{'-object'}->coord;
		splice(@{$self->{'_FEATS'}}, $i, 1);
		last;	
	}

	return $self;
}

sub _get_features {
    my $self = shift;
    return @{$self->{'_FEATS'}};
}

sub _initialise_features {
    my $self = shift;
    foreach my $f ($self->_get_features) {
	$f->up($self);
	$f->down(undef);		
	$f->left(undef);		
	$f->right(undef);
    }
    return $self;
}


#########################################
#########################################

1;

