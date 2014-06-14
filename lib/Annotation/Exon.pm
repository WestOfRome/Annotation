#!/usr/bin/perl

package Annotation::Exon;

use Annotation;
@ISA = qw(Annotation);

use Annotation::Contig;
use GlobalVars;

#########################################
#########################################

# autoloaded methods 

our $AUTOLOAD;

my %AUTOMETH = (
	start => undef,
	stop  => undef,
	strand => undef,
	intron => [0,0]
); 

sub AUTOLOAD {
    my $self = shift;    
   my $name = uc((split(/\:/, $AUTOLOAD))[-1]); # much faster 
    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

#########################################
#########################################

=head2

    New Exon object. Must specify: 

    START => coord,
    STOP => coord, 
    STRAND => ±1,
    INTRON => [score, score]

=cut

sub new {
    my $class = shift;        
    my $args = {@_};
    my $self  = $class->SUPER::new(@_);  
    
    foreach my $k ( map {uc($_)} keys %AUTOMETH) {
	if ($self->id == 0) {
	    $self->{$k} = $AUTOMETH{$k};
	} else {
	    $self->throw("$k undefined -- $args->{$k}") 
		unless exists $args->{$k} && defined $args->{$k};
	    $self->{$k} = $args->{$k};
	}
    }  

    $self->{'_CREATOR'} = (caller(1))[3];
    
    #this can occur in morph method. fix morph, then reinstate. 
    #$self->throw("Malformed exon: @_") if $self->{'START'} > $self->{'STOP'};	
    
    bless $self, $class;
}

#########################################
#########################################

=head2 morph(-method => 'expand|contract', -frame => '1|2|3', 
    -step => 3, -stop => '\*', -terminus => '3|5')

    This method is depracated! 

    Modify exon boundaries according to rule.
    e.g. expand 3' in frame 2 and steps of 3 until a '*' is found

    If you call directly you are probably doing something wrong.     

    NB: The responsibility for not expanding through a GAP object 
    NNNN or other is currently borne by calling method. ->morph 
    does not check intron scores to see if it can legally alter 
    boundary position. That is your job.

=cut 

# this method can be rewritten using _top_tail method
# probably slower but would require only pseudo, not codon obj.

sub morph {
    my $self = shift;
    my $args = {@_};
    my $length = $self->up->up->length;
    
    $args->{'-stop'} = '\*' unless exists $args->{'-stop'};
    $args->{'-method'} = 'expand' unless exists $args->{'-method'};
    $args->{'-terminus'} = 0 unless exists $args->{'-terminus'};
    $args->{'-step'} = $TRIPLET unless exists $args->{'-step'};	
    $args->{'-frame'} = undef unless exists $args->{'-frame'};
    
    # contract|shrink|decrease vs expand|increase|grow 
    
    $args->{'-firststep'} = $args->{'-step'}; # used for positioning 
    $args->{'-step'} *= -1 if $args->{'-method'} =~ /con|shr|dec/; 	
    $args->{'-terminus'} = 1 if $args->{'-terminus'} =~ /end|term|3/;
    $args->{'-terminus'} = -1 if $args->{'-terminus'} =~ /start|init|5/;

    ###############
    # SETUP
    ###############
    
    # make pseudo exon in reading frame of gene. 
    # use this to do iterations. if all goes well, update self.
    
    my $pseudo = $self->clone;
    
    # get frame of first and last nucleotides in exon 
    
    my ($fr1, $fr2) = (0,0); 
    if (defined $args->{'-frame'} && ref($args->{'-frame'}) eq 'ARRAY') {
	$fr1 = $args->{'-frame'}->[0];
	$fr2 = $args->{'-frame'}->[1];	
    } else {
	$fr1 = $self->frame(-nucleotide => 1); 
	$fr2 = $self->frame(-nucleotide => $self->length); 
    }
    
    # and calculate adjustments so we are workign in the 
    # open reading frame of the gene
    # necessary even for frame 0 (eg 1st) exons 
    # because we always mode in one codon

    my ($adj1, $adj2) = (0, 0);
    if ($fr1 == 0) {$adj1 = 3;}						
    elsif ($fr1 == 1) {$adj1 = 2;}		
    else {$adj1 = 1;}
    $pseudo->start(-adjust => $adj1, -R => 1);
    
    if ($fr2 == 0) {$adj2 = -1;}				
    elsif ($fr2 == 1) {$adj2 = -2;}		
    else {$adj2 = -3;}
    $pseudo->stop(-adjust => $adj2, -R => 1);		
    
    # annoying but necessary to deal with orfs at contig terminii
    # (we use adj - now x,y - again later)

    my ($x, $y);
    if ($self->strand == -1) {
	$x = $adj2;
	$y = $adj1;
    } else {
	$x = $adj1;
	$y = $adj2;
    }
    
    # make one codon exon to test next 'step'. a mine sweeper.  
    # can probably rewrite to just use $pseudo but no time...
    # codon should always be one step ahead of pseudo.
    
    my $codon = $pseudo->clone;
    
    ###############
    # EXON 5'PRIME 
    ###############
    
    my $five;
    
  FIVE: 
    goto THREE if $args->{'-terminus'} == 1;
    
    # place codon at start of exon ....
    
    $codon->start($pseudo->start(-R => 1), -R => 1);
    $codon->stop($codon->start);	
    
    # .... and get first codon (just inside exon)
    
    $codon->start(-adjust => -$args->{'-firststep'}, -R => 1);			
    $codon->stop(-adjust => -1, -R => 1);
    
    # get next codon until stopping condition is met
    # three stopping conditions:
    # 1 - find codon of interest
    # 2 - get to end of contig (expand only)
    # 3 - exon stop precedes start (contract only)
    
    until ($codon->sequence(-molecule => 'aa') =~ /$args->{'-stop'}/            # 1
	   || ($codon->start+$x < 1) || ($codon->stop+$y > $length)             # 2 
	   || (($pseudo->_terminal_dist(-query => 'distal', -method => 'start') # 3 
	   - $self->_terminal_dist(-method => 'stop')) < $args->{'-step'}) ) {
	#print $codon->sequence(-molecule => 'aa'), $codon->start, $codon->stop, $args->{'-stop'}, $args->{'-terminus'}; #DEVIN 

	# move codon (start and stop)
	
	$codon->start(-adjust => -$args->{'-step'}, -R => 1);
	$codon->stop(-adjust => -$args->{'-step'}, -R => 1);

	# adjust pseudo (start only)
	
	if ($args->{'-method'} =~ /con|shr|dec/) {
	    $pseudo->start(-adjust => -$args->{'-step'}, -R => 1);		
	} else {	
	    $pseudo->start(-adjust => -$args->{'-step'}, -R => 1) if	
		# this is not to do with stopping condition. ensures 
		# codon ahead of pseudo. important for first step:
		# codon starts behind but must get ahead of pseudo
		$pseudo->_terminal_dist(-method => 'start') - 
		$codon->_terminal_dist(-method => 'start') >
		abs($args->{'-step'});
	}
    }
	
    $five = $codon->sequence(-molecule => 'aa');
    
    ###############
    # EXON 3'PRIME 
    ###############
    
    my $three;
    
  THREE: 
    goto FINISH if $args->{'-terminus'} == -1;
    
    $codon->stop($pseudo->stop(-R => 1), -R => 1);
    $codon->start($codon->stop);
    $codon->start(-adjust => 1, -R => 1);
    $codon->stop(-adjust => $args->{'-firststep'}, -R => 1);

    until ($codon->sequence(-molecule => 'aa') =~ /$args->{'-stop'}/
	   || ($codon->start+$x < 1) || ($codon->stop+$y > $length)
	   || (($pseudo->_terminal_dist(-query => 'distal', -method => 'stop')
	   - $self->_terminal_dist(-method => 'start')) < $args->{'-step'}) ) {
	
	$codon->start(-adjust => $args->{'-step'}, -R => 1);
	$codon->stop(-adjust => $args->{'-step'}, -R => 1);
	
	if ($args->{'-method'} =~ /con|shr|dec/) {
	    $pseudo->stop(-adjust => $args->{'-step'}, -R => 1);
	} else {
	    $pseudo->stop(-adjust => $args->{'-step'}, -R => 1) if
		$pseudo->_terminal_dist(-method => 'stop') - 
		$codon->_terminal_dist(-method => 'stop') >
		abs($args->{'-step'});
	}
    }
    
    $three = $codon->sequence(-molecule => 'aa');			

    #################
    # wrap up 
    #################
    
  FINISH: 
    
    # restore frame of exon 
	
    $pseudo->start(-adjust => - $adj1, -R => 1);
    $pseudo->stop(-adjust => - $adj2, -R => 1);	
    
    # should account for frame above ...

    foreach my $meth (qw(start stop)) {
	$pseudo->$meth(-adjust => $args->{'-step'}, -R => 1) 
	    if $pseudo->$meth(-R => 1) < 0;
	$pseudo->$meth(-adjust => -$args->{'-step'}, -R => 1) 
	    if $pseudo->$meth(-R => 1) > $self->up->up->length;
    }

    $self->throw("Reading frame has not been honoured") unless 
	$pseudo->length%abs($args->{'-step'}) == 
	$self->length%abs($args->{'-step'});
    
    # ACCEPT CHANGES?
    
    # 1. if expand: has it reduced dist to contigs ends?
    # 2. if contract: accept if required STOP found (not sufficient
    # to just contract all the way to the contig edge). 
    # 3. start and stop have NOT 'crossed' 
    # -- checking for htis is now above in stopping conditions.  

    if ($args->{'-method'} =~ /ex|inc|gr/) { # EXPAND unless through GAP 

	unless ($pseudo->sequence =~ /N{10,}/) { # do not expand through gaps. last minute hack 
	    $self->start(-new => $pseudo->start(-R => 1), -R => 1)
		if $pseudo->_terminal_dist(-method => 'start')
		< $self->_terminal_dist(-method => 'start');	
	    
	    $self->stop(-new => $pseudo->stop(-R => 1), -R => 1)
		if $pseudo->_terminal_dist(-method => 'stop')
		< $self->_terminal_dist(-method => 'stop');
	    }
	
    } elsif ($pseudo->length > 0) { # CONTRACT unless lenght now negative 	
	
	unless ($args->{'-terminus'} == 1) { # 5' only 
	    $self->start(-R => 1, -new => $pseudo->start(-R => 1))
		if $five =~ /$args->{'-stop'}/;		
	}
	
	unless ($args->{'-terminus'} == -1) { #  3' only 		
	    $self->stop(-R => 1, -new => $pseudo->stop(-R => 1))
		if $three =~ /$args->{'-stop'}/;		
	}
    }
    
    # clean-up 

    $pseudo->DESTROY;
    $codon->DESTROY;
    return $self;
}

=head2 frame(-nucleotide => i, -distance => 1, -first => 1, -last => 1)

    Return the frame of the specified nt in the exon relative to the 
    start of the gene. If no nucleotide is specified it defaults to 1 
    (1st base) and therefore returns the frame of the exon itself.

    With  -distance returns not the frame but the distance from 
    the first nucleotide of the gene. Performing (D-1)%3 gives the frame. 

=cut

sub frame {
    my $self = shift;
    my $args = {@_};

    if ( exists $args->{'-nucleotide'} ) {
	$args->{'-nucleotide'} = $args->{'-nucleotide'};
    } elsif ( exists $args->{'-nuc'} ) {
	$args->{'-nucleotide'} = $args->{'-nuc'};
    } elsif ( exists $args->{'last'} ) {
	$args->{'-nucleotide'} = $self->length;
    #} elsif ( exists $args->{'-first'} ) {
	#$args->{'-nucleotide'} = 1;
    } else {
	$args->{'-nucleotide'} = 1;
    }

    $args->{'-distance'} = undef unless exists $args->{'-distance'};
    
    my $length = 0;	
    foreach my $ex ($self->up->stream) {
	if ($ex->_internal_id == $self->_internal_id) {
	    $length += $args->{'-nucleotide'};
	    if (defined $args->{'-distance'}) {
		return $length;
	    } else {
		last;
	    }
	} else {$length += $ex->length;}
    }
    
    #return 0 if $length==0; # DEVIN 31-5-8 : hack. check later. 
    $self->throw("Bad args: $args->{'-nucleotide'}") unless $length > 0;	
    return ($length-1)%3; # nt 1 is frame 0. same as for substr.  
}

=head2 reconstitute

    Recompose the exon by calling predict and returning minimal 
    set of open reading frames required to tile range. For a valid 
    exon - one large reading frame in (at least.. unfortunately) 
    one frame - the identical single-exon structure should be returned. 
    
    In scalar context constitutes a test. Return value of 1 means single ORF.

=cut 

sub reconstitute {
    my $self = shift;
    return $self->exonify(-minimum => 9, -overlap => 0, -frame => undef);
}

=head2 exonify(-minimum => 9, -overlap => 0, -frame => undef)

    Find all ORFs in the range specified by the exon that 
    are greater than -minimum. Returns an array of exon objs. 

    -overlap => 0 returns non-overlapping ORFs chosen by size. 
    -overlap => 1 returns ALL orfs > -min 

    NB: 
    We do not currenlty detect STOP codon creation across boundaries. 
    
=cut

sub exonify {
    my $self = shift;             
    my $args = {@_};
    
    $args->{'-minimum'} = 9 unless exists $args->{'-minimum'};
    $args->{'-overlap'} = 0 unless exists $args->{'-overlap'};
    $args->{'-frame'} = undef unless exists $args->{'-frame'}; # all frames 

    # some very light checking 
    
    #$self->throw("Below threshold $args->{'-minimum'} ".$self->length) 
    return undef unless $self->length >= $args->{'-minimum'};
    
    # transfer sequence to a pseudo-contig so we can use 'annotate'
    
    my $pseudo = Annotation::Contig
	->new(
	SEQUENCE => $self->sequence(-R => 1),
	SCAFFOLD => -$INFINITY,
	UP => $self->up->up->up # attach to genome 
	);
    
    # call the daddy method 

    $pseudo->predict(
	-strand => 1, # seq (above) automtically returns correct strand 
	-frame => $args->{'-frame'},
	-orf_min => $args->{'-minimum'},
	-term_min => $args->{'-minimum'},
	-exon_min => $args->{'-minimum'},		
	-optimise => 0 # do not require M at start 
	);

    # crude way to select ORF subset 

    $pseudo->purge(
	-overlap => $args->{'-overlap'}, # 0 = no overlaps
	-method => 'length',
	-order => 'down',
	-warn => 0,
	-force => 1
	);

    # end of Contig methods. 
    # peel off the "exon" from each called Orf (and pseudocontig)
    # and return it as array of fragments in caller coordinate space.  
    # caller can decide whether to keep or toss. 

    my @ret;
    foreach my $orf ($pseudo->stream) {	
	$self->throw("Exon not %$TRIPLET")
	    unless $orf->length%$TRIPLET == 0; 		
	$self->throw if scalar($orf->stream) > 1; 

	# terminal * 

	my $frag = $orf->down;
	my ($top, $tail) = $orf->_top_tail();	
	$frag->stop(-adjust => -3, -R => 1) if $CODONS{$tail} =~ /\*/;   # remove trailing * 
	# unless abs( $self->stop(-R => 1) - $frag->stop(-R =>1) ) < 3;

	# correct strand -- annotate runs on +1 only

	if ($self->strand == -1) {
	    $orf->_invert_coords; # must be pointing to $orf/$pseudo
	    $orf->strand(-1);
	    $frag->strand(-1);
	}	

	# correct coords + prep fopr return 
	
	$frag->start(-new => ($frag->start(-R => 1) + $self->start-1),-R => 1);		
	$frag->stop(-new => ($frag->stop(-R => 1) + $self->start-1),-R => 1);		

	$orf->remove(-object => $frag, -warn => 0); # silently remove sole exon
	push @ret, $frag;
    }

    map { $pseudo->remove(-object => $_, -log => 0, -warn => 0) } $pseudo->stream; 
    $pseudo->DESTROY;

    return map { $_->{_DEPRACATE}=1; $_ } @ret; # artificial origin is noted 
}

=head2 type()

    Exon types are arbitrary codes used internally: 
    0 = only exon
    1 = initial
    2 = internal
    3 = terminal 

=cut

sub type {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-order'} = 0 unless exists $args->{'-order'};
    $args->{'-introns'} = 1 unless exists $args->{'-introns'};
    
    if ($args->{'-order'}) {
	
    	if ($self->up->exons == 1) {
	    return 0;
    	} elsif ($self->up->down(-direction => 1) eq $self) {
	    return 1;
    	} elsif ($self->up->down(-direction => -1) eq $self) {
	    return 3;
    	} else {return 2;}
    	
    } elsif ($args->{'-introns'}) {
    
	my $l = $self->intron(-direction => 'left');	
	my $r = $self->intron(-direction => 'right');	
	
	if (($l == 0 || $l == $INFINITY) && 
	    ($r == 0 || $r == $INFINITY)) {
	    return 0;
	} elsif (($l == 0 || $l == $INFINITY) && $r > 0) {
	    return 1;	# M locked 
	} elsif (($r == 0 || $r == $INFINITY) && $l > 0) {
	    return 3;	# * locked 		
	} elsif ($l > 0 && $r > 0) {
	    return 2;
	} else {$self->throw();}
	
    } else {$self->throw();}	
}

=head2 intron(-direction => 'left|right', -object => exon_obj)

    2 functions. 

    If direction specified, the score of the left|right intron is returned.
    
    If an object is passed and it is a neighbouring object, the intron
    coords/sequence is determined and returned as an object.

=cut 

sub intron {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-direction'} = 'left' unless exists $args->{'-direction'};
    
    if (exists $args->{'-object'}) {
	$self->throw unless $self->isa(ref($args->{'-object'}));	
	my ($l,$r) = sort {$a->start <=> $b->start} ($self,$args->{'-object'});
	($l,$r) = ($r,$l) if $self->strand == -1;
	my $start = $l->stop(-R => 1);
	my $stop = $r->start(-R => 1);	
	
	my $intron = Annotation::Exon
	    ->new(
	    START => $start,
	    STOP => $start, # OK -- see below 
	    STRAND => $self->strand,
	    #INTRON => [$INFINITY,$INFINITY],
	    INTRON => [$l->intron(-direction => 'right'), $r->intron(-direction => 'left')],
	    UP => $self->up
	    );
	
	$intron->start(-R => 1, -new => $start);
	$intron->start(-R => 1, -adjust => +1);
	$intron->stop(-R => 1, -new => $stop);
	$intron->stop(-R => 1, -adjust => -1);
	
	return $intron;
    }
    
    my $val;
    if ($args->{'-direction'} eq 'left') {
	$self->{'INTRON'}->[0] = $args->{'-new'} if 
	    exists $args->{'-new'};
	$val = $self->{'INTRON'}->[0];
    } else {
	$self->{'INTRON'}->[1] = $args->{'-new'} if
	    exists $args->{'-new'};
	$val = $self->{'INTRON'}->[1];
    }
    
    return $val;	
}

#########################################
#########################################

=head2 sameframe 

    Returns true if first nt of 2 exons are in same 
    frame relative to the contig end. 

=cut 

sub sameframe {
    my $self = shift;
    my $args = {@_};
    
    $self->throw unless $self->isa( ref($args->{'-object'}) );
    $self->throw unless $self->up->up eq $args->{'-object'}->up->up;
    $self->throw unless $self->strand == $args->{'-object'}->strand;

    my $f1 = $self->_terminal_dist(-query => 'proximal', -method => 'start')%3;
    my $f2 = $args->{'-object'}->_terminal_dist(-query => 'proximal', -method => 'start')%3;
	
    return ( $f1==$f2 ? 1 : 0 );
}

=head2 _adjust(adjustment)
    
    Move exon location on sequence. 

=cut

sub _adjust {
    my $self = shift;
    my $args = {@_};
    my $adj = shift;

    $self->throw unless $adj =~ /^\-{0,1}\d+$/;
    
    $self->start( $self->start + $adj );
    $self->stop( $self->stop + $adj );

    return $self;
}

=head2 overlap(-fraction => 0)

    Return the number of overlapping bases between 2 exon annotations. 
    Or fraction if -fraction => 1.

=cut 

sub overlap {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-fraction'} = undef unless exists $args->{'-fraction'};
    
    $self->throw("$args->{'-object'} is not a ".ref($self)." object")
	unless ref($args->{'-object'}) eq ref($self);
    
    my ($x, $y) = sort {$a->start <=>
			    $b->start} ($self, $args->{'-object'});
    
    my $overlap;	
    if ($x->stop > $y->stop) {
    	$overlap = $y->length;    
    } elsif ($x->stop >= $y->start) {
    	$overlap = ($x->stop - $y->start)+1;
    } elsif ($x->stop < $y->start) {
    	$overlap = 0;
    }
    
    if ($args->{'-fraction'}) {
	my ($i,$j) = sort {$a->length <=> 
			       $b->length} ($self, $args->{'-object'});
	if ($i->length == 0) {
	    return 1;
	} else {
	    return $overlap/$i->length;
	}
    }
    
    return $overlap;
}

=head2 sequence(-molecule => 'aa', -relative => 1|0)

    Return the sequence of the object. 
    If relative specified will return the coding strand, 
    otherwise just the top strand. 
    If molecule => 'aa' will try to translate. 

=cut

sub sequence {
	my $self = shift;
    my $args = {@_};
	
	$args->{'-relative'} = 1 if 
		exists $args->{'-rel'} || exists $args->{'-R'}		
		|| $args->{'-molecule'} =~ /aa|prot|amino|pep/i;				
	$self->throw unless $self->up && $self->up->up;			
				
	my $seq = substr(
		     $self->up->up->sequence, 
		     $self->start-1, $self->length
		     );
		     
	if (exists $args->{'-relative'} && $self->strand == -1) {
	    $seq = join('', reverse(map {tr/ATGCN/TACGN/; $_;} split(//, $seq)));
	} 
	
	if ($args->{'-molecule'} =~ /aa|prot|amino/i) { 
	    my $prot;
	    for (my $x = 0; $x < length($seq); $x += $TRIPLET) {
        	my $nirnberg = substr($seq, $x, $TRIPLET);
        	if (exists $CODONS{$nirnberg}) {
		    $prot .= $CODONS{$nirnberg};
        	} else {$prot .= "X";}
	    }
	    $seq = $prot;
	}

	return $seq;	
}

=head2 length()
=cut 

sub length {
	my $self = shift;
	return(1+ $self->stop - $self->start);
}

=head2 genbank_tbl() 
=cut 

sub genbank_tbl {
    my $self = shift;
    my $args = {@_};
    
    $self->throw unless exists $args->{'-count'};
    $self->throw unless exists $args->{'-fh'};
    my $fh = $args->{'-fh'};
    
    my @bump = (3 x undef);

    print {$fh} $self->start(-R=>1),$self->stop(-R=>1),'exon';
    print {$fh} @bump,'number '.$args->{'-count'};
    
    return 1;
}

=head2 output(-fh => 'STDOUT')

    => type, start, stop, length, strand, left_score, right_score

=cut 

sub output {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-level'} = 2 unless exists $args->{'-level'};
    
    my $fh;
    if (exists $args->{'-fh'}) {
	$fh = $args->{'-fh'};
    } else {$fh = STDOUT;}
    
    my @r;
    if ($args->{'-level'} == 1) {
	push @r, "=>",$self->name,$self->start,$self->stop,$self->length,$self->strand;
    } elsif ($args->{'-level'} == 2) {			

	my @str;
	if ( ! $self->left ) {
	    @str = ($self->strand ==1 ? ('m7G-Cap', 'Donor') : ( 'Donor', 'm7G-Cap') )
	} elsif ( ! $self->right ) {
	    @str = ($self->strand ==1 ? ( 'Aceptor', 'PolyA-Tail') : ('PolyA-Tail', 'Aceptor') )
	} else {
	    @str = ($self->strand == -1 ? ('Donor', 'Aceptor') : ('Aceptor', 'Donor') )
	}

	push @r, "=>",$self->type,$self->start, $self->stop,		
	$self->strand,$self->length,
	# this ordering policy is so that the score/frame info 
	# matches the order of the coords. The user must check 
	# strand to understand the order of the exons/ boundaries
	# in the transcript but after that everything is consistent. 
	($self->strand ==1 ? $self->intron(-direction => 'left') : $self->intron(-direction => 'right') ),
	($self->strand ==1 ? $self->intron(-direction => 'right') : $self->intron(-direction => 'left') ),
	($self->strand ==1 ? $self->frame : $self->frame(-nucleotide => $self->length) ),
	($self->strand ==1 ? $self->frame(-nucleotide => $self->length) :  $self->frame),
	@str ;			
    }
    
    push @r, @{$args->{'-append'}} if $args->{'-append'};
    push @r, $self->_internal_id if $args->{'-debug'};
    unshift @r, @{$args->{'-prepend'}} if $args->{'-prepend'};
    
    # do the deed 
    
    print $fh @r;

    return $self;
}	

sub gff {
    my $self = shift;
    my $args = {@_};	

  $self->throw unless $args->{'-source'};

  my $fh = (exists $args->{'-fh'} ? $args->{'-fh'} : STDOUT);
    
    my %data = (
        ID => $self->_internal_id,
	Parent => $self->up->_internal_id,
	parent_name => $self->up->name,
	length => $self->length,
	'5prime_score' => $self->intron( -direction => 'left'),
	'3prime_score' => $self->intron( -direction => 'right'),
	'5prime_frame' => $self->frame(),
	'3prime_frame' => $self->frame(-nucleotide => $self->length),
        );

    push my @gff, 
    $self->up->id, 
    $args->{'-source'},#( map {s/\s+/_/; $_} $self->up->up->source),
    $args->{'-relationship'},
    $self->start,
    $self->stop,
    '.',
    ($self->strand > 0 ? '+' : '-' ),
    '.',
    join(';', map { $_."=".(defined $data{$_} ? $data{$_} : '.') } keys %data);
    
    print $fh @gff;
    if ( $args->{'-relationship'} eq 'exon' && $self->up->coding) {
	$gff[2] = 'CDS';
	$gff[7] = $self->frame;
	print $fh @gff;
    }
    return 1;
}

#########################################
#########################################

# override methods 

=head2 start(-relative => 1, -adjust => -8, -new => 123)

    Can use -new to specify new coord. 
    Or use -adjsut to modify existing coord. 
    -relative => 1 to operate in relative rather than absolute coors. 

=cut 

sub start {
    my $self = shift;
    my $args = {@_};	
	
    $args->{'-relative'} = 1 if exists $args->{'-rel'} || exists $args->{'-R'};	

    ####### 2014 / 6 / 11 
    #return $self->SUPER::start(@_) unless exists $args->{'-relative'};
    return $self->SUPER::start(@_) if $#_ <= 0; 
    #######

    my $new;	
    #if ($self->strand == -1) {
    if ($args->{'-relative'} && $self->strand == -1) { # added relative tag  2014 / 6 / 11 
	if (exists $args->{'-adjust'}) {
	    $new = $self->stop - $args->{'-adjust'};
	} elsif (exists $args->{'-new'}) {
	    $new = $args->{'-new'};
	} 
	$self->stop($new) if defined $new;
	return $self->stop;
    } else {
	if (exists $args->{'-adjust'}) {
	    $new = $self->start + $args->{'-adjust'};
	} elsif (exists $args->{'-new'}) {
	    $new = $args->{'-new'};
	} 
	$self->start($new) if defined $new;
	return $self->start;
    }	
}

=head2 stop(-relative => 1, -adjust => -8, -new => 123)

    Can use -new to specify new coord. 
    Or use -adjsut to modify existing coord. 
    -relative => 1 to operate in relative rather than absolute coors. 

=cut 

sub stop {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-relative'} = 1 if exists $args->{'-rel'} || exists $args->{'-R'};	

    ####### 2014 / 6 / 11 
    #return $self->SUPER::stop(@_) unless exists $args->{'-relative'};
    return $self->SUPER::stop(@_) if $#_ <= 0; 
    #######

    my $new;	
    #if ($self->strand == -1) {
    if ($args->{'-relative'} && $self->strand == -1) { # added relative tag  2014 / 6 / 11 
	if (exists $args->{'-adjust'}) {
	    $new = $self->start - $args->{'-adjust'};
	} elsif (exists $args->{'-new'}) {
	    $new = $args->{'-new'};
	} 
	$self->start($new) if defined $new;
	return $self->start;
    } else {
	if (exists $args->{'-adjust'}) {
	    $new = $self->stop + $args->{'-adjust'};
	} elsif (exists $args->{'-new'}) {
	    $new = $args->{'-new'};
	} 
	$self->stop($new) if defined $new;
	return $self->stop;
    }	
}

#########################################
#########################################

# private methods 

sub _top_tail {
    my $self = shift;
    my $dna = $self->sequence(-relative => 1);
    my $start = substr($dna, 0, 3);
    my $stop = substr($dna, length($dna)-3, 3);  
    return ($start, $stop);
}

sub _invert_coords {
    my $self = shift;
    my $start = $self->start;
    my $length = $self->up->up->length;
    
    $self->start(
	$length - 
	($self->stop - 1)
	);
    $self->stop(
	$length - 
	($start - 1)
	);
    
    return $self;
}

# use contig terminii as reference points 

=head2 _terminal_dist(-query => 'min|max|PROXIMAL|distal', 
    -method => 'START|stop')

    Returns the distance from the start|stop coord of the current
    exon to one of the two contig ends. 

    Proximal/distal refers to the relationship between the TSS and contig end. 

=cut

sub _terminal_dist {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-query'} = 'proximal' unless exists $args->{'-query'};
    $args->{'-method'} = 'start' unless exists $args->{'-method'}; 			
    my $method = $args->{'-method'};		
    
    $self->throw unless $self->up && $self->up->up;
    
    my ($min, $max, $ref);		
    if ($args->{'-query'} =~ /^m/) {
	
	($min, $max) = sort {$a <=> $b} (
	    $self->$method(-R => 1) - 1, 
	    $self->up->up->length - $self->$method(-R => 1)
	);
	
    } elsif (
	($args->{'-query'} =~ /^p/ && $self->strand == -1) ||
	($args->{'-query'} =~ /^d/ && $self->strand == 1)  ){
	
	if ($method eq 'start') {$ref = $self->up->up->length;}
	else {$ref = 1;}
	
    } elsif (
	($args->{'-query'} =~ /^p/ && $self->strand == 1)   ||
	($args->{'-query'} =~ /^d/ && $self->strand == -1)  ){
	
	if ($method eq 'start') {$ref = 1;}
	else {$ref = $self->up->up->length;}

    } else {$self->throw;}
    
    $min = abs($self->$method(-R => 1) - $ref)
	if $args->{'-query'} =~ /^[p|d]/;
    
    if ($args->{'-query'} =~ /^max/) {
	return $max;
    } else {return $min;}
}

1;
