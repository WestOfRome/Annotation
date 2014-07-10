package ExApp;

use GlobalVars;
use Annotation::Feature;
#use File::Temp qw(tempfile);
use LWP::Simple;

my %APPS;

#########################################
#########################################

=head2 blast(-query => 'file', -db => 'file', -program => 'blastp', 
    -formatdb => 0)

    Run BLAST on a seqeunce set and return an hash (keys are sequence ids) 
    of arrays, in which each element is a hit. The returned hits can be 
    processed to ORFs using _blast2orfs().

  NB: To run BLAST on a single gene and gather homology info use $orf->update().
    The blast() method is intended for largescale BLASTing and assumes additional 
    processing with one of:
    
    1. $genome->_blast2orfs() -- convert BLAST hits to NR ORF set 
    2. $genome->homology() -- calls blast and process return depending on DB type (TY etc)
    3. $orf->impute() -- does a variety of processing on returned @hits

    By default we do not format the DB though some methods (eg $genome->homology()) 
    do this prior to calling blast().

=cut 

sub blast { # avg 163ms/call
    my $self = shift;
    my $args = {@_};
	
    $args->{'-application'} = 'blast.blastall' unless exists $args->{'-application'};	
    $args->{'-program'} = 'blastp' unless exists $args->{'-program'};		
    $args->{'-qmol'} = 'aa' unless exists $args->{'-qmol'};		
    $args->{'-dbmol'} = 'aa' unless exists $args->{'-dbmol'};
    $args->{'-cpus'} = 4 unless exists $args->{'-cpus'};
    $args->{'-formatdb'} = 0 unless exists $args->{'-formatdb'};

    # make sure query OK 
	
    $args->{'-query'} = $self->_write_temp_seq(
	-seq => $args->{'-qseq'}, # an array ref 
	-molecule => $args->{'-qmol'}
	) unless exists $args->{'-query'}; # a quest fasta file 
    map {$self->throw unless -e $_ && -s $_} $args->{'-query'};

    # make sure DB file is good 

    $args->{'-db'} = $self->_write_temp_seq(
	-seq => $args->{'-dbseq'}, # an array ref 
	-molecule => $args->{'-dbmol'}
	) unless exists $args->{'-db'};	# a DB file || NCBI url 

    # do we need to format it? 
    # or is it NR?

    if ( -e $args->{'-db'}  ) { # db is a file 
	if ( $args->{'-formatdb'} == 1 ) {
	    my $formatdb = $self->_external_app(); # formatdb is default 
	    my $farg = '-p F' unless $args->{'-dbmol'} =~ /aa|prot|amino/;   
	    `$formatdb $farg -i $args->{'-db'}`;
	    $self->throw("Formatdb failed: ".`ls $args->{'-db'}`)
		unless (-e $args->{'-db'}.".nhr" || -e $args->{'-db'}.".phr");
	}
    } else {
	$args->{'-application'} = 'blast.blastcl3';
	$args->{'-params'} = ' -e 1 -I T -v 1 -b 1 -K 1 -a 10 ';
    }
    
    # run BLAST 
    
    $args->{'-params'} .= " -p $args->{'-program'} -i $args->{'-query'} ".
	"-d $args->{'-db'} -m 8 -a $args->{'-cpus'} ";# unless exists $args->{'-params'};

    my $binary = $self->_external_app(-application => $args->{'-application'});	
    print " $binary $args->{'-params'} 2> /dev/null " if $args->{'-verbose'};
    my $blr = `$binary $args->{'-params'} 2> /dev/null `;
    $self->_cleanupfile( $args->{'-query'} );
    return undef unless $blr;

    # parse   

    my $blast;
    foreach my $l (split/\n/, $blr) {
	my @r = split/\t/, $l;

	my ($start, $stop, $strand) = 
	    ($r[7] > $r[6] ? ($r[6], $r[7], 1) : ($r[7], $r[6], -1) );       		

	# returns simple HASH summary data structs
	# _process_BLAST can be called to return Orf objects instead

	push @{$blast->{$r[0]}}, {
	    HIT => $r[1],
	    EVALUE => $r[-2],
	    SCORE => $r[-1],
	    START => $start,
	    STOP => $stop,
	    STRAND => $strand
	};
    }			

    return $blast; # hash of arrays 
}

=head2 _blast2orfs(-queries => [orf objects], -blast => \%blasthash, 
    -evalue => 1e-10)
    
    Process BLAST results and return an array of ORF objects. These are
    attached to the sequences (contigs/genome) that were used as queries. 
    -homology specifies the evidence code and how 
    the reults of this BLAST will be interpreted. 
    
=cut 

sub _blast2orfs {
    my $self = shift;
    my $args = {@_};
    
    $self->throw("Required must supply array of query objects ".@_)
	unless defined $args->{'-queries'} && 
	ref($args->{'-queries'}) eq 'ARRAY';
    $self->throw("Required must supply BLAST report ".@_)		
	unless defined $args->{'-blast'} && 
	ref($args->{'-blast'}) eq 'HASH';
    $self->throw("Evalue cutoff required ".@_)		
	unless defined $args->{'-evalue'}; 
    
    # process the BLAST result using the tools at hand 
    
    my @new;
    foreach my $inter (@{$args->{'-queries'}}) {
	next unless my $blast = 
	    $args->{'-blast'}->{$inter->_internal_id};
	
	# since there may be multiple genes in the intergenic
	# will use a contig object to represent the hits 
	
	my $result = $inter->up->_pseudocontig; # clone contig 

	# use ORF objects to represnt hsps 
	# and add to the result object 
	
	foreach my $bl (@{$blast}) {		
	    next unless $bl->{'EVALUE'} <= $args->{'-evalue'};		

	    ########################################################
	    # parse to an HSP object. restore coords. split on STOP codons.  
	    ########################################################

	    my $hsp = Annotation::Orf
		->new(
		START => $bl->{'START'}+$inter->start-1, # restore contig coords 
		STOP => $bl->{'STOP'}+$inter->start-1,
		STRAND => $bl->{'STRAND'},
		UP => $result # attach to result object not contig. for now.
		);
	    
	    # use exonify to predict ORFs in this region 
	    # if HSP is good should be reconstituted unchanged. 
	    # if has STOPs will be returned as a set of exons.
	    my $hsp_range = $hsp->down; # get Exon object 
	    # call exonify to deal with STOP codons. Returns Exon array.  
	    # -frame => 0 because BLAST HSPs cannot have Frame-shifts.
	    next unless my @exons = $hsp_range->exonify(-frame => 0); # throw on error?
	    # exons are really created by BLAST so we update.
	    # used by other methods to know how to treat these exons. 
	    map { $_->_creator((caller(0))[3]) } @exons; 
	    # remove the original single range if @exons
	    map { $hsp->add(-object => $_ ) } @exons; # use hsp_range if no @exons?
	    $hsp->remove(-object => $hsp_range, -force => 0) unless ! @exons;
	    
	    ########################################################
	    # ensure the HSP object knows its SCORE 
	    ########################################################

	    # under revised implemenation only score 
	    # required for blastscore() below 

	    $hsp->accept('AA', $bl);	
	    
	    # need rich data for V1 cluster/merge methods below.
	    # not for V2. 
	    
	    foreach my $k (keys %HOMOLOGY) {
		if ( $bl->{'HIT'} =~ /^$HOMOLOGY{$k}/ ) {		    
		    $hsp->accept($k, $bl);
		    last;
		}
	    }
	    $result->add(-object => $hsp);
	}
	next unless $result->down; # nothing better than evalue 
	$result->index;

	########################################################
	# group HSPs into NR piles (same regions hits many seqs in DB). 
	########################################################
	
	# V1: gather HSPs based on hits to same gene 	
#	my $clusters = 
#	    $result->cluster(
#		-cluster => 'homology', # shared hits 
#		-param => 1,
#		-bases => 200, # have some constraint on hits 
#		-frame => 0    # do not require shared codons  
#	    );
	
	# V2: we pile up overlapping HSPs in same frame 
	# and compress to a single model. we assume, not prove homology. 
	# we essentially leave the problem of merging HSPs for later. 
	
	my $clusters = 
	    $result->cluster(
		-cluster => 'overlap', # overlapping 
		-param => $NONZERO,          # >= 80% of smaller  
		-frame => 1            # at least 1 codon in same frame 
	    );
	
	########################################################
	# choose best HSP (by score) from each cluster
	########################################################
	
	foreach my $hsp_pile ( grep { $#{$_} >0 }  values %{$clusters}) {
	    my ($main, @r) = sort {$b->score('blast') <=> $a->score('blast')} @{$hsp_pile}; # we know they are in frame 
	    $main->optimise(-feature => 'stop'); # ensure we get final STOP codon. ok? 
	    map { $result->remove(-object => $_, -log => 0, -warn => 0) } @r;
	    #$main->optimise(-feature => 'junction'); # not required under V2
	}

	########################################################	
	# purge remaining nonhomologous overlapping HSPs. expect few.  
	########################################################

	$result->purge(
	    -overlap => .5,
	    -method => 'blastscore',
	    -order => 'down',
	    -force => 1
	    );
	
	########################################################	
	# genes ready. remove from result obj & attach to contig. 
	########################################################

	foreach my $hit ($result->stream) {
	    $result->remove(-object => $hit, -force => 1);
	    $inter->up->add(-object => $hit); # add to correct contig (self=genome)
	    push @new, $hit;
	}
	$result->DESTROY;
    }
    
    return @new;
}

=head2 bl2seq(-evalue => 1e-10, -object => obj|file, -identity => 0|1 )

    Run bl2seq on a pair of sequences. 
    Returns overlap in nucleotides (ie number of bases in HSPs). 
    Or %identity (average over HSPs) with -identity switch. 

=cut 

sub bl2seq {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-evalue'} = 1e-10 unless exists $args->{'-evalue'};	
    $args->{'-application'} = 'blast.bl2seq' unless exists $args->{'-application'};	

    my $f1 = $self->_write_temp_seq();
    my $f2 = ( $self->isa(ref($args->{'-object'})) ? $args->{'-object'}->_write_temp_seq() : $args->{'-object'} );
    map {$self->throw unless -e $_ && -s $_ && -T $_} ($f1, $f2);

    # 

    $args->{'-params'} = "-p blastp -i $f1 -j $f2 -D 1" 
	unless exists $args->{'-params'}; 

    my $binary = $self->_external_app(-application => $args->{'-application'});
    my $blr = `$binary $args->{'-params'} 2> /dev/null`;
    map { $self->_cleanupfile($_) } ($f1,$f2);
    return 0 unless $blr;
    
    my $overlap = 0;	
    my @hsps;
    foreach my $l (split/\n/, $blr) {
	next if $l =~ /^\#/;
	my @r = split/\t/, $l;
	$overlap += ($r[3]*3) if $r[-2] <= $args->{'-evalue'};
	push @hsps, { LEN => $r[3], ID => $r[2] } if $r[-2] <= $args->{'-evalue'}; 
    }

    if ( $args->{'-identity'} ) {
	my ($tot, $id);
	map { $tot += $_->{LEN} } @hsps;
	map { $id += $_->{ID}*($_->{LEN}/$tot) } @hsps;
	return $id;
    }

    return $overlap ; # in nucleotides 
}

=head2 exonerate2(-object => $orf, -homolog => 'undef|species|id|file|object', 
    -model => 'global|local', -return => 'ALL|score|overlap')

    Use exonerate to align 2 proteins. They can be two objects,
    an object and a homolog or 2 objects and a homolog. In the latter 
    case we calculate an overlap score based on the fraction of the 
    homolog that is doubly-tiled by the 2 object sequences. 

    We return the alignment score and coverage for each relevant 
    attribute (self/object/compare) as ([score,cover]...[score,cover]).
    In the case of the first 2, score is the DP score and coverage
    the fraction of the homolog that is tiled. The last is the 
    overlap score (-1..1) and the fraction of the homolog covered
    by either self or object (ie the union of self and object). 
    The overlap score is positive for uniquely tiled regions and 
    negative for doubly tiled.

    -homolog : unless an argument is supplied, an homolog will
    be gathered. To do a direct comparison between 2 proteins,
    -homolog must be set to undef (or direct). 

    -model : choose between global and local alignment. Global is 
    default and required for accurate score calculation. In other 
    cases local recommended for decent speed or if gene model is 
    truncated. 

=cut 

sub exonerate2 {
    my $self = shift;
    my $args = {@_};

    # define defaults 

    $args->{'-object'} = undef unless exists $args->{'-object'};
    $args->{'-homolog'} = $args->{'-reference'} if exists $args->{'-reference'};
    $args->{'-return'} = 'all' unless exists $args->{'-return'};
    # 
    $args->{'-application'} = 'exonerate' unless exists $args->{'-application'};
    $args->{'-params'} = " --bestn 1 " unless exists  $args->{'-params'};
    $args->{'-extend'} = undef unless exists $args->{'-extend'};
    $args->{'-model'} = 'global' unless exists $args->{'-model'};
    # 
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};

    #

    $self->throw unless ! $args->{'-object'} || $self->isa( ref( $args->{'-object'} ) );

    ####################################################    
    # prep exonerate. tedious!  
    ####################################################
    
    $self->throw unless my $binary = 
	$self->_external_app(-application => $args->{'-application'});	  

    my $output = ' --showalignment 0 --showvulgar 0 --showtargetgff 0 ';
    my $output2 = ' --showalignment 1 --showvulgar 0 --showtargetgff 0 ';
    my $ryo = ' --ryo "%s %qab %qae %tab %tae\n" ';
    my $post = ' 2> /dev/null  | grep -v [a-z]  ';

    my $model = ' --model '.( 
	$args->{'-model'} =~ /g/i  
	? ' affine:global --exhaustive 1 '
	: ' affine:local --exhaustive 0 '
	);

    my $key = join('.', 'EXONERATE','PROTEIN', uc($args->{'-model'} =~ /g/i ? 'global' :'local') ); 
    # ***** $orf->score() relies on this. change with caution ******

    $args->{'-params'} .= " --gapextend $args->{'-extend'} " if $args->{'-extend'};

    ####################################################
    # Homolog special case: 
    # self + object but no homolog --> direct comparison 
    ####################################################

    if ( (exists $args->{'-homolog'} && ! defined $args->{'-homolog'}) 
	 || $args->{'-homolog'} eq 'direct' ) {
	$self->throw( $args->{'-homolog'} ) unless $self->isa( ref( $args->{'-object'} ) );
	$args->{'-homolog'} = $args->{'-object'};
	delete $args->{'-object'};
    }

    ####################################################
    # Homolog standard case: 
    # we want a file (homolog) to score the gene against
    # we can accept: object / file / geneName / species / undef 
    ####################################################

    my ($species, $homolog,$seqlen,$file,$protect);
    if ( ! exists $args->{'-homolog'} ) { # use method call 	
	if ( $args->{'-object'} ) {	
	    $homolog = $self->homolog( -object => $args->{'-object'} );		
	} else {
	    $homolog = $self->homolog( -fast => 1);		
	}
	if ( $homolog ) {
	    ($file) = $self->fetch(-id => $homolog, -molecule => 'aa', -file => 1);
	} #else {return undef;}
	return undef unless $file; # this is pretty questionable 
	
    } elsif ( exists $HOMOLOG{ $args->{'-homolog'} } ) { # use given species best hit 
       	($species, $homolog,$seqlen,$file) = 
	    $self->homolog(
		-species => [ $args->{'-homolog'} ],
		($args->{'-object'} ? (-object => $args->{'-object'}) : ())
	    );	
    } elsif ( $self->isa(ref($args->{'-homolog'})) ) { # use object 
	$file = $args->{'-homolog'}->_write_temp_seq;

    } elsif ( -e $args->{'-homolog'} ) { # use file 
	$file = $args->{'-homolog'};
	$protect=1;

    } elsif ( $args->{'-homolog'} =~ /^>/ ) { # string 
	$file = $self->_write_temp_file( -data => \$args->{'-homolog'} );

    } elsif ( $args->{'-homolog'} ) { # assume gene name and fetch from DB 
	($file) = $self->fetch(-id => $args->{'-homolog'}, -molecule => 'aa', -file => 1);

    } else { $self->throw; }

    # we have a homolog! 

    $self->oliver and 
	$self->throw("$args->{'-homolog'},$species, $homolog,$seqlen,$file") unless -e $file;

    ####################################################
    # initialize required attributes and data structures 
    ####################################################

    my @queries =  grep { defined } ($self, $args->{'-object'});
    map { $_->data($key => 0) } @queries; # initialize 
    
    # init scoring array 

    my %align;
    unless ( $args->{'-return'} eq 'score' ) {
	chomp($seqlen = `grep -v '>' $file | sed s/[[:space:]]//g | sed s/[[:space:]]+//g | wc -c `) 
	    unless $seqlen;
	$seqlen =~ s/\s+//; # required for score normalization 
	%align = map { $_->name => [ split//, (0 x $seqlen) ] } @queries;
    }

    ####################################################    
    # perform alignments 
    ####################################################

    foreach my $seq ( @queries ) {
	$seq->data( $key => 0 );
	
	my $qfile = $seq->_write_temp_seq;
	my $cmd = " $binary $args->{'-params'} $model -q $qfile -t $file $output $ryo $post ";

	print $cmd if $args->{'-verbose'};
	print ` $binary $args->{'-params'} $model -q $qfile -t $file $output2 $ryo  ` 
	    if $args->{'-verbose'} >= 2;

	my $score=0; # NB global default is -1. if run and no hit we set to 0. 
	foreach my $hsp ( split/\n/, ` $cmd ` ) {
	    print $hsp if $args->{'-verbose'};
	    $hsp =~ s/^\s+//;
	    my @r=split/\s+/, $hsp;
	    $score += ($r[0] =~ /^\d+$/ ? $r[0] : 0);
	    next if $args->{'-return'} eq 'score';
	    # for overlap... 
	    map { $align{ $seq->name }->[$_]++ } ($r[3]..$r[4]);
	    map { $align{ 'COMBINE' }->[$_]++ } ($r[3]..$r[4]);
	}
	$seq->data( $key => $score);
	$self->_cleanupfile($qfile);
    }
    $self->_cleanupfile($file) unless $protect;

    # return score if in standard use case ... 

    return $self->data( $key ) if $args->{'-return'} eq 'score';

    ####################################################
    # compute coverage stats if requested 
    ####################################################

    my %cover;
    foreach my $i ( 0..$#{$align{'COMBINE'}} ) {
	foreach my $track ( keys %align ) {
	    next unless $align{ $track }->[$i] >= 1; # ignore empty AAs 
	    if ( $track eq 'COMBINE' ) { 
		$cover{ 'SCORE' } += ( $align{ $track }->[$i] == 1 ? 1 : -1 ); # 1 or 2 
		$cover{ $align{ $track }->[$i] }++;
		$cover{ 'TOTAL' }++;
	    } else {
		$cover{ $track } += $align{ $track }->[$i];
	    }
	}
    }
    $cover{ 'TOTAL' }=1 unless $cover{ 'TOTAL' };
    # print map { "$_:$cover{$_}" } keys %cover;

    ####################################################
    # work out the return values .... score/overlap/all
    ####################################################

    # return the overlap score only.. 

    if ( $args->{'-return'} eq 'overlap') {
	return( $cover{ 'SCORE' }/$cover{ 'TOTAL' } );
    }

    my @ret = map { [$_->data($key), $cover{$_->name}/$seqlen] } @queries;
    
    # if multiple queries retun an additional overap score and fractoin tilinig statistic 
    # else just return (score, overlap)

    return ( $#queries>0 ? ( @ret, [$cover{ 'SCORE' }/$cover{ 'TOTAL' }, $cover{ 'TOTAL' }/$seqlen] ) : @{$ret[0]} ) ;
}

=head2 exonerate(-protein => self, -dna => self, -safe => 1, 
    -object => undef, -extend => 5000, -start => 1, -stop => len,
    -intron => -30)

    Use exonerate to align -protein to -dna and restrict the 
    search to the region defined by either (-start,-stop) or
    (-object,-extend). Return best (one) ORF object. -intron
    can be used to alter intron penalty. 

    NB: Returned objects _NOT_ added to contig by default. 

    The caller can be an orf or a contig provided all 
    requirements are ultimately met (dna,protein,range).
    See also wise().

    my $newgene = 
    $ctg->exonerate(-protein => 'TUB2', -start => 100, -stop => 1200);

=cut 

sub exonerate {
  my $self = shift;
  my $args = {@_};

  # 
  $args->{'-dna'} = (ref($self) =~ /Orf/i ? $self->up : $self ) unless $args->{'-dna'};
  $args->{'-object'} = (ref($self) =~ /Orf/i ? $self : undef ) unless exists $args->{'-obect'};
  $self->throw if $args->{'-object'} eq $args->{'-dna'};
  # 
  $args->{'-extend'} = 5000 unless exists $args->{'-extend'};
  $args->{'-start'} = 1 unless exists $args->{'-start'};
  $args->{'-stop'} = $args->{'-dna'}->length unless exists $args->{'-stop'};  
  # 
  $args->{'-protein'} = $args->{'-prot'} if exists $args->{'-prot'};
  $args->{'-protein'} = undef unless $args->{'-protein'};
  # 
  $args->{'-application'} = 'exonerate' unless exists $args->{'-application'};
  $args->{'-params'} = " --model protein2genome --exhaustive 0 --bestn 1 " 
      unless exists  $args->{'-params'};
  $args->{'-intron'} = -30 unless exists $args->{'-intron'};
  $self->throw unless $args->{'-intron'} <= 0 ;
  $args->{'-params'} .= " --intronpenalty $args->{'-intron'} ";
  # 
  $args->{'-safe'} = 1 unless exists $args->{'-safe'};
  $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
  
  # 

  $self->throw unless my $binary = 
      $self->_external_app(-application => $args->{'-application'});	  

  ###############################################
  # make a reduced contig to do localized search 
  ###############################################

  my ($pseudo,$adjust) = 
      $args->{'-dna'}->_pseudocontig(
	  ( $args->{'-object'} 
	    ? ( -object => $args->{'-object'}, -extend => $args->{'-extend'} ) 
	    : ( -start => $args->{'-start'}, -stop => $args->{'-stop'} ) 
	  )
      );
  $self->throw unless $pseudo;
  my $dna = $pseudo->_write_temp_seq(-molecule => 'dna');
  $self->throw unless -e $dna;
  $pseudo->DESTROY;

  ###############################################
  # proteins : this is identical to calling 'reference()' but ref() is 
  # an orf method and we amy not have 
  ###############################################
  
  my ($prot,$protect);
  if ( $args->{'-protein'} =~ /::/ ) {
      $prot = $args->{'-protein'}->_write_temp_seq(-molecule => 'aa');
  } elsif ( -e $args->{'-protein'} ) { 
      $prot = $args->{'-protein'};
      $protect=1;
  } elsif ( $args->{'-protein'} ) {
      ($prot) = $self->fetch( -id => $args->{'-protein'}, -molecule => 'aa', -file => 1 );
  } elsif ( $args->{'-object'} && $args->{'-object'}->data('_GENE') <= 1e-10 ) {
      my $best = $args->{'-object'}->data('GENE');
      ($prot) = $self->fetch( -id => $best, -molecule => 'aa', -file => 1 );
  } else { $self->throw; }
  $self->throw($prot) unless -e $prot;

  ###############################################
  # prep command 
  ###############################################

  my $outformat = "  --showalignment 0 --showcigar 0 --showvulgar 0 --showtargetgff 1 ";
  my $postprocess = " | grep -v '#' | cut -f 3-7  | grep -v exonerate | grep -v Hostname  ";
  my $cmd = "  $binary $args->{'-params'} $outformat -q $prot -t $dna $postprocess ";

  # debug

  my $outcmd = "  $binary $args->{'-params'} -q $prot -t $dna ".
      ($args->{'-verbose'} >= 2 ? " --showalignment 1 --showtargetgff 1 " : " $outformat $postprocess " );
  print "$cmd\n". ` $outcmd ` if $args->{'-verbose'}; 
 
  ###############################################
  # compose cmd, run, parse, fix coords, make object
  ###############################################

  my %hash;
  my $blockcount=0;
  foreach my $hitblock ( split/gene/, ` $cmd 2>/dev/null ` ) {
      my @lines = split/\n/, $hitblock;
      next if $#lines==0;
      $blockcount++;
      foreach my $line (@lines) {
	  $line =~ s/^\s+//;
	  my @r = split/\t/, $line;
	  if ( $r[0] eq 'exon' ) {
	      push @{$hash{$blockcount}{'EXONS'}}, [@r[1,2],$r[4]];
	  } elsif ( $r[0] eq 'similarity' ) {
	      $hash{$blockcount}{'SCR'} = $r[3];
	  } #else { $self->throw; }
      }
  }
  $self->_cleanupfile($prot) unless $protect;
  $self->_cleanupfile($dna);

  ###############################################
  # test all intervening genes in chosen species 
  ###############################################

  #my ($best) = sort { $hash{$b}{'SCR'} <=> $hash{$a}{'SCR'} } keys %hash;
  #return undef unless $best && $hash{$best}{'SCR'} > 0;
  #print ">$best", map { "$_:$hash{$_}{'SCR'}" } keys %hash if $args->{'-verbose'};

  print map { "$_:$hash{$_}{'SCR'}" } grep {defined} keys %hash if $args->{'-verbose'};
  
  my @orfs;
  foreach my $best ( sort { $hash{$b}{'SCR'} <=> $hash{$a}{'SCR'} } 
		     grep { $hash{$_}{'SCR'} > 0 } keys %hash ) {

      my @exons = sort { $a->[0] <=> $b->[0] } @{ $hash{$best}{'EXONS'} };
      my $orf = $args->{'-dna'}->_exons2orf( @exons );
      $orf->adjust( $adjust );
      #$orf->data('SCORE' => $hash{$best}{'SCR'});
      $orf->data('EXONERATE.DNA.LOCAL' => $hash{$best}{'SCR'} );
      
      $orf->output( -prepend => [ '>>>'.$orf->score ] ) if $args->{'-verbose'};

      # deal with unreported frameshifts etc 
      # we miss a 12AA extension to a YCL002W induced by a frameshift.
      # we are also vulnerable to STOP codon creation across boundaries. 
      # exonify() should be updated
      
      unless ( $orf->translatable ) {
	  foreach my $ex ( $orf->stream ) {
	      my $min = ( $ex->length < 60 ? 3 : 30);
	      if ( my @ex = grep {defined} $ex->exonify(-minimum => $min, -overlap => 0, -warn => 0) ) {
		  $orf->remove(-object => $ex, -warn => 0);
		  map { $orf->add(-object => $_) } @ex;
	      }
	  }
      }
      
      $orf->index;
      #$orf->DESTROY and return undef
      $self->oliver and $self->throw unless $orf->translatable; # Ty GAG . others? 
      
      # try snag terminal STOP 
      my ($top,$tail) = $orf->_top_tail;
      $orf->exons(-query => 'last')->stop(-R => 1, -adjust => +3) unless $CODONS{$tail} eq '*';
      
      $args->{'-dna'}->add( -object => $orf ) if $args->{'-safe'} == 0;  
      push @orfs, $orf;
  }

  return @orfs;
}

=head2 _locate_hmm 
=cut 

sub _locate_hmm {
    my $self = shift;
    my $anc = shift;

    $args->{'-db'} = $ENV{'YGOB_HMMER3_LIB'} unless $args->{'-db'};

    $self->throw("HMMER: $args->{'-db'}") 
	unless -e $args->{'-db'} && -d $args->{'-db'};

    my @sp = split/\//,$anc;
    
    my $hmm =  $args->{'-db'}.'/'.$sp[-1];
    $hmm .= '.h3m' unless $hmm =~ /\.h3m$/;
    $self->warn("No HMM file: $hmm") and return undef unless -e $hmm; # avoid horrible crashing feeling..
    return $hmm unless wantarray;

    # minlen is a cutoff below which we opt to ignore hits. 
    # we use 120nt (~1 domain) as default unless gene is short 
    # in which case we go smaller

    chomp(my @modelLen = split/\s+/, `head -5 $hmm | grep LENG`);
    my $modelLen = ( ($modelLen[1] && $modelLen[1] > 10) ? $modelLen[1] : 10)*3; # convert to nt
    my $minlen = int($modelLen < 360 ? $modelLen/3 : 120); # nucleotide lengths 
    
    return ($hmm,$minlen);
}

=head2 wise(-hmm =>, -dna =>, -object => , -extend => 5000, 
    -start => , -stop => , -model => '6:23L|21:93')

    Align -hmm to -dna using Genewise (-model) and return 
    the highest scoring positively-suppported Orf object in
    the region defined either by (-object,-extend) or 
    (-start,-stop).
    
    NB: Returned objects _NOT_ added to contig by default. 

    The caller can be an orf or a contig provided all 
    requirements are ultimately met (dna,hmm,range).
    See also exonerate().

    my $newgene = 
    $gene->exonerate(-hmm => 'Anc_2.14', -extend => 2500);

=cut

sub wise { #  avg 3.29s/call (for ~20Kb regions)
    my $self = shift;
    my $args = {@_}; 

    $args->{'-dna'} = (ref($self) =~ /Orf/i ? $self->up : $self ) unless $args->{'-dna'};
    $args->{'-object'} = (ref($self) =~ /Orf/i ? $self : undef ) unless exists $args->{'-obect'};
    $self->throw if $args->{'-object'} eq $args->{'-dna'};
    # 
    $args->{'-extend'} = 5000 unless exists $args->{'-extend'};
    $args->{'-start'} = 1 unless exists $args->{'-start'};
    $args->{'-stop'} = $args->{'-dna'}->length unless exists $args->{'-stop'};  
    # 
    $args->{'-hmm'} = undef unless exists $args->{'-hmm'};
    #$args->{'-model'} = '21:93' unless exists  $args->{'-model'};
    $args->{'-model'} = '6:23L' unless exists  $args->{'-model'};
    $args->{'-model'} =~ s/\://;
    # 
    $args->{'-application'} = 'wise2.genewise' unless exists $args->{'-application'};
    $args->{'-params'} = " -both -tabs -hmmer -alg $args->{'-model'} -gff -kbyte 10000 " 
	unless exists  $args->{'-params'}; #  -splice flat -intron tied
    # 
    $args->{'-safe'} = 1 unless exists $args->{'-safe'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    $args->{'-minlen'} = 18 unless exists  $args->{'-minlen'};

    ######################################################
    # 
    ######################################################
    
    $self->throw if $args->{'-object'} eq $args->{'-dna'};    
    $self->throw unless my $binary = 
	$self->_external_app(-application => $args->{'-application'});	  

    ######################################################    
    # make a reduced contig to do more localized search 
    ######################################################
    
    my ($pseudo,$adjust) = 
	$args->{'-dna'}->_pseudocontig(
	    ( $args->{'-object'} 
	      ? ( -object => $args->{'-object'}, -extend => $args->{'-extend'} ) 
	      : ( -start => $args->{'-start'}, -stop => $args->{'-stop'} ) 
	    )
	);
    $self->throw unless $pseudo;
    my $dna = $pseudo->_write_temp_seq(-molecule => 'dna');
    $self->throw unless -e $dna;
    $pseudo->DESTROY;

    ######################################################    
    # Get the HMM -- default to local otherwise query allele 
    ######################################################

    my $hmm = ( $args->{'-hmm'} ? $args->{'-hmm'} : $args->{'-object'}->ygob );
    $self->warn($hmm) and return undef unless $hmm;

    if ( $ENV{'YGOB_HMMER3_LIB'} ) {
	$hmm = $ENV{'YGOB_HMMER3_LIB'}.'/'.$hmm unless $hmm =~ /\//; 	
	$hmm .= '.h2m' unless  $hmm =~ /\.h2m$/; 
	$hmm =~ s/\.h3m//;
    } else {
	my $api = 'http://www.saccharomycessensustricto.org/SaccharomycesSensuStrictoResources//ygobhmms/'.$hmm;
	$api .= '.h3m' unless  $hmm =~ /\.h3m$/; 
	
	my ($fh, $file) = $self->_tempfile('XXXX');
	close($fh);
	system( " curl $api > $file 2> /dev/null " );
	my ($fh2, $file2) = $self->_tempfile('XXXX');
	close($fh2);
	system( " hm3.hmmconvert -2 $file > $file2 " );
	$hmm = $file2;
    }

    $self->warn($hmm) and return undef unless -e $hmm;
    
    ######################################################
    # compose cmd, run, parse, fix coords, make object 
    ######################################################
    
    my $cmd = " $binary $hmm $dna $args->{'-params'} ";
    print {STDERR} "$cmd\n". ` $cmd 2>/dev/null ` if $args->{'-verbose'} >= 2;

    my %hash;
    foreach my $r ( split/\n/, ` $cmd 2>/dev/null ` ) {
	$r =~ s/^\s+//;
	my @r=split/\t/, $r;
	if ( $r[2] eq 'match' ) {
	    $hash{$r[6].$r[8]}{'SCR'}=$r[5];
	} elsif ($r[1] eq 'GeneWise' ) {
	    push @{ $hash{$r[6].$r[8]}{'DATA'} }, 
	    [ (sort {$a <=> $b} @r[3..4]), $r[6] ] if $r[2] eq 'cds';
	}
    }
    $self->_cleanupfile($dna);
    return () unless %hash;

    ######################################################
    # process valid hits. must tolerate tandems etc.  
    ######################################################

    print map { "$_:$hash{$_}{'SCR'}" } grep {defined $hash{$_}{'SCR'}} keys %hash 
	if $args->{'-verbose'};
    
    my @orfs;
    foreach my $best ( sort { $hash{$b}{'SCR'} <=> $hash{$a}{'SCR'} } 
		       grep { $hash{$_}{'SCR'} > 0 } keys %hash ) {
	
	# 
	
	my @exons = sort { $a->[0] <=> $b->[0] } @{ $hash{$best}{'DATA'} };
	my $orf = $args->{'-dna'}->_exons2orf( @exons ); 
	$orf->adjust( $adjust );
	next unless $orf->length >= $args->{'-minlen'};

	#
	
	$orf->data('WISE' => $hash{$best}{'SCR'} );
	$orf->data('_WISE' => $args->{'-model'} );    
	$orf->oliver if $args->{'-verbose'};
	
	######################################################
	# deal with invalid ORFs (Fshifts etc) returned by Exonerate 
	######################################################
	
	unless ( $orf->translatable ) {
	    foreach my $ex ( $orf->stream ) {
		my $min = ( $ex->length < 60 ? 3 : 30);
		if ( my @ex = grep {defined} $ex->exonify(-minimum => $min, -overlap => 0, -warn => 0) ) {
		    $orf->remove(-object => $ex, -warn => 0);
		    map { $orf->add(-object => $_) } @ex;
		}
	    }
	}
	$orf->index;
	#$orf->DESTROY and return undef
	$self->throw unless $orf->translatable; # Ty GAG . others? 
	
	######################################################
	# try snag terminal STOP and finish up 
	######################################################
	
	my ($top,$tail) = $orf->_top_tail;
	$orf->exons(-query => 'last')->stop(-R => 1, -adjust => +3) unless $CODONS{$tail} eq '*';

	# 
	
	$args->{'-dna'}->add( -object => $orf ) if $args->{'-safe'} == 0;
	push @orfs, $orf;
    }

    return @orfs;
}

=head2 phyml(-object => [], -outgroup => '1|GeneID', -topology => 'undef|1|nwkformat', 
    -sowh => 0..N, -pvalue => .05, -speed => 1|2|3) 

    Return an hash containing nwk HKY tree and lnL.
    For SOWH test returns 2xHash and a result object. 
    Returns undef if we cannot run. 

    Defaults to orthogroup if no objects supplied. 
    Uses (clade) topology as search start if true/supplied. 
    Speed param defaults to medium (2): 
    1 => --pinv e --nclasses 4 --alpha e --search BEST --n_rand_starts 10 
    2 => --pinv e --nclasses 4 --alpha e --search NNI --n_rand_starts 1 [ ~ 20X faster ]
    3 => --search NNI --n_rand_starts 1 [ ~ 5X faster ]

    Performs the SOWH test with -sowh numbers of replicates if -sowh >= 1.
    Returns 1/0 for reject/accept the null (that free is no better than 
    supplied topolgy) with pvalue of -pvalue. 

    Goldman, N., Anderson, J. P., and Rodrigo, A. G. (2000) 
    Likelihood-based tests of topologies in phylogenetics. Systematic Biology 49: 652-670.
    http://www.ebi.ac.uk/goldman/tests/SOWHinstr.html

=cut 

sub phyml {
    my $self = shift;
    my $args = {@_}; 
    
    # generic params 
    $args->{'-application'} = 'phyml' unless exists $args->{'-application'};
    $args->{'-object'} = [ $self->orthogroup ] unless exists $args->{'-object'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};
    $args->{'-molecule'} = 'codon' unless exists $args->{'-molecule'};
    # phyml params 
    $args->{'-speed'} = 2 unless exists $args->{'-speed'};
    $args->{'-boots'} = 0 unless exists $args->{'-boots'}; 
    #$args->{'-params'} = " -d nt -b $args->{'-boots'} -m HKY85 -f m -t e --quiet " 
    $args->{'-params'} = ( $args->{'-molecule'} =~ /^[ap]/i 	
			   ? " 1 i 1 $args->{'-boots'} JTT e 4 e " 
			   : " 0 i 1 $args->{'-boots'} HKY e e 4 e " 
	) unless exists $args->{'-params'};	
    # tree params 
    $args->{'-topology'} = $self->up->up->tree unless exists $args->{'-topology'};
    $args->{'-outgroup'} = 1 unless exists $args->{'-outgroup'};
    $args->{'-min'} = 90 unless exists $args->{'-min'};
    # SOWH params 
    $args->{'-sowh'} = 0 unless exists $args->{'-sowh'};
    $args->{'-pvalue'} = .05 unless exists $args->{'-pvalue'};

    my $fherr = \*STDOUT;
    
    ######################################### 
    # deal with params 
    ######################################### 
    
    $self->throw unless $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->throw unless $self->isa(ref($_)) } @{$args->{'-object'}};
    return undef unless $#{$args->{'-object'}} >= 1;

    ######################################### 
    # get phyml set speed-based run params 
    # -- currently disabled until we can upgrade phyml 
    ######################################### 

    $self->throw unless my $binary = 
	$self->_external_app(-application => $args->{'-application'});

    if ( $args->{'-speed'} == 1 ) {
	#$args->{'-params'} .= " --pinv e --nclasses 4 --alpha e --search BEST --n_rand_starts 10 ";
    } elsif ( $args->{'-speed'} == 2 ) {
	#$args->{'-params'} .= " --pinv e --nclasses 4 --alpha e --search NNI ";	
    } else {
	#$args->{'-params'} .= " --pinv 0 --nclasses 1 --search NNI ";	
    }

    ######################################### 
    # get outgroup 
    ######################################### 

    my ($id,$dna,$aa);
    if ( $args->{'-outgroup'} ) {
	$id = ( 
	    length($args->{'-outgroup'}) > 1 
	    ? $args->{'-outgroup'}
	    : $self->outgroup(-object => $args->{'-object'}).''
	    );
	print {$fherr} "Outgroup!\n".$self->output(-string=>1) 
	    and return undef unless length($id) > 1;
	($dna,$aa) = $self->fetch(-id => $id, -relabel => 'Outg' );
	print "Outgroup: $id" if $args->{'-verbose'};
    }

    ######################################### 
    # require topology file ? 
    # add outgroup if required.
    ######################################### 

    my ( $topfile);
    if ( $args->{'-topology'} ) {	
	$args->{'-topology'} = $self->up->up->tree unless  $args->{'-topology'} =~ /^\(/;
	$self->throw unless $args->{'-topology'}  =~ /^\(/;

	# sketchy! we have trifurcating ordered tree by default. hack. 
	unless ( ! $args->{'-outgroup'} || $args->{'-topology'} =~ /^>Outg/ || ! $id ) { 
	    $args->{'-topology'} =~ s/\,/\,\(/;
	    $args->{'-topology'} =~ s/\;$/\)\;/;
	    $args->{'-topology'} =~ s/^\(/\(Outg,/;
	} 

	print $args->{'-topology'} if $args->{'-verbose'};
	$topfile = $self->_write_temp_file( -data => \$args->{'-topology'} );
    }

    ###################################################
    # get an alignment
    ###################################################
    
    my $align = $self->align(	
	-format => 'phylip', 
	-molecule => $args->{'-molecule'},
	-object => $args->{'-object'},
	($args->{'-outgroup'} ? (-rawseq => [$dna,$aa]) : () ), # fasta formatted string 
	-gaps => 0 # must be ungapped for SOWH 
	);
    print {$fherr} "Align!" and return undef unless -s $align > 100;
    my @info = split/\s+/, `head -1 $align`; # needed for SOWH 
    print {$fherr} "Gaps!" and return undef unless $info[2] >= $args->{'-min'};

    ######################################### 
    # basic tree drawing 
    ######################################### 

    unless ( $args->{'-sowh'} ) {
	my $cmdline = 
	    #"$binary ".($topfile ? " --inputtree $topfile " : "")." $args->{'-params'} -i $align &>/dev/null ";
	    "$binary $align  $args->{'-params'} ".($topfile ? " $topfile " : "BIONJ")." n y &>/dev/null ";
	print( $cmdline )  if $args->{'-verbose'};
	system( $cmdline );    
	
	my ($treefile, $statsfile) = map { $align.$_ } qw(_phyml_tree.txt _phyml_stat.txt);
	print {$fherr} "Tree/phyml 1!" and return undef unless -s $treefile > 50;
	my $ml_param = _parsePhymlStatsFile( $statsfile );
	chomp(my $ml_tree = `cat $treefile `);
	return { TREE => $ml_tree, LNL => $ml_param->{'LNL'} };
    }

    ###################################################
    # SOWH topology test 
    ###################################################
    
    # if we are testing topology fit to species tree, 
    # we run once with forced top and once w/ out, then compare to null. 

    $self->throw unless my $sqbin = 
	$self->_external_app(-application => 'seq-gen');
    my ($sqfh, $sqfile) = $self->_tempfile("annot_phyml_XXXXXX");
    
    # 

    my %res;    
    my @res;
    my $null = 'lr';

    foreach my $rm ('lr', 'tlr') { # lr is constrained to the null topology 

	###################################################
	# real data 
	###################################################

	#my $cmdline = " $binary -o $rm ".(
	#    $rm eq $null ? " --inputtree $topfile " : undef
 	#    )." $args->{'-params'} -i $align &> /dev/null ";	
	my $cmdline =  "$binary $align  $args->{'-params'} ".($topfile ? " $topfile " : "BIONJ")." ".
	    ($rm eq 'lr' ? 'n' : 'y')." y &>/dev/null ";
	print( $cmdline )  if $args->{'-verbose'};
	system( $cmdline );
	# 
	my ($treefile, $statsfile) = map { $align.$_ } qw(_phyml_tree.txt _phyml_stat.txt); #_stats.txt
	print {$fherr} "Tree/phyml 2!" and return undef unless -s $treefile > 50;
	my $ml_param = _parsePhymlStatsFile( $statsfile );

	chomp(my $ml_tree = `cat $treefile `);
	push @res, { TREE => $ml_tree , LNL => $ml_param->{'LNL'} };

	###################################################
	# seq-gen 
	# we create simulations once only 
	# but use for both lr and tlr 
	###################################################

	my $rep_count;
	if ( $rm eq $null ) {
	    my %seqgen = (
		'l' => $info[2],
		'n' => $args->{'-sowh'}*2, # we need enough to choose from 
		'm' => 'HKY',
		'g' => 4,
		'a' => $ml_param->{'ALPHA'},
		'i' => $ml_param->{'PINV'},
		't' => $ml_param->{'TSTV'},
		'f' => $ml_param->{'BASEF'}
		);
	    
	    my $sq_cmd = " $sqbin ".( join(' ', map { " -$_$seqgen{$_} " } keys %seqgen) ).
		" < $treefile 2> /dev/null " ;
	    
	    # seq-gen inserts weird special characters into files sometimes
	    # so we have to catch output, remove junk and send rest to a file
	    #my $sqstring = join('', grep {m/[\w\s_\-\.]/i} split//, ` $sq_cmd `);

	    print( $sq_cmd )  if $args->{'-verbose'};
	    my ($exp,@reps) = split/Outg/, ` $sq_cmd `;
	    my ($exp_seq,$exp_len)= grep {/\w/} split/\s+/,$exp;
	    
	    #print $exp_seq, scalar(@{$args->{'-object'}})+2,
	    #$exp_len, $info[2],
	    #scalar(@reps), $args->{'-sowh'};

	    print {$fherr} "Sims/seq-gen !" and return undef unless 
		$exp_seq == scalar(@{$args->{'-object'}})+2 && 
		$exp_len == $info[2] && 
		scalar(@reps) >= $args->{'-sowh'};

	  REP: foreach my $rep ( @reps ) {
	      $rep = 'Outg'.$rep; # restore
	      my @lines = split/\n/, $rep;
	      pop(@lines) unless scalar(@lines)==$exp_seq; # ditch trailing junk 
	      my $newrep;
	      foreach my $i (0..$#lines) {
		  my ($label,$seq)=split/\s+/,$lines[$i];
		  my $clean = join('', grep {m/[ATGCN]/i} split//,$seq);
		  next REP unless length($clean)==$exp_len;
		  $newrep.="$label\t$clean\n";
	      }
	      print {$sqfh} $exp.$newrep;
	      $rep_count++;
	      last if $rep_count==$args->{'-sowh'};
	  }	    
	    close $sqfh;
	    print {$fherr} "Sims/seq-gen: $rep_count !" and return undef 
		unless $rep_count==$args->{'-sowh'};
	    #`cp $sqfile /Users/devin/examine`;
	}

	###################################################	
	# run phyml for simulated sequences. 
	# we allow reoptimization of params (ie full SOWH).
	###################################################
	
	#my $cmdline = " $binary -o $rm ".(
	#    $rm eq $null ? " --inputtree $topfile " : undef
 	#    )." -n $args->{'-sowh'} $args->{'-params'} -i $sqfile &>/dev/null ";
	my $simparams = $args->{'-params'};
	$simparams =~ s/([01]\s+i\s+)\d+/\1$args->{'-sowh'}/ || die($args->{'-params'});
	my $cmdline =  "$binary $sqfile  $simparams ".($topfile ? " $topfile " : "BIONJ")." ".
	    ($rm eq 'lr' ? 'n' : 'y')." y &>/dev/null ";
	print( '>'.$cmdline ) if $args->{'-verbose'};
	system( $cmdline );	
	
	my ($treefile, $statsfile) = map { $sqfile.$_ } qw(_phyml_tree.txt _phyml_stat.txt); # _stats.txt
	print {$fherr} "Tree/phyml 3! :  $treefile ".(-s  $treefile) and return undef unless -s $treefile > 50;
	my @recs = _parsePhymlStatsFile( $statsfile );
	map { $res{$_}{$rm}=$recs[$_] } 0..$#recs;	
    }   
    
    $self->_cleanupfile( $topfile );
    $self->_cleanupfile( $sqfile );

    # tally results 

    my $real = ($res[1]->{LNL} - $res[0]->{LNL}); # tlr - lr == free - constrained 
    my @delta = sort {$a <=> $b}  map { ($res{$_}{'tlr'}->{LNL} - $res{$_}{'lr'}->{LNL}) } grep {!/REAL/} keys %res;
    my ($mean,$sd) = _calcMeanSD(@delta);

    my $best=0; # this is a counter for how many rands we best 
    for my $i (0..$#delta) { # these are sorted...
	$best=$i if $real > $delta[$i];
	#print $i, $best, $delta[$i], $real;
	last if $delta[$i] > $real;
    }   
    my $emp_pval = 1-(($best+1)/scalar(@delta));

    # @res = FIX(species,null), FREE(ml,alt)    

    my %ret = (
	FIX => sprintf("%.3f", $res[0]->{LNL}),
	FREE => sprintf("%.3f", $res[1]->{LNL}),
	DELTA => sprintf("%.3f", $real),
	MEAN => sprintf("%.3f", $mean),
	SD => sprintf("%.3f", $sd),
	N => scalar(@delta),
	PVAL => $emp_pval,
	RESULT => ( $emp_pval <= $args->{'-pvalue'} ? 'Alternative' : 'Species')
	);

    return \%ret;
}

sub _parsePhymlStatsFile {
    my $file = shift;
    return undef unless -e $file;
    #print $file;

    local $/ = 'Montpellier';
    open(my $fh, $file);
    #my ($junk,@chunk) = <$fh>;
    my (@chunk) = <$fh>;
    close $fh;

    my @records;
    if ( $chunk[0] =~ /Likelihood\s+Discrete\s+Number/) { # different version of phyml 
	foreach my $line ( split/\n/, $chunk[0] ) {
	    next unless $line =~ /^\s+\#\d+/;
	    my @r=split/\s+/, $line;
	    my %hash = (
		LNL => $r[3],
		ALPHA => $r[6],
		PINV => $r[7],
		TSTV => $r[8]
		);
	    push @records, \%hash; 
	}
    } else {
    foreach my $chnk ( @chunk ) {
	my %hash;
	foreach (split/\n/, $chnk) {
	    my @r=split/\s+/;
	    if ( /transversion/ ) {
		$hash{'TSTV'} = $r[-1];
	    } elsif (/likelihood/i) {
		$hash{'LNL'} = $r[-1];
	    } elsif (/Gamma\sshape\sparameter/) {
		$hash{'ALPHA'} = $r[-1];
	    } elsif (/invariant/) {
		$hash{'PINV'} = $r[-1];
	    } elsif (/f\([ATGC]\)\=/) {
		$hash{'BASEF'} .= ( $hash{'BASEF'} ? ",$r[-1]" : $r[-1] );
	    }
	}
	push @records, \%hash;
    }
    }

    return ( ! wantarray && $#records == 0 ? pop(@records) : @records );
}

sub _calcMeanSD {
    # return undef if $#_ <= 0;
    my ($mean, $sd);
    map { $mean += $_ } @_;
    $mean /= ($#_+1);
    map { $sd += ($mean-$_)**2 } @_;
    return (  sprintf("%.3f", $mean),  sprintf("%.3f", sqrt($sd/($#_+1))  ));
}

=head2 trnascan(-score => 20)

    Run trnascan on the calling object. 
    Adds both ORF objects (representing tRNAs with exons etc) and 
    Feature objects to the calling Contig object. 

    Features are used to limit gene prediction later. 

    Score specifies the COVE score cutoff.

=cut 

sub trnascan {
    my $self = shift;
    my $args = {@_};

    # this is a copy of the blast routine above.  
    # but returns feature objects rather than simple hashes
    
    $args->{'-application'} = 'tRNAscan-SE' unless exists $args->{'-application'};
    $args->{'-params'} = ' -q ' unless exists $args->{'-params'};
    $args->{'-score'} = 20 unless exists $args->{'-score'};
    
    foreach my $f ($self->_get_features) {
	$f->DESTROY if $f->feature =~ /TRNA/;
    }
    
    # make temp file and run

    my $file = $self->_write_temp_seq();
    $self->throw unless -e $file && -s $file > 20;

    my $binary = 
	$self->_external_app(-application => $args->{'-application'});	
    my $mc = `$binary $args->{'-params'} $file`;		
    $self->_cleanupfile($file);
    return undef unless $mc;
    my ($header, $hits) = split/\-+\n/, $mc;

    # 
    
    foreach my $x (split/\n/, $hits) {
	my @r = split/\s+/, $x;
	next unless $r[-1] >= $args->{'-score'};
	
	# make ORF object for tRNA
	
	my $strand = 1;
	$strand = -1 if $r[2] > $r[3];		
	my $orf = Annotation::Orf
	    ->new(
	    START => 1,
	    STOP => 3,
	    STRAND => $strand,
	    UP => $self
	    );
	
	$orf->remove(-object => $orf->down, -warn => 0, -log => 0);
	$self->add(-object => $orf);
	$orf->data('RNA', $r[-1]);
	$orf->data('GENE', $r[4].":".$r[5]);
	
	# get exons 
	
	foreach my $off (2, 6) {
	    my ($i, $j) = @r[$off, $off+1];
	    last unless $i > 0;
	    
	    # get exon ojbect in absolute coords 	
	    
	    my ($x,$y) = sort {$a <=> $b} ($i,$j);
	    
	    if (my $first = $orf->down) {
		my $z = $first->stop;
		$first->stop($x);
		$x = $y;
		$y = $z;
	    }
	    
	    my $exon = Annotation::Exon
		->new(
		START => $x,
		STOP => $y,
		STRAND => $strand,
		INTRON => [$INFINITY, $INFINITY]
		);
	    
	    $orf->add(-object => $exon);
	    
	    # get features in relative coords 	
	    
	    if ($i > $j) {
            	$i = $self->length - $i + 1;
            	$j = $self->length - $j + 1;
	    } 
			
	    my $f1 = Annotation::Feature
		->new(
		FEATURE => 'TRNA_1',
		SCORE => $INFINITY,
		COORD => $i,
		STRAND => 0
		);
	    
	    my $f2 = Annotation::Feature
		->new(
		FEATURE => 'TRNA_2',
		SCORE => $INFINITY,
		COORD => $j,
		STRAND => 0
		);
	    
	    $f2->link($f1);
	    $f1->link($f2);
	    $self->_add_feature(-object => $f1);
	    $self->_add_feature(-object => $f2);				
	}
	
	$orf->index;
    }
    
    return $self;
}


=head2 hmmer(-hmm => 'hmm_file' -model => 'what the hmm represents',
    -score => 4.5, -application => 'hmmfs|hmmls')

    Run HMMER on the calling Genome object and places models of
    DNA features such as introns. Adds feature objects to the genome. 

    Use hmmfs for introns and hmmls for telomeres. 

    For proteins use hmmer3. 

=cut 

sub hmmer {
    my $self = shift;
    my $args = {@_};
    
    # this is a copy of the blast routine above.  
    # but returns feature objects rather than simple hashes

    $args->{'-application'} = 'hmmer.hmmfs' unless exists $args->{'-application'};
    $args->{'-score'} = 4.5 unless exists $args->{'-score'};
    $args->{'-params'} = '-c ' unless exists $args->{'-params'};
    my $model = uc($args->{'-model'});

    $self->throw("Bad model file: $args->{'-hmm'}")
	unless exists $args->{'-hmm'} && -e $args->{'-hmm'};
    $self->throw("Model must be specified: $args->{'-model'}")
	unless exists $FEATURES{$model};
    
    # make temp file (no fasta method for contigs)
    # and make an hash to map ids to contigs
    
    my @contigs = (ref($self) =~ /Contig/ ? ($self) : $self->stream);
    my %hash = map {$_->_internal_id => $_} @contigs;
    my $file = $self->_write_temp_seq(-seq => [@contigs]);

    # run and parse 
    
    my $binary = $self->_external_app(-application => $args->{'-application'});	
    my $mc = `$binary $args->{'-params'} $args->{'-hmm'} $file`;

    $self->_cleanupfile($file);
    return undef unless $mc;

    my ($header, $mt1, $params, $mt2, $out) = split/(\s\-)+\n/, $mc;        

    # hmmls and hmmfs have different ouput formats
    
    my $hits;
    if ($args->{'-application'} =~ /fs/) {
	my ($subheader, $mt3, $hi) = split/(\-)+\n/, $out;
	$hits = $hi;
    } else {$hits = $out;}
    
    # pretty easy from here -- one result per line..
    
    my %ls;
    foreach my $l (split/\n/, $hits) {
	my @r = split/\s+/, $l;
	next unless $r[0] >= $args->{'-score'};
	
	# hmmls and hmmfs have different ouput formats		
	
	if ($args->{'-application'} =~ /fs/) {			
	    next unless $r[3] == 1 && # ignore incompletes ?
		$r[4] == $FEATURES{$model}->{'LENGTH'}; 
		} else { # convert hmmls -> hmmfs 
		    $l =~ /f\:\s*(\d+)\s*t\:\s*(\d+)\s*Target\:\s*(\d+)/;
		    $r[1] = $1;
		    $r[2] = $2;
		    $r[5] = $3;
		}
	
	# coords are absolute in hmmer but inverted to indicate a 
	# revese strand hit. convert to relative coords so 
	# they are ready for use in _link_features.

	#################################################
	# feature objects know their strand but, based on 
	# coords, always think they are on +1 strand
	#################################################
		
	my ($gt, $ag, $str) = ($r[1], $r[2], 1);		
	if ($r[1] > $r[2]) {
	    $str = -1;
            $gt = $hash{$r[5]}->length - $gt + 1;
            $ag = $hash{$r[5]}->length - $ag + 1;
	}
	
	# make a adjustments -- model may be longer than feature
	
	$gt += $FEATURES{$model}->{'START'};
	$ag += $FEATURES{$model}->{'STOP'};
	
	# may have multiple models for each feature 
	# eg INTRON.SCER, INTRON.EGOS 
	
	my $type;
	if ($model =~ /\./) {
	    ($model, $type) = split/\./, $model; # \//
	}
	
	# return individual feature object for model start and
	# stop -- free to use combinations later. 
	
	my $f1 = Annotation::Feature
	    ->new(
	    FEATURE => $model."_1".$type,
	    SCORE => $r[0],
	    COORD => $gt,
	    STRAND => $str
	    );			
	my $f2 = Annotation::Feature
	    ->new(
	    FEATURE => $model."_2".$type,
	    SCORE => $r[0],
	    COORD => $ag,
	    STRAND => $str
	    );	
	
	
	if ($args->{'-application'} =~ /fs/) {			
	    $f2->link($f1);
	    $f1->link($f2);		
	    $hash{$r[5]}->_add_feature(-object => $f1);
	    $hash{$r[5]}->_add_feature(-object => $f2);	
	} else {push @{$ls{$r[5]}{$str}}, ($f1, $f2);}
    }	
    
    # if using select single best hit combination 
    
    if (%ls) {
	foreach my $c ($self->stream) {
	    my $cid = $c->_internal_id;
	    foreach my $str (keys %{$ls{$cid}}) {
		my @r = sort {$a->coord <=> $b->coord} @{$ls{$cid}{$str}};
		$ls{$cid}{'LEN'}{$str} = $r[-1]->coord - $r[0]->coord;
		$ls{$cid}{$str} = \@r;			
	    }
	    my ($l,$s) = sort {$ls{$cid}{'LEN'}{$b} <=> $ls{$cid}{'LEN'}{$a}} keys %{$ls{$cid}};			
	    my $i = shift(@{$ls{$cid}{$l}});
	    my $j = pop(@{$ls{$cid}{$l}});
	    next unless ($i && $j) && (ref($i) eq ref($j));				
	    $i->link($j);
	    $j->link($i);
	    $c->_add_feature(-object => $i);
	    $c->_add_feature(-object => $j);			
	}		
    }
    
    return $self;
}

=head2 hmmer3( -db => hmmerfile -application => hmmerscan|hmmersearch|phmmer
    -evalue => 1, -drop => 1e20, -max_hits => 30) 

    Use hmmer3 to detect homology to protein DB. Caller is sually an Orf.
    Returns an array of hits similar to the blast() method. These must be 
    further processed for most applications. 

  NB: To run HMMER3 on a single gene and gather homology info use $orf->update().

    See hmmer() for placing models in DNA seqeunces. 

=cut

sub hmmer3 {
    my $self = shift;
    my $args = {@_};
    
    # hmmscan -- Search a sequence against a profile HMM database
    # hmmsearch -- Search a profile HMM against a sequence database
    # .... HMM FSA 
    # phmmer -- seach a fsa against a FSA DB 
    # .... query DB 

    $args->{'-application'} = 'hm3.hmmscan' unless exists $args->{'-application'};
    $args->{'-db'} = ($args->{'-application'} =~ /phmmer/i ? $ENV{YGOB_AA_DB} : $ENV{YGOB_HMMER3_DB}) 
	unless exists $args->{'-db'};
    #
    $args->{'-evalue'} = 1 unless exists $args->{'-evalue'};
    $args->{'-drop'} = 1e20 unless exists $args->{'-drop'};
    # 
    $args->{'-params'} = '' unless exists $args->{'-params'};
    $args->{'-bias'} = $INFINITY unless exists $args->{'-bias'};
    $args->{'-max_hits'} = 30 unless exists $args->{'-max_hits'};
    
    # prepare binary 
    
    my $binary = 
	$self->_external_app(-application => $args->{'-application'});	
    
    # prep FSA: caller can be contig or orf etc 

    my @seq = ( ref($self) =~ /Orf/i ? $self : $self->stream);
    my $fasta = 
	$self->_write_temp_seq(-remove_stop => 1, -seq => \@seq);  

    # prep hmmfile. in case of phmmer may be a second FSA 

    my @order;
    if ( $args->{'-application'} =~ /phmmer/i ) {
	$self->throw unless -e $args->{'-db'} && -f $args->{'-db'} && 
	    -T $args->{'-db'} && -s $args->{'-db'};
	@order = ( $fasta, $args->{'-db'} );
    } else {
	$self->throw unless -e $args->{'-db'} && -T $args->{'-db'} && 
	    -s $args->{'-db'};
	@order = ( $args->{'-db'}, $fasta );	
    }

    # prep outfile 
    
    my ($tblFH,$tblout) = 
	$self->_tempfile("annot_hmmer3_XXXXXX"); 
    close( $tblFH );
    
    # 

    print("$binary --tblout $tblout $args->{'-params'} @order &> /dev/null")
	if $args->{'-verbose'};
    system("$binary --tblout $tblout $args->{'-params'} @order &> /dev/null");
    $self->_cleanupfile($fasta);

    my $string =  `cat $tblout`;
    print $string and $self->throw if $string =~ /Error\:/;
    
    # 

    my $oldE;
    my @hits;
    foreach my $line (split/\n/, $string) {
	next if $line =~ /^\#/;
	my ($target, $ac1, $query, $ac2, $evalue, $score, $bias, @others)=split/\s+/, $line;

	if (  $args->{'-verbose'} ) {
	    print $line;
	    print $target, $query, $evalue, $score, $bias;
	}

	last if $evalue > $args->{'-evalue'}; # max eval to record 
	last if ($oldE && (( $evalue / ($oldE > 0 ? $oldE : $INFINITY )) > $args->{'-drop'})); # min drop 
	$oldE = $evalue;
	next if $bias > $args->{'-bias'};

	push @hits, {
	    QUERY => $query,
	    HIT => $target,  
	    EVALUE => ($evalue eq '0' ? (0.0)*1 : $evalue), # force 0 to number  
	    SCORE => $score,
	    # lies ... want to fake BLAST-like return hash 
	    START => 1,
	    STOP => 1000, # ugh 
	    STRAND => 1
	};
 
	last if $#hits >= $args->{'-max_hits'}; # max hits to keep 
    }
    $self->_cleanupfile($tblout);

    return @hits;
}

=head2 hcnf(-evalue => 1e-10, -minimum => 3, -coverage => .7)

    Run BLAST on a genome and identify High Copy Number Families
    using single-linkage clustering. Uses only sequences classified
    as 'HYPO', 'NOVEL', 'INTER', 'REPEAT'. 

    Core specifies the minimum number of sequences in cluster to qualify as HCNF. 

=cut

sub hcnf {
	my $self = shift;
	my $args = {@_};	

	$args->{'-evalue'} = 1e-10 unless exists $args->{'-evalue'};		
	$args->{'-coverage'} = 0.7 unless exists $args->{'-coverage'}; 
	$args->{'-minimum'} = 3 unless exists $args->{'-minimum'};

	$self->throw("-evalue set too low $args->{'-evalue'}")
	    unless $args->{'-evalue'} < 1;
	
	# prep a databse and make an index 

	my @r = $self->orfs(
	    -restrict => ['HYPO', 'NOVEL', 'INTER'], # 'REPEAT']
	    );
	return $self if $#r == -1;	
	my $db = $self->_write_temp_seq(-seq => \@r);
	my %index = map { $_->_internal_id => $_ } @r;
	
	# blast 
	
	my $blast = $self->blast(
	    -program => 'blastp',
	    -query => $db,
	    -db => $db
	    );  

	# make a hit matrix 

	my %data;	
	foreach my $query (keys %{$blast}) {
	  H:foreach my $hsp (@{$blast->{$query}}) {
	      next H if $query == $hsp->{HIT}; # self hit. 
	      next H unless $hsp->{EVALUE} <= $args->{'-evalue'};
	      next H unless $hsp->{STRAND} == 1;
	      next H unless (($hsp->{STOP}-$hsp->{START}+1)*3)/($index{$query}->length) 
		  > $args->{'-coverage'};
	      next H if exists $data{$query}{$hsp->{HIT}}; # top HSP foreach gene only 
	      #print $query, $hsp->{HIT}, (($hsp->{STOP}-$hsp->{START}+1)*3), ($index{$query}->length);
	      $data{$query}{$hsp->{HIT}} = $hsp->{EVALUE};		
	  }
	}

	# single linkage clustering

	my %cluster;
	foreach my $Q (keys %data) { # _internal_id 
	    $cluster{$Q} = [$index{$Q}] unless exists $cluster{$Q};

	    foreach my $H (keys %{$data{$Q}}) {
		$cluster{$H} = [$index{$H}] unless exists $cluster{$H};
		next if $cluster{$Q} eq $cluster{$H};

		my ($hi, $lo) = sort { $#{$cluster{$b}} <=> $#{$cluster{$a}} } ($Q, $H);

		foreach my $obj (grep {/\w/} @{$cluster{$lo}}) {
		    push @{$cluster{$hi}}, $obj;
		    $cluster{$obj->_internal_id} = $cluster{$hi};
		}
	    }
	}

	# throw out redundant clusters we do not want and add genes to families 
	
	my %seen;
	my @unique = grep { @{$_} >= $args->{'-minimum'} } grep { !$seen{$_}++ } values %cluster;

	# do we see the expected pattern of a high density of 
	# mutual hits ? 
	

	# 

	foreach my $i (0..$#unique) {
	    my $cid = 'T'.lc($self->organism).($i+1);
	    print $cid, scalar(@{$unique[$i]}) if $args->{'-verbose'};	    
	    foreach my $orf (@{$unique[$i]}) {
		$orf->data('HCNF' => );
		$orf->data('_HCNF' => $cid);
	    }
	}

	return $self;
}

=head2 fetch(-id => $id, -molecule => 'both|aa|nt', -file => 1
    -db => databaseFile, -relabel => undef)

    Return a fasta sequence(s) from a DNA database file. 
    By default we return [$nt,$aa] but the -molecule argument 
    can be used to select only one and the -file argument can 
    be used to return files instead. 

    -relabel sets the id to a defined value such as 'Outg'.

=cut 

sub fetch {
    my $self = shift;
    my $args = {@_};

    $args->{'-db'} = $ENV{'YGOB_DNA_DB'} unless exists $args->{'-db'};
    $args->{'-relabel'} = undef unless exists $args->{'-relabel'};
    $args->{'-file'} = undef unless exists $args->{'-file'};
    $args->{'-molecule'} = 'both' unless exists $args->{'-molecule'};

    $self->throw unless $args->{'-id'} && $args->{'-db'} =~ /\w/;
    $self->throw(  $args->{'-db'} ) unless -e  $args->{'-db'};
    my $caller = (caller(1))[3];

    my ($id,$nt,@r,$aa);
    if ( $args->{'-id'} =~ /^\d+$/ ) {
	#return undef;
	$id = $args->{'-id'};
	$aa = $self->_fetch_ncbi($args->{'-id'});
    } else {
	$self->throw($args->{'-id'}) unless $args->{'-id'} =~ /\d/;
	my $app = $self->_external_app(-application => "blast.fastacmd" ) || die;

	# chomp(my $seq = `grep -A 1 $id $db`);
	chomp(my $seq = `$app -d $args->{'-db'} -s $args->{'-id'}`);    

	$self->warn(" $caller : $args->{'-db'} / $args->{'-id'} : $seq ") and return undef unless 
	    $seq && $seq =~ /^\>([^\n]+)\n[GPAVLIMCFYWHKRQNEDSTX]*/;
	($id, @r) = split/\n/, $seq;
	$nt = join('', @r);
	my @codons = map { substr($nt, $_, 3) } ( grep { $_%3==0 } (0..length($nt))  );
	$aa = join('', map { $CODONS{$_} } @codons);	
    }

    # 
    
    my $idline = ($args->{'-relabel'} ? ">".$args->{'-relabel'} : $id)."\n";
    my $aa_rec = $idline.$aa;
    my $nt_rec = $idline.$nt;

    # 

    my @ok = ( $args->{'-molecule'} =~ /^b/i 
	       ? ($nt_rec, $aa_rec) 
	       : ( $args->{'-molecule'} =~ /^[ap]/i ? $aa_rec : $nt_rec )
	);

    # 

    my @ret = ( $args->{'-file'} ? ( map { $self->_write_temp_file( -data => \$_ ) } @ok ) : @ok );

    return (wantarray || $#ret > 0 ? @ret : $ret[0] );
}

sub _fetch_ncbi {
    my $self = shift;
    my $id = shift;

    my $content = get("$eutils$id"); # GBSeq_sequence
    $content =~ /GBSeq_sequence\>([^\<]+)/;
    my $protein = uc($1);
    return $protein;
}


=head2 align(-application => 'fsa|clustalw', -molecule => 'dna|aa', -format => 'hash|phylip|clustal|fasta|etc')

    Default is protein assisted DNA align. DNA or AA only align by request. 

=cut 

sub align {
    my $self = shift;
    my $args = {@_};

    $args->{'-application'} = 'fsa' unless exists $args->{'-application'};
    $args->{'-molecule'} = 'codon' unless exists $args->{'-molecule'};
    $args->{'-format'} = 'clustal' unless exists  $args->{'-format'};
    $args->{'-gaps'} = 1 unless exists  $args->{'-gaps'};
    $args->{'-object'} = [ $self->orthogroup ] unless exists $args->{'-object'};
    $args->{'-rawseq'} = undef unless exists $args->{'-rawseq'};

    my $range = "-extract=".join('..', @{$args->{'-range'}}) if $args->{'-range'};

    my $readseqfmt = ( $args->{'-format'} eq 'hash' ? 'fasta' : $args->{'-format'} ); 

    $self->throw unless $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->throw unless $self->isa(ref($_)) } @{$args->{'-object'}};
    $self->throw if $args->{'-gaps'}==0 && $args->{'-application'} ne 'fsa';

    my %uniq;
    my @obj = grep { !$uniq{$_->unique_id}++} ($self, @{$args->{'-object'}});
    my $aa = $self->_write_temp_seq(-molecule => 'aa', -seq => \@obj, -id => 'organism');
    my $dna = $self->_write_temp_seq(-molecule => 'dna', -seq => \@obj, -id => 'organism');

    # this shit is ugly.. but needed to sneak in an outgroup 

    if ( $args->{'-rawseq'} ) {
	open(AA, ">>$aa");
	print AA $args->{'-rawseq'}->[1];
	close(AA);
	open(DNA, ">>$dna");
	print DNA $args->{'-rawseq'}->[0];
	close(DNA);
    }

    my ($TMP_HANDLE1, $dna_aln) = 
	$self->_tempfile("annot_align_dnaXXXXXX");
    my ($TMP_HANDLE2, $aa_aln) = 
	$self->_tempfile("annot_align_alnXXXXXX");
    my ($TMP_HANDLE3, $fmt_aln) = 
	$self->_tempfile("annot_align_fmtXXXXXX");

    # we need all these things for clustalw 

    my $clustalFile = ( $args->{'-molecule'} eq 'dna' ? $dna : $aa );
    $clustalFile =~ s/\.\w+$//;
    $clustalFile .= '.aln';
    my $prefix = ( $args->{'-application'} eq 'fsa' ? undef : ' -INFILE=' );

    #  get executables 

    my $fsa = $self->_external_app(-application => $args->{'-application'} );
    my $revtrans = $self->_external_app(-application => 'revtrans');

    # run aligners 
    
    my $align;
    if ( $args->{'-molecule'} eq 'aa' ) {
	#print `cat  $prefix$aa`;
	system("$fsa $prefix$aa > $aa_aln 2> /dev/null");	
	$align = ( $args->{'-application'} eq 'fsa' ? $aa_aln : $clustalFile );
    } elsif ( $args->{'-molecule'} eq 'dna' ) {
	system("$fsa $prefix$dna > $dna_aln 2> /dev/null");	
	$align = ( $args->{'-application'} eq 'fsa' ? $dna_aln : $clustalFile );
    } else {
	system("$fsa $prefix$aa  > $aa_aln 2> /dev/null ");
	$align = ( $args->{'-application'} eq 'fsa' ? $aa_aln : $clustalFile );
	$self->warn( $prefix.$aa ) unless -e $align && -s $align;
	system("$revtrans $dna $align -O fasta -match name > $dna_aln "); # autodetects format 
	$align = $dna_aln;
    }
    close($TMP_HANDLE1); # open 3 files 
    close($TMP_HANDLE2); # .... 
    close($TMP_HANDLE3); # close 3 files 

    # remove gaps 
    # only works for fasta 

    if ( $args->{'-gaps'}==0 ) {
	my %hash;
	open(my $fh, $align);
	{
	    local $/ = ">";
	    <$fh>;
	    while ( my $s = <$fh> ) {
		chomp($s);
		my ($id, @r)=split/\n/, $s;
		$id =~ /^(\w+)/ || die($s);
		$hash{ $1 } = [split//, join('', @r)];
	    }
	}
	close $fh;

	my ($x) = keys %hash;
      POS:for (my $i = $#{$hash{$x}}; $i >= 0; $i--) {
	    foreach my $k (keys %hash) {
		if ( $hash{$k}->[$i] eq '-' ) {
		    map { splice(@{$_}, $i, 1) } values %hash;
		    next POS;
		}
	    }
	}
	return \%hash if $args->{'-format'} eq 'hash';

	my ($TMP_HANDLE9, $gap_aln) = $self->_tempfile("annot_align_gapXXXXXX");	    
	map { print $TMP_HANDLE9 ">$_\n".join('', @{$hash{$_}}) } keys %hash;
	close $TMP_HANDLE9;
	$align = $gap_aln;
    }

    # reformat the align if required. 
    # readseq will detect format fasta / clustal 

    my $rseq = '/Users/devin/Code/bin/readseq.jar';
    $self->throw unless -e $rseq;
    system("java -jar $rseq $align -format=$readseqfmt $range -p > $fmt_aln 2> /dev/null "); # stupid java

    map { $self->_cleanupfile($_) } ($aa, $dna, $aa_aln, $dna_aln); # unlink all files except fmt_aln
    $self->_cleanupfile( $gap_aln ) if $gap_aln ;

    # 

    if ( $args->{'-format'} eq 'hash'  ) {
	my %hash;
	open(my $fh, $fmt_aln);
	{
	    local $/ = ">";
	    <$fh>;
	    while ( my $s = <$fh> ) {
		chomp($s);
		my ($id, @r)=split/\n/, $s;
		$id =~ /^(\w+)/ || die($s);
		$hash{ $1 } = [split//, join('', @r)];
	    }
	}
	close $fh;
	$self->_cleanupfile( $fmt_aln );
	return \%hash;
    } elsif ( $args->{'-format'} =~ /3.2/ ) { # stupid readseq 
	open(my $fh, $fmt_aln); # open 2 files 
	my ($TMP_HANDLE4, $ret_aln) = $self->_tempfile("annot_align_retXXXXXX");	    
	{
	    local $\ = "";	
	    my $head = <$fh>;
	    $head =~ s/I//;
	    print $TMP_HANDLE4 $head;
	    print $TMP_HANDLE4 $_ while <$fh>;
	}
	close($TMP_HANDLE4); # close 2 files 
	close($fh);
	$self->_cleanupfile( $fmt_aln );
	return $ret_aln;
    } else {
	return $fmt_aln;
    }
}

=head2 paml( -object => [ORFS], -method => 'codeml|yn00', -range => [],
    -length_ratio => .25, -min_codons => 30, -application => 'fsa')
    
=cut 

sub paml {
    my $self = shift;
    my $args = {@_};

    ########################################
    # some cursory checking 
    ########################################

    $self->throw unless $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->throw unless $self->isa(ref($_)) } @{$args->{'-object'}};
    $args->{'-length_ratio'} = undef unless exists $args->{'-length_ratio'};
    $args->{'-min_codons'} = 30 unless exists $args->{'-min_codons'};
    $args->{'-method'} = ( $#{$args->{'-object'}} == 0 ? 'yn00' : 'codeml' ) 
	unless exists  $args->{'-method'};

    ########################################
    # get app and 4 files ready 
    ########################################

    my $paml = $self->_external_app(-application => 'paml.'.$args->{'-method'} );

    my ($TMP_HANDLE, $pamlFile) = 
	$self->_tempfile("annot_paml_XXXXXX");

    ########################################
    # align and validate 
    ########################################

    my $alnFile = $self->align(
	-molecule => 'codon', 
	-format => 'phylip3.2', 
	-object => $args->{'-object'}, 
	-range => $args->{'-range'}, 
	-application => ($args->{'-application'} || 'fsa')
	);

    if ( defined $args->{'-length_ratio'} ) {
	chomp(my $alnLen = `head -1 $alnFile | cut -f 3 -d ' ' `);
	my $ratio = $alnLen / $args->{'-object'}->[0]->length; 
	#print $alnLen, $ratio;
	return () if $ratio <= $args->{'-length_ratio'} || 
	    $alnLen/3 <= $args->{'-min_codons'};
    }

    ########################################
    # run PAML 
    ########################################

    my $treeFile = $self->_writeNewickTree( @{$args->{'-object'}} ) 
	if $args->{'-method'} eq 'codeml';
    my $ctrlFile = &_writePamlControlFile($args->{'-method'}, $alnFile, $pamlFile, $treeFile); # !! not a tmp file !! 
    my @files = ($ctrlFile, $alnFile, ($args->{'-method'} eq 'codeml' ? $treeFile : ()));

    # run if we have what we need 
    # immediatley start to clean up 

    system("$paml &> /dev/null ") if (-s $alnFile && -s $ctrlFile);
    map { $self->_cleanupfile($_) } @files; # delete 3 + close 1 
    close($TMP_HANDLE); # we must still delete this 1 later 

    ########################################
    # parse output if it exists 
    ########################################

    my (%ka, %ks);
    my ($ka, $ks) = (undef, undef);
    if ( -s $pamlFile ) {
	if ( $args->{'-method'} eq 'codeml' ) {
	    ($ka) = map {s/\s+//g; $_} `grep 'tree length for dN:' $pamlFile | cut -d ':' -f 2`;
	    ($ks) = map {s/\s+//g; $_} `grep 'tree length for dS:' $pamlFile | cut -d ':' -f 2`;
	} else {
	    my @orgs = map { $_->organism }  @{$args->{'-object'}};
	    my @pw = `grep '\+\-' $pamlFile | grep -v SE`; # assumes one PW comparison
	    foreach my $line ( @pw ) {
		my @r =  split/\s+/, $line;
		($ka, $ks) = ($r[8], $r[11]);
		my $org1 = $orgs[ $r[1]-1 ];
		my $org2 = $orgs[ $r[2]-1 ];
		$ks{$org1}{$org2}=$ks;
		$ks{$org2}{$org1}=$ks;
		$ka{$org1}{$org2}=$ka;
		$ka{$org2}{$org1}=$ka;
		# print $ka, $ks, $org1, $org2, $ks{$org2}{$org1}, $ka{$org2}{$org1};
	    }
	    
	    ########################################
	    # this is a terrible hack. DEVIN. 
	    # To be replaced by proper Orthogroup class. 
	    $self->data( KA_MATRIX => \%ka ) if $self->orthogroup;
	    $self->data( KS_MATRIX => \%ks ) if $self->orthogroup;
	    ########################################

	}
    } else { $self->warn($self->name." $paml $! "); } # complain 

    $self->_cleanupfile($pamlFile); # unlink the last file 
    return ($ka, $ks);
}

sub _writeNewickTree {
    my $self = shift;
    my $tree = $self->up->up->tree; # string 
    $self->throw unless $tree;

    foreach my $x ( $self, @_ ) { # @_ checked by caller 
	my $org = $x->up->up->organism;
	my $iid = $x->name;
	#print $tree, $org, $iid;
	#$tree =~ s/$org/$iid/ || $self->throw("$org, $iid, $tree");
    }
    $tree .= ';' unless $tree =~ /\;$/;

    my ($TMP_HANDLE, $treeFile) = 
	$self->_tempfile("annot_tree_XXXXXX");
    print $TMP_HANDLE $tree;
    close $TMP_HANDLE;
    return $treeFile;
}

sub _writePamlControlFile {
    my ($x, $in, $out, $tree) = @_;
    my $retfile = $x.".ctl";

    open(PAML, ">$retfile");
    if ( $x eq 'yn00' ) {
	print PAML "
      seqfile = $in * sequence data file name
      outfile = $out          * main result file
      verbose = 0  * 1: detailed output (list sequences), 0: concise output

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*       ndata = 1
";
    } else {
	print PAML "
      seqfile = $in * sequence data filename
     treefile = $tree      * tree structure file name
      outfile = $out   * main result file name

        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

*        ndata = 10
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = dat/mtArt.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 0
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0   * set to 3 for non-neutral branch site 
      	      	   * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1  * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0  * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 4  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.

";
    }
    close(PAML);

    return( $retfile );
}

=head2 lda(-intergenic => 300, -minimum => 10)

    Run linear discriminant analysis on putative genes ('HYPO', 'NOVEL', 'HSP')
    and use codon usage freqs (based on 'REAL' ORF objects and intergenics) 
    to classify.

    Uses only intergenics above specified size. 

=cut

sub lda {
    my $self = shift;
    my $args = {@_};
    
    #$self->throw("Check parsing");
    
    $args->{'-minimum'} = 100 unless exists $args->{'-minimum'};
    $args->{'-intergenic'} = 300 unless exists $args->{'-intergenic'};
    $args->{'-application'} = 'R' unless exists $args->{'-application'};
    $args->{'-parameters'} = '--vanilla --slave' unless exists $args->{'-parameters'};

    # prepare trainign data 

    map { $_->evaluate } $self->orfs;
    my @r = $self->orfs(-restrict => ['REAL']);
    my @n = $self->orfs(-intergenic => $args->{'-intergenic'});			
    $self->warn("Insufficient data for LDA: $#n, $#r") and return $self
	unless ($#r > $args->{'-minimum'} && $#n > $args->{'-minimum'});

    # make training sets 
    
    my ($train,$group, $test);
    foreach my $x (@r, @n) {
	next unless my $freqs = $x->codons;
	next unless ref($freqs) eq 'ARRAY' && $#{$freqs} == 60;
	$train .= join("\t", @{$freqs})."\n";
	$group .= $x->assign."\n";
    }
    my $TRAIN = $self->_write_temp_file(-data => \$train);
    my $GROUP = $self->_write_temp_file(-data => \$group);

    # lets work on the test set. hash to manage refs
    
    my %hash;
    my $count;
    foreach my $x ($self->orfs(-restrict => ['HYPO', 'NOVEL', 'HSP'])) {
	next unless my $freqs = $x->codons;
	next unless ref($freqs) eq 'ARRAY' && $#{$freqs} == 60;
	$test .= join("\t", @{$freqs})."\n";
	$hash{++$count} = $x;
    }	
    my $TEST = $self->_write_temp_file(-data => \$test);
    
    # R code 
    
    my $command = "
library(MASS)
list<-scan(\'$GROUP\', what=\'character\')
group<-factor(list)
train<-read.table(\'$TRAIN\', header=FALSE)
matrix<-as.matrix(train)
ldafac<-lda(matrix, group, tol=0)
test<-read.table(\'$TEST\', header=FALSE)
q=predict(ldafac, test, dimen=1)
write.table(q$posterior, quote=F, col.names=F, row.names=F)";

    my $R = $self->_write_temp_file(-data => \$command);
    
    # run and retrieve 
    
    my $binary = $self->_external_app(-application => $args->{'-application'});
    $self->throw("Run error: $binary $args->{'-parameters'} <$R") unless 
        my $r = `$binary $args->{'-parameters'} <$R`;
    map { $self->_cleanupfile($_) } ($TRAIN, $GROUP, $TEST, $R);
    return undef unless $r;

    $self->throw unless scalar(keys %hash) == scalar(split/\n/,$r);

    my $count;
    foreach my $line (split/\n/, $r) {
	my @r = split/\s+/, $line;
	my $pcoding = sprintf("%.4f",$r[2]*1);
	$hash{++$count}->data('LDA' => $pcoding);
	$hash{$count}->evaluate;
    }
	
    return $self;
}

sub _write_temp_file {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-directory'} = '/tmp/' unless exists $args->{'-directory'};
    $args->{'-unlink'} = 1 unless exists $args->{'-unlink'};
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    
    $self->throw("no data provided: $args->{'-data'}") unless 
	exists $args->{'-data'} && ref($args->{'-data'}) eq 'SCALAR';

    # 

    my @meth = split/::/, (caller(1))[3]; # / 
    my $filename = join("_",("annot",  substr($meth[-1], 0, 6) ,"XXXXXX"));
    my ($TMP_HANDLE, $TMP_FILE) = $self->_tempfile($filename);
    print ++$totalfilecount, &_open_file_count(), $self->_method, $TMP_FILE if $args->{'-verbose'};
    print $TMP_HANDLE ${$args->{'-data'}};    
    close $TMP_HANDLE;

    return ($TMP_FILE);
}

sub _open_file_count {
    return  scalar(split/\n/, `lsof | grep $$  `) - 9;
    # binary (1) : /opt/local/bin/perl5.12
    # lib (3) : /usr/lib/libutil1.0.dylib : /usr/lib/dyld : /private/var/db/dyld/dyld_shared_cache_x86_64
    # tty (3) : 3 x /dev/ttys004 (STDIN/OUT/ERR)  
    # dir (1) : working dir 
    # script (1) : 
}

=head2 _write_temp_seq(-molecule => 'aa|dna')

=cut 

sub _write_temp_seq {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-directory'} = '/tmp/' unless exists $args->{'-directory'};
    $args->{'-unlink'} = 1 unless exists $args->{'-unlink'};
    $args->{'-molecule'} = (ref($self) =~ /Orf/i ? 'aa' : 'dna' ) 
	unless exists $args->{'-molecule'};
    $args->{'-name'} = undef unless exists  $args->{'-name'};
    $args->{'-id'} = 'internal' unless exists $args->{'-id'};
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};

    push @{$args->{'-seq'}}, $self unless exists $args->{'-seq'};
    map {$self->throw("") unless /::/ } @{$args->{'-seq'}}; # / :\/ #  
	     
    # 

    my @meth = split/::/, (caller(1))[3]; # / 
    my $filename = join("_",("annot", substr($meth[-1], 0, 6) ,"XXXXXX"));
    my ($TMP_HANDLE, $TMP_FILE) = $self->_tempfile($filename);
    print ++$totalfilecount,&_open_file_count(), $self->_method, $TMP_FILE if $args->{'-verbose'};

    # 

    foreach my $seq (@{$args->{'-seq'}}) {
	# tRNA-scan and others choke on long lines (somewhere around 500Kb) 
	# so write out in classic fasta format. 
	# we trim terminal * for HMMER. 
	
	my $id = $seq->_internal_id;
	if ( $args->{'-id'} eq 'name') {
	    $id = $seq->name;
	} elsif ( $args->{'-id'} eq 'organism' ) {
	    my $obj = $seq->up;
	    $obj = $obj->up until (ref($obj) =~ /Genome/i);
	    $id = $obj->organism;
	}

	print $TMP_HANDLE ">$id";
	
	my @seq = split//, ( 
	    ref($seq) =~ /ORF/i
	    ? $seq->sequence(-nostop => 1, -molecule => $args->{'-molecule'})
	    : $seq->sequence # contig 
	);   
	
	my $line_length=60;
	until ($#seq==-1) {
	    my $lim = (scalar(@seq) > $line_length ? $line_length : scalar(@seq));
	    my $tseq = join('', splice(@seq, 0, $lim));
	    print $TMP_HANDLE $tseq;	    
	    # print $#seq, $lim, $tseq; 
	}
    }	 
    close $TMP_HANDLE;
    
    return ($TMP_FILE);
}

sub _tempfile {
    my $self=shift;
    my $filename=shift;
    
    # we ignore any other params..

    unless ( $TEMP_DIR_STAMP ) {
	$TEMP_DIR_STAMP='ANNA'.time;
	mkdir("/tmp/$TEMP_DIR_STAMP/") ||  $self->throw("/tmp/$TEMP_DIR_STAMP/");
    }

    if ( $filename =~ /(X+)$/) {
	my $n = length($1);
	my $val = int(rand(9)+1);
	$val .= int(rand(9)+1) until length($val)>=$n;
	$filename =~ s/X{$n}$/$val/;
    }
    
    my $filename = "/tmp/$TEMP_DIR_STAMP/$filename";

    #$self->throw("$filename exists") if -e $filename;
    $filename.=int(rand(9)+1) while -e $filename;
    open(my $fh, ">$filename") ||  $self->throw("$filename open() fail");
    $self->throw("$filename not created") unless -e $filename;
    
    return($fh,$filename);
}

sub _cleanupfile {
    my $self = shift;
    my $file = shift;
    die($file) unless -e $file;
    unlink($file);
    return 1;
}

sub _external_app {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-application'} = 'blast.formatdb' unless exists $args->{'-application'};
    return $APPS{ $args->{'-application'} } if exists $APPS{ $args->{'-application'} };

    if ( $args->{'-application'} =~ /\w+\.(\w+)/ ) {          
        my $alt = $1;                                
        chomp(my $binary = `which $alt`);            
        if ( $binary && -e $binary && -x $binary ) {               
            $args->{'-application'} = $alt;          
        }                                            
    }

    # find app and ensure executable 
    
    $self->throw("Cannot locate $args->{'-application'} ($binary) - $!")
	unless chomp(my $binary = `which $args->{'-application'}`);
    $self->throw("Binary not executable: $binary - $!") 
        unless -e $binary && -x $binary;

    $APPS{ $args->{'-application'} } = $binary;

    return $binary;
}

1;
