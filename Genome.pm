#!/usr/bin/perl

package Annotation::Genome;

use Annotation;
use GlobalVars;

use Storable qw(retrieve nstore);
use Roman;

$Storable::canonical=1;

@ISA = qw(Annotation);

use Annotation::Contig;
use Annotation::Orf;

#########################################
#########################################

# autoloaded methods 

our $AUTOLOAD;

my %AUTOMETH = (
    organism => undef,
    strain  => undef,
    sample => undef,
    coverage => undef,
    source => undef,
    wgd => undef
); 

sub AUTOLOAD {
    my $self = shift;    
    my $name = uc((split(/\:/, $AUTOLOAD))[-1]); # much faster 
    #$name =~ s/.*://;   # strip fully-qualified portion
    $self->{$name} = shift if @_; # set 
    return($self->{$name});       # get 
}                   

my $STORE = 0;

#########################################
#########################################
 
=head2 new()

    ORGANISM => 'Kluyveromyces polysporus',
    STRAIN => 'DSMZ 70294',
    SAMPLE => 1, # in case multiple individuals from same strain are sequenced
    COVERAGE => 7.8, [ previously quality ]
    SOURCE => 'Scannell et al. PNAS 2007'

    These attributes are set automatically: 

    ANNOTATION => user,
    DATE => date,
    _VERSION => software vresion number,

=cut 

sub new {
    my $class = shift;        
    my $args = {@_};

    my $base = Annotation::Contig->new(ID => 0);
    my $self  = $class->SUPER::new(@_, DOWN => $base);  
    $base->up($self);
    
    foreach my $k (keys %AUTOMETH) {
	$k = uc($k);
	if ($self->id == 0) {
	    $self->{$k} = $AUTOMETH{$k};
	} else {
	    $self->throw("$k undefined -- $k") 
		unless exists $args->{$k} && defined $args->{$k};
	    $self->{$k} = $args->{$k};
	}
    }    
    
    $self->{'ANNOTATION'} = $ENV{'USER'} unless exists $self->{'ANNOTATION'};
    $self->{'HISTORY'} = [] unless exists $self->{'HISTORY'}; # {DATE=>, FILE=>, SCRIPT=>}
    # some of these will be depracated 
    chomp(my $date = `date`);
    $self->{'DATE'} = $date unless exists $self->{'DATE'};
    $self->{'_VERSION'} = $VERSION unless exists $self->{'_VERSION'};
    $self->{'_STAMP'} = $TIME unless exists $self->{'_STAMP'};
    
    bless $self, $class;
}

#########################################
#########################################

=head2 annotate( -library => HmmLibPath )
    
    Run feature finding, open reading frame prediction, homology-based
    prediction and reconciliation methods using generic params.

=cut 

sub annotate {
    my $self = shift;        
    my $args = {@_};

    $args->{'-library'} = $ENV{'ANNA_HMMER2_LIB'} unless exists $args->{'-library'};
    $args->{'-restart'} = undef unless exists $args->{'-restart'};
    $args->{'-skip'} = 0 unless exists  $args->{'-skip'};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};
    $args->{'-overlap'} = 0.8 unless exists $args->{'-overlap'};

    goto uc($args->{'-restart'}) if $args->{'-restart'};

    ############################################
    # 1. basic prediction and basic BLASTing 
    ############################################

  CORE:     

    # place non-ORF feature models. 
    # for protein search -log2(db_size/expected_hits) is a good heuristic
    # all numbers here are heuristics. 

    $self->features(
        -intron => [$args->{'-library'}.'/INTRON.hmm', 4.5], 
        #-telomere => [$ENV{'ANNA_HMMER2_LIB'}.'/TEL.hmm', 10],
	#-ars => [$hmm.'/ARS.hmm', -1],   
	#-centromere => [$hmm.'/CEN.hmm', -1]
	);
    
    # find all orfs above context specific size 
    
    foreach my $contig ($self->stream) {
        $contig->predict(
            -orf_min => 300,       # this is standard. find short / introns later. 
            -term_min => 60,       # scer 70-75% coding. gene at break >70%
            -exon_min => 30,       # better check this but should be ok
            -feat_min => 30,       # noncoding. telomeres etc 
            -intron_min => 30,     # 50 more accurate?
            -intron_max => 1500,   # scer max 1050. mode 100. median 200. 	    
	    # these should be moved elsewhere 
            -trnascan => 20, # COVE score. 
	    #-snoscan  => 14  # COVE score. 
            );
	
	# seek homology/YGOB info by running BLAST and HMMER. 
	# these data required to compute synteny scores below. 
	# we cannot run exonerate2 at this point as it uses
	# LOSS to identify pillars and homologs. 
	# we just populate the homology arrays. 
	# update() runs evaluate() by default 

	map { $_->update( -store => 1, -verbose => 2, -exonerate => undef, -hsp => undef ) } $contig->orfs; 
    }

    $self->backup( -tag => 'core' );
    
    ############################################
    # 2. build synteny (LOSS) matrices. 
    # required to use impute() / pillar() / homolog() etc. etc
    # these in turn allow homology based gene prediction etc 
    ############################################

  HOMOLOGY:

    $self->_synteny( -nulldist => 2 );    

    # 0. add homology information foreach gene. 
    # 1. compute a synteny matrix. 
    # 2. use syn mat to perform more carefully gene identification / calling 
    # 3. identify YGOB pillar so we can call homolog() etc
    # convert to an iterative version like syntenic_orthologs()?

    foreach my $orf ( $self->orfs ) {
	next if $orf->rank < -1;
	next unless $orf->id >= $args->{'-skip'};

	# this is nasty. but running BLAST/HMMER again is worse. 
	# correct by implementing lightweight HIT object that permits 
	# holding adequate hits for inference (20?). 

	my $blast = $orf->data('_STORE_BLAST');
	$orf->data( '_STORE_BLAST' => undef ); ## ugh.
	$orf->impute( -blast => $blast ) if $blast; 
	undef $blast;

	# with benefit of synteny we now choose YGOB pillar, best homolog etc.
	# can run local/global laignment scores etc. these can in turn be used 
	# to compose the fragment test. 

	$orf->pillar( -force => 1 );
        next unless my $homolog = 
	    $orf->homolog( -verbose => 1, -model => 'global' ); # -fast not possible 

	$orf->exonerate2( -model => 'global', -return => 'score' );
	$orf->exonerate2( -model => 'local', -return => 'score' );
	$orf->fragment();

	# becuase we now have access to proper homologs, we can 
	# reoptimise gene prediction using WISE/exonerate (reoptimise).
	# in safe mode reoptimise() does not alter the caller but 
	# returns a new model. we keep both for VARIANTS section below. 
	# the new models are added to the genome automatically and
	# reoptimize calls update() so the new models also have
	# BLAST info. the old model does not have update called on it
	# hence, we run exonerate2() above. 

	$orf->reoptimise( -reference => $homolog, -safe => 1, -verbose => $args->{'-verbose'} );
    }
    $self->backup( -tag => 'homology' );

    ############################################
    # 3. purge and merge
    ############################################

  VARIANTS:

    foreach my $contig ( $self->stream ) {
	$contig->purge2( -overlap => $args->{'-overlap'}, -merge => 1 );
    }
    
    $self->backup( -tag => 'variants' );

    return $self;
}


=head2 dumpOrthogroups(-all => 0, -clade => 'Stricto')

    Walk along reference genome and output identifying info for 
    homologs etc. 

=cut 

sub dumpOrthogroups {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-all'} = undef unless exists $args->{'-all'}; # print genes not in orthogroups 
    $args->{'-clade'} = 'Stricto' unless exists $args->{'-clade'};

    $self->throw unless $self->bound;

    my @org = map { /(\w)(\w+)/; uc($1).lc($2); } ($self->organism, $self->bound);

    my @header = (
	'#Name','Orthogroup','Gene_Predicted','YGOB_Predicted','SGD','SGD_Status',
	(map { $_.'_Name' } @org),
	(map { $_.'_Gene' } @org),
	(map { $_.'_YGOB' } @org),
	(map { $_.'_Length' } @org),
	(map { $_.'_Fragments' } @org),
	(map { $_.'_Synteny' } @org)
	);

    my $log = $self->organism.'.allgenes';
    open($fh, ">$log") || die;
    print $fh @header; 

    foreach my $orf ( map { $_->stream } $self->stream ) {
	
	next unless ( (my @og = $orf->orthogroup) || $args->{'-all'} );	
	unshift(@og, $orf);
	push @og, ( (undef) x scalar($self->bound) ) unless $#og == scalar($self->bound);

	my ($sgd_name,$ygob_name) = $orf->identify;
	my $sgd = ( $orf->_debug ? $orf->_debug : undef );
	my ($official, $class)=( $sgd  ?  ($sgd->gene, $sgd->data('CLASS')) : ('NA','NA') );

	print $fh 
	    $orf->name, ($orf->ogid || 'NA'),
	    ($sgd_name || 'NA'), ($ygob_name || 'NA'), 
	    ($official || 'NA'), ($class || 'NA'),
	    (map { ( defined($_) ? ($_->name || die) : 'NA') } @og), 
	    (map { ( defined($_) ? ($_->gene || 'NA') : 'NA') } @og), 
	    (map { ( defined($_) ? ($_->ygob || 'NA') : 'NA') } @og),
	    (map { ( defined($_) ? ($_->length || die) : 'NA') } @og),
	    (map { ( defined($_) ? ($_->exons || die) : 'NA') } @og),
	    (map { ( defined($_) ? ($_->loss || 0) : 'NA') } @og);
    }
    close $fh;

    return $self;
}

=head2 dumpOrthologousCoding()

    We impose a relatively high synteny burden here. We are not short coding regions
    and we do not want shitty data confusing results. 

=cut 

sub dumpOrthologousCoding {
    my $self = shift;
    my $args = {@_};

    $args->{'-min'} = 30 unless exists  $args->{'-min'}; # alpha factor 
    $args->{'-max'} = 15000 unless exists  $args->{'-max'}; # IRA1/2 
    $args->{'-gap'} = .25 unless exists  $args->{'-gap'};
    $args->{'-synteny'} = 2 unless exists $args->{'-synteny'};    
    $args->{'-spanning'} = 1 unless exists $args->{'-spanning'};
    $args->{'-sensustricto'} = undef unless exists  $args->{'-sensustricto'};
    $args->{'-clade'} = 'Stricto' unless exists $args->{'-clade'};
    $args->{'-fasta'} = 1 unless exists $args->{'-fasta'};

    $self->throw unless $self->bound;

    # 

    my %seen;
    my %count;
    my (@ok);
    foreach my $c ( $self->stream ) {
      ORF: foreach my $o ( grep {$_->orthogroup} 
			   grep {$_->coding || $_->assign =~ /PSEUDO/} $c->stream ) {
	  next if exists $seen{ $o->name };
	  my ($sgd_name, $ygob_name) = $o->identify;	  
	  $seen{$ygob_name}++; # record any high quality OG 

	  # do we index by the scer gene or the ygob locus? 
	  # ygob preferred but do scer for S3. 

	  if ( $o->ohnolog =~ /Orf/ && ! $args->{'-sensustricto'} ) { # DEVIN 
	      my $ohno = $o->ohnolog; 
	      $seen{ $ohno->name }++;

	      next unless $ohno->orthogroup;
	      next unless $ohno->ygob;
	      map { next ORF unless $_->ygob eq $o->ygob }  ($ohno,$ohno->orthogroup);
	      map { next ORF unless 	      
			$_->querysynteny(
			    -spanning => $args->{'-spanning'},
			    -difference => 7,
			    -distance => 3,
			    -restrict => ['YGOB']
			) >= $args->{'-synteny'}
	      } ($ohno, $ohno->orthogroup);
	      push @ok, [ $o, $ohno];
	  } else {
	      $seen{ $sgd_name }++;
	      push @ok, [ $o ];	      
	  }

      }
    }


    #############################################
    # 
    my @org = map { /(\w)(\w+)/; uc($1).lc($2); } ($self->organism, $self->bound);
    
    my @header = (
	'#Name','Orthogroup','SGD','SGD_Status','Gene_Predicted','YGOB_Predicted','CopyNumber', 'Quality',
	(map { $_.'_Name' } @org),
	(map { $_.'_Gene' } @org),
	(map { $_.'_YGOB' } @org),
	(map { $_.'_Length' } @org),
	(map { $_.'_Fragments' } @org),
	(map { $_.'_Synteny' } @org)
	);
    
    my $log = $args->{'-clade'}.".orthologs";
    open(my $fh, ">$log") || die;
    print $fh @header;

    #############################################

    # lets cycle through the candidate OGs 
    # check sequence etc. 

  OK: foreach my $og ( @ok ) {
      my @og = map { $_,$_->orthogroup} @{$og};

      my ($sgd_name,$ygob_name)= $og->[0]->identify;
      
      if (  $args->{'-sensustricto'} ) {
	  unless ( $sgd_name ) {
	      map { $_->output(-prepend => ['NOSGD'] ) } @og;
	      $sgd_name = 'NOSGD';
	  }
	  unless (  $seen{ $sgd_name } == 1 || $sgd_name eq 'NOSGD' ) {
	      $og->[0]->output(-prepend => ['NOUNIQ', $sgd_name ] );
	      #next;
	  }
      } else {
	  unless ( $seen{ $og->[0]->ygob } == 1 ) { 
	      $og->[0]->output( );
	      next;
	  }
      }

      $count{'5.Unique'}++;

      # valildate sequences 

      foreach my $o ( @og ) {
	  $o->output and $self->warn if $o->length >= $args->{'-max'};
	  next OK unless $o->length >= $args->{'-min'} &&
	      $o->length <= $args->{'-max'};
	  my %comp;
	  map { $comp{$_}++ } split//, $o->sequence(-molecule => 'dna');
	  next OK if $comp{'N'}/$o->length >= $args->{'-gap'};
      }
      $count{'6.Sequence'}++;

      # Are they high quality 

      my ($m,$sd) = _calcMeanSD( map {$_->length} @og );
      my @int = sort {$a <=> $b} map {$_->data('INTRONS')} @og;
      my @stop = sort {$a <=> $b} map {$_->data('STOP')} @og;
      my @exons = sort {$a <=> $b} map {$_->exons} @og;
      my @ass = sort {$a cmp $b} map {$_->assign} @og; 

      my $tag='highqual';
      if ( ($sd/$m) >= .05 ) {
	  $tag = 'length';
      } elsif ( $int[0]!=$int[-1] || $stop[0]!=$stop[-1] || $exons[0] != $exons[-1] ) {
	  $tag = 'struture';
      } elsif ( $ass[0] ne $ass[-1] ) {
	  $tag = 'assignment';
      }
      $count{'7.HighQual'}++ if $tag eq 'highqual';

      $count{'8.Ohnolog'}++ if $#{$og} == 1;

      # 
      
      my $copynumber = ( $#{$og} == 0 ? 1 : 2 );
      my $name = join('_', 'OG'.(++$counter), $sgd_name, $ygob_name.':'.$copynumber, $tag).'.fsa';
      my $sgd = ($og->[0]->_debug =~ /\:\:/ ? $og->[0]->_debug : undef );
      my ($official, $class)=( $sgd  ?  ($sgd->gene, $sgd->data('CLASS')) : ('NA','NA') );

      print $fh 
	  $name, $counter, 
	  ($official || 'NA'), ($class || 'NA'),
	  ($sgd_name || 'NA'), ($ygob_name || 'NA'), 
	  $copynumber, $tag, 
	  (map {$_->name || die } grep {defined} @og), 
	  (map {$_->gene || 'NA' } grep {defined} @og), 
	  (map { $_->ygob || 'NA' } grep {defined} @og),
	  (map {$_->length || die } grep {defined} @og),
	  (map {$_->exons || die } grep {defined} @og),
	  (map {$_->querysynteny(-distance => 3, -difference => 7, -spanning => 0) || 0 } grep {defined} @og);

      # 

      map { $_->ogid($counter) unless $_->ogid > 1 } @og;

      # output 

      next unless $args->{'-fasta'};

      print STDERR $name if -e $name;
      open(my $fsa, '>'.$name);
      foreach my $o ( @og ) {
	  print $fsa $o->sequence(-format => 'fasta', -decorate => 1, -molecule => 'dna');
      }
      close $fsa;
  }
    close $fh;
    
    map { print $_, $count{$_} } sort keys %count if $args->{'-verbose'};

    return $self;
}

=head2 dumpOrthologousIntergenics(-ignorerepeats => 0, 
    -min => 30, -max => 3000, -gap => .5)

=cut 

sub dumpOrthologousIntergenics {
    my $self = shift;
    my $args = {@_};

    $args->{'-min'} = 30 unless exists  $args->{'-min'};
    $args->{'-max'} = 3000 unless exists  $args->{'-max'};
    $args->{'-gap'} = .50 unless exists  $args->{'-gap'};

    $args->{'-ignorerepeats'} = undef unless $args->{'-ignorerepeats'};
    $args->{'-clade'} = 'Stricto' unless exists $args->{'-clade'};
    $args->{'-exclude'} = undef unless exists $args->{'-exclude'};

    $self->throw unless my @species = grep { uc($_) ne uc($args->{'-exclude'}) } $self->bound;

    my $log = $args->{'-clade'}.'.intergenics';
    unless ( $args->{'-noalign'} ) { open(LOG, ">$log") || die };

    # we give ohnologs unique ids so we can distinguish.
    # all non-ohnologs are considered 1-of-1.

    my %ohnoid;
    foreach my $o ( grep { $_->ohnolog } $self->orfs ) {
	next if exists $ohnoid{ $o->name };
	$ohnoid{ $o->name } = 1;
	$ohnoid{ $o->ohnolog->name } = 2;
    }

    # 

    my %seen;
    my %count;
    foreach my $c ( $self->stream ) {
      ORF: foreach my $o ( $c->stream ) {
	  $count{'0.FEAT'}++;
	  next unless $o->assign ne 'GAP';

	  $count{'1.ORF'}++;
	  
	  next unless $o->orthogroup;
	  $count{'2.OG'}++;
	  
	  # get a bounding pair 
	  
	  my ($left,$right);
	  if ( $left = $o->neighbour(-direction => 'left' ) ) {
	      if ( $args->{'-ignorerepeats'} ) {
		  $left = $left->neighbour(-direction => 'left' ) while 
		      ($left && $left->evidence =~ /TY|LTR/ && ! $left->ogid );
		  next ORF unless $left;
	      }
	      next unless $left->orthogroup;
	      next if $left->overlap(-object => $o);
	      $right = $o;
	  } else { next; }
	  
	  next unless scalar( $left->orthogroup ) == scalar($self->bound);
	  next unless scalar( $right->orthogroup ) == scalar($self->bound);
	  $count{'3.BOUND'}++;
	  
	  # 
	  
	  foreach my $sp ( @species ) {
	      next ORF unless my ($Lsp,$Rsp) = $right->$sp->neighbours; # right gene. neighbours of ortho 
	      next ORF unless my $left_sp = $left->$sp; # left gene. ortho only. 
	      my @match = grep { $_ eq $left_sp } ($Lsp, $Rsp); # is left_sp == [LR]sp ? 		
	      next ORF unless $#match == 0;
	  }
	  $count{'5.SYNORTH'}++;
	  
	  # check for overlaps 
	  
	  foreach my $sp ( @species ) {
	      my $left_sp = $left->$sp; 
	      my $right_sp = $right->$sp;
	      next ORF if $left_sp->overlap(-object => $right_sp);		
	  }
	  $count{'6.NOOLAP'}++;	
	  
	  # check for orientation  
	  
	  foreach my $sp ( @species ) {
	      my $left_sp = $left->$sp;
	      my $right_sp = $right->$sp;		
	      my $fac = ( $left_sp->id < $right_sp->id ? 1 : -1);
	      next ORF unless $left_sp->strand == $left->strand*$fac && $right_sp->strand == $right->strand*$fac;
	  }	
	  $count{'7.ORIENT'}++;
	  
	  #################################
	  
	  # get base names for neighboring genes 

#	  my %sgc;
#	  map { $sgc{$_->data('SGD')}++ } ($o,$o->orthogroup);
#	  my ($sgd_name) = grep {/Y/} sort { $sgc{$b} <=> $sgc{$a} } keys %sgc;
	  
	  my %rname;
	  foreach my $g ($right, $right->orthogroup) {
	      if ( $args->{'-clade'} eq 'Stricto' ) { 
		  $rname{ ( $g->data('GENE') =~ /Y[A-P][LR]\d{3}/ ? $g->data('GENE') : $g->data('SGD') ) }++;
	      } else {
		  $rname{ ( $g->ygob ? $g->ygob : $g->data('GENE') ) }++;
	      }
	  }
	  my ($rname) = sort { $rname{$b} <=> $rname{$a} } keys %rname;
	  my %lname;
	  foreach my $g ($left, $left->orthogroup) {
	      if ( $args->{'-clade'} eq 'Stricto' ) { 
		  $lname{ ( $g->data('GENE') =~ /Y[A-P][LR]\d{3}/ ? $g->data('GENE') : $g->data('SGD') ) }++;
	      } else {
		  $lname{ ( $g->ygob ? $g->ygob : $g->data('GENE') ) }++;
	      }
	  }
	  my ($lname) = sort { $lname{$b} <=> $lname{$a} } keys %lname;

	  next unless $rname && $lname;
	  next unless ($rname =~ /Anc|\:/ && $lname =~ /Anc|\:/) || $args->{'-clade'} eq 'Stricto';
	  
	  # get relationhsip of intergenic to neighbouring genes  
	  # get neighbour copy-number 
	  
	  my $rdir = ( $right->strand == 1 ? 'u' : 'd');
	  my $ldir = ( $left->strand == 1 ? 'd' : 'u');	
	  my $r_reg = ($rdir eq 'u' ? 'upstream' : 'downstream');
	  my $l_reg = ($ldir eq 'u' ? 'upstream' : 'downstream');	  
	  my $lnum = ( $left->ohnolog ? $ohnoid{ $left->name }.'of2' : '1of1');
	  my $rnum = ( $right->ohnolog ? $ohnoid{ $right->name }.'of2' : '1of1');
	  
	  # make bounding genes names. 

	  my $bound1 = $lname.$ldir.($args->{'-clade'} eq 'Stricto' ? undef : $lnum); 
	  my $bound2 = $rname.$rdir.($args->{'-clade'} eq 'Stricto' ? undef : $lnum); 
	  my $temp = $bound1.'_'.$bound2;
	  
	  # decide on the output orientation. 
	  # we want diff clades to name files reproducibly.
	  # eg Anc_1.123 _ Anc_1.124 NOT Anc_1.124 _ Anc_1.123
	  # we must also deal with tRNA naming system. 
	  
	  my $ori;
	  if ( $temp =~ /Anc_(\d+)\.(\d+).+_Anc_(\d+)\.(\d+)/ ) {
	      if ( $1 < $3 ) {
		  $ori = 1;
	      } elsif ( $1 > $3 ) {
		  $ori = -1;
	      } elsif ( $4 < $2 ) {
		  $ori = -1;
	      } else {
		  $ori = 1;
	      } 
	  } elsif ( $temp =~ /(\w{3}\:\w{3}).+_(\w{3}\:\w{3})/ ) {
	      my ($x,$y) = sort ($1,$2);
	      $ori = ( $temp =~ /^$x/ ? 1 : -1);
	  } elsif ( $temp =~ /\w{3}\:\w{3}.+_Anc_\d+\.\d+/ ) {
	      $ori = -1;
	  } elsif ( $temp =~ /Anc_\d+\.\d+.+_\w{3}\:\w{3}/ ) {
	      $ori = 1;
	  } else { 
	      #$ori = 1; # SENSUSTRICTO 
	      $self->warn($temp); next;
	  }
	  
	  # make the name and output info 
	  
	  my (@log,$intername);
	  if (  $ori == 1  ) {
	      $intername = $bound1.'_'.$bound2;
	      push @log, $left->name, $left->strand, $l_reg, $lname, $lnum,
	      $right->name, $right->strand, $r_reg, $rname, $rnum;
	  } else {
	      $intername = $bound2.'_'.$bound1;
	      push @log, $right->name, $right->strand, $r_reg, $rname, $rnum, 
	      $left->name, $left->strand, $l_reg, $lname, $lnum;
	  }
	  $self->warn("$intername : $seen{$intername}") if exists $seen{$intername};
	  $seen{$intername}++;
	  
	  # write to logfile 
	  
	  print LOG $intername, @log;
	  print $intername, @log;
	  
	  #######################################
	  # make output. 
	  #######################################
	  
	  # prep the data..

	  my %data;
	  foreach my $sp ( $self->organism, @species ) {

	      my ($l_sp, $r_sp) = ( $sp eq $self->organism ? ($left,$right) : ($left->$sp,$right->$sp) );
	      my $order = ( $l_sp->id < $r_sp->id ? 1 : -1 ); # same as ref or diff
	      
	      # get seq 
	      
	      my $seq = $l_sp->intergenic(-object => $r_sp, -strand => $ori*$order);
	      next unless length( $seq ) >= $args->{'-min'} && length( $seq ) <= $args->{'-max'};

	      # qc 
	      
	      my %comp;
	      map { $comp{$_}++; $comp{'TOT'}++ } split//, $seq;
	      next if $comp{'N'}/$comp{'TOT'} > $args->{'-gap'};

	      # 
	      
	      $sp =~ /(\w)(\w+)/;
	      $sp = uc($1).lc($2);
	      my $string = ">$sp:$intername [".
		  join('/', (  $ori*$order == 1 
			       ? (1, $l_sp->name, $l_sp->strand, $r_sp->name, $r_sp->strand) 
			       : (-1, $r_sp->name, $r_sp->strand*-1, $l_sp->name, $l_sp->strand*-1)
		       ))."]";
	      
	      $data{$string} = $seq;
	  }
	  next unless scalar( keys %data ) == scalar(@species)+1;	  
	  $count{'8.SEQ_OK'}++;
	  next if $args->{'-noalign'};

	  open (my $interFH, ">$intername".'.fsa') || $self->throw("Could not open: $intername .fsa"); 
	  map { print $interFH "$_\n".$data{$_}->sequence } keys %data;
	  close($interFH);
      }
    }
    close LOG;

    map { print $_,$count{$_} } sort keys %count if $args->{'-verbose'};
    
    return $self;
}


=head2 dumpDataForYGOB(-norepeats => undef) 

    prepare 3 fiels for YGOB:
    _data.tab : unique_id, strand, chr, id
    _orth.tab : unique_id, anc ortholog is exists 
    _ohno.tab : unique_id, ohnolog_unique_id

=cut

sub dumpDataForYGOB {
    my $self = shift;
    my $args = {@_};

    $args->{'-ignorerepeats'} = undef unless  $args->{'-ignorerepeats'};

    $self->throw( $self->organism ) unless my @species = $self->bound;
    
    foreach my $genome ( $self, map { $self->bound($_) } @species ) {
	
	my $record = scalar( $genome->orfs );
	my $og =  scalar(grep { $_->ogid } $genome->orfs);
	$genome->throw("No Orthogroups") unless $og;
	#$genome->warn($genome->organism." few Orthogroups ( $og / $record ): May malfucntion") unless $og/$record >= .5;
	
	# 
	
	my %hash;
	open(FSA, ">".$genome->organism."_ygob.aa");
	open(my $fh, ">".lc($genome->organism)."_data.tab");
	open(my $ohno, ">".lc($genome->organism)."_ohno.tab");
	
	my %uniq;
	foreach my $o  ( grep { $_->assign ne 'GAP' } $genome->orfs ) {
	    my $c = $o->up;
	    
	    next unless $o->assign ne 'REPEAT' || ! $args->{'-ignorerepeats'};
	    
	    # data file 
	    my $id = ( $genome->organism eq 'Scer' ? ($o->data('SGD') || $c->id.'.'.$o->id) : $c->id.'.'.$o->id );
	    $id = $o->data('GENE') if $o->assign =~ /RNA/;
	    
	    my $string = $o->output(
		-string => 1, 
		-append => [$o->ygob, ($o->ohnolog ? $o->ohnolog->name : undef), 
			    $o->data('SGD'), $o->data('SCER'), $o->ontology(-attribute => 'process')]
		);
	    $string =~ s/\t/\; /g;
	    print $fh $o->name, $o->strand, $c->id, $id , $string;
	    
	    # ohnologs 	    
	    print $ohno $o->name, $o->ohnolog->name and $uniq{$o->name}++ 
		if ( $o->ohnolog && ! exists $uniq{$o->ohnolog->name});
	    
	    # orthologs 
	    # this is now handled below 
	    
	    next if $o->assign =~ /RNA/; # these are not handled for homology etc 
	    print FSA $o->sequence(-molecule => 'aa', -format => 'fasta', -decorate => 1);	    
	}
	close $fh;
	close $orth;
	close $ohno;
	close FSA;
    }
    
    # we now implement the orthologs policy from dumpOrthoGroups;

    my %fh;
    foreach my $x ( $self, map { $self->bound($_) } @species ) {
	open(my $fh, ">".lc($x->organism)."_orth.tab");	
	$fh{$x->organism}=$fh;
    }

    foreach my $orf ( grep { $_->orthogroup } map { $_->stream } $self->stream ) {
	my %id;
	map { $id{ ( $_->ygob ? $_->ygob : $_->data('GENE') ) }++ } ($orf, $orf->orthogroup);
	my @order = sort { $id{$b} <=> $id{$a} } keys %id;
	
	# we want to represent all orthogroups but how do we include those that 
	# are not well mapped to YGOB oethogrousp?
	# we take a policy where if a majority of genes map to a column, then we 
	# accept that. otherwise we do not map into columns, but map the genes in the OG
	# to each other (via the ref). 
	
	if ( $order[0] =~ /Anc/ && ( scalar(keys %id)==1 || $id{$order[0]}>scalar($self, $self->bound)/2 ) ) {
	    foreach my $g ($orf,  $orf->orthogroup) {
		my $fh = $fh{$g->organism}; 
		print $fh $g->name, $order[0];
	    }
	} else {
	    map { my $fh = $fh{$_->organism}; print $fh $_->name, $orf->name; } $orf->orthogroup;
	}
    }
    
    #

    return $self;
}

#####################################
# subroutines : simple queries / accessors 
#####################################

=head2 n50(-cutoff => 95) 

    Use cutoff param to obtain N95 etc. N50 is default. 

=cut 

sub n50 {
  my $self = shift;        
    my $args = {@_};

    $args->{'-cutoff'} = 50 unless exists $args->{'-cutoff'};

    my ($T, $RT) = (0,0); 
    map { $T += $_->length } $self->stream;
    foreach my $ctg ( sort { $b->length <=> $a->length } $self->stream ) {
	$RT += $ctg->length;
	return $ctg->length if $RT >= $T*($args->{'-cutoff'}/100);
    }
}

=head2 g50(-cutoff => 95) 

    Calculate the N50 (or 95 etc) equivalent for number of genes on a contig. 

=cut 

sub g50 {
    my $self = shift;        
    my $args = {@_};

    $args->{'-cutoff'} = 50 unless exists $args->{'-cutoff'};

    my @lengths = map { scalar(grep {$_->assign =~ /REAL|PSEUDO|HSP|NOVEL|HYPO|RNA/} $_->stream) } $self->stream;

    my ($T, $RT) = (0,0);     
    map { $T +=  $_} @lengths;

    foreach my $genes ( sort { $b <=> $a } @lengths ) {
	$RT += $genes;
	return $genes if $RT >= $T*($args->{'-cutoff'}/100);
    }
}

=head2 ontology(/path/to/GeneOntology.sto)

    Gene Ontology organized by Ancestral Gene. 

    The file should be of form: 
    % = {
    Anc_x.y => {
        P => ['arbitrary', 'string', 'values', ...],
        F => ['string', 'values'],
        C => ['string', 'values'],
    }
    Anc_x.z => {
        P => ['string', 'values'],
        F => ['string', 'values'],
        C => ['string', 'values'],
    }
    .....
    }

    The file is loaded immediatley and data acessible to ORF objects.

=cut 

sub ontology {
    my $self = shift;        
    my $args = {@_};
    if ( @_ ) { 
	return undef unless my $go = retrieve(shift);
	$self->{ONTOLOGY} = $go;
    }
    return $self->{ONTOLOGY};
}

=head2 find(-contig => i, -orf=> i)
       find(-internal => i)
       find(-scer => i)
       find('Spar_03.209')
 
    Returns the corresponding Contig or ORF object. 

=cut

sub find {
	my $self = shift;
	my $args = {@_};

	my $org = $self->organism;
	my ($regex) = map { qr/$_/ } ( $org.'_(\d+)\.(\d+)' );

	my ($key,$val);
	if (scalar(@_) == 1 ) {
	    my @decomp = &Annotation::Orf::_decompose_gene_name( $_[0] ) unless $_[0] =~ /^(\d+)$/;

	    if ( $_[0] =~ $regex ) { # a gene name in this organism 
		$args->{'-contig'} = $1;
		$args->{'-orf'} = $2;
	    } elsif ( $_[0] =~ /^(\d+)$/ ) { # a contig in this oirgnaism 
		$args->{'-contig'} = $1;
	    } elsif ( @decomp ) {
		$key = $decomp[0];
		$val = $_[0];
	    }
	    $key='SCER' if $key eq 'SGD';

	} elsif ( scalar( grep {!/contig|orf|intern/} keys %{$args} ) == 1 ) { # we have _supplied_ key 
	    ($key) = ( map {s/\-//; uc($_)}  grep {!/contig|orf|intern/} keys %{$args} );
	    $self->throw($key) unless exists $HOMOLOGY{$key} || $key eq 'GENE';
	    $val = $args->{ '-'.lc($key) };
	}

	$self->throw("Must specify an orf, contig, internal id, homology! ".join(' ',@_) ) unless 
	    $args->{'-contig'} =~ /\d/ || $args->{'-internal'} =~ /\d/ || ($key && $val);

	# homology mode 

	if ( $key ) {
	    return grep { $_->data($key) eq $val } $self->orfs;
 	} 

	# object id search 

	foreach my $c ($self->stream) {
	    if ($args->{'-internal'}) { # internal id search 
		return $c if $c->_internal_id == $args->{'-internal'};
	    } elsif (! $args->{'-orf'}) {# contig search 
		return $c if $c->id == $args->{'-contig'};
		next;
	    } else { # orf search 
		next unless $c->id == $args->{'-contig'};
	    }
                
	    foreach my $o ($c->stream) {
		if ($args->{'-internal'}) {
		    return $o if $o->_internal_id == $args->{'-internal'};
		} elsif ($args->{'-orf'}) {
		    return $o if $o->id == $args->{'-orf'};
		} else {$self->throw;}
	    }
	}
	
	$self->warn("Could not find: @_");
	return $self;	
}

=head2 load(-rename => 1, -degeneracy => .2, -minimum => 100)

    Load a fasta file of contigs into a genome object. 

    -rename tells software to rename contigs with a sane naming system. 
       => use it unless you have a good reason not to. 

    -degeneracy throws out contigs where the fraction of non standard 
    (ie degenerate) bases is > degeneracy. N is not treated as degenerate
    since it is handled well by code. All other letters are treated as 
    unknown (and NOT capable of generating a STOP codon) when annotating. 

    -minimum toss contigs below thi size. 

=cut

sub load {    
    my $self = shift;
    my $args = {@_};
    
    # some basic checking 
    
    $args->{'-debug'} = undef unless exists $args->{'-debug'};
    $args->{'-rename'} = 1 unless exists $args->{'-rename'};
    $args->{'-degeneracy'} = .2 unless exists $args->{'-degeneracy'};	
    $args->{'-minimum'} = 100 unless exists $args->{'-minimum'};	
    $args->{'-verbose'} = 0 unless  $args->{'-verbose'};
    my $flag = '***';

    # all ok ?
    
    $self->throw("Bad sequence file -- $args->{'-file'}")
	unless $self->_valid_file(@_);

    # 
    chomp(my $date = `date`);
    $self->history({
	DATE => $date,
	FILE => $args->{'-file'},
	SCRIPT => $0,
	STAMP => $TIME,
	STEP => '__START__'
	});

    #

    local $/ = ">";
    open(IN, $args->{'-file'});
    <IN>; # burn one

  SEQ: while (my $nl = <IN>) {
      chomp($nl);
      my ($id, @r) = split/\n/, $nl;
      my $seq = join('', @r);      
      my $length = length($seq);

      my $kmercov = -1; 
      ($id,$kmercov) = split/\s+/, $id; # common format >scaffold67 45.3 
	  
      unless ($length >= $args->{'-minimum'}) {
	  my $string = join("\t",  ("Length", $id, $length, $flag) );
	  $self->log(-text => $string);
	  print {STDERR} $string if $args->{'-verbose'};
	  next SEQ;
      }
      
      # check composition
      
      my %comp = (GAP => 0, X => 0, A => 0, T => 0, C => 0, G => 0);
      $comp{GAP} = (($seq =~ s/\s//g) || 0);
      $comp{X} = (($seq =~ s/X/N/g) || 0);
      map {$comp{$_}++} split//, $seq;
      
      my $err = join('; ', map {$_."=".$comp{$_}} keys %comp); 
      map { $self->throw("Non-alpha\t$id\t$length\t$err") if /[^A-Z]/ } 
      (grep {!/GAP/} keys %comp);
      
      unless ($comp{N}/$length < $args->{'-degeneracy'}) {
	  my $string = join("\t", "Degenerate", $id, $length, $flag);
	  $self->log(-text => $string);
	  print {STDERR} $string if $args->{'-verbose'};
	  next SEQ;
      }  
  
      # rename contig?
      
      if ($args->{'-rename'}) {		
	  $id = $args->{'-rename'};
	  $args->{'-rename'}++;
      } else { $self->throw("Numbers only in contig name.\nName/$id/$length/$err") if $id =~ /\D/; } 

      # truncate sequence for debugging purposes ?
      
      if ( $args->{'-debug'} ) {
	  my $start = (ref($args->{'-debug'}) eq 'ARRAY' ? $args->{'-debug'}->[0] : 1);
	  my $stop = (ref($args->{'-debug'}) eq 'ARRAY' ? $args->{'-debug'}->[1] : $args->{'-debug'});
	  $self->throw unless $stop >= $start;
	  $self->throw unless ($start >= 1 && $start <= $length); 
	  $self->throw unless ($stop >= 1 && $stop <= $length);
	  $seq = substr($seq, $start-1, ($stop-$start)+1);
      }

      # make data structure 
      
      my $contig = Annotation::Contig
          ->new(
	  SEQUENCE => $seq,
	  SCAFFOLD => $id,
	  KMERCOVERAGE => $kmercov,  # DEVIN 20110303 
	  ID => $id
	  );
      
      # add GAP features on the fly 
      # add a 'JOIN' feature. add to both strands.  
      
      my $offset=0;
      foreach my $match ($seq =~ /N{1,}/g) {	      
	  my $coord=index($seq, $match, $offset);
	  # print $id, $offset, $coord, length($match);
	  my $j1 = Annotation::Feature
	      ->new(
	      COORD => $coord+1,
	      SCORE => $INFINITY,
	      FEATURE => 'NNNNN_1',
	      STRAND => 0,
	      UP => $contig					
	      );
	  my $j2 = Annotation::Feature
	      ->new(
	      COORD => $coord+length($match),
	      SCORE => $INFINITY,
	      FEATURE => 'NNNNN_2',
	      STRAND => 0,
	      UP => $contig
	      );
	  
	  $j2->link($j1);
	  $j1->link($j2);
	  $contig->_add_feature(-object => $j1);	
	  $contig->_add_feature(-object => $j2);

	  $offset = $j2->coord;
      }

      print STDERR "Adding", $contig->id, $contig->length,$contig->kmercoverage, $err
	  if $args->{'-verbose'}; 
      $self->add(-object => $contig);       
  }
    
    # did we get anything for our trouble?
    
    $self->throw("\nNo valid contigs loaded.\n") 
	unless $self->down;

    # finalize naming by contig size 

    if ( $args->{'-rename'} ) {
	my @contigs = sort { $b->length <=> $a->length } $self->stream;
	foreach my $i (0..$#contigs) {
	    $contigs[$i]->id( $i+1 );
	}    
    }

    return $self;
}

=head2 bound('Species')

Returns all species names for the bound genomes
or, with a species name argument, returns the object.

=cut 

sub bound {
    my $self = shift;
    my $arg = shift;
    if ( $arg ) {	
	return $self->{'BOUND'}->{$arg};
    } else {
	return keys %{$self->{'BOUND'}};
    }
}

=head2 iterate

    Cycle through all bound genomes.

=cut

sub iterate { return($_[0], (map { $_[0]->bound($_) } $_[0]->bound));}

=head2 _bind()

    Bind genomes into a gound genome object. 
    Called only by  syntenic_orthologs.

=cut
    
sub _bind {
    my $self = shift;

    $self->throw if $self->bound;
    map { $self->throw unless $self->isa(ref($_)); } @_;
    map { $self->throw if $_ eq $self } @_;
    map { $self->throw if $_->bound } @_;

    map { $self->{'BOUND'}->{$_->organism}=$_ } @_;
    $self->throw unless $self->bound;

    return $self;
}

#########################################
# subroutines : Core methods 
#########################################

=head2 syntenic_orthologs(-object => [genome1, genome2..genomeN],
    -greedy => 1, -verbose => 0|1, -debug=>Anc_23.132)
    
    Group genes from different species into orthology groups based on 
    homology and synteny. 

=cut 

sub syntenic_orthologs {
    my $self = shift;
    my $args = {@_};

    # set defaults 

    $args->{'-debug'} = undef unless exists  $args->{'-debug'};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};
    $args->{'-greedy'} = 1 unless exists $args->{'-greedy'};
    $args->{'-homology'} = 'YGOB' unless exists $args->{'-homology'};
    
    # check params 

    $self->throw("YGOB only implemented: $args->{'-homology'}") unless $args->{'-homology'} eq 'YGOB';
    $self->throw unless exists $args->{'-object'} && ref($args->{'-object'}) eq 'ARRAY';
    map { $self->isa(ref($_)) } @{$args->{'-object'}};
    $self->_bind( @{$args->{'-object'}} );

    my @genomes = ($self,  @{$args->{'-object'}});

    # set hardcoded vars 
    # log10(Eval) / SyntenyScr / DPlimit / LengthVar 

    my @thresholds = (
	[-400,50,1,0.05],
	[-300,40,2,0.1],
	[-200,30,3,0.15],
	[-100,20,4,0.2],
	[-50,20,4,0.2], # 
	[-10,10,5,0.25],
	[-5,10,5,0.25] # 
	);

    my $fh = STDOUT;
    my $fherr = STDERR;

    my %synteny;
    if ( -e 'syntenyhash.sto') {
	%synteny = %{retrieve('syntenyhash.sto')};
	#print scalar( keys %synteny );
    }

    ##################################################
    # ALGORITHM: 
    # We do banded syntenic ortholog search, resolving the 
    # easiest cases first and then using these as extra 
    # constrains to resolve harder cases. This is facilitated 
    # by switching back and forth between an "homology first,
    # synteny second" phase and a "synteny first, homology second"
    # phase. Once both phases are complete, we decrease stringency
    # and iterate phases. 
    ##################################################
    
    # Select thresholds for this iteration .... 
    
  THRESHOLDS: for my $thresh ( 0..$#thresholds ) {
      my ( $ygobmax, $synmin, $dpoglimit, $maxlenvar ) = @{ $thresholds[$thresh] };
      my $syndrop = $synmin/2; 

    FMSIZE: for my $fmsize ( ((0..5),10,20,50,100,100000) ) {
	print {$fh} ">>>$thresh/$fmsize" if  $args->{'-verbose'} >= 2;
	
	##################################################
	# PHASE 1. GROUP BY HOMOLOGY FIRST THEN SYNTENY
	##################################################
	
	# make an index of all the genes in the game
	
	my %hash;
      INDEX: foreach my $g (@genomes) {
	  foreach my $o ( grep {! $_->ogid} $g->orfs(-noncoding => 1) ) {
	      my $group = ( $o->data( $args->{'-homology'} ) ? $o->data( $args->{'-homology'} ) : $o->gene );
	      my $regex = $HOMOLOGY{ $args->{'-homology'} };
	      push @{$hash{$group}{$g->organism}}, $o if $group =~ /$regex|\:/;
	  }
      }
	map { delete $hash{$_} unless scalar(keys %{$hash{$_}})==scalar(@genomes) } keys %hash;
	
	# iterate over homology groups, split into synteny groups and evaluate 
	
      FAMILY: foreach my $anc (sort grep {/\w/} keys %hash) {
	  next unless $anc eq $args->{'-debug'} || ! $args->{'-debug'};
	  print {$fherr} $thresh,$fmsize,$anc, (map { "$_:".scalar(@{$hash{$anc}{$_}}) } keys %{$hash{$anc}})
	      if $args->{'-verbose'} >= 2;
	  
	  # choose reference species. 
	  # we have grouped genes by homology. within each family we will 
	  # group by synteny. below we will check that group memebership 
	  # is both strong and unambiguous.
	  # map { print $_,@{$hash{$anc}{$_}} } keys %{$hash{$anc}};
	  
	  map { $self->throw if $_->ogid } map {@{$_}} values %{ $hash{$anc} };
	  map { next FAMILY unless $#{$hash{$anc}{$_}} <= $fmsize } keys %{$hash{$anc}};
	  my ($reference) = sort { $#{$hash{$anc}{$a}} <=> $#{$hash{$anc}{$b}} } keys %{$hash{$anc}};
	  
	  # compute synteny between all pairs of candidates and genes 
	  # all genes above the synteny drop are included in the 
	  # candidate array. we later check that successful additions are 
	  # unique. geneA not added to OG1 because syn50 but drop10 and geneB
	  # added because syn40 and drop 20. we do not want to make an OG with geneB. 
	  
	  my @cands = @{$hash{$anc}{$reference}}; # candidate OGs 
	  my %clusters = map { $_->name => [$_]} @cands;	  
	SYNTENY: foreach my $gene ( map { @{$hash{$anc}{$_}} } grep {!/$reference/} keys %{$hash{$anc}} ) {
	    foreach my $cand ( @cands ) {
		$synteny{$gene->name}{$cand->name} = $gene->synteny_conserved(-object => $cand)
		    unless exists $synteny{$gene->name}{$cand->name};
	    }
	    
	    # add all genes over syndrop 
	    
	    my $last;
	    foreach my $cluster ( sort { $synteny{$gene->name}{$b} <=> $synteny{$gene->name}{$a} } map {$_->name} @cands ) {	
		my $synscr = $synteny{$gene->name}{$cluster};
		my $delta = ($last ? $synscr - $last : 0);
		$last = $synscr;
		last if $delta > $syndrop;
		push @{ $clusters{$cluster} }, $gene;
	    }
	}
	  # map { print $_,@{$clusters{$_}} } keys %clusters;
	  
	  # check each candidate cluster and see if we can make 
	  # high confidence 1:1 OGs
	  
	CANDIDATE: foreach my $cluster ( keys %clusters ) {
	    
	    # 1)) is the synteny data both compelling and unambiguous? 
	    
	    my @synset = grep { $_->organism eq $reference } @{ $clusters{$cluster} };
	    $self->throw( "$reference, $#set" ) unless $#synset == 0;
	    
	    foreach my $sp ( grep { $_ ne $reference } map {$_->organism} @genomes ) {
		next CANDIDATE unless my @synsort = 
		    sort { $synteny{$b->name}{$cluster} <=> $synteny{$a->name}{$cluster}  }
		grep { $_->organism eq $sp }  @{ $clusters{$cluster} };
		my @synscr = map { ($synteny{$_->name}{$cluster} || 0) } @synsort;
		next CANDIDATE unless ($#synsort == 0 || ($synscr[0] >= $synmin && ($synscr[0]-$synscr[1] >= $syndrop)) );
		push @synset, $synsort[0];
	    }
	    $self->throw unless $#synset == $#genomes;
	    
	    # 2)) now check homology data. does it meet criteria? 
	    # we also exclude variable length loci 
	    
	    my @og = grep { $_->logscore('ygob') <= $ygobmax } @synset;
	    next CANDIDATE unless $#og == $#genomes;
	    my ($m2,$s2) = _calcMeanSD(map {$_->length} @synset);
	    next CANDIDATE if $s2/$m2 >= $maxlenvar;
	    
	    # 3)) make the OG 	   

	    my ($o) = grep { $_->organism eq $self->organism } @og;
	    $o->_define_orthogroup( 
		-object => [grep { $_ ne $o } @og],
		);
	    nstore(\%synteny, 'syntenyhash.sto') if $o->ogid%100==0;
	    
	    if ( $args->{'-verbose'} >= 2 ) {
		my ($m0,$s0) = _calcMeanSD( map {$_->logscore('ygob')} $o->_orthogroup );
		my ($m1,$s1) = _calcMeanSD( map { $synteny{$_->name}{$cluster} } grep {$_->name ne $cluster} $o->_orthogroup  );
		print {$fh} $thresh, $fmsize, $anc, $cluster, substr($o->gene, 0, 7), 
		$ygobmax, $m0,$s0, 
		"$synmin ($syndrop)", $m1,$s1, 
		$maxlenvar, $m2,$s2, $o->ogid; 
	    }

	    # 4)) remove used genes from %hash so we can iterate. 
	    # if any species goes to 0, we remove the Anc since no OG can form. 
	    
	  REMOVE: foreach my $rm ( $o->_orthogroup ) {
	      for (my $i = $#{$hash{$anc}{$rm->organism}}; $i >= 0; $i-- ) {
		  splice(@{$hash{$anc}{$rm->organism}}, $i, 1) and next REMOVE if $hash{$anc}{$rm->organism}->[$i] eq $rm;
	      }
	      $rm->output(-fh => $fherr);
	      $self->throw;
	  }
	    map { delete $hash{$anc} and next FAMILY unless $#{$hash{$anc}{$_}} >= 0} map { $_->organism } @genomes;
	}	  
      } # YGOB FAMILY 
	
	##################################################
	# PHASE 2. GROUP BY SYNTENY FIRST THEN HOMOLOGY 
	##################################################
	
	# we have examined all families at the current thresholds. 
	# we now check to see whether there are any cases that can  
	# be resolved in a high confidence manner by DP in the inter-OG
	# regions (controlled by dpoglimit);
      
      INTER: foreach my $o ( grep { $_->ogid } $self->orfs(-noncoding => 1) ) {
	  next INTER unless my $l = $o->neighbour(-direction => 'left', -orthogroup => 1);
	  
	  my %genes = %{ $o->bracket(-object => $l, -trna => 1, -clean => 1) };
	  next INTER unless scalar(keys %genes)==scalar(@genomes);	  
	  map { next INTER unless scalar(@{$_}) <= $dpoglimit } values %genes;
	  
	  $self->throw unless my ($dp,$og) = 
	      $self->dpalign( -hash => \%genes, -match => 'local', -trna => 100 );
	  
	  foreach my $cand ( @{$og} ) {
	      # impose length criteria or examine PW scores ? 
	      my ($o) = grep { $_->organism eq $self->organism } @{$cand};
	      my @og = grep {$_ ne $o} @{$cand};
	      $o->_define_orthogroup( @og );	    
	      # 
	      if ( $args->{'-verbose'} >= 2 ) {
		  my ($m0,$s0) = _calcMeanSD( map {$_->logscore('ygob')} $o->_orthogroup );
		  my ($m1,$s1) = _calcMeanSD( map { $synteny{$_->name}{$o} } $o->orthogroup  );
		  my ($m2,$s2) = _calcMeanSD(map {$_->length} $o->_orthogroup );
		  my ($sgd,$gob) = $o->identify; 
		  print {$fh} $thresh, "$fmsize", $gob, $o->name, substr($o->gene, 0, 7), 
		  "DP/$dpoglimit", $m0,$s0, 
		  "DP/$dpoglimit", $m1,$s1, 
		  "DP/$dpoglimit", $m2,$s2, $o->ogid; 
	      }
	  }
      } # INTER 	

    } # FMSIZE 
      print {$fh} ">$thresh", scalar($self->orthogroups) if $args->{'-verbose'};
  } # THRESHHOLDS

    print {$fh} ">Final", scalar($self->orthogroups) if $args->{'-verbose'};    
    return $self unless $args->{'-greedy'};

  GREEDY: 

    ##################################################
    ##################################################
    # Algorithm proper is over but there might be some OGs left 
    # if (1) we accept non-bounded cases, (2) do not enforce 
    # limits on #'s of genes in intergenics and (3) allow
    # repeats to form OGs. The principal advantage of the
    # last is that it might create smaller bounded regions for 
    # subsequent OG searches. 
    my %success=(
	'NOLIMIT' => 0, 
	'UNBOUND' => 0, 
	'REPEAT' => 0
	);
    ##################################################

    my $key;
    my $iteration=0;
    my $old_count = 0;
    my $og_count = scalar( $self->orthogroups );
    until ( $old_count == $og_count ) {
	
	##################################################
	# 1) relax limit on # of genes in bounded region 
	##################################################

	$key = 'NOLIMIT';	
	foreach my $o ( grep { $_->ogid } $self->orfs(-noncoding => 1) ) {
	    next unless my $l = $o->neighbour(-direction => 'left', -orthogroup => 1);
	    my %genes = %{ $o->bracket(-object => $l, -trna => 1, -clean => 1) };
	    next unless scalar(keys %genes)==scalar(@genomes);
	    
	    $self->throw unless my ($dp,$og) = 
		$self->dpalign( -hash => \%genes, -match => 'local', -trna => 100 );
	    
	    foreach my $cand ( @{$og} ) {
		my ($o) = grep { $_->organism eq $self->organism } @{$cand};
		my @og = grep {$_ ne $o} @{$cand};
		$o->_define_orthogroup( @og );
		
		if ( $args->{'-verbose'} >= 2 ) {
		    my ($m0,$s0) = _calcMeanSD( map {$_->logscore('ygob')} $o->_orthogroup );
		    my ($m1,$s1) = _calcMeanSD( map { $synteny{$_->name}{$o} } $o->orthogroup  );
		    my ($m2,$s2) = _calcMeanSD(map {$_->length} $o->_orthogroup );
		    my ($sgd,$gob) = $o->identify; 
		    print {$fh} ++$iteration, $key.':'.++$success{$key}, $gob, $o->name, substr($o->gene, 0, 7), 
		    "DP", $m0,$s0, 
		    "DP", $m1,$s1, 
		    "DP", $m2,$s2, $o->ogid; 
		}		
	    }
	}
	print {$fh} ">$key", scalar($self->orthogroups) if $args->{'-verbose'};

	##################################################
	# 2) look for orthologs that are at unbounded locations (contig ends) 
	##################################################

	$key = 'UNBOUND';    	
	my %seen;
	foreach my $o ( grep { $_->assign !~ /GAP|FEATURE|REPEAT/ } grep { ! $_->ogid } $self->orfs(-noncoding => 1) ) {
	    my %genes = %{ $o->gather(-trna => 1, -clean => 1) };
	    next unless scalar(keys %genes)==scalar(@genomes);
	    
	    # required since gather() grabs genes to left and right. 
	    map { next if exists $seen{$o->name}; $seen{$o->name}++  } $genes{uc($o->organism)};
	    
	    $self->throw unless my ($dp,$og) = 
		$self->dpalign( -hash => \%genes, -match => 'local', -trna => 100 );
	    
	    foreach my $cand ( @{$og} ) {
		my ($o) = grep { $_->organism eq $self->organism } @{$cand};
		my @og = grep {$_ ne $o} @{$cand};
		$o->_define_orthogroup( @og );

		if ( $args->{'-verbose'} >= 2 ) {
		    my ($m0,$s0) = _calcMeanSD( map {$_->logscore('ygob')} $o->_orthogroup );
		    my ($m1,$s1) = _calcMeanSD( map { $synteny{$_->name}{$o} } $o->orthogroup  );
		    my ($m2,$s2) = _calcMeanSD(map {$_->length} $o->_orthogroup );
		    my ($sgd,$gob) = $o->identify; 
		    print {$fh} ++$iteration, $key.':'.++$success{$key}, $gob, $o->name, substr($o->gene, 0, 7), 
		    "DP", $m0,$s0, 
		    "DP", $m1,$s1, 
		    "DP", $m2,$s2, $o->ogid; 
		}
	    }
	  }
	print {$fh} ">$key", scalar($self->orthogroups) if $args->{'-verbose'};
	
	##################################################
	# 3) allow repeats 
	##################################################	

	$key = 'REPEAT';	
	foreach my $o ( grep { $_->ogid } $self->orfs(-noncoding => 1) ) {
	    next unless my $l = $o->neighbour(-direction => 'left', -orthogroup => 1);
	    my %genes = %{ $o->bracket(-object => $l, -trna => 1, -remove => [qw(GAP FEATURE)]) };	    
	    next unless scalar(keys %genes)==scalar(@genomes);
	    
	    $self->throw unless my ($dp,$og) = 
		$self->dpalign( -hash => \%genes, -match => 'local', -trna => 100 );

	    foreach my $cand ( @{$og} ) {
		my ($o) = grep { $_->organism eq $self->organism } @{$cand};
		my @og = grep {$_ ne $o} @{$cand};
		$o->_define_orthogroup( @og );
		
		if ( $args->{'-verbose'} >= 2 ) {
		    my ($m0,$s0) = _calcMeanSD( map {$_->logscore('ygob')} $o->_orthogroup );
		    my ($m1,$s1) = _calcMeanSD( map { $synteny{$_->name}{$o} } $o->orthogroup  );
		    my ($m2,$s2) = _calcMeanSD(map {$_->length} $o->_orthogroup );
		    print {$fh} ++$iteration, $key.':'.++$success{$key}, $gob, $o->name, substr($o->gene, 0, 7), 
		    "DP", $m0,$s0, 
		    "DP", $m1,$s1, 
		    "DP", $m2,$s2, $o->ogid; 
		}	
	    }
	}	
	
	$old_count = $og_count;
	$og_count = scalar( $self->orthogroups );
    }
    print {$fh} ">$key", scalar($self->orthogroups) if $args->{'-verbose'};
    
    $self->quality( -nulldist => 1, -sample => 1000 );

    return $self;
}

=head2 dpalign(-hash => \%hash, -global => 1, -order => [Scer, Spar, Smik etc],
    -match => 'local|global|ygob|sgd|ect', -mismatch => -INF, -inversion => -100, 
    -gap => 0, -trna => 100)

    Accepts an hash of arrays (Speices => [G1, G2, G3]) and returns
    an array of hashes where each element is a prospective orthogroup:
    [ {Scer=>,Sbay=>..}, {}, {} ] based on DP determined gene order. 

    -order : species are added to the alignment in this order (not tree).
    -global : do NW or SW alignemnt. Currentyl NW only. 
    -match : choose the scoring mechanism. By default we use exonerate local DP 
    but with -score can choose to use precompute dhomologies such as 'YGOB'. 
    This changes the scale of the scoring scheme. 
    -mismatch : -INF 
    -gap : gap penalty should be zero unless there is a very good reason.. 
    -inversion : we allow segments to be inverted but apply penalty to reverse 
    orientation. affine. -INF to disable (fwd aligned guaranteed better score). 
    -trna : supply a score to be used for matching tRNAs. 
    
    Inversion Penalty Scoring 
    Homologous genes on opposite strands are penalized (-inversion * 0.5). Thus, 
    with standard -inversion penalty of 100 (1), they are penalized 50 (0.5). 
    Simultaneously, segments that are inverted prior to aligning carry a base 
    penalty of -inversion * 0.75. Thus, at least 2 more genes must be inverted 
    reltaive to the reference than colinear ot favour inverting the intervening 
    region:  (-inversion * 0.5) * 2 > (-inversion * 0.75)
    
=cut 

sub dpalign {
    my $self = shift;
    my $args = {@_};

    # control DP algorithm 
    $args->{'-global'} = 1 unless exists $args->{'-global'};
    $args->{'-order'} = [ reverse(split/\s+/, (map { s/[\(\)\;\,]+/ /g; $_ } $self->tree)[0] ) ] 
	unless exists $args->{'-order'};
    $args->{'-order'} = [ map {uc($_)} @{$args->{'-order'}} ];
    $args->{'-reference'} = 0 unless exists $args->{'-reference'};

    # control scoring mechanisms
    $args->{'-match'} = 'local' unless exists $args->{'-match'};
    $args->{'-gap'} = 0 unless exists $args->{'-gap'}; # open + extend 
    $args->{'-mismatch'} = -$INFINITY unless exists $args->{'-mismatch'}; # matches are calculated below 
    $args->{'-inversion'} = -1 * ($args->{'-match'} =~ /local|global/i ? 100 : 1 ) 
	unless exists $args->{'-inversion'};
    $args->{'-trna'} = 100 unless exists $args->{'-trna'};

    # 
    $args->{'-verbose'}=0 unless exists $args->{'-verbose'};

    # check arguments 

    $self->throw("Smith Waterman not implemented!") unless $args->{'-global'};
    $self->throw unless ref($args->{'-hash'}) eq 'HASH'; 
    #map { print ".$_\t$#{$args->{'-hash'}->{$_}}" } keys %{ $args->{'-hash'} };
    #map { $self->throw($_) unless exists $args->{'-hash'}->{uc($_)} } ($self->organism, $self->bound);
    #map { return undef if $#{$args->{'-hash'}->{uc($_)}}==-1 }  ($self->organism, $self->bound);
    
    my $fh = STDOUT;
    #######################

    # make some vars. 
    # using "$ref => $_" below allows us to pass in fake orf arrays to 
    # align to. used for ohnolog detection 

    my %hash=%{$args->{'-hash'}};
    my $ref = shift( @{$args->{'-order'}});     
    my @ref = map { {$ref => $_, REF => $_} } @{ $hash{$ref} }; # uc($_->organism)
    my @keys = ($ref);

    # iterate over the order array 
    # and consider Fwd/Rev aligns for each new species 

    while ( my $add = shift @{$args->{'-order'}} ) {
	print "\n>>>$add" if $args->{'-verbose'};
 
	# establish vars 

	my @fwd = map { {uc($add) => $_, REF => $_} } @{ $hash{$add} }; # uc($_->organism)
	my @rev = reverse( @fwd );	
	my %dirs = (
	    FWD => { ARRAY => \@fwd, MATRIX => undef, SCORE => undef, FACTOR => -.5 },
	    REV => { ARRAY => \@rev, MATRIX => undef, SCORE => undef, FACTOR => .5 }
	    );
	delete $dirs{REV} if $args->{'-inversion'} == -$INFINITY;
	
	my $score_mat = _make_dp_scoring_matrix( 
	    \@ref, \@fwd, 
	    $args->{'-match'}, 
	    $args->{'-mismatch'}, 
	    $args->{'-trna'}
	    );
	
	# do the DP 
	
	foreach my $dir ( keys %dirs ) {    
	    my @add = @{$dirs{$dir}->{ARRAY}};
	    print $fh undef, undef,(map {substr($_->{REF}->data('SGD'),0,7)} @add) 
		if $args->{'-verbose'};
	    
	    # complete matrix collecting scores and traceback 
	    
	    my $dp_mat = _init_dp_matrix( \@ref, \@add, $args->{'-global'} );
	    foreach my $row ( 0 .. $#{$dp_mat} ) {
		foreach my $col (  0 .. $#{$dp_mat->[$row]} ) {
		    next if $row==0 && $col==0;
		    my $row_obj = $ref[ $row-1 ]->{REF} if $row>0;
		    my $col_obj = $add[ $col-1 ]->{REF} if $col>0;

		    # we penalize matches that are on the opposite strand 
		    # but there is no reward for being same strand 
		    
		    my $strand_score = ( 
			($row>0 && $col>0) 
			? $row_obj->strand *  $col_obj->strand * $dirs{$dir}->{FACTOR} * $args->{'-inversion'}
			: 0
			);			
		    $strand_score = 0 if $strand_score > 0;
		    
		    # conventional DP scoring choices 
		    
		    my %sort_hash;
		    $sort_hash{'DIAG'} = 
			($dp_mat->[$row-1]->[$col-1]->{'SCORE'} + 
			 $score_mat->{ $row_obj->name }->{ $col_obj->name } + 
			$strand_score) 
			if ($row>0 && $col>0);
		    $sort_hash{'UP'} =
			($dp_mat->[$row-1]->[$col]->{'SCORE'} + $args->{'-gap'}) if $row > 0;
		    $sort_hash{'LEFT'} = 
			($dp_mat->[$row]->[$col-1]->{'SCORE'} + $args->{'-gap'}) if $col > 0;
		    $sort_hash{'STOP'} = 0 if $args->{'-global'}==0;
		    
		    # 
		    
		    my ($best) = sort { $sort_hash{$b} <=> $sort_hash{$a} } keys %sort_hash;
		    $dp_mat->[$row]->[$col] = {SCORE => $sort_hash{$best}, PATH => $best};
		    $dirs{$dir}->{SCORE} = $sort_hash{$best};
		}
		
		if ( $args->{'-verbose'} ) {
		    foreach my $var ( 'SCORE', ($args->{'-verbose'} >=2 ? 'PATH' : ()) ) {
			print $fh $var, ($row>0 ? substr($ref[ $row-1 ]->{REF}->data('SGD'),0,7) : $dir), 
			(map {$_->{$var}} @{$dp_mat->[$row]});
		    }
		}
	    }
	    print $fh "" if $args->{'-verbose'};	    
	    $dirs{$dir}->{MATRIX} = $dp_mat;
	}
	
	# choose between fwd and reverse orientations by applying inversion penalty 
	
	$dirs{'REV'}->{SCORE} += $args->{'-inversion'}*(0.75);
	my ($bestDir,$other) = sort { $dirs{$b}->{SCORE} <=> $dirs{$a}->{SCORE} } keys %dirs;
	my $key = ($bestDir eq 'REV' ? '-' : '').uc($add);
	push @keys, $key;
	print $fh $bestDir, $dirs{$bestDir}->{SCORE}, "($dirs{$other}->{SCORE})" if $args->{'-verbose'};
	my $dp_mat = $dirs{$bestDir}->{MATRIX};
	my @add = @{$dirs{$bestDir}->{ARRAY}};
	
	# perform the traceback and redefine the ref so we can iterate 
	
	my @path;
	my ($i,$j)=($#{$dp_mat}, $#{$dp_mat->[0]}); 
	until ( $i == 0 && $j ==0 ) {
	    my $data = $dp_mat->[$i]->[$j];
	    unshift @path, [$i, $j,$data->{PATH}];
	    print @{$path[0]} if $args->{'-verbose'};
	    if ( $data->{PATH} eq 'UP' ) {
		$i--;
	    } elsif ( $data->{PATH} eq 'LEFT' ) {
		$j--;
	    } elsif ( $data->{PATH} eq 'DIAG' ) {
		$i--; 
		$j--;
	    } elsif ( $data->{PATH} eq 'INIT' ) {
		last;
	    } else { $self->throw("not implemented: $data->{PATH} $i,$j"); }	    
	}
	
	for my $i ( 0..$#path ) {	    
	    if ($path[$i]->[2] eq 'LEFT' ) {
		splice(@ref, $i, 0, undef);
	    } elsif ($path[$i]->[2] eq 'UP') {
		splice(@add, $i, 0, undef);
	    } 
	}

	if ( $args->{'-verbose'} ) {
	    print $fh map { $_->[2] } @path;
	    print $fh map { ( $_ ? $_->{REF}->name : $_) } @ref;
	    print $fh map { ( $_ ? $_->{REF}->name : $_) } @add;
	    print $fh map { ( $_ ? $_->{REF}->data('SGD') : $_) } @ref;
	    print $fh map { ( $_ ? $_->{REF}->data('SGD') : $_) } @add;
	}

	for my $i ( 0..$#path ) { ####### NNB WILL NOT WORK FOR TREE-BASED ALIGN #######
	    $ref[$i]->{$key}=$add[$i]->{$add} if $add[$i]->{$add};  # transfer new species data
	    my ($long) = sort { $b->length <=> $a->length } grep {$_} values %{$ref[$i]}; 
	    if ( $args->{'-reference'} ) {
		my ($ans) = grep { $_->organism eq $args->{'-reference'} } grep {$_} values %{$ref[$i]};
		$long = $ans if $ans;
	    }
	    $ref[$i]->{'REF'} = $long; # define a new 'REF'
	}
	print $fh join("\t", map { ( $_ ? $_->{REF}->data('SGD') : $_) } @ref)." *\n" if $args->{'-verbose'};
    }

    # ensure matrix propely defined before return 
    # all speceis all positions and no ref 

    my @og;
    foreach my $pos (@ref) {
	delete $pos->{REF};
	push @og, [values %{$pos}] if (grep { $pos->{$_} } @keys) == scalar(@keys);
	map { $pos->{$_}=undef unless exists $pos->{$_} } @keys;
    }

    # output? 
    
    if ( $args->{'-verbose'} ) {
	print $fh ">";
	foreach my $sp ( @keys ) {
	    print $fh ( $sp, map { 
		( $_->{$sp}
		  ? substr($_->{$sp}->data('GENE'),0,7).'_'.
		  ($_->{$sp}->length*($_->{$sp}->strand ? $_->{$sp}->strand : 1)*( $sp =~ /^\-/ ? -1 : 1 ) )
		  : (undef,undef) ) } @ref );
	}
    }
    
    return (wantarray ? (\@ref, \@og) : \@ref );
}

sub _init_dp_matrix {
    my ($r1,$r2,$score) = @_;
    die unless ref($r1) eq 'ARRAY' && ref($r2) eq 'ARRAY';
    return undef unless $#{$r1}>-2 && $#{$r2}>-2;

    my @mat;
    for my $i (0..scalar(@{$r1})) {
	for my $j (0..scalar(@{$r2})) {
	    $mat[$i]->[$j]={SCORE => 0, PATH => 'INIT'};
	}
    }

    return \@mat;
}


sub _make_dp_scoring_matrix {
    my ($r1,$r2,$score,$mm,$trna) = @_;
    die unless ref($r1) eq 'ARRAY' && ref($r2) eq 'ARRAY';
    return undef unless $#{$r1}>-1 && $#{$r2}>-1;
    die unless $score && $mm;

    my %mat;
    foreach my $i ( map {$_->{REF}} @{$r1} ) {
	foreach my $j ( map {$_->{REF}} @{$r2} ) {
	    if ( $i->assign eq 'TRNA' || $j->assign eq 'TRNA' ) {
		$mat{ $i->name }{ $j->name } = ( $trna && $i->gene eq $j->gene ? $trna : $mm );
	    } elsif ( $i->coding && $j->coding && $score =~ /local|global/i ) {
		my $ex = $i->exonerate2(-object => $j, -homolog => undef, -return => 'score', -model => lc($score) );
		$mat{ $i->name }{ $j->name } = ( $ex > 0 ? $ex : $mm );
	    } else {
		$mat{ $i->name }{ $j->name } = 
		    ( $i->data($score) && $i->data($score) eq $j->data($score) ? 1 : $mm );
		#print STDOUT $i->name, $j->name, $i->data($score) eq $j->data($score), $mat{ $i->name }{ $j->name };
	    }
	    $mat{ $j->name }{ $i->name } = $mat{ $i->name }{ $j->name };
	}
    }    

    return \%mat;
}

=head2 ygob()

    Method to identify the most likely YGOB pillar for all loci.
    Uses information from multiple genes and YGOB loci.
    
    Incomplete. Not available. Use $orf->pillar();

=cut 

sub ygob {
    my $self = shift;
    my $args = {@_};
    $self->throw;

    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};

    my $key = 'PILLAR';
    
    foreach my $c ( $self->stream ) {
	
	# 1. first pass: gather and order YGOB pillar candidates 
	
	my @order;
	foreach my $o ( $c->orfs ) { 	    
	    my $hash = $o->lookup( -query => $o->gene ) if $o->gene;
	    my %count = ( $hash->{ID} => $hash) if $hash; # seed with 'GENE' 
	    
	    # iterate ..
	    foreach my $sp ( grep {!/MITO|SPOM|SGD/} keys %HOMOLOGY ) { # allow YGOB 
		next unless my $hit = $o->data( $sp );
		#next unless $o->data( '_'.$sp) <= 1e-5;
		next unless my $loss = 
		    $o->synteny( -lod => 1, -distance => 3, -restrict => [$sp] );
		my $hash = $o->lookup( -query => $hit );
		if ( $hash && $hash->{'ID'} ) {
		    $count{ $hash->{'ID'} } = $hash unless exists $count{ $hash->{'ID'} };	
		    $count{ $hash->{'ID'} }->{'LOSS'} += $loss;
		}
	    }
	    my ($first,@others) = sort { $b->{LOSS} <=> $a->{LOSS} } values %count;
	    $first->{DELTA} = sprintf("%.3f", $first->{LOSS} - (@others ? $others[0]->{LOSS} : 0));
	    push @order, [$o, $first, @others];
	    
	    $o->data( $key => $first->{'ID'} );
	    $o->data( '_'.$key => $first->{'DELTA'} );
	    
	    print {STDERR} $o->name, $o->data( $key ), $o->data( '_'.$key )
		if $args->{'-verbose'} >= 1;
	}

	next;
	
	# this appears unnecessary as long as we dope in the 'GENE' 
	# hits that proposed by impute. also means that we can call it as an 
	# ORF method. 
	####################################################################################
	# 2. second pass make decisions. use good pillars to frame ambiguous ones. 
	
	my ($left,$right);
	for my $i ( 0..$#order ) {
	    my ($o,$first,@others)= @{$order[$i]};

	    if ( $first->{DELTA} >20 ) {
		#print $o->name, $first->{DELTA}, $first->{ID}, @{$first->{SGD}};
		$left = $i;
		next;
	    } elsif ( @others ) {
		$right = $i+1;
		$right++ until ( $right >= $#order || $order[$right]->[1]->{DELTA} > 20 );
		$right='' if $right == $#order;
	    } else { next; }	   

	    # compute neighbour compatibility scores for all proposed locatoins 

	    foreach my $pillar (grep {defined} ( $first,@others ) ) {  # get score by summing over all species 
		foreach my $sp ( grep { exists $pillar->{$_} && defined $pillar->{$_} } grep {!/MITO|YGOB/} keys %HOMOLOGY ) {
		    
		    my @biscr; # 2 possibilites for orientation so we take best 
		    foreach my $gene ( grep {defined} @{ $pillar->{$sp} } ) {			
			
			my $biscr=0; # two possible orientatinos again (for each of left + right..)
			if ( my $l = $order[$left]->[1]->{$sp} ) {
			    my ($hi,$lo) = sort {$b <=> $a} grep {defined} map { $o->_loss($gene, $_, ($i-$left)) }  (grep {defined} @{$l});
			    $biscr+= $hi;
			}		
			if ( my $r = $order[$right]->[1]->{$sp} ) {
			    # left and right pillars should be locked wrt each other but too much hassle.. 
			    my ($hi,$lo) = sort {$b <=> $a} grep {defined} map { $o->_loss($gene, $_, ($right-$i)) }  (grep {defined} @{$r});
			    $biscr+= $hi;
			}
			push @biscr, $biscr;
		    }
		    
		    # choose best orientation for this species and sum ...
		    my ($hi,$lo) = sort { $b <=> $a } @biscr;
		    $pillar->{'_SCR'} += $hi;
		} #HOMOLOGY 
	    }

	    # 
	    
	    print ">".$o->name;
	    map { print $_->{_SCR},$_->{ID},$_->{LOSS},@{$_->{SGD}} } sort {$b->{_SCR} <=> $a->{_SCR} } ($first,@others);
	    $o->oliver;

	}
    }

    # 

    return $self;
}

=head2 subtelomeres( -window => 20, -max => 50, -sd => -1.5 ) 
    
    Identify subtelomeric boundaries based on syntenic decay
    at contig ends and label all genes as LEFT_SUBTEL,
    INTERNAL,RIGHT_SUBTEL. Uses BLAST homolgy information.
    
    Briefly, we compare contig terminii to a sample of 
    contig-internal -window length regions across genome.
    We extend putative subtelomeric regions inward from terminii
    until synteny within -sd of internal value or upto -max genes.
    -max also used to define chr internal regions so should be set 
    moderately high. -sd => {-1.5 .. -2.5} usually places boundary 
    within 2-3 genes of true boundary (as inferred from YGOB). 

    As justification for the present approach, a plot of syn scores 
    along the chr ( eg plot( cumsum( x$synlogodds ) ) ) shows clear 
    plateaus at start and end. Typically 20ish genes with syn ~20 
    whereas internal genes have ~200;

  NB: 
    1. We do not currently use copy-number profiles though these might
    be used to develop a set of telomere-positive regions that would
    allow learning of subtel synteny parameters that would allow 
    use of a machine-learning classifier (with internal regions 
    params learned from regions chosen by terminal-dist).
    2. We do not controld the total number of subtelomeres in the 
    genome. In the event that the edges of all contigs or all of 
    small contigs tend to exists in low synteny arrested, they will 
    be labelled as SUBTEL.

=cut 

sub subtelomeres {
    my $self = shift;
    my $args = {@_};

    $args->{'-window'} = 20 unless exists $args->{'-window'};
    $args->{'-sd'} = -1.5 unless exists $args->{'-sd'};
    $args->{'-max'} = 50 unless exists $args->{'-max'};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};

    # define vars and initiailize 
    
    my $fh = *STDOUT;
    my $key = 'SUBTEL';
    my $syndist = 5;
    
    map { $_->data( $key => 0 ) } $self->orfs;
    map { $_->data( '_'.$key => 'INTERNAL' ) } $self->orfs;

    # select chr internal regions and compue synteny in windows

    my ($old,@internal);
    foreach my $o ( grep { $_->telomere >= 50 } $self->orfs ) {
	next if $old && $o->distance( -object => $old) < $args->{'-window'}/2; # allow some overlap..
	$old = $o;
	my @scores = map { $_->synteny(-logodds => 1, -distance => $syndist) } grep { $_->rank <= 1 }
	$o->context(-distance => int( $args->{'-window'}/2 ), -self => 1, -variant => -1);
	push @internal, [ _calcMeanSD( @scores ) ];
    }    

    # compute simple stats 
    
    my @mean = map { $_->[0] } @internal;
    my ($m,$sd,$n) = _calcMeanSD( @mean );
    $self->throw unless $n >= 20;
    $self->throw unless $sd > 0;

    # exmaine each chr arms, moving inwards from terminii 
    
    foreach my $chr ( $self->stream ) {
	$self->warn("Skipping") and next if $chr->length < $syndist;

	# pick subtel genes 
	
	my $max = ( scalar( $chr->stream ) < $args->{'-max'} ? scalar( $chr->stream ) : $args->{'-max'} );	
	my @orfs = sort {$a->id <=> $b->id} $chr->stream;
	my %arms = (
	    'LEFT' => [ @orfs[0..$max] ],
	    'RIGHT' => [ (reverse(@orfs))[0..$max] ]
	    );
	
	# move in from terminii + compute window scores 

	foreach my $dir ( keys %arms ) {
	    my $runtot;
	    foreach my $o ( @{ $arms{$dir} } ) {
		my @scores = map { $_->synteny(-logodds => 1, -distance => $syndist) } grep { $_->rank <= 1 }
		$o->context(-distance => int( $args->{'-window'}/2 ), -self => 1, -variant => -1);
		$self->throw unless @scores;
		#next unless scalar(@scores) >= $args->{'-window'}/2;
		my ($t_m, $t_sd, $t_n) = _calcMeanSD( @scores );
		my $zscr = ($t_m-$m)/$sd;
		$o->data( $key => $zscr );
		$o->data( '_'.$key => ( $zscr < $args->{'-sd'} ? $key.'_'.$dir : 'INTERNAL') );
		# 
		$runtot += $zscr;
		print {$fh} $o->data( '_'.$key ), $chr->id, $o->id, 
		$o->telomere+0, $o->telomere(-bp => 1),
		sprintf( "%.1f", $o->synteny(-logodds => 1, -distance => $syndist)), $o->data('SGD'),
		($m,$sd,$n), ($t_m, $t_sd, $t_n), $o->data( $key ), $runtot if $args->{'-verbose'} ;		
	    }

	    # try to call acutal telomeres? 
	    # not at the moment...

	    my $subtel=0;
	    if ( $subtel < -$INFINITY ) {
		my $subtel = Annotation::Orf
		    ->new(
		    START => ($dir eq 'left' ? 1 : $chr->length),
		    STOP => ($dir eq 'left' ? 1 : $chr->length),
		    STRAND => 0
		    );
		$chr->add( -object => $subtel );		
	    }
	}
    }

    return $self;
}


=head2 _synteny( -nulldist => 2|1, -distance => 20, -limit => 50 )
    
    Build log-odds scoring table for synteny metric. Consider 
    +/- distance genes around caller and compute probabilities for
    homologs that are within limit of caller homolog. Larger 
    treated as non-syntenic. Two alternative randomizations are 
    available: (1) genes singly randomly relocated in genome. This
    approximates mislabelling by BLAST or transposition. The 
    synteny scores relative to neighbouring genes are correlated. 
    (2) randomly permute the entire genome. This approximates the
    true random chance of finding neighbouring genes in proximity
    in another genome assuming no common ancestry. In practice there
    seems to be little difference between (1) and (2).

    NB: 
    (1) We do not correct for number of positions actually contributing to 
    score (ie neighbours x homologs) which may vary depending on genomic
    location (ie contig end) and gene conservation. 
    (2) At present we return the score but do not compute an E-value. 
    
=cut 

sub _synteny {
    my $self = shift;
    my $args = {@_};
    
    # modes
    $args->{'-object'} = undef unless exists $args->{'-object'};
    $args->{'-nulldist'} = undef unless exists $args->{'-nulldist'};    
    # params 
    $args->{'-distance'} = 20 unless exists $args->{'-distance'};
    $args->{'-limit'} = 50 unless exists $args->{'-limit'};
    # generic 
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};

    $self->throw unless
	scalar(grep {$args->{$_}} qw(-nulldist -object)) == 1;
    
    #

    my $key = 'SYNTENY';
    my $fherr = *STDERR;
    
    # 
    
    if ( my $orf = $args->{'-object'} ) {
	$self->throw("$self,$orf") unless $self->down->down->isa(ref($orf));

	my $mat = 
	    $orf->_compose_synteny_delta_matrix(@_); # [pos]->{species} == delta 
	
	if ( $args->{'-verbose'} >= 2 ) {
	    print {$fherr} $orf->data('SGD'), map { $args->{'-distance'}-$_ } 0..$#{$mat};
	    foreach my $sp ( @{$args->{'-restrict'}} ) {
		print {$fherr} $sp, $orf->data($sp), map { $mat->[$_]->{$sp} } 0..$#{$mat};
	    }
	}
	
	# 
	
	$self->throw unless my $hash = $self->{ $key };
	$self->throw if $args->{'-limit'} > $hash->{'LIMIT'};
	$self->throw if $args->{'-distance'} > $hash->{'DISTANCE'};

	# recenter the new matrix on the null matrices 

	my $adj = $hash->{'DISTANCE'} - $args->{'-distance'};	
	
	# 
	
	my ($m0,$m1,$ratiosum,$positions)=([],[],0,0);
	foreach my $i (0..$#{$mat}) {
	    next if $i == $#{$mat}/2;
	    if ( defined $mat->[$i] ) {
		$positions++; # not presently used 
	    } else { next; }

	    foreach my $sp ( @{$args->{'-restrict'}} ) {
		my $val = $mat->[$i]->{$sp}; 
		# see  _compose_synteny_delta_matrix() for setting of undef 
		next unless defined $val; # if undef mean we have no homolog. 0 OK (good..)
		$val = $args->{'-limit'} if abs($val) >= $args->{'-limit'}; # deal with infinitiy 
		$m0->[$i]->{ $sp } = $self->_access_synteny_null_dist('NULL_DIST', $i+$adj, $sp, $val);
		$m1->[$i]->{ $sp } = $self->_access_synteny_null_dist('SYN_DIST', $i+$adj, $sp, $val);
		$ratiosum += log($m1->[$i]->{ $sp }/$m0->[$i]->{ $sp }); # /log(10) ?
		#
		print {$fherr} $i, $sp, $val,
		$m0->[$i]->{ $sp }, $m1->[$i]->{ $sp }, $m1->[$i]->{ $sp }/$m0->[$i]->{ $sp }, 
		log($m1->[$i]->{ $sp }/$m0->[$i]->{ $sp }), $ratiosum, sprintf("%.3f", $ratiosum/$positions), 
		$orf->data($sp) if $args->{'-verbose'};
	    } 
	}
	return $ratiosum;
    }

    #############################
    # Null Dist 
    
    # initilize data structure 
    # we ignore genes more than $args->{'-limit'} from the query

    my ($realmat,$randmat);
    foreach my $i ( 0..$args->{'-distance'}*2 ) {
	foreach my $sp ( keys %HOMOLOGY ) {
	    %{$realmat->[$i]->{$sp}} = map { $_ => 1 } ( -$args->{'-limit'} .. $args->{'-limit'} );
	    %{$randmat->[$i]->{$sp}} = map { $_ => 1 } ( -$args->{'-limit'} .. $args->{'-limit'} );
	    $realmat->[$i]->{$sp}->{'TOTAL'} = $args->{'-limit'};
	    $randmat->[$i]->{$sp}->{'TOTAL'} = $args->{'-limit'};
	}
    }

    # index orfs for fast access 

    my $index=0;
    my %index = map { ++$index => $_ } $self->orfs;
    my %store = map { $_ => $index{$_}->_data } keys %index;

    # 

    if ( $args->{'-nulldist'} == 2 ) {
	
	foreach my $mode ( 0,1 ) {
	    
	    my $currmat;
	    if ( $mode == 1 ) { # randomize 
		my @orfs = $self->orfs;
		foreach my $o ( $self->orfs ) {
		    my $swap = splice(@orfs, int(rand( $#orfs )), 1);
		    # $self->swap( -object => $swap );
		    my ($d1,$d2) = ($o->_data, $swap->_data);
		    $o->_data( $d2 );
		    $swap->_data( $d1 );
		}
		$self->throw unless $#keys == -1;
		$currmat = $randmat;
	    } else { $currmat = $realmat; }
	    
	    # compute synteny data 
	    
	    foreach my $o ( $self->orfs ) {
		my $mat = 
		    $o->_compose_synteny_delta_matrix( -distance => $args->{'-distance'}, -restrict => undef);
		foreach my $i (  0 .. $args->{'-distance'}*2 ) {
		    foreach my $sp ( keys %HOMOLOGY ) {
			$currmat->[$i]->{$sp}->{ 'TOTAL' }++;
			$currmat->[$i]->{$sp}->{ $mat->[$i]->{$sp} }++;
			#print {$fherr} $mode, $o->name, $i, $sp, $mat->[$i]->{$sp}, $currmat->[$i]->{$sp}->{ $mat->[$i]->{$sp} };
		    }
		}
	    }	    

	    # restore associations between genes and data
	    
	    if ( $mode == 1 ) {
		map { $index{$_}->_data( $store{$_} ) } keys %index;
	    }
	}

    } elsif ( $args->{'-nulldist'} == 1 ) { # change location of foregound Orf only 
	
	foreach my $o ( $self->orfs ) {
	    my $mat = 
		$o->_compose_synteny_delta_matrix( -distance => $args->{'-distance'}, -restrict => undef);

	    my $store = $o->_data;
	    my $sample = int(rand($index))+1;
	    $o->_data( $index{ $sample }->_data );
	    my $swapmat = $o->_compose_synteny_delta_matrix( -distance => $args->{'-distance'} , -restrict => undef);
	    $o->_data( $store );
	    
	    foreach my $i (  0 .. $args->{'-distance'}*2 ) {
		foreach my $sp ( keys %HOMOLOGY ) {
		    $realmat->[$i]->{$sp}->{ 'TOTAL' }++;
		    $realmat->[$i]->{$sp}->{ $mat->[$i]->{$sp} }++;
		    # 
		    $randmat->[$i]->{$sp}->{ 'TOTAL' }++;
		    $randmat->[$i]->{$sp}->{ $swapmat->[$i]->{$sp} }++;
		    #print {$fherr} $o->name, $i, $sp, $swapmat->[$i]->{$sp}, $randmat->[$i]->{$sp}->{ $swapmat->[$i]->{$sp} };
		}
	    }
	}
    }

    # make randomized data and store originals 
    
    $self->{$key} = {
	'LIMIT' => $args->{'-limit'},
	'DISTANCE' => $args->{'-distance'},
	'NULL_DIST' => $randmat,
	'SYN_DIST' => $realmat
    };
    
    return $self;
}

sub _access_synteny_null_dist {
    my $self = shift;
    my ($key, $pos, $species, $value) = @_;

    ##########################################################
    # test for availablilty of sysnteny data / LOSS metric 
    # this is used by pillar() and (via pillar) by homolog();
    if ( ! @_ ) {
	return ( exists $self->{SYNTENY} ? $self : undef); # TRUE/FALSE. NOT 0.
    }
    ##########################################################
    
    $self->throw unless $self->{SYNTENY};

    $self->throw($key) unless my $data = $self->{SYNTENY}->{$key};
    $self->throw unless ref($data) eq 'ARRAY' || $key !~ /_DIST/;
    return $data unless defined $pos;

    $self->throw unless my $dist = $data->[$pos]->{$species};
    $self->throw unless ref($dist) eq 'HASH'; 
    return $dist unless defined $value;
    
    # my ($max) = sort { $b <=> $a } keys %{$dist}; # takes too much time 
    $value = $self->{SYNTENY}->{'LIMIT'} if $value > $self->{SYNTENY}->{'LIMIT'};
    $self->throw($value) unless $dist->{$value} > 0;

    $self->throw unless exists $dist->{'TOTAL'} && $dist->{'TOTAL'} > 0;    
    return ( wantarray ? ($dist->{$value},$dist->{'TOTAL'}) : ($dist->{$value}/$dist->{'TOTAL'}) );
}

=head2 quality( -object => $orf, -nulldist => #samples, -equalize => 1|0) 

    Compute quality of Orthogroups based on population statistics. 
    We compute distributions for %ID, synteny and length variation among
    members of an orthogropu and then compare query to obtain Z-score for
    each of 3 criteria. 

    With no arguments compares all OGs to the existing null distribution. 

    -object : compute quality of specified OG 
    -nulldist : build the null distribution based on -nulldist OGs. 
    -equalize : weight the 3 quality metrics equally

=cut 

sub quality {
    my $self = shift;
    my $args = {@_};

    $args->{'-object'} = undef unless exists $args->{'-object'};
    $args->{'-unit'} = 20; # unless exists $args->{'-unit'};
    $args->{'-nulldist'} = undef unless exists $args->{'-nulldist'};
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    $args->{'-equalize'} = 1 unless exists $args->{'-equalize'};
    $args->{'-sowh'} = 20 unless exists  $args->{'-sowh'};

    $self->throw if $args->{'-object'} && $args->{'-nulldist'};
    my $key = 'QUALITY';

    ##############################    
    # P2-ICEB-ULUR-AHUC-ALAK
    
    my $count=0;
    if ( $args->{'-nulldist'} ) {

	my @ogs = grep { $_->assign !~ /RNA/ } $self->orthogroups();
	$self->throw if $#ogs < $args->{'-nulldist'};
	
	my %hash;
	foreach my $o ( sort { $a->_internal_id <=> $b->_internal_id } @ogs ) {
	    
	    # for %ID we index by species pair 
	    
	    foreach my $i ( sort {$a->organism cmp $b->organism} $o->_orthogroup ) {
		foreach my $j ( sort {$a->organism cmp $b->organism} $o->_orthogroup ) {
		    next unless $i->organism lt $j->organism;
		    my $tag = $i->organism.':'.$j->organism;
		    push @{$hash{'PID'}{ $tag }}, $i->bl2seq(-object => $j, -identity => 1);		    
		}
	    }

	    # synteny 

	    foreach my $i ( sort {$a->organism cmp $b->organism} $o->_orthogroup ) {
		foreach my $j ( sort {$a->organism cmp $b->organism} $o->_orthogroup ) {
		    next unless $i->organism lt $j->organism;
		    my $tag = $i->organism.':'.$j->organism;
		    push @{$hash{'SYN'}{ $tag }}, $i->synteny_conserved(-object => $j);
		}
	    }
	    
	    #########################################
	    # depracated : index by telomere distance 	    
	    if ( 1==0 ) {
		foreach my $i ( grep { $_->up->length >1e5} $o->_orthogroup ) {
		    my $tdist = $i->telomereDistance;
		    my $index = int( ($tdist + $args->{'-unit'} ) / $args->{'-unit'} );
		    $index = 10 if $index > 10;
		    push @{ $hash{'SYN'}{ $i->organism.$index } }, $i->synteny( -hyper => );
		}
	    }
	    #########################################

	    # no indexing for length 
	    
	    my ($m,$sd) = _calcMeanSD( map {$_->length} $o->_orthogroup );
	    my $loglen = sprintf("%.1f", log($m)/log(10));	
	    push @{$hash{ 'LEN' }{ $loglen }}, $sd;	    
	    last if ++$count >= $args->{'-nulldist'};
	}
	
	# 

	delete $self->{$key};
	foreach my $i ( keys %hash ) {
	    foreach my $j ( sort keys %{$hash{$i}} ) {
		print $i,$j,$#{$hash{$i}{$j}}, _calcMeanSD( @{$hash{ $i }{ $j } }) 
		    if $args->{'-verbose'};
		$self->warn("LOW COUNT: $i,$j,$#{$hash{$i}{$j}}") if $#{$hash{$i}{$j}} <= 20;
		$self->{$key}->{ $i }->{ $j } = 
		    [ _calcMeanSD( @{$hash{ $i }{ $j } }), scalar( @{$hash{ $i }{ $j } } ) ];
	    }
	}

    } elsif ( $args->{'-object'} ) {

	$self->throw unless exists $self->{$key} &&
	    exists $self->{$key}->{'PID'};
	$self->throw unless $self->down->down->isa(ref( $args->{'-object'} ));	
	my ($qual,$o) = ($self->{$key}, $args->{'-object'});

	my ($total,$len_z,$pid_z,$syn_z);

	# 
	
	my @org = map { $_->organism } $o->_orthogroup;
	my @tags = map { [sort ($org[$_],$org[$_-1])] } 0..$#org;
	foreach my $pair ( @tags ) {
	    my ($i) = grep { $_->organism eq $pair->[0] } $o->_orthogroup; 
	    my ($j) = grep { $_->organism eq $pair->[1] } $o->_orthogroup; 
	    my $pid = $i->bl2seq(-object => $j, -identity => 1);		
	    
	    my $tag = join(':', @{$pair});
	    $self->throw($tag) unless exists $qual->{'PID'}->{$tag};
	    my ($mean,$sd, $n) = @{$qual->{'PID'}->{$tag}};
	    
	    if ( $sd > 0 ) {
		my $zscr = ($pid - $mean ) / $sd;
		$pid_z += ( $zscr / ($args->{'-equalize'} ? scalar(@org) : 1) );
		#print $i->name , $j->name, $tag, $pid, $mean, $sd, $zscr, $pid_z;
	    }

	    # syntney ...

	    my $syn = $i->synteny_conserved(-object => $j);

	    $self->throw($tag) unless exists $qual->{'SYN'}->{$tag};
	    my ($mean,$sd, $n) = @{$qual->{'SYN'}->{$tag}};
	    
	    if ( $sd > 0 ) {
		my $zscr = ($pid - $mean ) / $sd;
		$syn_z += ( $zscr / ($args->{'-equalize'} ? scalar(@org) : 1) );
		#print $i->name , $j->name, $tag, $pid, $mean, $sd, $zscr, $pid_z;
	    }	    
	}

	#########################################
	# depracated : does not provide relevant info 
	if ( 0 == 1 ) {
	    foreach my $i ( $o->_orthogroup ) {
		my $binom = $i->binomialsynteny;
		# 
		my $tdist = $i->telomereDistance;
		my $index = int( ($tdist + $args->{'-unit'} ) / $args->{'-unit'} );
		$index = 10 if $index > 10;
		my $tag = $i->organism.$index;
		my ($mean,$sd, $n) = @{$qual->{'SYN'}->{ $tag }};
		
		if ( $sd > 0 ) {
		    my $zscr = ($binom - $mean ) / $sd;
		    $syn_z += ( $zscr / ($args->{'-equalize'} ? scalar(@org) : 1) );
		    #print $i->name , $tdist, $tag, $binom, $mean, $sd, $zscr, $syn_z;  
		}
	    }
	}
	#########################################

	# 

	if ( $o->_orthogroup ) {
	    my ($m_og,$sd_og) = _calcMeanSD( map {$_->length} $o->_orthogroup );
	    my $tag = sprintf("%.1f", log($m_og)/log(10));	
	    my ($mean, $sd, $n) = @{$qual->{'LEN'}->{ $tag }};	
	    if ( $sd > 0 ) {
		my $zscr = ($sd_og - $mean ) / $sd;
		$len_z -= $zscr; # unlike others, lower is better than higher 
		#print {STDERR} $o->identify, $tag, $m_og, $mean, $sd, $sd_og, $zscr, $len_z;  
	    }
	}

	my $total = $len_z + $pid_z + $syn_z;
	my @res = map { sprintf("%.1f", $_) } ( $total, $pid_z, $syn_z, $len_z );
	return ( wantarray ? @res : $res[0] ) ;
	
    } else {
	    
	$self->quality( -nulldist => 1000 ) unless exists $self->{$key};
	
	my $check=0;
	foreach my $og ( grep {$_->coding}  $self->orthogroups ) {

	    my @qual = $self->quality( -object => $og, @_ );
	    my $res = ( $qual[0]/3 < -1 || $qual[1] < -2 || $qual[2] < -2 || $qual[3] < -2 ? 1 : 0);
	    $check++ if $res;

	    #$og->kaks;
	    my $sowh = $og->phyml( -sowh =>  $args->{'-sowh'} ) if $args->{'-sowh'};
	    my @sowh = map { $sowh->{$_} } qw(DELTA PVAL RESULT);
	    $og->summarize( @qual, $res, @sowh ); # ++$count, $check,
	}    
    }
    
    return $self;
}


sub ohnologComparative2 {
    my $self = shift;
    my $args = {@_};

    $self->throw unless $self->bound;
    $args->{'-refresh'} = undef unless exists $args->{'-refresh'};
    $args->{'-orthogroup'} = 1 unless exists $args->{'-orthogroup'};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};
    
    my $fh = STDERR;
    my $count= scalar($self->bound)+1;

    print scalar( grep {$_->ohnolog} $self->orfs);    
    print scalar( grep {$_->ogid} $self->orfs);
    #map { print $_->species, scalar( grep {$_->ohnolog} map {$_->stream} $_->stream ) } $self->iterate;

    ################################
    # prep work 
    ################################

    # for ohnologs  ...
    map { $_->ohnologs(-window => 7, -drop => 5, -penalty => 20) } $self->iterate if $args->{'-refresh'};
    my %og = map { $_->ogid => $_ } grep {$_->ogid} map {$_->stream} $self->stream;  

    # for orthologs ...
    my %kaks;
    foreach my $k (qw(KA KS KAKS)) {
	my ($dn,$dn_sd) = &_calcMeanSD( grep {$_>0 && $_ <10} map {$_->data($k)} $self->orfs);
	$kaks{$k}{'MEAN'}=$dn;
	$kaks{$k}{'SD'}=$dn_sd;
	print $k, $dn,$dn_sd;
    }

    ################################
    # 
    ################################

    open($fho, ">pre.ohnologs");
    map { print $fho $_->name, $_->ohnolog->name, $_->ygob, $_->sgd, $_->ogid, $_->identify } 
    grep { $_->ohnolog } $self->orfs;
    close($fho);

    return $self;
}

=head2 ohnologComparative

    Compare ohnologs across species and leverage 
    information to capture all likely ohnos across all. 

    Possibilities:
    1. OG1 ohnos are consistent, OG2 has no extra 
    2. OG1 ohnos are consistent, OG2 has extra
    3. OG1 ohnos are not consistent, OG2 has extra 

    We use _dp_align to align sister regions across species
    and leveage as much information as possible to identify 
    related orthologs and paralogs. 

=cut 

sub ohnologComparative {
    my $self = shift;
    my $args = {@_};

    $self->throw unless $self->bound;
    $args->{'-refresh'} = undef unless exists $args->{'-refresh'};
    $args->{'-orthogroup'} = 1 unless exists $args->{'-orthogroup'};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};
    
    my $fh = STDERR;
    my $count= scalar($self->bound)+1;

    print scalar( grep {$_->ohnolog} $self->orfs);    
    print scalar( grep {$_->ogid} $self->orfs);
    #map { print $_->species, scalar( grep {$_->ohnolog} map {$_->stream} $_->stream ) } $self->iterate;

    ################################
    # prep work 
    ################################

    # for ohnologs  ...
    map { $_->ohnologs(-window => 7, -drop => 5, -penalty => 20) } $self->iterate if $args->{'-refresh'};
    my %og = map { $_->ogid => $_ } grep {$_->ogid} map {$_->stream} $self->stream;  

    # for orthologs ...
    my %kaks;
    foreach my $k (qw(KA KS KAKS)) {
	my ($dn,$dn_sd) = &_calcMeanSD( grep {$_>0 && $_ <10} map {$_->data($k)} $self->orfs);
	$kaks{$k}{'MEAN'}=$dn;
	$kaks{$k}{'SD'}=$dn_sd;
	print $k, $dn,$dn_sd;
    }

    ################################
    # 
    ################################

    open($fho, ">pre.ohnologs");
    map { print $fho $_->name, $_->ohnolog->name, $_->ygob, $_->sgd, $_->ogid, $_->identify } 
    grep { $_->ohnolog } $self->orfs;
    close($fho);

    # go through each orthogroup and see where we can "fill"
    # an OG with ohnologs using info borrowed from others. 

    my %seen;
    foreach my $o ( sort {$b->length <=> $a->length} grep {$_->orthogroup} map {$_->stream} $self->stream ) {
	next if exists $seen{$o->ogid};
	$seen{$o->ogid}=1;
	my $ohno = $o->_ohnolog_consistency;	
	
	################################
	# Special Case 1
	# Total consistency -- nothing to do
	################################	
	
	next if scalar(keys %{$ohno})==1 && (keys %{$ohno})[0] != 0; # perfect ohno or singleton
	my ($ogid, @other) = sort { $#{$ohno->{$b}} <=> $#{$ohno->{$a}} } grep {$_>0} keys %{$ohno};	
	next unless ! $args->{'-debug'} || $ogid == $args->{'-debug'};

	#next unless $#{$ohno->{$ogid}}==0 && $#{$ohno->{'-1'}}==3;	
	print (">\nALL\t".$o->name, $o->identify, (map { $_.'='.scalar(@{$ohno->{$_}}) } sort {$a <=> $b} keys %{$ohno}))
	    if $args->{'-verbose'} >= 2;

	################################
	# Special Case 2
	# No OG info so cannot create Ohnologs 
	# might be able to create OGs if adequate synteny 
	################################
	
	if ( ! $ogid ) {
	    next unless $args->{'-orthogroup'} == 1;
	    # any combination of 0/-1 is OK
	    # can define a locus based on 0s
	    # if consistent search here for each -1
	    # if we get a full set can try make OG

	    # simplified version (no -1 search): 
	    if ( scalar(@{$ohno->{'0'}}) == $count ) {
		my @x=map {$_->ohnolog} @{$ohno->{'0'}};
		my @scores;
		for my $i (0..($#x-1)) {
		    for my $j (($i+1)..$#x) {
			push @scores, $x[$i]->synteny_conserved(-object => $x[$j]) || 0;
		    }
		}
		
		print &_calcMeanSD(@scores);;
		#map { $_->output(-ohno => 1, -og => 1) } @x;
		
		my ($m,$sd) = &_calcMeanSD(@scores);
		if ( $m > 30 && $sd < 5) {
		    my ($ref) = grep {$_->organism eq $self->organism} @x;
		    $ref->_define_orthogroup(-object => [grep {$_ ne $ref} @x], -kaks => 1);
		    foreach my $k (qw(KA KS KAKS)) {
			if ( ($ref->data($k)-$kaks{$k}{'MEAN'}) > 2*$kaks{$k}{'SD'}) {
			    $ref->_dissolve_orthogroup();
			}
		    }
		}		
	    }
	    next;
	}

	################################
	# Special Case 3
	# >1 OG. Simplify or skip. 
	# need to code a mechanism for breaking ohnologs. 
	################################

	if ( (@other || @{$ohno->{'0'}}) || # require ohnologs to be broken. do later.  
 	     scalar( @{$ohno->{$ogid}} ) == 1) { # temp.. weak evidence 
	    next;
	}
	
	################################
	# only one candidate OG 
	# this is a relatively low bar. 
	# 1) We allow breaking ohnologs that are not in OGs on the recip side only. 
	# 2) The condition for breaking is that the weight of evidence for the recip
	# OG is greater than all other possibilities and only one OG is implicated. 
	# eg 123 => [scer,sbay,spar], 0 => [smik], -1 => [skud] . we break the smik link. 
	# 3) we make 123 OGs for smik and skud. 
	# this is not really a good approach. we should be weighing the evidence for the 
	# alternative ohnoogs scenarios with some kind of weight based on the evidence 
	# from other species. should probably be iterative like the OG creation method. 
	################################
	
	my $og = $og{$ogid};
	my $recip = $og->_ohnolog_consistency;
	my ($rid, @recip) = sort { $#{$recip->{$b}} <=> $#{$recip->{$a}} } grep {$_>0} keys %{$recip};	
	print ("RECIP\t".$og->name, $og->identify, (map { $_.'='.scalar(@{$recip->{$_}}) } sort {$a <=> $b} keys %{$recip}))
	    if $args->{'-verbose'} >= 2;
	next unless @{$recip->{$rid}}>=$bound/2 && ! @recip; 	

	##########################
	# break ohnologs that are not in OGs (0 -> -1)
	##########################	

	map { $_->ohnolog(undef) } @{$recip->{'0'}}; # we break these 
	
	##########################
	# make new ohnologs 
	##########################	

	for my $noh ( grep {!$_->ohnolog} $o->_orthogroup ) {
	    my $meth = uc($noh->organism);
	    my $pair = $og->$meth || $og;
	    map { $_->throw if $_->ohnolog } ($noh,$pair);
	    $noh->ohnolog( $pair );
	}
	$seen{$og->ogid}=1;
	
	# report ...
	my $ohno = $o->_ohnolog_consistency;
	print ("NEW\t",$o->name, $o->identify, (map { $_.'='.scalar(@{$ohno->{$_}}) } sort {$a <=> $b} keys %{$ohno}))
	    if $args->{'-verbose'};
    }
    
    ##########################
    # DEBUGGING 
    ##########################	

    foreach my $o ( sort {$b->length <=> $a->length} grep {$_->orthogroup} map {$_->stream} $self->stream ) {
	my $ohno = $o->_ohnolog_consistency;
	next if scalar(keys %{$ohno})==1 && (keys %{$ohno})[0] != 0; # perfect ohno or singleton
	#print ++$yy,(map { $_.'='.scalar(@{$ohno->{$_}}) } sort {$a <=> $b} keys %{$ohno});
	#map { $_->output(-ohno => 1, -og => 1) } grep {defined} map {$_->ohnolog} map { @{$_} } values %{$ohno};
    }
    
    print scalar( grep {$_->ogid} $self->orfs);    
    print scalar( grep {$_->ohnolog} $self->orfs);    
    
    open($fho, ">post.ohnologs");
    map { print $fho $_->name, $_->ohnolog->name, $_->ygob, $_->sgd, $_->ogid, $_->identify } 
    grep { $_->ohnolog } $self->orfs;
    close($fho);
    exit;

    #map { print $_->species, scalar( grep {$_->ohnolog} map {$_->stream} $_->stream ) } $self->iterate;    
    return $self;
}

=head2 ohnologs(-window => 7, -max => auto|i, -drop => 5, 
    -penalty => 20, -debug => 'Anc')
    
    Call ohnolog pairs. Not fully parametrized and should really 
    be comparing to a proper null background model but performs 
    pretty well. When ohnologs are pooled via orthogroups the 
    numbers reach 98% agreement with YGOB. 
    
    -window : defines the region around 

=cut 

sub ohnologs {
    my $self = shift;
    my $args = {@_};

    $args->{'-debug'} = undef unless exists $args->{'-debug'};
    $args->{'-window'} = 7 unless exists $args->{'-window'};
    $args->{'-max'} = sprintf("%d", $args->{'-window'}*2)+1 unless exists $args->{'-max'};
    $args->{'-drop'} = 5 unless exists $args->{'-drop'};
    $args->{'-penalty'} = 20  unless exists $args->{'-penalty'};
    
    #######################################
    # Initialize: 
    # Remove all pre-existing ohnolog data 
    # Index protein coding genes by YGOB locus id 
    #######################################     
    
    my $attr = '_ohno_temp_var'.int(rand(100));
    my $ohno = scalar( grep { $_->ohnolog } $self->orfs );
    
    my %hash;
    foreach my $orf ( grep { $_->ygob } map { $_->stream } $self->stream ) {
        $orf->data($attr => undef);
        if ( my $oh = $orf->ohnolog ) {
            $oh->ohnolog(undef); # undefs both 
	}
        push @{ $hash{ ( $orf->assign =~ /RNA/ ? $orf->data('GENE') :  $orf->ygob ) } }, $orf;
    }
    
    #######################################
    # Iterate through YGOB families and process each family independently.
    # We take the YGOB designations on trust. 
    #######################################

    my $family_c;
    my %count;   
    foreach my $anc ( grep { $#{$hash{$_}} >= 1 } keys %hash) {
        next unless ( ! $args->{'-debug'} || $anc eq $args->{'-debug'} );	
        $count{'0.2copy'}++;

	#######################################
        # Iterate through familiy members. Nominate candidate ohnologs or quit: 
	# 1. exclude tandem duplicates / near neighbours 
	# 2. compute synteny vs pre-WGD ancestor (YGOB) and store
	# 3. sort by synteny and accept top two as candidate pair
	#######################################
	
        my $old;
        foreach my $cand ( sort {$a->name <=> $b->name} @{$hash{$anc}} ) {
            # ignore tandems 
            if ( $old && $old->up == $cand->up ) {
                my @inter = grep { $_->assign ne 'GAP' }  $old->intervening($cand);
                #print map {$_->name} @inter;
                $old=$cand and next if $#inter <= 2;
            }

            # asses synteny 
            my $newsyn = $cand->querysynteny(
                -spanning => 0, 
                -difference => 10, 
                -distance => 5,
                -restrict => ['YGOB'] #, 'KLAC', 'SKLU', 'KWAL']
                );
            $cand->data($attr => $newsyn);
            $old = $cand;
	    
	    # output
            if ($args->{'-debug'}) {
		$cand->oliver(-append => [$newsyn]);
		$cand->context();
	    }
	}

	# sort by synteny

        my ($one,$two,$other) = map { $_->data($attr => undef); $_ }  sort { 
            $b->data($attr) <=> $a->data($attr)
        } grep { $_->data($attr) } @{$hash{$anc}};
	
        next unless $one && $two; # Terminate families that have no candidate pair.
        $count{'1.synteny'}++;

	#######################################
	# Terminate candidates that fail a more stringent test for local duplication. 
	#######################################
        
        next if $one->up->id == $two->up->id && 
            #( $one->distance(-object => $two) < 50 || 
            $one->distance(-object => $two, -bases => 1) < 5e4; 
        #);        
        $count{'2.nontandem'}++;

	#######################################
        # Collect syntey data for designated candidate pair. 
	# Collect for 3rd choice also-- a control of sorts.
	# Data collected are indices of neighbouring gene on the ancestral chr.
	# Cand1 (Anc_4.7) :  1,4,6,<7>,9,11,15
	# Cand2 (Anc_4.7) :  2,3,5,<7>,9,10,12
	# Cand3 (Anc_4.7) :  ,,,<7>,,19,
	#######################################
       
        
        $one->ygob =~ /(\d+)\.(\d+)/;
        my ($chr,$pos) = ($1,$2);
	
	my $onecount;
        my @one =  map { $_->ygob =~ /\.(\d+)/; $1 } grep { $_->ygob =~ /_$chr\./ } 
	map { $onecount++; $_ } grep { $_->ygob ne $one->ygob }
        ( map { $one->traverse(-direction => $_, -distance => $args->{'-window'}) } ('left','right') );
	
	my $twocount;        
	my @two = map { $_->ygob =~ /\.(\d+)/; $1 } grep { $_->ygob =~ /_$chr\./ } 
	map { $twocount++; $_ } grep { $_->ygob ne $two->ygob }        
        ( map { $two->traverse(-direction => $_, -distance => $args->{'-window'}) } ('left','right') );
	
	my $othercount=0;
	my @other = map { $_->ygob =~ /\.(\d+)/; $1 } grep { $_->ygob =~ /_$chr\./ } 	
	map { $othercount++; $_ } grep { $_->ygob ne $other->ygob }
	( map { $other->traverse(-direction => $_, -distance => $args->{'-window'}) } ('left','right') ) 
	    if $other;
	
        next unless @one && @two; # cannot make ohnos without synteny 
	die unless scalar( grep {defined} @one)>0;
        $count{'3.syntenydata'}++;

	#######################################
	# Exclude false positives due to breaking of a contig in the middle of a gene:
	# Cand1 (Anc_4.7) :  1,4,6,<7>,,,
	# Cand2 (Anc_4.7) :  ,,,<7>,9,10,12	
	# We toleate candidates that either :
	# 1. the genes align well 
	# 2. are both clearly internal to a contig. 
	#######################################
	
        my @o_sort = sort {$a <=> $b} @one;
        my @t_sort = sort {$a <=> $b} @two;
	if ( ( @o_sort[-1] < $pos && $pos < @t_sort[0] ) || ( @t_sort[-1] < $pos && $pos < @o_sort[0] ) ) {
            # if they align might be ok..
            my $align = $one->overlap(-evalue => 1e-5, -object => $two, -compare => 'seq');
            
            # it they are clearly internal..
            my @ol = $one->traverse(-direction => 'left', -distance => 5);
            my @or = $one->traverse(-direction => 'right', -distance => 5);
            my @tl = $two->traverse(-direction => 'left', -distance => 5);
            my @tr = $two->traverse(-direction => 'right', -distance => 5);
            my $internal=1 if ($#ol >= 2 && $#or >= 2) || ( $#tl >= 2 && $#tr >= 2 ); # only require 1.. 
            
            unless ( $align >= 0.5 || $internal ) {
                $one->output(-fh => $fh,  -append => [$two->name]);
                $two->output(-fh => $fh, -append => [$one->name]);
                next;
            }
        }
        $count{'4.nonterminii'}++;

	#######################################
	# If there is a third gene in the family - presumptive non-ohnolog - 
	# we use this to establish a null by comparing delta between second 
	# candidate gene and the extra gene. 
	# Otherwise we default to a predefined value.
	#######################################
	
	my ($mat,$scr,$nullstatistic)=(undef,0,$args->{'-penalty'});
	
        if ( @other ) {
	    $nullstatistic=0;
            for my $i (0..$#other) { map { $mat->[$i]->[$_] = abs( $other[$i] - $two[$_] ) }  (0..$#two); } 
            $nullstatistic += $scr while ( $scr = &_findMatrixMin($mat, $args->{'-penalty'}) );
            $nullstatistic += $args->{'-penalty'} * 
		( $othercount < $twocount ? $othercount-scalar(@other) : $twocount-scalar(@two) );
            $nullstatistic /= ( $othercount < $twocount ? $othercount : $twocount );
        }
	
	#######################################
        # Compute a test statistic and sanity check. 
	# The test statistic is a function of the extent of interleaving: 
	# 1. construct a distance matrix for each possible pair of genes
	# across the two sister regions. D is based on ancestral gene index. 
	# 2. iteratively choose neighbours to minimize the distance between 
	# them. constrain to pairs that are consistent with chromosomal order.
	# 3. Apply a penalty based on the lack of information encoded in the
	# shorter region. 
	#######################################

	my ($mat,$scr,$teststatistic)=(undef,0,0);	
        for my $i (0..$#one) { map { $mat->[$i]->[$_] = abs( $one[$i] - $two[$_] ) }  (0..$#two); } 
        $teststatistic += $scr while ( $scr = &_findMatrixMin($mat, $args->{'-penalty'}) ); 
        $teststatistic += $args->{'-penalty'} * 
	    ( $onecount < $twocount ? $onecount-scalar(@one) : $twocount-scalar(@two) );
        $teststatistic /= ( $onecount < $twocount ? $onecount : $twocount );
	
	# 
	
        print ++$family_c, $anc, $#one, $#two, $#other, 
	sprintf("%.2f", $teststatistic), sprintf("%.2f",$nullstatistic), 
	$onecount,$twocount,$othercount, $args->{'-penalty'}, $args->{'-max'};
        next unless $teststatistic <= $args->{'-max'};
        $count{'5.score'}++;
	
	#######################################
        # close it out... apply 2 primary criteria 
	#######################################

        next unless ($nullstatistic - $teststatistic) >= $args->{'-drop'};
        $count{'6.drop'}++;

        $one->ohnolog( $two );
        #$two->ohnolog( $one );
    }
    close $fh;

    map { print $_,$count{$_} } sort keys %count if $args->{'-verbose'};
    exit;
    
    map { $_->data($attr => 'delete') } $self->orfs;
    return $self;
}

#######################################
#######################################

sub ohnologs2 {
    my $self = shift;
    my $args = {@_};

    $args->{'-debug'} = undef unless exists $args->{'-debug'};
    $args->{'-window'} = 7 unless exists $args->{'-window'};
    $args->{'-verbose'} = 0 unless exists $args->{'-verbose'};

    $args->{'-synteny'} = 0 unless exists $args->{'-synteny'};
    $args->{'-homology'} = 'YGOB' unless exists $args->{'-homology'};

    $args->{'-max'} = sprintf("%d", $args->{'-window'}*2)+1 unless exists $args->{'-max'};
    $args->{'-drop'} = 5 unless exists $args->{'-drop'};
    $args->{'-penalty'} = 20  unless exists $args->{'-penalty'};

    #######################################
    # make some vars for the next phase 
    
    my $fh = STDOUT;
    my $attr = '_ohno_temp_var'.int(rand(100));
    my $ohno = scalar( grep { $_->ohnolog } $self->orfs );    
  
    # 

    my $keyX = 'YGOB';    
    my $key0 = 'ANC';
    my $key1 = uc($self->organism).'1';
    my $key2 = uc($self->organism).'2';
    
    # 
    
    my %matrix = (
	'OHNO' => 4,
	'GAP' => 0,
	'CROSS' => 2,
	'SAME' => 0
	);
    
    my %anc;
    foreach my $orf ( grep { #$_->assign =~ /RNA/ || 
	$_->ygob } map { $_->stream } $self->stream ) {
	$orf->data($attr => undef);
	if ( my $oh = $orf->ohnolog ) {
	    $oh->ohnolog(undef);
	}
	push @{ $anc{ ( $orf->assign =~ /RNA/ ? $orf->data('GENE') :  $orf->ygob ) } }, $orf;
    }

    #######################################
    # cycle through all Anc families 

    my %count;   
    foreach my $anc ( grep { $#{$anc{$_}} >= 1 } keys %anc) {
	next unless ( ! $args->{'-debug'} || $anc eq $args->{'-debug'} );	
	$anc =~ /Anc_(\d+)\.(\d+)/ || $self->throw;
	my ($chr,$pos) = ($1,$2);
	
	###################
	# cycle through homologs and get synteny 
	###################
	
	my $old;
	foreach my $cand ( sort {$a->name <=> $b->name} @{$anc{$anc}} ) {
	    # randomly choose 1 of 2 tandems 
	    if ( $old && $old->up == $cand->up ) {
		my @inter = grep { $_->assign ne 'GAP' }  $old->intervening($cand);
		$old=$cand and next if $#inter <= 2;
	    }
	    
	    # asses synteny 	    
	    $cand->data($attr => $cand->hypergob);
	    $old = $cand;
	    
	    if ( $args->{'-debug'} ) {
		$cand->oliver(-append => [ $cand->data($attr) ]);
		map { $_->oliver } $cand->context( -all => 1 , -distance => 6 ); 
		print;
	    }
	}
	
	######################################	
	# basic synteny qualification -- we set the bar LOW
	
	my ($x, @init) = map { $_->data($attr => undef); $_ }  
	sort { $b->data($attr) <=> $a->data($attr) } 
	grep { $_->assign ne 'REPEAT' }
	grep { $_->data($attr) > $args->{'-synteny'} } @{$anc{$anc}};
	next unless $x;
	
	# _filter_tandems will choose one member from each tandem array.
	# tandem arrays are defined by number of interventing genes. 
	# not a good method but performs OK for yeast. 	
	# 2 x windowsize prevents same genes aligning to each other
	# in alignment phase.
	# in scer we miss YBR024W and YBR037C ohnologs. not sure how to recover. 
	
	my @cand = $x->_filter_tandems(-object => \@init, -distance => $args->{'-window'}*2);
	next unless $cand[0] && $cand[1];

	######################################	
	# We use a metric which rewards interleaving of 
	# genes on opposite regions (exact matches are also rewarded)
	# and thus can recognize duplicate origin for segments
	# that have no remaining duplicates. 
	######################################		

	print ">>$anc [$key0, $key1, $key2]" if $args->{'-verbose'};
	
	my %sisters;
	CAND: for my $i (0..($#cand-1)) {
	    my @i_array = $cand[$i]->context(-distance => $args->{'-window'}, -self => 1);

	    # define order and valid range. only Ancs in reference are 
	    # admitted for scoring so we constrain ultimate range/score without 
	    # without having to edit the actual gene order array i_array. 
	    
	    my (%hash,$i_rt);	
	    my @i_sort = sort {$a <=> $b} grep {defined} map { $_->ygob =~ /\.(\d+)/; $1 } grep { $_->ygob =~ /_$chr\./ } @i_array;
	    my @i_prune = &_prune_from_ends(\@i_sort, $args->{'-window'}+1);	   
	    map { $i_rt += ($i_prune[$_]>$i_prune[$_-1] ? 1 : -1) } 1..$#i_prune; # determine order on chr
	    
	    # TRACK 1!
	    $hash{ $key1 } = [ $i_rt > 0 ? @i_array : reverse @i_array ];
	    
	    for my $j ( ($i+1)..$#cand) {
		my @j_array = $cand[$j]->context(-distance => $args->{'-window'}, -self => 1);
		
		my $j_rt;
		my @j_sort = sort {$a <=> $b} grep {defined} map { $_->ygob =~ /\.(\d+)/; $1 } grep { $_->ygob =~ /_$chr\./ } @j_array;
		my @j_prune = &_prune_from_ends(\@j_sort, $args->{'-window'}+1);
		map { $j_rt += ($j_prune[$_]>$j_prune[$_-1] ? 1 : -1) } 1..$#j_prune; # determine order on chr

		# TRACK 2!
		$hash{ $key2 } = [ $j_rt > 0 ? @j_array : reverse @j_array ];

		# TRACK 0! #####################
		# use the min and max indices on the relevant chr to make the 
		# ancestral gene order. should prune ends...
		
		my ($max) = sort {$b <=> $a} ($i_prune[-1],$j_prune[-1]);
		my ($min) = sort {$a <=> $b} ($i_prune[0],$j_prune[0]);		

		my $genome = $self->clone;
		$genome->organism( $key0 );
		my $contig = ref($self->down)->new(SEQUENCE => 1e5 x 'ATGC', ID => 1);
		$genome->add(-object => $contig);
		#print $genome->organism, $min, $max;
		
		for my $pos ( $min..$max ) {
		    my $homol = 'Anc_'.$chr.'.'.$pos;
		    my $orf = ref($cand[$i])
			->new(
			START => $pos*10,
			STOP => ($pos*10)+3,
			STRAND => 1,
			UP => undef
			);
		    $orf->data( $keyX => $homol );
		    $contig->add(-object => $orf);
		    push @{$hash{ $key0 }}, $orf;
		}
		$contig->index;
		
		################################		
		# align two regions to the ancestor 
		
		my $align = $self->dpalign(
		    # alignment 
		    -hash => \%hash, 
		    -order => [$key0, $key1, $key2 ], #$key0 must be first 
		    -reference => $key0,
		    -global => 1,
		    # scoring 
		    -match => $keyX,
		    -mismatch => -1000,
		    -gap => 0,
		    -inversion => -1,
		    -trna => 1,
		    # 
		    -verbose => 0
		    );
		
		# we scrub anything that is not relevant to scoring.
		# they were kept to aid for alignment.
		
		my (@clean);
		foreach my $q ( grep { $_->{$key0} } @{$align} ) {				
		    foreach my $k ( grep {/\-/} keys %{$q} ) {
			my ($jnk,$k2) = split/\-/,$k;
			$q->{$k2}=$q->{$k};
			delete $q->{$k};
		    }
		    #$score{'GAP'}++ unless $q->{$key1} || $q->{$key2};
		    push @clean, $q; # if $q->{$key1} || $q->{$key2};
		}
		
		my ($ccc,$ddd);
		my ($q,$qq) = grep { $clean[$_]->{$key1} || $clean[$_]->{$key2} } 0..$#clean;
		next CAND unless defined $q && defined $qq;
		until ( $q==0 && ($qq-$q <= $args->{'-window'}) ) {
		    splice(@clean, 0, ($q+1));
		    ($q,$qq) = grep { $clean[$_]->{$key1} || $clean[$_]->{$key2} } 0..$#clean;
		    next CAND unless defined $q && defined $qq;
		    $self->throw if ++$ccc >= 1e4;
		}
		my @rev = reverse @clean;
		my ($q,$qq) = grep { $rev[$_]->{$key1} || $rev[$_]->{$key2} } 0..$#rev;
		next CAND unless defined $q && defined $qq;
		until ( $q==0 && ($qq-$q <= $args->{'-window'}) ) {
		    splice(@rev, 0, ($q+1));
		    ($q,$qq) = grep { $rev[$_]->{$key1} || $rev[$_]->{$key2} } 0..$#rev;
		    next CAND unless defined $q && defined $qq;
		    $self->throw if ++$ddd >= 1e4;
		}
		my @clean = reverse @rev;

		################################		
		# score the alignment 
		# the nice way to do this is have a log odds score 
		# of observing the alignment under two models: 
		# 1) a wgd occurred
		# 2) no wgd occurred
		# the second is kind of tricky for 2 regions: 
		# 1. there are many different mutational paths. we should weight and sum over. 
		# we need weights for duplications/deletions of given lengths. 
		# 2. for insertions the chance of three going in close in the genome 
		# and in the right order is very small. probs will get tiny fast.
		# the first is relatively easy. we just count deletions and assume max size. 
		# (is this the problem that Gavin solved?)
		#####
		# ... do we have the right null? 
		# Given that a wgd occurred, what are odds that this pair of regions, are sisters? 
		# This is actually closer to comparing the first choice sister region pair
		# to the second choice sister pair and asking can we infer the right pair with 
		# confidence. If no second choice, what is the null?
		# This is exactly what we were testing in old implementation but with 
		# fairly rough and ready scoring and no real alignmnet step. 
		# Under this Conception we would not try to score the "no wgd model"
		# - we can just assume it - and would instead compute all pairs and later 
		# compare sister_choice_1 to sister_choice_2. 
		# This sure gets easier if we can ignore the non-wgd model. 
		# We would score: 
		# 1) duplicates
		# 2) order consistent with the anc order [no circularity because YGOB built from other genomes]
		# 3) lack of gaps relative to anc
		#####
		# references... 
		# my PNAS paper. last SI.  
		# gavin paper 
		# felsenstein 
		#####

		my %score = (
		    'OHNO' => 0,
		    'GAP' => 0,
		    'CROSS' => 0,
		    'SAME' => 0
		    );	

		foreach my $k ( 0..$#clean ) {
		    my ($row,$old) = ($clean[$k], $clean[$k-1]);
		    #$self->throw unless $row->{$key0}->data( $keyX ) =~ /Anc_\d+\.(\d+)/;
		    #my $catch = $1;
		    #$self->throw unless $old->{$key0}->data( $keyX ) =~ /Anc_\d+\.(\d+)/;
		    #my $delta = abs($1 - $catch);
		    #my ($o_score, $g_score)=(0,0);    
		    #$o_score = $matrix{'OHNO'} if (defined $row->{$key1} && defined $row->{$key2} ? 1 : 0);
		    #$g_score = $matrix{'GAP'}*($delta-1) if 
		    #(defined $row->{$key1} && defined $old->{$key2} || defined $row->{$key2} && defined $old->{$key1});
		    #$score += $o_score + $g_score;		    		   

		    if ( $k==0 ) {
			$score{'OHNO'}++ if $row->{$key1} && $row->{$key2};
		    } else {
			if ( $row->{$key1} && $row->{$key2} ) { 			
			    $score{'OHNO'}++;
			} elsif ( ($row->{$key1} && $old->{$key1}) || ($row->{$key2} && $old->{$key2}) ) {
			    $score{'SAME'}++;
			} elsif ( ($row->{$key1} && $old->{$key2}) || ($row->{$key2} && $old->{$key1}) ) {
			    $score{'CROSS'}++;
			} else {}#$self->throw;}
		    }
		    
		    print {$fh} ($k,
				 (map { "$_:$score{$_}:".$score{$_}*$matrix{$_} } keys %matrix), '}', 
				 (map { "$_:".($row->{$_} =~ /\:\:/ ? $row->{$_}->data($keyX) : 'NA')} ($key1,$key2) )
		    ) if $args->{'-verbose'};			
		}
		
		my $score=0;
		map { $score+=$score{$_}*$matrix{$_} } grep {!/gap/i} keys %matrix; 
		$sisters{ $i.'.'.$j }={
		    G1 => $cand[$i],
		    G2 => $cand[$j],
		    L1 => scalar(@i_sort),
		    L2 => scalar(@j_sort),
		    SCR => $score
		};
	    } #j
	} #i
	next unless %sisters;
	
	########################
	# how do we do the statistics? 
	########################
	# do we need to do 2 tests?
	# are we better than second best option (ie can be reliably distinguished)?
	# are we better than null?
	# how do we exclude cases where the gene truly does not have an ohnolog, 
	# dut does have a duplicate elsewhere in the genome. 
	# in these cases there exists a well aligned region that we are not 
	# considering. how do we compare candidate to regions w/o ohnolog?
	
	my ($best, $second) = sort { $sisters{$b}->{SCR} <=> $sisters{$a}->{SCR} } keys %sisters;	
	my $bar = ($second ? $sisters{$second}->{SCR} : 8);
	#next unless $sister{$best} > $bar;
	print 
	    $anc, $sisters{$best}->{SCR}, $bar, 
	    $sisters{$best}->{G1}->name, 
	    $sisters{$best}->{G2}->name if $args->{'-verbose'};
	$sisters{$best}->{G1}->ohnolog( $sisters{$best}->{G2} ); # we set reciprocals automatically
    }
    close $fh;
    
    # 
    
    map { print $_,$count{$_} } sort keys %count if $args->{'-verbose'};

    map { $_->data($attr => 'delete') } $self->orfs;
    return $self;
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

=head2 report()

    Write various reports.

=cut 

sub report {
    my $self = shift;
    my $args = {@_};

    my @species = qw(Scer Spar Smik Skud Sbay);

    #open(my $ty, ">".$TIME."_repeats.fsa"); # fasta repeat files. for BLASTIng.  
    #open(my $hgt, ">".$TIME."_hgt.tab");    # HGT candidates 
    #open(my $deep, ">".$TIME."_nosgd.tab"); # no SGD hit. no YGOB hit. interesting? 

    my $ob = $self->ontology;
    
    ######################################
    # tRNAs
    ######################################
    
    my (%trna, %ygob);
    foreach my $g ( $self->iterate ) {
	foreach my $c ( $g->stream ) {
	    foreach my $o ( $c->stream ) {
		
		if ( $o->assign =~ /RNA/ && $o->gene =~ /\:/) {
		    push @{$trna{$o->gene}{$g->organism}}, $o->name;

		} elsif ( $o->assign =~ /REPEAT/ ) {
		    #print {$ty} $o->sequence(-molecule => 'dna', -format => 'fasta', -decorate => 1) ;
		    
		} elsif ( $o->evidence eq 'NCBI' && $o->data('NCBI') <= 1e-5 ) {
		    #$o->oliver(-fh => $hgt, -append =>[$o->data('NCBI')] );

		} elsif ( $o->assign eq 'REAL' && $o->evalue('ygob') <= 1e-5 ) {
		    push @{$ygob{$o->ygob}{$g->organism}}, $o;
		    push @{$names{$o->ygob}{$g->organism}}, ($o->ohnolog ? '*' : undef).$o->name;
		    
		} elsif ( $o->assign eq 'REAL' && $o->evidence ne 'YGOB' && 
			  $o->sgd !~ /^Y[A-P][LR]\d{3}/ && $o->gene !~ /^Y[A-P][LR]\d{3}/ ) {
		    #$o->oliver(-fh => $deep, -append => [$o->logscore('gene'), $o->gene] );
		}            
	    }
	}
    }
    
    ######################################
    # tRNAs
    ######################################

    my $name = $TIME."_tRNA.tab";
    print STDERR $name;
    open(my $fh, ">$name" );

    foreach my $tr ( keys %trna ) {
	print {$fh} ">$tr [".
	    join("|", ( map { scalar(@{$trna{$tr}{$_}}) || 0 } @species ) )
	    ."]"; 
	map { print {$fh} $_, (sort @{$trna{$tr}{$_}}) } @species;
    }
    
    ######################################
    # YGOB
    ######################################

    my $name = $TIME."_ygob.tab";
    print STDERR $name;
    open(my $fh, ">$name" );

    foreach my $yg ( grep {/\w/} keys %names ) {
	my @go = @{ $ob->{$yg}->{'P'} };
	my $go = join(';', @go);
	print {$fh} ">$yg [".
	    join("|", ( map { scalar(@{$names{$yg}{$_}}) || 0 } @species ) )
	    ."] $go"; 
	map { print {$fh} $_, (sort @{$names{$yg}{$_}}) } @species;
    }

    ######################################
    # Ohnologs
    ######################################

    my $name = $TIME."_ohnologs.tab";
    print STDERR $name;
    open(my $fh, ">$name" );
    print $fh qw(#count SGD YGOB Gene ExonsSGD Exons ExonsDelta Critique SGD);

    foreach my $yg ( keys %ygob ) {	
	next unless my @x = grep { $_->ohnolog } map { @{$ygob{$yg}{$_}} } @species;
	
	my @go = @{ $ob->{$yg}->{'P'} };
	my $go = join(';', @go);
	print {$fh} ">$yg [".
	    join("|", ( map { scalar(@{$ygob{$yg}{$_}}) || 0 } @species ) )
	    ."] $go"; 
        
	foreach my $sp ( @species ) {
	    print {$fh} $sp, 
	    map { ($_->ohnolog ? '*' : undef).$_->name.'/'.$_->gene } (sort { $a->up->id <=> $b->up->id } @{$ygob{$yg}{$sp}});
	}
    }

    ######################################
    # INTRONS
    ######################################
    
    my $name = $TIME."_introns.tab";
    print STDERR $name;
    open(my $fh, ">$name" );
    print $fh qw(#count SGD YGOB Gene ExonsSGD Exons ExonsDelta Critique SGD);

    foreach my $o ( grep {$_->orthogroup || $_->_debug} $self->orfs ) {    
        my @int = grep { $_->data('INTRONS')>0 } ($o,$o->orthogroup);
        my $sgd = $o->_debug;
        my $sgd_ex = ($sgd ? $sgd->stream : 1);
        next unless @int || $sgd_ex>1;
	
        my @intx = sort {$a <=> $b} map { $_->stream+0 } ($o,$o->orthogroup);
        my $m = $intx[ ($#intx%2==0 ? $#intx : ($#intx+1))/2 ];
        
	my $exons='NA';
	if( my $sgdx = $o->_debug ) {
	    $exons = join('/', (map {$_->length} $sgdx->stream));
	}
	
	
	print {$fh} ++$icount, ($o->identify), $sgd_ex, $m, ($sgd_ex-$m), sprintf("%.2f",$o->critique), $exons, 
	( map { join('/', $_->name, (map {$_->length} $_->stream) ) } ($o, $o->orthogroup) ),
	$o->_genericlinks;
    }
    
    return $self;
}


=head2 merge(-synteny => 2, -distance => 200, -align => .25)
    
    Goes through genome and looks for contigs that can be merged
    into super-contigs on basis of annotated genes and synteny
    information. Does all the hardwork of adjusting coords,
    fixing gene structures at the edges to be joined and removing
    old Contig objects from Genome object.

    It uses 5 tests to determine that a pair of HSPs are suitable
    to guide contig merger: 
    1) Same gene at each contig end
    2) No other hits to gene elsewhere in genome
    3) Synteny extends from gene
    4) Orientations make sense 
    5) The length of the 2 fragments < length of database hit

    -synteny : multifunctional. Loosely, it sets the number of genes
    either side of the break point that are considered in determining
    whether synteny is conserved across the break. 
    
    -distance : max distance of HSPs from the edge of their contigs.
    
    -align : max fraction of the shorter HSP that can align to the 
    longer when checking for homology between the 2 HSPS. A very
    liberal e-value threshold is used, making this test pretty 
    conservative. 

=cut 

sub scaffold {
    my $self = shift;
    my $args = {@_};

    $args->{'-restrict'} = [ grep {!/MITO|SPOM/} keys %HOMOLOGY] unless exists $args->{'-restrict'};
    $args->{'-intervening'} = 3 unless exists $args->{'-intervening'};
    $args->{'-overlap'} = 0.5 unless exists $args->{'-overlap'};

    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    
    ###############################################
    # 
    ###############################################

    my $dist = 7;
    my @context = (-self => 1, -repeat => -1, -trna => -1, -variant => 0);
    my $attr = '_SCAF_TEMP_VAR';
    my $attr2 = '_SCAF_TEMP_GENE';

    goto DP if -e 'juncFile.sto';

    #my @junctions;
    #use Tie::RefHash;
    my %hash;
    #tie %hash, 'Tie::RefHash::Nestable';
    
    ###############################################
    # do a little prep-work.can eliminate many comparisons
    # by comparing chr only ..
    ###############################################

    my @contigs=( grep { $_->orfs >= 2 } $self->stream );
    foreach my $c ( @contigs ) {
	foreach my $d (-1,+1) {
	    my $o = $c->down(-direction => $d);
	    map { $o->chromosome(-species=> $_, -distance => $dist, -self => 1) } @{$args->{'-restrict'}};
	}
    }

    ###############################################
    # do all by all comparison of contigs ends 
    ###############################################
    
    # we are interested in 2 criteria: 
    # 1. some definite evidence that A-B are connected. shared homolog. 
    # 2. synteny across the breaksite. 

  C1: for (my $i = 0; $i < $#contigs; $i++) {
      foreach my $ori1 (-1,1) { # orientation: current or flip 
	  my $t1 = $contigs[$i]->down(-direction => -1*$ori1); # we use RHS of contig 
	  #$t1 = $t1->$dir until ! $t1 || $t1->rank < 2; # no LDA -- require homology 
	  #next if $t1->rank <= 0; # REPEATS and tRNAs
	  #next if $t1->_terminal_dist2(-direction => $dir) > $args->{'-distance'};
	  #next unless my @t1 = $t1->chromosome(-species=> $key, -distance => $dist, -self => 1);
	  
	  #next unless $t1->ygob =~ /442/;

	C2: for (my $j = $i+1; $j <= $#contigs; $j++) {
	    foreach my $ori2 (-1,1) { 
		my $t2 = $contigs[$j]->down(-direction => $ori2); # we use LHS : 
		
		###############################################
		# do coarse filter. 
		###############################################

		my %sc;
		foreach my $sp ( @{$args->{'-restrict'}} ) {
		    # compare indices of genes with pre-computed CHR
		    next unless $t1->chromosome(-species => $sp) eq $t2->chromosome(-species => $sp);
		    next unless my @t1 = $t1->chromosome(-species=> $sp, -distance => $dist, -self => 1);
		    next unless my @t2 = $t2->chromosome(-species=> $sp, -distance => $dist, -self => 1);
		    my $delta = abs($t1[1]-$t2[1]);
		    $sc{ $sp } = [$t1[0], $t1[1], $t2[1]] if ( $t1[0] eq $t2[0] && $delta <= 3*$dist);
		}
		next unless %sc;

		###############################################
		# do fine filter 
		###############################################
		
		# we need pseudocontigs here to do before/after comparisons. 
		
		my $cc1 = $contigs[$i]->clone();
		map { $_->transfer(-to => $cc1, -warn => 0) } map { $_->clone(-data => 1) } 
		$t1->traverse(-distance => 3*$dist, -direction => (-1*$ori1 == 1 ? 'right' : 'left') );
		$cc1->index;
		my $ct1 = $cc1->down(-direction => $ori1);

		my $cc2 = $contigs[$j]->clone();
		map { $_->transfer(-to => $cc2, -warn => 0) } map { $_->clone(-data => 1) } 
		$t2->traverse(-distance => 3*$dist, -direction => ($ori2 == 1 ? 'right' : 'left') );
		$cc2->index;
		my $ct2 = $cc2->down(-direction => $ori2);

		# get sysnteny values _before_ merger

		my @syn1 = map { $_->synteny(-loss => 1, -distance => $dist, -restrict => $args->{'-restrict'}) } 
		$ct1->context( -distance => $dist, @context); 
		my @syn2 = map { $_->synteny(-loss => 1, -distance => $dist, -restrict => $args->{'-restrict'}) } 
		$ct2->context( -distance => $dist, @context );

		# this is the important bit #############
		# make fused contig to retest synteny 
		# we always fuse conitg1-3' <=> 5'-contig3 
		# we therefore flip contig1 if required and 
		# provide orientation guidance to fuse()
		# NB fuse expects an instruction: 1=current, -1=invert 
		$cc1->_invert_contig if $ori1 == -1;
		my ($fuse,$gap) = $cc1->fuse( -object => $cc2, -orientation => $ori2);
		$fuse->index;
		#########################################

		# recalc synteny _after_ merger so we can compare 
		# shame that we do not have statistics.
		
		my @syn3 = map { $_->synteny(-loss => 1, -distance => $dist, -restrict => $args->{'-restrict'}) } 
		$gap->context( -distance => $dist, @context );

		# 'score' below this is the key number 
		# we do a heuristic filter here but it is really 
		# used much later. 
		
		my ($s1,$s2,$s3);
		map { $s1+=$_ } @syn1;
		map { $s2+=$_ } @syn2;
		map { $s3+=$_ } @syn3;
		my $score = sprintf("%.1f",$s3 - $s1 - $s2);
		my ($min) = sort { $a <=> $b} ($s1,$s2,$s3);
		$fuse->DESTROY and next unless $score >= $min/2;
		
		# pump out some data...
		
		if ( $args->{'-verbose'} ) {
		    print "\n>>>".join('/',map { $contigs[$_]->id } ($i,$j)).": $score [$s1,$s2,$s3]";
		    if ( $args->{'-verbose'} >=2 ) {
			$t1->show( -window => $dist*2, -species=>[] );
			print undef, '>>', map {sprintf("%.1f",$_)} @syn1;
			$t2->show( -window => $dist*2, -species=>[] );
			print undef, '>>', map {sprintf("%.1f",$_)} @syn2;
			$gap->show( -window => $dist*2, -species=>[] );
			print undef, '>>', map {sprintf("%.1f",$_)} @syn3;
		    }
		}

		###############################################
		# choose bounding genes and test gene order 
		###############################################

		# we consider each species to find one with consistent 
		# gene order etc. we will use this for all further tests.

	      SPECIES: foreach my $sp ( grep { $sc{$_} }  @{$args->{'-restrict'}} ) {
		  next unless my ($chr,$lv,$rv)=@{$sc{$sp}}; ## how?
		  
		  # identify bounding gene limits  
		  
		  my ($min,$max) = sort {$a <=> $b} ($lv,$rv);
		  $min-=.5 if $min =~ /\./;
		  $max+=.5 if $max =~ /\./;
		  
		  # get orientation and bounding gene on each side. 
		  
		  my $ori=0;
		  my @bound;
		  foreach my $dir ('left','right') {

		      my @genes; 
		      foreach my $orf ( grep { $_->data($sp) } 
					$gap->context( -distance => $dist, -direction => $dir, @context ) ) {
			  my  ($lsp,$lchr,$lid) = 
			      &Annotation::Orf::_decompose_gene_name($orf->data($sp));
			  next unless $lchr eq $chr;
			  $orf->data($attr => $lid);
			  if ( $dir eq 'left' ) {
			      unshift @genes, $orf;
			  } else { push @genes, $orf; }
		      }
		      next SPECIES unless $#genes>=1;
		      push @bound,$genes[0]; # point gene
		      
		      # do "straight line" test  
		      # closest genes to gap are 0,1,2..
		      
		      my $cmp=1;
		      $cmp++ until ( ! exists $genes[$cmp] || $genes[0]->data($attr) != $genes[$cmp]->data($attr) );
		      next unless exists $genes[$cmp];
		      $ori += ( $genes[$cmp]->data($attr) > $genes[0]->data($attr) ? 1 : -1);
		      
		      # verbose
		      
		      #print "+$sp",$dir,$bound[-1]->name, $bound[-1]->data($sp), $cmp,
		      #$genes[$cmp]->name, $genes[$cmp]->data($sp), $ori if $args->{'-verbose'};		      
		  }

		  push @{$sc{$sp}}, @bound if $ori==0;
	      }

		###############################################
		# choose best species for detailed examination 
		###############################################	      		
		
		foreach my $sp ( keys %sc ) {
		    delete $sc{$sp} and next unless $#{$sc{$sp}}==4;
		    my  ($lsp,$lchr,$lid) = 
			&Annotation::Orf::_decompose_gene_name($sc{$sp}->[-2]->data($sp));
		    my  ($rsp,$rchr,$rid) = 
			&Annotation::Orf::_decompose_gene_name($sc{$sp}->[-1]->data($sp));
		    $self->throw unless $lchr && $rchr && $lchr eq $rchr;
		    push @{$sc{$sp}}, ( (sort {$a <=> $b} ($lid,$rid)), abs($lid-$rid)) 
			unless $lchr =~ /\D/;
		}
		$fuse->DESTROY and next unless keys %sc;
		
		###############################################
		# test all intervening genes in chosen species 
		###############################################

		my ($best) = sort { $sc{$a}->[-1] <=> $sc{$b}->[-1] } keys %sc;		
		$fuse->DESTROY and next if $sc{$best}->[-1] >= 10; ### HARD CODED
	
		my ($ch,$l1,$l2,$o1,$o2,$min,$max,$delta)=@{ $sc{$best} };
		my ($left,$right)= sort { $a->start <=> $b->start } ($o1,$o2); 
		my ($locus,$adjust) = $fuse->
		    _pseudocontig( 
			-start => $left->start, 
			-stop => $right->stop, 
		    );

		for my $i ( $min..$max ) {
		    my $gene = ($best eq 'YGOB' ? 'Anc' : $best).'_'.$ch.'.'.$i;
		    # 
		    my @cand = grep { $_->length > 60 } grep {defined} 	 ### HARD CODED 
		    (
		     $best eq 'YGOB' ? 
		     $locus->wise( -hmm => $gene, -verbose => 0, -safe => 1) :
		     $locus->exonerate( -protein => $gene, -verbose => 0, -safe => 1) 
		    );
		    # 
		    #map { $_->update(-initialize => ['EXONERATE.DNA.LOCAL']) } 
		    map { $_->transfer(-from => $locus, -to => $fuse, -warn => 0) }
		    map { $_->data($attr2 => $gene);$_} 
		    map { $_->adjust($adjust) } @cand if @cand;
		}
		
		###############################################
		# require gene model to span junction 
		###############################################
		
		my @span;
		my %frags;
		foreach my $x (grep {$_ ne $gap } $gap->context(-distance => $dist, @context) ) {
		    if ( $x->start <= $gap->start && $x->stop >= $gap->stop ) {
			push @span, $x->data($attr2);
			$x->oliver( -prepend => ['CROSS'] ) if $args->{'-verbose'};
		    } elsif ( $x->start <= $gap->start ) {
			push @{$frags{$x->data($best)}{'left'}},$x;
		    } else {
			push @{$frags{$x->data($best)}{'right'}},$x;
		    }
		}

		# 

		foreach my $hg ( keys %frags ) {
		    next unless exists $frags{$hg}{'left'} && $frags{$hg}{'right'};
		    my ($left) = sort { 
			$a->distance(-object => $gap, -bases => 1) <=> 
			    $b->distance(-object => $gap, -bases => 1)
		    } @{ $frags{$hg}{'left'} };
		    my ($right) = sort { 
			$a->distance(-object => $gap, -bases => 1) <=> 
			    $b->distance(-object => $gap, -bases => 1)
		    } @{ $frags{$hg}{'right'} };

		    #

		    next unless $left && $right;
		    next unless $left->strand == $right->strand;
		    next if $left->distance( -object => $right, -bases => 1) > 2000;  ### HARD CODED
		    #next if $left->distance( -object => $right) > $args->{'-intervening'};
		    
		    next unless my $homolog = ( 
			$hg !~ /^Anc_/ ? $hg : $left->homolog(			    
			    -pillar => $hg, 
			    -other => $right, 
			    -model => 'local'
			)
			); 
		    
		    my $olap = $left->exonerate2( 
			-object => $right, 
			-homolog => $homolog, 
			-model => 'local', 
			-return => 'overlap' 
			);
		    
		    map { $_->oliver(-prepend => ["FRAGS"], -append => [$olap,$hg,$homolog]) } ($left, $right) 
			if $args->{'-verbose'};
		    push @span, $hg if $olap >= $args->{'-overlap'};
		}

		###############################################
		# build scoring matrix and clean up 
		###############################################

		my $support='-none-';
		if ( @span ) {
		    my %count;
		    map { $count{$_}++ } @span;
		    ($support) = sort {$b <=> $a} @span; 
		    $hash{$t1->_internal_id.'.'.$t2->_internal_id} = 
		    {
			'SX' => [$s1,$s2,$s3],
			'SCR' => $score,
			'HOM' => $support
		    };
		}
		print $support, $contigs[$i]->id*$ori1, $contigs[$j]->id*$ori2, $t1->name, $t2->name, $score
		    if $args->{'-verbose'};
		
		$fuse->DESTROY; # fuse = cc1.cc2 
	    } # O2
	} # C2
      } # O1
  } # C1
    print if $args->{'-verbose'};

  DP:
    if (%hash) {
	nstore(\%hash, 'juncFile.sto');
    } else {
	%hash = %{ retrieve('juncFile.sto') };
    }

    ###############################################
    # select best pairs from matrix and merge contigs
    ###############################################

    my @fix;
    my %seen;
    foreach my $id (sort {$hash{$b}->{'SCR'} <=> $hash{$a}->{'SCR'}} keys %hash) {

	# has terminus been used already?

	my ($iid1,$iid2) = split/\./, $id;
	next if exists $seen{ $iid1 } || exists $seen{ $iid2 };
	map { $seen{ $_ }++ } ($iid1,$iid2);

	# get objects 

	$self->throw unless my $o1 = $self->find(-internal => $iid1);
	$self->throw unless my $o2 = $self->find(-internal => $iid2);
	$self->throw unless my $c1 = $o1->up;
	$self->throw unless my $c2 = $o2->up;

	# orientations 

	if ( $c1->down(-direction => 1) eq $o1 ) {
	    $c1->_invert_contig;
	} elsif ( $c1->down(-direction => -1) eq $o1 ) {
	} else { $self->throw;}
	if ( $c2->down(-direction => -1) eq $o2 ) {
	    $c2->_invert_contig;
	} elsif ( $c2->down(-direction => 1) eq $o2 ) {
	} else { $self->throw;}

	#
 
	if ($args->{'-verbose'}) {
	    $o1->show( -window => $dist*2, -species=>[] );
	    $o2->show( -window => $dist*2, -species=>[] );
	}

	print ++$DDD, $c1->id, $c2->id, $o1->name, $o2->name, $hash{$id}->{'SCR'};

	my ($fuse,$gap) = $c1->fuse(-object => $c2);
	$fuse->index;
	push @fix, [$gap, $hash{$id}->{'HOM'}];

	if ( $args->{'-verbose'}) {
	    $gap->show( -window => $dist*2, -species=>[] );
	    $gap->oliver(-prepend => [$hash{$id}->{'HOM'}], -append => [$hash{$id}->{'SCR'}, $gap->description] );
	    print;
	}
    }

    ###############################################
    # fix genes at contig borders 
    ###############################################

    foreach my $fix ( @fix ) {
	my ($gap,$gene)=@{$fix};
	my ($locus,$adjust) = $gap->up->
	    _pseudocontig( 
		-object => $gap,
		-extend => 12e3 # got really weird WISE errors here. results depends on input length. 
	    );
	
	next unless my @cand = grep { $_->length > 60 } grep {defined} 	 ### HARD CODED 
	(
	 $gene =~ /^Anc_/ ? 
	 $locus->wise( -hmm => $gene, -verbose => 0, -safe => 1) :
	 $locus->exonerate( -protein => $gene, -verbose => 0, -safe => 1) 
	);
	
	foreach my $cand ( @cand ) {
	    $cand->adjust($adjust);
	    $cand->transfer(-from => $locus, -to => $gap->up, -warn => 0);
	    $cand->update();
	    $cand->oliver if $args->{'-verbose'};
	}

	# use purge2() method to clean up overlaps ect. 
	# purge automatically runs merge() 
	# so fragments are automatically gathered 

	$gap->up->purge2(); 
    }

    return $self;
}

sub merge {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-synteny'} = 2 unless exists $args->{'-synteny'}; # this is used for LOTS 
    $args->{'-distance'} = 200 unless exists $args->{'-distance'};
    $args->{'-align'} = .25 unless exists $args->{'-align'};
    $args->{'-homology'} = 2 unless exists $args->{'-homology'};

    # we want to be able to use hash keys as objects (or maybe objects as keys..)
    
    use Tie::RefHash;
    undef my %hash;
    tie %hash, 'Tie::RefHash::Nestable';

    ###############################################
    # index contigs
    ###############################################

   my %chk;
  C1: foreach my $c1 ($self->stream) {
      next unless $c1->orfs(-restrict => ['REAL', 'HSP'])
          >= $args->{'-synteny'}+1; 
      
      foreach my $dir ('left', 'right') {
          my $t1 = $c1->down(-direction => $dir);
          my $s1 = $t1->synteny(
              -distance => $args->{'-synteny'}, 
              -difference => 5,
              -spanning => 0,
              );
          
          next unless $s1 >= $args->{'-synteny'};       
          next unless $t1->assign =~ /REAL|HSP/;
          next unless $t1->_terminal_dist <= $args->{'-distance'};
          
	C2: foreach my $c2 ($self->stream) {
	    next C1 if $c2 eq $c1 || 
		exists $chk{$c1->_internal_id}{$c2->_internal_id};
	    $chk{$c2->_internal_id}{$c1->_internal_id}=1;
	    next unless $c2->orfs(-restrict => ['REAL', 'HSP']) 
		>= $args->{'-synteny'}+1;
	    
	    foreach my $dir ('left', 'right') {
		my $t2 = $c2->down(-direction => $dir);	
		my $s2 = $t2->synteny(
		    -distance => $args->{'-synteny'}, 
		    -difference => 5,
		    -spanning => 0,
		    );

		next unless $s2 >= $args->{'-synteny'};
		next unless $t2->assign =~ /REAL|HSP/;	
		next unless $t2->_terminal_dist <= $args->{'-distance'};
		
		###############################################
		# detailed comparison of candidates + add to matrix 
		###############################################

		$count{'CAND'}++;

		# 1) homology ? 
		
		next unless 
		    $t1->homology(-object => $t2) >= $args->{'-homology'} ||
		    ($t1->ygob && ( $t1->ygob eq $t2->ygob ) );
		$count{'HOM'}++;

		next if $t1->overlap(
		    -object => $t2,
		    -compare => 'seq') > $args->{'-align'};

		$count{'OLAP'}++;

		# 2) orientation ?	
		
		my ($d1, $o1) = ('left', 1);	
		if (! $t1->left) {
		    $d1 = 'right';
		    $o1 = -1;
		} 					
		my ($d2, $o2) = ('left', 1);	
		if (! $t2->left) {
		    $d2 = 'right';
		    $o2 = -1;
		} 		
		next unless $o1*$o2*$t1->strand*$t2->strand == -1;
		$count{'ORI'}++;

		# 3) spanning synteny ?
		
		my @g1 = grep { $_ ne $t1 } $t1->traverse(
		    -distance => $args->{'-synteny'}, -direction => $d1);
		my @g2 = grep { $_ ne $t2 } $t2->traverse(
		    -distance => $args->{'-synteny'}, -direction => $d2);
		
		my (%spent, $score);
	      G1: foreach my $i (reverse(@g1)) {
		G2: foreach my $j (reverse(@g2)) {			
		    next G2 unless $i->synteny(
			-object => $j,
			-difference => $args->{'-synteny'}*2,
			);		    
		    $score++ unless exists $spent{$j->_internal_id};
		    $spent{$j->_internal_id}=1;
		    next G1; 
		}
	      }
		$count{'SYN'}++;

		# 4) length checks: A + B < HIT 
		# 5) genomic uniqueness: number of hits to X in genome == 2  

		# build scoring matrix 
		
		$hash{$t1}{$t2} = {
		    O1 => $o1,
		    O2 => $o2,
		    SCR => $score
		};
	    }
	}
      }
  }

    map { print $_, $count{$_} } keys %count;

    ###############################################
    # select best pairs from matrix and merge 
    ###############################################

    open(my $fh, ">".$self->organism.'.joins');

    foreach my $t1 (sort {$b->up->orfs('REAL') <=> 
			      $a->up->orfs('REAL')} keys %hash) {
	my ($t2, @excess) = sort {$hash{$t1}{$b}->{SCR} <=> 
				      $hash{$t1}{$a}->{SCR}} keys %{$hash{$t1}};
	next if @excess; # 
	next unless $t1 && $t2;
	next unless $t1->up && $t2->up;

	########################################
	$t1->up->output(-fh => $fh);
	print $fh "-";
	$t2->up->output(-fh => $fh);
	print $fh "-";
	########################################
	
	my $o1 = $hash{$t1}{$t2}->{O1};
	my $o2 = $hash{$t1}{$t2}->{O2};
	my $c1 = $t1->up;
	my $c2 = $t2->up;

	# orient the contigs for the merger 

	unless ($o1 == 1) {
	    $c1->sequence($c1->revcomp);
	    foreach my $o ($c1->stream) {
		$o->_invert_coords;
		$o->strand($o->strand*-1);
		foreach my $ex ($o->stream) {
		    $ex->strand($o->strand);
		}
		$self->throw unless $o->translatable;
	    }
	}
		
	unless ($o2 == -1) {
	    $c2->sequence($c2->revcomp);
	    foreach my $o ($c2->stream) {
		$o->_invert_coords;
		$o->strand($o->strand*-1);
		foreach my $ex ($o->stream) {
		    $ex->strand($o->strand);
		}
		$self->throw unless $o->translatable;
	    }	
	}

	# merge two contigs 

	my $offset = $c1->length + 100;
	my $gap = 'N' x 100;
	my $new = $c1->sequence.$gap.$c2->sequence;
	$c1->sequence($new);
		
	foreach my $o ($c2->stream) {
	    $c2->remove(-object => $o, -force => 1);
	    foreach my $ex ($o->stream) {	
		$ex->start($ex->start + $offset);
		$ex->stop($ex->stop + $offset);	
	    }
	    $c1->add(-object => $o);
	    $self->throw unless $o->translatable;
	}
	$c1->index;	
	$self->remove(-object => $c2, -force => 1);		
	$c2->DESTROY;	

	# place a gap object - this is a bit cumbersome..

	my $gap = Annotation::Orf
	    ->new(
	    START => $offset-99,
	    STOP => $offset,
	    STRAND => 0
	    );

	my $ex = $gap->exons(-query => 'first');
	$ex->introns(-direction => 'left', -new => $INFINITY);
	$ex->introns(-direction => 'right', -new => $INFINITY);

	$c1->add(-object => $gap);
	$gap->data('NNNN' => $INFINITY);
	$gap->evaluate(-structure => 0, -validate => 0);
	$c1->index;	
	
	# the gene that started it all... 
	
	$t1->merge(-object => $t2, -force => 1);
	$c1->index;		

	# store the contig id history 
	
	my $scaf_history = $c1->id.'+'.$c2->id;
	$scaf_history =~ s/$c1->id/$c1->data/ if $c1->data;
	$scaf_history =~ s/$c2->id/$c2->data/ if $c2->data;
	$c1->data($scaf_history);
	
	# 

	$c1->output(-fh => $fh);
	print $fh '*';
    }
    
    close($fh);
    return $self;
}


#########################################
# subroutines : gene search 
#########################################

=head2 recover(-intergenic => 90, -aa => [aa_db.fsa, 1e-3])
    
    Similar to homology() method but runs on intergenics not
    ORFs. The purpose is to recover genes with good homology     
    that were missed by the initial search for ORFs because
    the minimum length limit was too retrictive. 
        
    -intergenic sets a new limit on the size of interenic
    to be BLASTed. 

=cut

sub recover {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-intergenic'} = 90 unless exists $args->{'-intergenic'};	
    $args->{'-rank'} = 1.5 unless exists $args->{'-rank'};	

    foreach my $h (grep { ref($args->{$_}) eq 'ARRAY' }  keys %{$args}) {	
	my $homology = uc($h);
	$homology =~ s/^-//;
	my ($db, $e) = @{$args->{$h}};
	
	# make sure database ok
	
	if (!$db) {
	    $self->warn("No db supplied: $db ($homology)");next;
	} elsif (ref($db) eq ref($self)) {
	    $db = $self->_write_temp_seq(-seq => [$db->orfs]);
	} 
	$self->throw unless $self->_valid_file(-file => $db);
	
	map { $_->evaluate; } $self->orfs;

	# make intergenics
	
	return $self unless my @inters = 
	    # grep { $_->left->assign ne 'GAP' || $_->right->assign ne 'GAP' }  
	    $self->orfs(
		-intergenic => $args->{'-intergenic'}, # specify min size 
		-rank => $args->{'-rank'}, # w/ defult require homology or better 
		-noncoding => 1            # take account of tRNAs
	    );

	print 'Gaps',scalar(@inters) if $args->{'-verbose'};

	# run BLAST 
	
	next unless my $blast = 
	    $self->blast(
		-program => 'blastx',
		-qseq => \@inters,
		-qmol => 'dna',
		-db => $db,
		-formatdb => 1,
		-cpus => 4
	    );
	
	# process BLAST report. get multiple hits from each intergenic. 
	# returns new NR multi-exonic gene models 
	
	my @new = $self->_blast2orfs(
	    -blast => $blast,			
	    -queries => \@inters,
	    -evalue => $e # we use the evalue here to prune results 
	    );
	
	foreach my $hit ( @new ) {	   
	    # get the pillar homolog... first get HMMER data...
	    $hit->update(-blast => undef, -exonerate => undef, -hsp => undef); # get revised BLAST etc 
	    $hit->pillar();
	    my $hom = $hit->homolog(); # -fast not possible 
	    # try fix the structure...
	    $hit->reoptimise( -reference => $hom ) if $hom;
	    # full blast....
	    $hit->update(); # get revised BLAST etc 
	    $hit->oliver(-creator => 1) if $args->{'-verbose'};
	}

	# ensure each contig is properly organized 
	
	map { $_->index } $self->stream;  
    }

    return $self;
}

=head2 seek() 

    In this method we actively look for ('seek') all YGOB Anc genes that are
    not represented by at least 1 syntenic famliy member in the genomic locations 
    where they are expected to reside (eg Anc_2.114). There are 2 run modes that 
    differ in how the candidate mising genes and locations in which to seek them
    are compiled. 

    Runmode 1: 
    We proceed in 3 steps. First, we compile a list of all (potentially)
    missing genes. Next we gather and search the syntenic locations in which 
    they are expected to reside. Finally, we compare to existing models to deal
    with complications arising from mis-assignmnet (ie present but labelled 
    as another gene) or from tandems etc.

    Example: 
    In the case of Anc_2.114, find Anc_2.112 .. Anc_2.116 then excise regions 
    around all the genes left and right of these (2.111 -> 2.117). We 
    annotate the 10Kb region using a query dependent minimum ORF size 
    equal to 2/3 lenght of the query. In this case 34*2. We seek to 
    recover the gene multiple times and then compare to all the genes in 
    the region to see if we alrady detected it but had failed to correctly
    identify it. 

    In general, we find a match in the region for small genes. Where the query 
    is long, it tends to be legitametly missing. But may be elsewhere. 

    Runmode 2: 

=cut 

sub seek {
    my $self = shift;
    my $args = {@_};

    return $self->_seek(@_) if $args->{'-orthogroups'} || $args->{'-og'};

    $self->throw("Hardcoded-- fix that!") if exists $args->{'-synteny'};

    $args->{'-db'} = $ENV{'YGOB_HMMER3_DB'} unless $args->{'-db'};
    # 
    $args->{'-window'} = 5 unless exists $args->{'-window'};    # number of orfs around target genes 
    $args->{'-extend'} = 5000 unless exists $args->{'-extend'}; # bp around each orf 
    $args->{'-length'} = 165 unless exists $args->{'-length'};  # models longer than this get weaker scrutiny 
    # 
    $args->{'-debug'} = undef unless exists  $args->{'-debug'};
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    $args->{'-fast'} = 1 unless exists $args->{'-fast'};
    
    $self->throw if $args->{'-slow'} || ! $args->{'-fast'};

    ##############################################
    # collect all required files 
    ##############################################

    $self->throw("$args->{'-db'}") unless -e $args->{'-db'};
    my $dbdir = $args->{'-db'};
    $dbdir =~ s/ygob\.hmm$/YGOB\//; 
    $self->throw unless -e $dbdir;
    my @hmms =  <$dbdir/*h3m>;
    $self->throw unless $#hmms >= 1000;

    ##############################################
    # prep some vars 
    ##############################################

    my $attr = '_TEMP_SYN_HOM';
    my $attr2 = '_TEMP_NEW_GENE';
    foreach my $o ( $self->orfs ) {
	$o->data($attr => undef);
	$o->data($attr2 => undef);
    }
    
    my $expect = ($self->wgd == 1 ? 2 : 1);

    ##############################################
    # collect all Anc models and recorded hits 
    ##############################################

    my %max;
    foreach my $model (@hmms) {
	$model =~ /Anc_(\d+)\.(\d+)/;
	$max{ $1 } = $2 if $2 > $max{$1}; 
    }
    #map { print "$_\t$max{$_}"; } keys %max;

    my %ygob;
    foreach my $o ( grep { $_->evalue('ygob') < 1e-5 } $self->orfs ) {
	$o->ygob  =~ /Anc_(\d+)\.(\d+)/;
	push @{$ygob{ $o->ygob }}, $o;
	$o->data($attr => $o->syntenic+0); # flag syntenic homologs for later
	# +0 required to get teh scalar valued return. 
	# syntenic has hard coded synteny test values
	# -- should be replaced with an adaptive set that
	# is based on N50 G50 etc
	#$o->synteny(-hyper => 1,-distance => 5,-restrict => ['YGOB']) >= 5
	#$o->synteny(-logodds => 1,-distance => 5) >= 25;
    }
 
    ##############################################
    # walk along anc chromosomes + identify ancs with 0(/1) syn hom. 
    ##############################################

    my %ranges;
    my $counter;
    foreach my $chr (sort {$a <=> $b} keys %max ) { # CHR 1.. 8
	for my $pos (1..$max{$chr}) {
	    my $name = 'Anc_'.$chr.'.'.$pos;
	    next unless ! $args->{'-debug'} || $name eq $args->{'-debug'};
	  
	    # how many believable hits? we make assumptions based on 
	    # pre/post WGD status of the organism 
	    
	    my $hitcount = # must be syntenic and non-fragments. esp low N50 genomes.
		scalar( grep { $_->partial > -1 } grep {$_->data($attr)==1} @{$ygob{$name}} );
	    next if $hitcount >= $expect;

	    # look up HMM info and accept longer hits without question 
	    # but we want to be sure about short. <165AA : ~700 HMMs / 15%. 	  
	    # this should be replaced by a more thoughtful subroutine... 

	    next unless my $pillar = $YGOB->access($name);
	    next unless my $expLen = $pillar->{'MEAN'}/2; # mean length in pillar / 2
	    next if ($expLen > $args->{'-length'} && $hitcount >= 1 ); # we get the long ones OK 
	    
	    # look up genes that hit neighbours of target as path to the 
	    # expected syntenic locations. include some buffer. 
	    
	    #print ">$pos";
	    my %regions;
	    my $from = ( $pos-$args->{'-window'} < 1 ? 1 : $pos-$args->{'-window'});
	    my $to = ( $pos+$args->{'-window'} > $max{$chr} ? $max{$chr} : $pos+$args->{'-window'});	
	    foreach my $q (grep {$_ != $pos} ( $from .. $to ) ) { 
		#print 'Anc_'.$chr.".$q", (map {$_->name} @{$ygob{ 'Anc_'.$chr.".$q" }});
		map { push @{$regions{$_->up->id}}, (grep {defined} ($_->left, $_, $_->right)) } 
		@{$ygob{ 'Anc_'.$chr.".$q" }};
	    }

	    # parse collections of genes into genomic regions based on syntenic proximity

	    foreach my $chr ( keys %regions ) {
		my @sort = sort {$a->start <=> $b->start} @{$regions{$chr}};
		$counter++;
		foreach my $x (@sort) {
		    my $last = $ranges{$counter}{'ORFS'}->[-1];
		    next if $last eq $x;
		    $counter++ if ( $last && $x->distance(-object => $last) > 5 );
		    push @{$ranges{$counter}{'ORFS'}}, $x;
		    $ranges{$counter}{'X'}=$hitcount;
		    $ranges{$counter}{'ID'}=$name;
		}
	    }	  
	}	
    }

    ##############################################
    # process ranges:
    # remove those that already have a syntenic homolog 
    # compute ranges for other regions 
    ##############################################

    print STDERR scalar(keys %ranges);
  RANGE:foreach my $i (sort {$a <=> $b}  keys %ranges ) {
      my @sort = sort {$a->start <=> $b->start} @{$ranges{$i}{ORFS}};
      my @test = ($#sort==0 ? @sort : ( $sort[0],$sort[0]->intervening($sort[-1]),$sort[-1]) );

      foreach my $int ( @test ) { 
	  delete $ranges{$i} and next RANGE if 
	      (grep { $int eq $_ } grep { $_->data($attr)==1 } @{$ygob{ $ranges{$i}{'ID'} }}); # 
      }
      
      $ranges{$i}{'START'} = $ranges{$i}{ORFS}->[0]->start;
      $ranges{$i}{'STOP'} = $ranges{$i}{ORFS}->[-1]->stop;
      $ranges{$i}{'CONTIG'} = $ranges{$i}{ORFS}->[0]->up;
  }
    print STDERR scalar(keys %ranges);

    ##############################################
    # This is the major step. 
    # 1. search for missing Anc genes using $contig->search()
    # 2. group based on synteny 
    # 3. validate that the proposed YGOB is correct 
    ##############################################

    my $oldhmm; # just for output 
    foreach my $i (sort {$a <=> $b}  keys %ranges ) {
	my $contig = $ranges{$i}{CONTIG};
	
	##############################################
	# use WISE to search in each location 
	# much faster than current version of search() 
	# -> very large number of pointless HMMER calls 
	# we optinally run search at a later stage on a 
	# much tighter region. 
	##############################################	

	my @cand;
	$contig->index;	
	@cand = $contig->search(
	    -start => $ranges{$i}{'START'},
	    -stop => $ranges{$i}{'STOP'},
	    -extend => $args->{'-extend'},
	    -hmm => $ranges{$i}{'ID'},
	    -safe => 0, # if off, we add to calling contig 
	    -fast => $args->{'-fast'}
	    );
	$contig->index;
	
	##############################################
	# output?
	##############################################
	
	if ( $args->{'-verbose'} >=1 ) {
	    print "\n>>>$ranges{$i}{ID} [$ranges{$i}{'HMMLEN'} bp, ".scalar(@{$ygob{$ranges{$i}{ID}}})."/$ranges{$i}{'X'}]"
		if $oldhmm ne $ranges{$i}{ID};
	    print '>'.$ranges{$i}{'START'}, $ranges{$i}{'STOP'}, 
	    (map { $_->name.'|'.$_->ygob } grep {$_->up} map {$ranges{$i}{ORFS}->[$_]} (0,-1)), (@cand+0);
	    $oldhmm = $ranges{$i}{ID};
	}
	next unless @cand;
	
	##############################################
	# do they overlap existing models? or are they novel? 
	##############################################

	# gather everything to compare old and new 
	# cluster using commoncodon method 

	map {$_->data($attr2 => 1) } @cand; # label candidates as such 
	
	my $clusters = 
	    $contig->cluster( 
		#-object => [ grep {$_->rank > 0} @cand ], # added to stream. no pseudos etc 
		-param => $NONZERO, 
		-frame => 1, 
		-start => $ranges{$i}{'START'},
		-stop => $ranges{$i}{'STOP'}
	    );

	##############################################
	# 1. reduce the complexity by choosing 1 new model and
	# one old model at each locus (=cluster).
	# 2. compare old and new if both exist. 
	# 3. validate the YGOB identity 
	##############################################
	
	my %uniq;
	foreach my $cl ( grep {!$uniq{$_}++ } grep {$#{$_}>=0} values %{$clusters} ) {
	    my @old = grep { ! $_->data($attr2) } @{$cl}; 
	    my @new = grep { $_->data($attr2) } @{$cl}; 
	    
	    ##############################################
	    # what are possibilities? new:old
	    # 1:1, 0:0, 1:0, 0:1, N:N. 0:0 does not exist. 
	    # convert complex cases to simple [01]:[01] scneario 
	    ##############################################

	    $contig->simplify( -object => \@old ) if @old; # these should not arise?  
	    $contig->simplify( -object => \@new ) if @new; # best to get this sorted out..
	    my ($old) = grep {$_->up} grep {defined} @old;
	    my ($new) = grep {$_->up} grep {defined} @new;

	    ##############################################
	    # use update to get the required data only: BLAST, synteny 
	    ############################################## 

	    if ( $new ) {
		$new->update( -hmmer => undef, -exonerate => undef, -hsp => undef ); # run BLAST/synteny only 
		$new->pillar( -force => 1 ); # sets PILLAR 
		if ( my $hom = $new->homolog ) { # sets HOMOLOG 
		    # exonerate better than WISE? search() useful? 
		    # if there is an old gene we do not bother running search()
		    # if there is not, then we run with lenient params 
		    $new->reoptimise( -hmm => undef, -orfmin => ($old ? undef : $minlen), -reference => $hom );
		}
	    }

	    ###############################################
	    # 2 choices. ORF model. YGOB label. 
	    # for former we choose least work. 
	    # for latter we compare either NEW to OLD, or NEW to NULL.  
	    ##############################################
	    
	    my ($tag,$rem,$keep,$delta);
	    if ( $old && $new && $old->overlap(-object => $new) ) {
		# we test the overlap condition since sometimes running 
		# reoptimise above improves the model and removes the 
		# overlap. in this case we can just accept. 
		$tag = 'old';		

		($keep, $delta) = # select the best model. choose() will fast return if identical 
		    $old->choose(-object => $new, -reference => 1); # consensus reference 
		$rem = ( $keep eq $old || $delta == 0 || $delta == $INFINITY ? $new : $old ); # which?

	    } elsif ( $new ) { # simple case 
		$tag = 'new';

		$keep = $new; # be consistent with above 		
		$rem = $new->clone(); # _correct_ YGOB hit? clone to do unbiased HMMER
		$rem->update( -blast => undef, -exonerate => undef, -hsp => undef ); # run HMMER 
		
	    } else {next;} # we do not care 

	    ##############################################
	    # use synteny to choose the 'true' YGOB homolog. 
	    # compare the forced YGOB hit and the "free" YGOB hit w/ synteny.
	    ##############################################
	    
	    unless ( $keep->ygob eq $rem->ygob ) {
		#map { print $_->name, $_->hypergob, $_->ygob, $_->score('ygob'), $_->evalue('ygob'),$_->data($attr2) } ($keep,$rem);
		my ($gob) = sort {$b->hypergob <=> $a->hypergob} ($keep, $rem);
		my $hit = {
		    HIT => $gob->ygob,
		    EVALUE => $gob->evalue('ygob'),
		    SCORE => $gob->score('ygob')
		}; # gene was present just not recognized 
		$keep->accept( 'YGOB' => $hit ); # impose preferred YGOB data 
	    }

	    ##############################################
	    # finally process the model to be kept and the one to be DESTROY()ed 
	    ##############################################

	    # be careful here. in one condition defined above 
	    # we want to keep the new AND the old. when removing
	    # remove only $rem. be specific.
	    $contig->remove( -object => $rem, -warn => 0 ) 
		if $old && $tag eq 'old'; # **NB** clones not added to contig()
	    $rem->DESTROY;
	    
	    # finally, we fill out missing exonerate data since we know we are keeping 
	    $keep->update( -blast => undef, -hmmer => undef ) if $keep eq $new; # run exonerate/HSP (+ syn again)
	    $keep->oliver(-prepend => [$tag],-append => [$keep->syntenic, $keep->data($attr2)]); 
	    $keep->data($attr2 => undef);
	    $count{$tag}++;
	} # clusters 
    } # range 

    map { print STDERR "$_:$count{$_}" } keys %count if $args->{'-verbose'};

    ###########################################
    # do some cleanup 
    ###########################################

    foreach my $o ( $self->orfs ) {
	$o->data($attr => 'delete');
	$o->data($attr2 => 'delete');
    }
    map { $_->index } $self->stream;

    $self->summarize if $args->{'-verbose'} >= 2;

    return $self;
}

=head2 _seek(-min, -max) 
    
    Use seek() with runmode 2. 
    
=cut 

sub _seek {
    my $self = shift;
    my $args = {@_};

    $args->{'-min'} = 120 unless exists $args->{'-min'};
    $args->{'-max'} = 12000 unless exists $args->{'-max'};
    # 
    $args->{'-extend'} = 200 unless exists $args->{'-extend'};
    # 
    $args->{'-debug'} = undef unless exists  $args->{'-debug'};
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    $args->{'-fast'} = 1 unless exists $args->{'-fast'};
    
    $self->throw if $args->{'-slow'} || ! $args->{'-fast'};
    $self->throw unless $args->{'-max'} && $args->{'-min'};
    #$self->throw unless $args->{'-hmm'};

    ##############################################
    # prep some vars 
    ##############################################
    
    my @keys = (uc($self->organism), $self->bound);
    print scalar( grep { $_->ogid } $self->orfs );# if $args->{'-verbose'};
    
    ##############################################
    # we cycle through every gene-gene pair in the ref 
    # and test whether anything is missing in any species
    ##############################################
    
    my @success;
    foreach my $c ( $self->stream ) {
	$c->index;

      ORF:foreach my $o ( grep { $_->orthogroup } $c->stream ) {		
	  #next unless $o->name eq $args->{'-debug'};
	  next unless my $l = $o->neighbour(-direction => 'left', -orthogroup => 1);
	  #last unless my $r = $o->neighbour(-direction => 'right', -orthogroup => 1);

	  ##############################################	  
	  # gather all relevant genes 
	  ##############################################
	  
	  my $gather = $o->bracket( 
	      -object => $l, 
	      -max => $args->{'-max'}, 
	      -min => $args->{'-min'},
	      -og => -1
	      );
	  my %hash = %{ $gather };
	  
	  ##############################################	  
	  # can/should we do something? 
	  ##############################################
	  
	  next unless scalar(keys %hash) == scalar(@keys);
	  my $sum=0;
	  map { $sum+=scalar(@{$hash{$_}}) } keys %hash;
	  next if $sum==0;
	  next if $sum > 2*(scalar keys %hash); # too complex 
	  #print map { "$_:".scalar(@{$hash{$_}}) } keys %hash;

	  ##############################################
	  # compose matrix of missing genes 
	  ##############################################
	  
	  my %anc;
	  foreach my $org ( keys %hash ) {
	      foreach my $o ( @{$hash{$org}} ) {
		  if ( exists $anc{$o->ygob}{$org} ) {
		      my ($low) = sort {$a->logscore('ygob') <=> $b->logscore('ygob')} ($o, $anc{$o->ygob}{$org});
		      $anc{$o->ygob}{$org}=$low;
		  } else {
		      $anc{$o->ygob}{$org}=$o;
		  }
	      }
	  }
	  
	  ##############################################
	  # search for missing genes in relevant locus 
	  ##############################################
	  
	  foreach my $anc ( grep {$_}  keys %anc ) {
	      next unless scalar( %{$anc{$anc}} ) >= 2;
	      next if scalar( keys %{$anc{$anc}} ) == scalar(@keys);	      
	      print $anc, keys %{$anc{$anc}} if $args->{'-verbose'};
	      
	      #########################################
	      # search in each species 
	      #########################################
	      
	      foreach my $sp ( map {lc($_)} grep { ! exists $anc{$anc}{uc($_)} } @keys ) {
		  my ($start, $stop, $contig);
		  if ( $sp eq lc($self->organism) ) {
		      $start = $l->stop;
		      $stop = $o->start;
		      $contig = $o->up;
		  } else {
		      my $ortho = $o->$sp;
		      my $lortho = $l->$sp;		    
		      ($start, $stop) = ($ortho->index < $lortho->index ? 
					 ($ortho->stop, $lortho->start) : 
					 ($lortho->stop, $ortho->start) );
		      $contig = $ortho->up;
		  }		    
		  my $d = $stop - $start;
		  next unless $d <= $args->{'-max'} && $d >  $args->{'-min'};
		  #print $anc, $sp, $d, values %{$args};

		  #########################################
		  # region is elligible. run search and report top result 
		  #########################################
		  
		  $contig->index;	  
		  next unless my ($best) = 
		      $contig->search(
			  -start => $start, 
			  -stop => $stop, 
			  -extend => $args->{'-extend'},
			  # 
			  -hmm => $anc, # hmm name 
			  -top => 1,    # only return single best hit 
			  # 
			  -fast => 1,   # accelerated WISE search 
			  -safe => 0
		      );

		  # does this gene model already exist? 
		  # we do not try to do all the work that seek() does. 

		  $contig->index;		  
		  if ( $best->variants ) {
		      $contig->remove( -object => $best, -warn => 0 );
		      $best->DESTROY;
		      next;
		  }
		  
		  $anc{$anc}{uc($sp)}=$best;
		  $best->output( -prepend => ['NEW']) if $args->{'-verbose'};
	      }#sp 
	      next unless scalar( grep {/\w/} keys %{$anc{$anc}} ) == scalar(@keys);
	      
	      # we completed orthogroup. build data structure. 
	      
	      $self->throw unless my $ref = $anc{$anc}{uc($self->organism)};
	      $ref->_define_orthogroup(-object => [map { $anc{$anc}{uc($_)} } @keys[1..$#keys]], -verbose => 1);
	  }#anc
      }#$o -- focus gene 
    }
    
    #########################################
    # go back and do some clean up 
    #########################################
    
    foreach my $g ( $self, map { $self->bound($_) } $self->bound ) {
	map { $_->index } $g->stream; 
    }
    
    print scalar( grep { $_->ogid } $self->orfs );
    
    return $self;    
}

#########################################
# subroutines : wrappers 
#########################################

=head2 homology(-aa => [aa_db.fsa, 1e-5], -ty => [ty_db.fsa, 1e-5], etc)

    This method uses nonstandard argument structure and is closely tied
    to the %EVIDENCE hash defined in GlobalVars.pm

    In brief, each argument (-aa, -ltr etc) must correspond to a piece 
    of evidence in the %EVIDENCE hash and the specifc task of this 
    method is to allow gathering of evidence from seqeunce databases.

    e.g. -aa corresponds to the entry AA in the %EVIDENCE hash 
    and is associated with the *inference* 'REAL' -- a real 
    protein coding gene. 

    when $genome->homology(-aa => [aa_db.fsa, 1e-5]) is called
    every ORF in $genome is blasted against aa_db.fsa (presumably
    a db of known yeast genes) and best-hits <= 1e-5 are added to 
    the evidence stack for that ORF with the code AA. 

    This can be done for many different databases and when
    $orf->evaluate is called later, all the evidence will be 
    considered and an inference made (based on the evidence
    hierarchy %EVIDENCE hash) on what kind of sequence this
    should be assigned as. 

    As an example, $genome->homology(-aa => [aa_db.fsa, 1e-5])
    may cause an ORF to be assigned as a REAL (protein coding gene)
    but a subsequent call $genome->homology(-ty => [ty_db.fsa, 1e-5])
    may result in evidence for TY homology and lead to the inference 
    'REPEAT'. 

    The evidence code 'NCBI' (eg -ncbi => ['nr', 1-e10]) causes slightly
    differnet behaviour: The BLAST query is sent to NCBI nr DB using
    blastcl3. Most importantly however, not all genes are used as 
    queries. Only those with no AA hit (local prtein DB) and no 
    YGOBHMM (local YGOB HMMER DB) hit are sent to NCBI. NCBI is slow.
    1e-10 ($e) is the threshold for deciding what should be sent to NCBI 
    AND the threshold for what should be believed. Eg with 1e-50 
    we adopt a BLAST many, believe few policy. 
    To limit number of queries, -ncbi is always triggered last.
    
=cut

sub homology {
    my $self = shift;
    my $args = {@_};

    # NCBI uses a TOE (time of execution) system to execute searches.
    # the first submission gets TOE = now and each subsequent sub 
    # gets a 60s penalty. the 10th sub has a TOE of 540s. 
    # if the queue time is <60s it is better to submit 
    # individual queries and avoid the TOE penalty for multiple
    # simultaneous submissions. 
    
    my $local_batch_size = 1000;
    my $ncbi_batch_size = 1;
    my $min_ncbi_length = 150; # HGTs tend to be small but lower too annoying 

    # NCBI last if present 

    my @order = grep {!/ncbi/} keys %{$args};
    push @order, grep {/ncbi/} keys %{$args};

    # iterate 

    foreach my $h (@order) {	
	my $homology = uc($h);
	$homology =~ s/^-//;

	# look at args 
	
	next unless ref($args->{$h}) eq 'ARRAY'; # allow for non-DB arguments 
	my ($db, $e, $prog) = @{$args->{$h}};				
	$prog = 'blastp' unless $prog;
	$self->_set_eval(-evidence => $homology, -set => $e) if defined $e;
        my $qmol = ($prog =~ /^t/ || $prog =~ /p$/ ? 'aa' : 'dna');
        my $dbmol = ($prog =~ /[px]$/ ? 'aa' : 'dna');
        
	# make sure target database ok
	
	$self->throw("No db supplied: $db ($homology)") unless $db;
	$db = $self->_write_temp_seq(				
	    -seq => [ $self->orfs ],
	    -molecule => $dbmol
	    ) if ref($db) eq ref($self);

	unless ( $db eq 'nr' ) { 
	    $self->_valid_file(-file => $db);
	    my $formatdb = $self->_external_app(); # formatdb is default 
	    my $farg = '-p F' unless $args->{'-dbmol'} =~ /aa|prot|amino/;   
	    `$formatdb $farg -i $db`;
	    $self->throw("Formatdb failed: ".`ls $db*`)
		unless (-e $db.".nhr" || -e $db.".phr");
	}
	
	# define query sequences 

	my @orfs = (
	    $homology eq 'NCBI' ? 
	    grep { $_->length > $min_ncbi_length } 
	    grep { $_->evalue('aa') > $e && $_->evalue('ygob') > $e && $_->evalue('ty') > $e } 
	    $self->orfs : 
	    $self->orfs # everythign must have BLAST for synteny methods to work 
	    );

	return $self if $#orfs == -1;
	my $submit = ( $homology eq 'NCBI' ? $ncbi_batch_size : $local_batch_size );
	my $old_time = time;

	# start executing BLAST jobs 

	for (my $i = 0; $i <= $#orfs; $i += $submit ) {  # takes pressure off blastcl3 
	    my $j = ($i+$submit-1 > $#orfs ? $#orfs : $i+$submit-1);

	    # make bacth and time queries 
	    
	    my $time = time;
	    my $querytime = sprintf("%.2f", ($time-$old_time)/($j-$i+1) );
	    print STDERR $#orfs, $i, $j, ($time-$old_time).'s', $querytime.'s', 
	    int($querytime*($#orfs-$j)/60).'m' if $i%10==0;
	    $old_time = $time;

	    # run BLAST 
	    
	    next unless my $blast = 
		$self->blast(
		    -program => $prog,
		    -qseq => [@orfs[$i..$j]], 
		    -qmol => $qmol,
		    -db => $db, # can be an url or ncbi DB or file 
		    -dbmol => $dbmol,
		);
		
	    # add data to ORF objects 
	    # no need to call _process_BLAST becasue only taking 
	    # topp hit each time 
	    
	  ORF:foreach my $orf ($self->orfs) {
	      next unless exists $blast->{$orf->_internal_id};	     
	      my @bl = sort {$a->{EVALUE} <=> $b->{EVALUE}} @{ $blast->{$orf->_internal_id} };

	      if ($homology eq 'AA') {
		  
		  $orf->homology(-blast => \@bl); # populate the homology structure 
		  
	      } elsif ($homology eq 'NCBI') { # process the BLAST array to remove unwanted hits 		  
		  if ( $self->organism =~ /Sbay|Scer|Smik|Spar|Skud|Sarb/ ) {
		      $bl[0]->{'HIT'} =~ /gi\|(\d+)\|/;
		      my $content = get("$eutils$1");
		      #$content =~ /ORGANISM\s+([^\n]+)\n/;
		      $content =~ /GBSeq_organism\>([^\<]+)/;# Saccharomyces cerevisiae S288c</GBSeq_organism>
		      my $org = $1;
		      while ( $org =~ /Saccharomyces/i ) {
			  pop(@bl);
			  next ORF unless @bl;
			  $bl[0]->{'HIT'} =~ /gi\|(\d+)\|/;
			  my $content = get("$eutils$1");
			  #$content =~ /ORGANISM\s+([^\n]+)\n/;
			  $content =~ /GBSeq_organism\>([^\<]+)/;# Saccharomyces cerevisiae S288c</GBSeq_organism>
			  $org = $1;
		      }
		  }
	      } # LTR / TY -- no specific actions 	      
	      next unless @bl && $#bl >= 0; # git anything?

	      # impute() does 2 things: 
	      # 1. it generically sets KEY, _KEY, __KEY values 
	      # 2. depending on the homology class, it determines a "best" hit
	      # by also using synteny. This is stored on the (_)(_)GENE attributes
	      # and is considered the most likley identity for the query. 
	      # In general, TY/LTR cannot set GENE, NCBI can set it only 
	      # if there is no better hit in the local DB (AA) and, in the common case,
	      #  extensive use of synteny is used to choose the best GENE. 
	      
	      my $new = $orf->impute(
		  -blast => \@bl,
		  -homology => $homology,
		  -evalue => $e
		  ); 
	      
	      # NCBI specific stuff. 
	      
	      if ( $new && $homology eq 'NCBI' ) {
		  $orf->evaluate;
		  $orf->output(-prepend => ['NCBI']);
	      }
	  }	
	}
    }

    map { $_->evaluate } $self->orfs; 
    
    return $self;
}

sub hmmer_ygob {
    my $self = shift;        
    my $args = {@_};
    $args->{'-verbose'} = 1 unless exists $args->{'-verbose'};

    my $orfs = scalar( $self->orfs );
    my ($totalorfs, $inittime) = (0,time);

    foreach my $orf ( $self->orfs ) {
	$orf->update( -blast => undef, -exonerate => undef, -synteny => undef ); # run hmmer only 
	$totalorfs++;
	my $totaltime = time - $inittime;
	my $unittime = sprintf("%.2f", $totaltime/$totalorfs);

	print $orfs, $totalorfs, $totaltime.'s', $unittime.'s',
	$orf->name, ($orf->ygob || 'NOYGOBHIT') , 
	sprintf("%.2f", ($unittime/60)*($orfs-$totalorfs) ).'m' 
	    if $args->{'-verbose'};
    }

    return $self;
}

=head2 learn

    Eventual replacement for the features method. We learn models
    from genome then use contigs->features to find all instacnes. 

=cut 

sub learn {
    my $self = shift;
    $self->throw;
}

=head2 features(-intron => ['INTRON.hmm', 4.5], -telomere => ['TEL.hmm', 10])

    Wrapper around a contig method to place HMMs. 

    This method uses the same syntax as homology() but applied to  
    sequence features, such as introns and telomere repeats, that are 
    well described by HMMs but maybe not by BLAST. 

=cut

sub features {
    my $self = shift;
    map { $_->features(@_) } $self->stream;
    return $self;
}

=head2 analysis(-lda => .95, -hcnf => 1e-10, -kaks => .05)
    
    This is a high-level method for running generic 
    evidence gathergin methods. 

    e.g. the lda() method performs linear discriminant analysis
    using codon usage frequencies and requires only a posterior
    probability to interpret the results. it can therefore be 
    dispatched using the analysis() method with other 
    similar analyses such as kaks(). 

    Like the homology() method the argument names must correspond
    to entries in the %EVIDNCE hash in GlobalVars.pm to be used
    properly to infer sequence types. The values are the evidence 
    threshold to be believed. 

=cut

sub analysis {
    my $self = shift;
    my $args = {@_};	
    
    # requires reliable orfs and non-orfs
    
    foreach my $contig ($self->stream) {
	foreach my $orf ($contig->stream) {
	    $orf->evaluate;
	}
    }
    
    # call analysis methods in turn 
    
    foreach my $evidence (keys %{$args}) {	
	my $method = $evidence;		
	$method =~ s/^-//;
	$self->_set_eval(-evidence => uc($method), -set => $args->{$evidence})
	    if defined $args->{$evidence} && $args->{$evidence} =~ /\d/;# defaults
	$self->throw($method) unless $self->can($method);
	$self->$method(@_);
    }
    
    return $self;
}

=head2 kaks()

    Run either codeml or yn00 on the calling ORF and 
    populate KA, KS and KAKS evidence parameters. 

=cut 

sub kaks {
    my $self = shift;             
    map { $_->kaks(@_) } $self->orfs;
    return $self;
}

=head2 atg(-align => 10, -verbose => 1)

    This is a wrapper around two orthogroup methods. 
    We use alignedStartCodon to test whether we believe we have the correct 
    start codon. If it returns false, we use atg_og to leverage information 
    acorss genomes to find a "best consensus" start codon. 

=cut 

sub atg {
    my $self = shift;
    my $args = {@_};
  
    $args->{'-verbose'} = undef unless exists $args->{'-verbose'};
    
    # 

    my $count;
  OG:foreach my $o ( grep { $_->assign =~ /REAL|HSP|NOVEL/ } grep {$_->orthogroup} $self->orfs ) {
      #next unless $o->name eq 'Scer_2.100';
      
      # skip or try get M ?

      map {$_->index; next unless $_->exons(-query => 'first')->length >= $TRIPLET*10 } ($o,$o->orthogroup);
      my $mal= $o->alignedStartCodon();
      
      my ($m,$sd) = _calcMeanSD( map { $_->length } ($o, $o->orthogroup) ); # this is just fluff 
      print ">".(++$count), $o->name, $o->gene, $m, $sd, $mal if $args->{'-verbose'};

      $o->atg_og(-verbose => $args->{'-verbose'}, -reference => undef) unless $mal;
  }
    
    return $self;
}

=head2 orfs(-intergenic => undef)

    Get an ordered stream of ORF objects for whole genome. 
    Dispatches the hard work to the Contig method of the 
    same name -- see there for arguments. 

=cut 

sub orfs {
    my $self = shift;
    return map { $_->orfs(@_) } $self->stream;
}

#########################################
#  subroutines : accounting 
#########################################

=head2 _validate_data_structures

    -repair : call reoptimise() on genes that do not validate. 
    -forgive : tolerate (and fix) genes that think they are in 
    OGs but are not. 
    
=cut 

sub _validate_data_structures {
    my $self = shift;
    my $args = {@_};

    $args->{'-repair'} = undef unless exists $args->{'-repair'};
    $args->{'-forgive'} = undef unless exists $args->{'-forgive'};

    # 

    my %index;
    if ( $args->{'-repair'} ) {
	foreach my $x ( grep { $_->orthogroup} $self->orfs ) {
	    my ($sgd,$anc) = $x->identify;
	    map { $index{$_->_internal_id } = $sgd } ($x,$x->orthogroup);
	}
    }

    # 
    
    my @genomes = ($self, (map { $self->bound($_) } $self->bound));;

    my %hash;
    foreach my $g ( @genomes ) {
	foreach my $c ($g->stream) {
	    foreach my $o ($c->stream) {

		# complete this method.. 
		
		# $o->validate;
		
		# is there a reciprocal link to the parent object? 
		# MUST be symmetrical. 
		# check exons

		print $g->organism.$c->id and $o->throw unless $o->up;
		map { print $o->name and $o->throw unless $_->up eq $o } $o->stream;
		
		# object structure ok ?

		unless ( $o->rank < -1 || $o->translatable ) {
		    if ( $args->{'-repair'} ) {
			my $refer =  (
			    exists $index{$o->_internal_id} && defined $index{$o->_internal_id} ? 
			    $index{$o->_internal_id} : 'GENE');
			$o->reoptimise(-reference => $refer);
			$o->evaluate(-validate => 1, -structure => 1); 
		    } else { 
			$o->output;
			print $o->aa;
			$self->throw;
		    }
		}

		# 
		
		if ( my $ohno = $o->ohnolog ) {
		    if ( $ohno->ohnolog ne $o ) {
			$o->output;
			$ohno->output;
			if ( ! $ohno->ohnolog ) {
			    $ohno->ohnolog($o);
			    print;
			} else {
			    $ohno->ohnolog->output;
			    $self->throw;
			}
		    }
		}
		
		# we now look for all possile objects that could somehow 
		# be missed by the stream method. 

		map { $o->oliver and $_->throw("LEFT:".$o->name) unless $_->up; } grep {defined} $o->left;
		map { $o->oliver and $_->throw("RIGHT:".$o->name) unless $_->up; } grep {defined} $o->right;
		map { $o->oliver and $_->throw("OG:".$o->name) unless $_->up; } grep {defined} $o->orthogroup;
		map { $o->oliver and $_->throw("OHNO:".$o->name) unless $_->up; } grep {defined} $o->ohnolog;

		# this is the setup for testing OG system. 
		
		if ( $o->ogid ) { # we detect with the _OG system .... 
		    $hash{ $o->organism }{ $o->_internal_id } = $o;
		}
	    }
	} 
    }
    
    # 
    
    my $ogid; # .... and compare to the ->orthogroup system 
    foreach my $c ( $self->stream ) {
	foreach my $o ( grep {$_->orthogroup} $c->stream ) {
	    $o->oliver and $o->throw unless exists $hash{ $o->organism }{ $o->_internal_id };
	    delete $hash{ $o->organism }{ $o->_internal_id };
	    
	    foreach my $og ( $o->orthogroup ) {
		$o->oliver and $o->throw unless exists $hash{ $og->organism }{ $og->_internal_id };
		delete $hash{ $og->organism }{ $og->_internal_id };
	    }
	}
    }

    # 
    
    if ( %hash ) {
	print "OG Problems...";
	my $dangling=0;
	foreach my $sp ( keys %hash ) {
	    print ">$sp ",scalar( keys %{$hash{$sp}} );
	    map { $_->oliver(-append => [$_->ohnolog]) } sort {$a->name cmp $b->name } values %{$hash{$sp}};
	    if ( $args->{'-forgive'} ) {
		map { $_->ogid(undef) } values %{$hash{$sp}};
	    } else {
		$dangling=1 if values %{$hash{$sp}};
	    }
	}
	$self->throw if $dangling;
    } else { print "looks OK..."; }
    
    return $self;
}

=head2 gff2genome()

    Read in a gff file (only tested for SGD format) and 
    port annotations to a cloned chromosomes. These are 
    asymmetrically attached to the existing genome object
    -- the new orfs/contigs can access their organism etc 
    but the cloned contigs are not returned when stream
    is called on the genome object. 

    NB: We use the CDS field ONLY and assume the following: 

    1. orf_classification\=([^;\s]+) 
    2. Parent\=([^;\s]+) -> we use this to compose genes from exons

=cut 

sub gff2genome {
    my $self = shift;
    my $args = {@_}; 

    $args->{'-source'} = 'SGD' unless exists $args->{'-source'};
    $args->{'-restrict'} = undef unless exists $args->{'-restrict'};
    $args->{'-dubious'} = undef unless exists $args->{'-dubious'};
    $self->throw($args->{'-gff'}) unless $args->{'-gff'} && -e $args->{'-gff'};

    my $fh = *STDERR;
    my $attr = 'CLASS';

    #########################################
    # Read in GFF 
    #########################################

    my %chr;
    open(GFF, $args->{'-gff'}) || $self->throw( $args->{'-gff'} );
    while ( my $line = <GFF> ) {
	next if $line =~ /^\#/;
	my @r=split/\t/, $line;	
	next unless $r[2] eq 'CDS';

	########################	
	# deal with chr format 
	########################

	$r[0] =~ s/chr//;
	my $chr = ( $r[0] =~ /[XIV]/ ? arabic($r[0]) : $r[0]);
	next if $chr =~ /\D/; #next if $r[0] =~ /micron|mito/i;

	########################
	# apply filters
	########################

	next if defined $args->{'-restrict'} && $args->{'-restrict'} != $chr;
	next unless $args->{'-dubious'} || $r[8] !~ /orf_classification\=Dubious/i;

	########################
	# compose genes from exons/etc 
	########################

	$r[8] =~ /Parent\=([^;\s]+)/ || die; 
	my $par = $1;
	$par =~ s/_mRNA//i;

	########################
	# load up a data structure 
	########################

	push @{$chr{$chr}{$par}{'EXONS'}}, 
	{
	    START => $r[3],
	    STOP => $r[4],
	    STRAND => ($r[6].($r[6] =~ /1/ ? undef : '1'))*1,
	    INTRON => [$INFINITY, $INFINITY]
	};

	$chr{$chr}{$par}{$attr} = ( $r[8] =~ /orf_classification\=([^;\s]+)/ ? $1 : 'OtherTE');
    }
    close GFF;

    #########################################
    # 
    #########################################

    my $newGenome = ref($self)->new
	(
	 ORGANISM => $self->organism,
	 STRAIN => $self->strain,
	 SAMPLE => $self->sample,
	 # 
	 COVERAGE => 'Unknown',
	 SOURCE => $args->{'-source'},
	 ANNOTATION => $args->{'-source'},
	 # 
	 WGD => $self->wgd
	);

    #########################################
    # make pseudo-objects for chrs and ORFS 
    #########################################
    
    foreach my $chr ( keys %chr ) {
	$self->throw unless my $ctg = $self->find(-contig => $chr);
	$self->throw unless my $clone = $ctg->clone;
	$newGenome->add( -object => $clone );

	# 

	foreach my $gene (keys %{$chr{$chr}}) {
	    my $orf = Annotation::Orf
		->new(
		START => 1,
		STOP => 1,
		STRAND => 1
		);	    
	    my $fex = $orf->down;
	    foreach my $ex ( @{$chr{$chr}{$gene}{'EXONS'}} ) {
		my $obj = Annotation::Exon->new( %{$ex} );
		$orf->add(-object => $obj);
	    }
	    $orf->remove(-object => $fex);
	    $fex->DESTROY;
	    $orf->{'STRAND'} = $orf->down->strand; # fix orf str
	    $orf->index;

	    # set some vars. most notably we provide the nonstandad CLASS attribute  
	    
	    $orf->data($attr => ($chr{$chr}{$gene}{$attr} || 'unknown') );
	    $orf->data('GENE' => $gene );	    
	    $orf->_creator( $args->{'-source'} );
	    $clone->add(-object => $orf); 
	}
	$clone->index;       
    }   

    return $newGenome;
}
 
=head2 compare(-gff => file, -object => genome, -restrict => contigid,
     -link => 1, -view => undef, -dubious => 0)

    Read (SGD style) GFF or genome obbject and compare to current annotation. 
    By Default we output counts of error types but with -guide we 
    initiate updates using exonerate. 

    -link establishes a reciprocal link between paired annotations 
    from each genome. eg the SGD TUB2 annotation and self-produced 
    TUB2 annotation. They are accessible by $sgd = $self->_debug()
    Note that this ~ doubles the size of the storable file.

=cut 

sub compare {
    my $self = shift;
    my $args = {@_};    

    $self->throw if $args->{'-guide'}; # old option
    #$args->{'-restrict'} = undef unless exists $args->{'-restrict'};
    #$args->{'-dubious'} = undef unless exists $args->{'-dubious'};
    $args->{'-view'} = undef unless exists $args->{'-view'};
    $args->{'-link'} = 1 unless exists $args->{'-link'};

    $self->throw if $args->{'-gff'} && $args->{'-object'};

    my $fh = *STDERR;
    my $attr = '_COMPLEX_OLAP_COUNT';

    ######################################### 
    # read in GFF or genome object 
    ######################################### 

    my $gff = ($args->{'-object'} ? $args->{'-object'} : $self->gff2genome( @_ ) );
    $self->throw unless $gff && $self->isa(ref($gff));
    $self->throw unless scalar($self->stream) == scalar($gff->stream);

    # might be alt version of same genome so must make 
    # _internal_id's unique
    
    if ( $gff->_internal_id eq $self->_internal_id ) {
	map { $_->{'_ID'} = &{$Annotation::unique_id} } map { $_->stream} $gff->stream;
	$args->{'-link'} = undef; # we do not link ortho models 
    }

    ######################################### 
    # use commoncodon to group gene models 
    ######################################### 

    my %labels;    
    foreach my $ctg ( $self->stream ) {
	$self->throw unless 
	    my $clone = $gff->find( -contig => $ctg->id );

	$ctg->index;

	# cluster all the orfs 

	my $clusters = $ctg->cluster( 
	    -object => [ $args->{'-object'} ? $clone->orfs : $clone->stream ], # add SGD
	    -frame => 1,                   # use commoncodon 
	    -param => $NONZERO,            # any overlap is OK 
	    -accelerate => 10              # how many neighbour genes to consider 
	    );

	
	# examine each cluster 

	foreach my $cl ( values %{$clusters} ) {

	    my @g1 = grep { $_->up->up eq $self } @{ $cl };
	    my @g2 = grep { $_->up->up ne $self } @{ $cl };
	    my ($me,$sgd)=($g1[0],$g2[0]);

	    my $tag='junk';
	    my @tag=();	    
	    if ( grep { $_->assign eq 'REPEAT' || $_->gene eq 'YNL054W-B' 
			    || $_->data('SGD') eq 'YNL054W-B'} @{$cl} ) {
		# this should be cleaned up ...
		$tag='repeat';
		push @tag, @g1; #@g2;

	    } elsif ( ($#g1 >0 && $#g2>=0) || ($#g2 >0 && $#g1>=0) ) { # 2:1 / 1:2 

		$tag = ( $#g1 == $#g2 ? 'complx' : 'strct'); # assume complex locus but correct 
		$me->data( $attr => scalar(@g1).'/'.scalar(@g2) );
		push @tag, @g1 ; #, @g2;

	    } elsif ( $#g1 >0 || $#g2 >0 ) {
		$tag = 'rdnt';
		push @tag, @g1,@g2;

	    } elsif ( $#g1 == -1) {
		$tag = 'missed';
		push @tag, @g2;

	    } elsif ( $#g2 == -1 ) {
		if ($me->rank < 3 && $me->rank >= -1 && $me->assign !~ /REP/) {
		    $tag='extra';
		    push @tag, @g1;
		}
	    } elsif ( $#g1==0 && $#g2==0 ) { # 1:1 
		if ( $me && abs($sgd->start - $me->start) <=3 && 
		     abs($sgd->stop - $me->stop) <=3 ) {		
		    $tag = 'good';		
		} elsif ($sgd->overlap( -object => $me, -contig => 0 ) >= 0.7) {
		    $tag = 'ok';
		} else {
		    $tag='poor';
		}
		push @tag, @g1;	
		$me->_debug($sgd) if $args->{'-link'};
	    } else { $self->throw("Unprogrammed case:  $#g1, $#g2 ");}
	    push @{$labels{$tag}}, @tag;
	} # clusters 
    }# chr


    ######################################### 
    # output
    ######################################### 
    
    foreach my $lab ( sort { $#{$labels{$b}} <=> $#{$labels{$a}} } keys %labels ) {
	my %cx;
	if ( $lab eq 'missed') { 
	    map { $cx{ $_->data('CLASS') || 'NoClass' }++ } @{$labels{$lab}};
	    print {$fh} $lab, (map { "$_:$cx{$_}" } sort keys %cx); 
	} elsif ($lab eq 'cmplx' || $lab eq 'strct' ) {
	    map { $cx{ $_->data($attr) }++ } grep { $_->data($attr) } @{$labels{$lab}};
	    print {$fh} $lab, (map { "$_:$cx{$_}" } sort keys %cx); 
	} else {
	    print {$fh} $lab, 
	    (( $#{$labels{$lab}} > 15 ? @{$labels{$lab}}+0 : 
	       join(' | ', map {$_->data('GENE') || $_->data('SGD')} @{$labels{$lab}} ) ) || 0);
	}
    }
    
    ######################################### 
    # more output?
    ######################################### 
    
    if ( $args->{'-view'} ) {
	my $xyz;
	foreach my $orf ( sort { $a->name <=> $b->name } @{$labels{$args->{'-view'}}} ) { 
	    $orf->evaluate(-force => 1);
	    $orf->glyph if $orf->data($attr);
	    #$orf->show(-species => []) if $orf->data($attr);
	    
	    $orf->output(
		-prepend => [
		     #$orf->score('local'),$orf->score('global'),$orf->data('HSP'),
		     #sprintf("%.1f",$orf->synteny(-loss => 1)),
		     #$orf->length,
		     ++$xyz,
		     $orf->data($attr)
		],
		-append => [
		     substr($orf->data('CLASS'), 0, 6), 
		     #($orf->gene ? $orf->lookup(-query => $orf->gene, -species=>'YGOB') : undef),
		],
		-recurse => 1
		);
	}
    }

    ######################################### 
    # get out of here
    ######################################### 
    
    map { $_->data( $attr => 'delete' ) } map { $_->stream } map { $_->stream } ($gff,$self);
    return $self; # unless $args->{'-guide'};
}

sub _guide_mode {
    
    ################################################################
    # accounting complete. missing -> search_region, inexact -> reoprimise
    ################################################################

    my $fixed;
    foreach my $atg ( map { @{$count{$_}} }  qw(atg stop both) ) {
	next unless my ($ctg,$sgd) = ($atg->up,$atg->data($sgdstoretag));
	next unless $sgd->translatable;  # ignore frame-shifted genes 
	my $e1 = ($sgd->length - $atg->length);
	#$atg->oliver(-append => [$e1]);
	#$sgd->oliver(-append => ['SGD']);
	$atg->reoptimise(-reference => $sgd, -extend =>  $sgd->length);
	my $e2 = ($sgd->length - $atg->length);
	#$atg->oliver(-append => [$e2]);
	$fixed++ if $e2==0;
   }

    # At this point every gene on chr3 is perfect except 6 genes we do not find. 
    # YCR108C | YCR050C | YCR024C-A | YCL001W-B | YCL048W-A | YCL058C
    # most have a good - if not wholly exonerating - explanation. 
    # eg one entirely covered by another larger protein.. 

    my $found;
    foreach my $sgd ( @{$count{'missed'}} ) {
	my $cl = $sgd->clone;
	my $chr = $self->find(-contig => $cl->up->id);
	$cl->transfer(-from => $cl->up, -to => $chr, -warn => 0);
	my $new = $chr->exonerate(
	    -protein => $cl, 
	    -object => $cl, # need transfer above to get coords .. 
	    -extend => 500, 
	    -intron => ($sgd->exons > 1 ? -10 : -30) 
	    );
	if ( $new && $new->length == $sgd->length )  {
	    $new->update;
	    $chr->add(-object => $new);
	    $new->_debug($sgd);
	    $found++;
	}
	$cl->DESTROY;
    }
    map { $_->index } $self->stream;

    print $fixed, $found;
    # OK. thats all of them....  
    ################################################################

    return $self;
}


sub _allObjectStream {
    my $self = shift;
    return( 
	$self, 
	$self->stream,
	(map { $_->stream } $self->stream),
	(map {$_->stream} map { $_->stream } $self->stream)
	);
}


################################################################
# output/logging routines 
################################################################

sub store { my $self = shift; return $self->backup(@_); }

=head2 backup(-bind => [any extra objects], -path => , -tag => '')

    Backup object to file using Storable. We try to save object 
    at all costs and will try work around broken args, filehandles etc. 

    Filename is chosen by _outfilenamehandle(), typically of form 
    /full/path/Spec_time.tag.suffix
    
    We automatically close and delete any open filehandles. 
    
=cut

sub backup {
    my $self = shift;
    my $args = {@_};

    $args->{'-bind'} = [] unless exists $args->{'-bind'};

    ########################################    
    # usually we throw for bad arguments but here we want to save
    # the storable at all costs. And alert user to problems. 
    ########################################   

    my $file = $self->_outfilenamehandle( 
	@_, -suffix => ($#{$args->{'-bind'}} >= 0) ? 'bgo' : 'sto');

    ########################################   
    # coerce -bind to an array reference, $bundle  
    ########################################   

    unless ( ref($args->{'-bind'}) eq 'ARRAY' ) {
	my $cargo = $args->{'-bind'};
	push @{ $args->{'-bind'} }, $cargo;
	$self->warn("Unexpected argument to -bind: $args->{'-bind'}");
	$file .= ".badcargo";
    }

    my $bundle = 
	( $#{$args->{'-bind'}} > -1 ? [$self, @{$args->{'-bind'}}] : [$self] );

    ########################################
    # check bundle for live code refs (filehandles) 
    # and close. they will prevent storing. 
    ########################################
    
    foreach my $genome ( @{$bundle} ) {
	if ( $self->isa(ref($genome)) ) {
	    next unless my $fh = $genome->{'_LOG'};
	    close($fh);
	    delete $genome->{'_LOG'};
	}
    }
    $bundle = $bundle->[0] if $#{$bundle} == 0; # ugh. 

    ########################################
    # write network format binary file 
    ########################################

    $self->throw("nstore error: $file / $!") unless
	nstore($bundle, $file);
    print STDERR "Genome stored to $file";
    
    return $self;
}

=head2 _outfilenamehandle(-path => `pwd`, -id => organism, -time => $TIME,
    -tag => , -suffix => undef, -fh => undef, -append => undef) 

    Return a filename (default) or filehandle based on supplied params:

    /full/path/Spec_time.tag.suffix

=cut 

sub _outfilenamehandle {
    my $self = shift;
    my $args = {@_};
    
    $args->{'-id'} = $self->organism unless exists  $args->{'-id'};
    $args->{'-time'} = $TIME unless exists $args->{'-time'};
    $args->{'-tag'} = ++$STORE unless exists $args->{'-tag'};
    $args->{'-suffix'} = undef unless exists $args->{'-suffix'};
    # 
    $args->{'-path'} = `pwd` unless exists $args->{'-path'};
    chomp( $args->{'-path'} );
    # 
    $args->{'-fh'} = undef unless exists $args->{'-fh'};
    $args->{'-append'} = undef unless exists $args->{'-append'};
    
    
    ########################################
    # construct a file name 
    ########################################
    
    my $file = $args->{'-id'}.'_'.$args->{'-time'}.'.'.$args->{'-tag'};
    $file .= '.'.$args->{'-suffix'} if $args->{'-suffix'};
    $file = $args->{'-path'}.'/'.$file if
	( -d $args->{'-path'} && -w $args->{'-path'} && $file !~ /\// );
    
    ########################################
    # check 
    ########################################
    
    if ( -e $file ) {
	$self->warn("$file exists!");
	$file.='.avoid-overwrite'.int(rand(10000000));
    }
    
    system("touch $file");
    $self->warn("$file: bad touch!") unless -e $file && -w $file;
    
    if ( $args->{'-fh'} ) {
       open(my $fh, ($args->{'-append'} ? '>>' : '>').$file) || 
	   $self->throw("Open: $file");
       return $fh;
    }
    
    return $file;
}

=head2 history(newhash|integer) 

    If hash adds the hash to the annotation history. 
    If integer returns the requested history item-- 
    0 corresponds to present, 1 to previous logged step etc. 
    With no arguments prints the entire annotation history. 

=cut 

sub history {
    my $self = shift;
    my $args = {@_};
    my $new = shift;
    
    if ( ref($new) eq 'HASH' ) {
	$self->throw() unless exists $new->{FILE} && exists $new->{DATE}; 
	push @{$self->{HISTORY}}, $new;
    } elsif ( $new =~ /^\d+$/ ) {
	$self->warn("No step $new in history!") and return undef if 
	    ($new > $#{$self->{HISTORY}} || $new < 0);
	# gah. should have used unshift instead of push above.. 
	return $self->{HISTORY}->[ ( $#{$self->{HISTORY}} - $new) ];
    } else {
	map { print $_->{STAMP},$_->{DATE},$_->{SCRIPT},$_->{STEP},$_->{FILE},$_->{NOTE} } 
	@{ $self->{HISTORY} };
    }

    return $self;
}

=head2 log(-object => $obj, -text => , -fh => , -echo => '1|fh')
=cut 

sub log {
    my $self = shift;
    my $args = {@_};

    $args->{'-echo'} = undef unless $args->{'-echo'};
    $args->{'-fh'} = $self->{'_LOG'} unless exists $args->{'-fh'};
    $self->throw unless $args->{'-object'} || $args->{'-text'};

    if ( $args->{'-object'} ) {
	$self->throw unless $args->{'-object'} =~ /::/ && $args->{'-object'}->can('output');
    }

    my @fh; 
    if ( $args->{'-fh'} ) {
	push @fh, $args->{'-fh'};
    } else {
	my $stamp = $self->history(0)->{STAMP}; # current run id 
	#print $stamp;
	open(my $fh, '>>'.$self->organism.'.'.$stamp.'.log') || $self->throw;
	$self->{'_LOG'} = $fh;
	$args->{'-fh'} = $fh;
	push @fh, $fh;
    }

    if ( $args->{'-echo'} ) {
	push @fh, ( $args->{'-echo'} == 1 ? \*STDERR{IO} : $args->{'-echo'} );
    }

    foreach my $fh ( @fh ) {
	print {$fh} ">$args->{'-text'}" if $args->{'-text'};
	$args->{'-fh'} = $fh;
	$args->{'-object'}->output(%{$args}) if $args->{'-object'};
    }

    return $self;
}

=head2 summarize(-nosynteny => 0|1)

    Print out genome statistics as described below. 
    -nosyntney leads to a significant speedup but
    omits 'Synteny' and 'Ancestral Orthogroup Size'.

SPECIES Zbis 
CONTIGS 499
N50     54532

GAP     NNNN:606
TRNA    RNA:269
PSEUDO  STOP:29
REPEAT  TY:1
REAL    SYNT:34 NCBI:5  AA:151  YGOBHMM:5010
HSP     HSP:57
NOVEL   LDA:57  LENGTH:23

Synteny (# syntenic neighbours):
1       -1:13   0:15    1:121   2:387
2       -1:17   0:10    1:115   2:527
3       -1:20   0:17    1:112   2:432
4       -1:13   0:6     1:97    2:273
5       -1:13   0:18    1:159   2:529
6       -1:9    0:7     1:94    2:271
7       -1:11   0:8     1:100   2:429
8       -1:24   0:22    1:177   2:656

Ancestral Orthogroup Size:
0:120
1:4114
2:360
3:73
4:19
5:3
6:5
7:2
9:3
13:1
14:2

=cut

sub summarize {
    my $self = shift;
    my $args = {@_};

    $args->{'-nosynteny'} = undef unless exists $args->{'-nosynteny'};
    $args->{'-fast'} = undef unless exists $args->{'-fast'};
    $args->{'-fh'} = undef unless exists $args->{'-fh'};

    my $fh;
    if (  $args->{'-fh'} ) {
	$fh =  $args->{'-fh'};
    } else {
	$fh = *STDOUT;
    }

    print $fh '';
    print $fh $self->organism;
    print $fh 'Contigs', scalar( $self->stream );
    print $fh 'N50', $self->n50;
    print $fh 'Feats', scalar($self->orfs);
    print $fh 'Ortho', scalar( grep { $_->ogid }  $self->orfs);
    print $fh 'Ohno', scalar( grep { $_->ohnolog } $self->orfs);
    print $fh 'G50', $self->g50;
    print $fh '';

    my %hash;
    my %syn;
    my %ygob;
    foreach my $contig ( $self->stream ) {
	foreach my $orf ( $contig->stream ) {
	    $orf->evaluate unless $args->{'-fast'};
	    $hash{ $orf->assign }{ $orf->evidence }++;
	    next if $args->{'-nosynteny'};

	    $syn{ $orf->data('SYNT') }++;

	    if ( $orf->ygob =~ /_(\d+)\.(\d+)/ ) {
		my ($anc,$ord) = ($1, $2);

		my $syn2 = $orf->querysynteny(
		    -distance => 2,
		    -difference => 5,
		    -spanning => 1,
		    -restrict => ['YGOB']
		    );
		my $syn1 = $orf->querysynteny(
		    -distance => 2,
		    -difference => 5,
		    -spanning => 0,
		    -restrict => ['YGOB']
		    );

		my $syn=0;
		if ( $syn2 ) {
		    $syn=2;
		} elsif ( $syn1 ) {
		    $syn=1;
		} 

		$ygob{ $anc }->[$ord]->{ $syn }++;
	    }
	}
    }
    
    # inference by evidence 

    my %order;
    foreach my $q ( sort {$EVIDENCE{$a}->{'RANK'} <=>
                             $EVIDENCE{$b}->{'RANK'}} keys %EVIDENCE ) {
	$order{ $EVIDENCE{$q}->{'INFER'} } = $EVIDENCE{$q}->{'RANK'};
    }

    foreach my $k (sort {$order{$a} <=> $order{$b}} keys %hash) {                             
	print $fh $k, map { "$_:$hash{$k}{$_}" } keys %{$hash{$k}};
    }
    print $fh '';

    return $self if $args->{'-nosynteny'};

    # synteny by YGOB chromosome 

    print $fh 'Synteny (# syntenic neighbours):';

    my %count;
    foreach my $chr (sort {$a <=> $b} keys %ygob ) {
	my %maxsyn;
	foreach my $i ( 1..$#{$ygob{$chr}} ) {
	    # get distribution of hits to models 
	    my $lcount=0;
	    map { $lcount+=$ygob{ $chr }->[ $i ]->{$_} } keys %{ $ygob{$chr}->[$i] } if $ygob{ $chr }->[ $i ];
	    $count{$lcount}++;
	    # things we missed completely 
	    $ygob{ $chr }->[ $i ]->{ -1 }++ unless $ygob{ $chr }->[ $i ];
	    # we just record the maximum 
	    my ( $max ) = sort { $b <=> $a } keys %{ $ygob{ $chr }->[ $i ] };
	    $maxsyn{ $max }++; 
	}
	print $fh $chr, map { "$_:$maxsyn{$_}" } sort {$a <=> $b} keys %maxsyn ;
    }

    print $fh '';
    print $fh 'Ancestral Orthogroup Size:';
    map { print $fh "$_:$count{$_}" } sort {$a <=> $b} keys %count;
    print $fh '';

    return $self;
}

#########################################
# subroutines : output 
#########################################

=head2 gff(-file => file.name, -exons => 1)

    Write gff file. Evidence is written an tag=value pairs.
    Use -exon switch to enable exon annotation output.
    Contigs and higher levels are not currently written. 

=cut 

sub gff {
    my $self = shift;
    my $args = {@_};

    $args->{'-source'} = 'ScannellZill2011G3' unless exists $args->{'-source'};
    
    my $fh;
    if (exists $args->{'-file'}) {
	$args->{'-file'} = '>'.$args->{'-file'} unless $args->{'-file'} =~ /^\>/;
        open($fh, $args->{'-file'}) || die($args->{'-file'});
    } else {
        $fh = STDOUT;
    }
    $args->{'-fh'} = $fh;

    print {$fh} "##gff-version 3";
    print {$fh} "#Produced by Annotation.pm",$AUTHOR,$YEAR;
    print {$fh} "#Details ".join(' ',map { lc($_).($self->$_ =~ /\s/ ? '="'.$self->$_.'"' : '='.$self->$_ ) } qw(ORGANISM STRAIN SAMPLE COVERAGE SOURCE WGD DATE ID));
    print {$fh} "#SourceTag ".$args->{'-source'};     
    print {$fh} "#Label ".$args->{'-label'};
    print {$fh} "#FreeText ".$args->{'-notes'};     
    
    foreach my $contig ($self->stream) {
	print {$fh} join(" ","##sequence-region",$contig->id,'1',$contig->length);       
        foreach my $orf ($contig->stream) {
            $orf->gff(%{$args});
        }
    }

    return $self;
}

#########################################
# subroutines : aliases
#########################################

sub species { my $self=shift; return $self->organism(@_); }
sub orthogroups { my $self = shift; return grep { $_->ogid } $self->orfs(-noncoding => 1); }

#########################################
# subroutines : private/experimental methods 
#########################################

sub _valid_file {
    my $self = shift;
    my $args = {@_};

	$args->{'-format'} = 'fasta' unless exists $args->{'-format'};	
    $self->throw("Must provide a sequence source (file or directory)")
        unless exists $args->{'-file'} && $args->{'-file'};	
		
    return 0 unless 
        -e $args->{'-file'} &&  # exists 
        -r $args->{'-file'} &&  # readable 
        -s $args->{'-file'} &&  # non-zero size
        -T $args->{'-file'};    # ASCII Text  

	if ($args->{'-format'} eq 'fasta') {
    	my $header = `head -2 $args->{'-file'}`;
    	return 1 if $header =~ /^>/;
    }
    
    return 0;
}

sub _findMatrixMin {
    my $m = shift;
    my $pen = (shift || 20)-1;

    die unless $pen;
    die unless ref($m) eq 'ARRAY'; 
    die unless ref($m->[0]) eq 'ARRAY'; 
    
    my $imax = $#{$m};
    my $jmax = $#{ $m->[0] };
    my ($ix, $jx, $min) = (undef, undef, 10000);

    return undef if $imax == -1 || $jmax == -1;

    for my $i (0..$imax) {
	for my $j (0..$jmax) { 
	    if ($m->[$i]->[$j] < $min) {
		$min = $m->[$i]->[$j];
		$ix = $i;
		$jx = $j; 
	    }
	}
    }
    
    # remove rows and columns 
    splice( @{$m}, $ix, 1);
    map { splice( @{$_}, $ij, 1) } @{ $m }; 
    
    #print {NOWHERE} 'M', $imax, $jmax, $ix, $jx, $min, $#{$m}, $#{ $m->[0] };

    $min = $pen if $min > $pen;
    return( $min+1 );
}

sub _calcMeanSD {
    &Annotation::throw("@_") if $#_ < 0;
    my ($mean, $sd);
    map { $mean += $_ } @_;
    $mean /= ($#_+1);
    map { $sd += ($mean-$_)**2 } @_;
    my @sort = sort { $a <=> $b } @_;
    return ( 
	sprintf("%.3f", $mean), 
	sprintf("%.3f", sqrt($sd/($#_+1))), 
	$#_+1,
	$sort[-1],
	$sort[0]
	);
}

#########################################
# UGLY JUNK 
#########################################


1;

