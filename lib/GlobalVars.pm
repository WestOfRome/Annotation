package GlobalVars;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(
		%FEATURES	
		%HOMOLOGY	
%SPECIES
		%EVIDENCE
@EVIDENCE
@REGEX
	    %CODONS
	    %CODONS_CANDIDA
	    %DNA
	    $TRIPLET
$START_CODON
$STOP_CODON
	    $INFINITY
$NONZERO
$AUTHOR
$YEAR
	    $TIME
		$VERSION
$eutils
		);

#########################################
#########################################

# used by all objects 

$AUTHOR = 'Devin Scannell';
$YEAR = 2012;
$VERSION = 2.3;
$TIME = time;

%CODONS = (
		NNN => 'X',
	   ###################### T ######################
	   TTT => "F", TTC => "F", TTA => "L", TTG => "L",      # LF   # T
	   TCT => "S", TCC => "S", TCA => "S", TCG => "S", TCN => "S", # C
	   TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",      # Y*   # A
	   TGT => "C", TGC => "C", TGA => "*", TGG => "W",      # CW   # G
	   ###################### C ######################
	   CTT => "L", CTC => "L", CTA => "L", CTG => "L", CTN => "L", # T
	   CCT => "P", CCC => "P", CCA => "P", CCG => "P", CCN => "P", # C
	   CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",      # HQ   # A
	   CGT => "R", CGC => "R", CGA => "R", CGG => "R", CGN => "R", # G
	   ###################### A ######################
	   ATT => "I", ATC => "I", ATA => "I", ATG => "M",      # IM   # T
	   ACT => "T", ACC => "T", ACA => "T", ACG => "T", ACN => "T", # C
	   AAT => "N", AAC => "N", AAA => "K", AAG => "K",      # NK   # A
	   AGT => "S", AGC => "S", AGA => "R", AGG => "R",      # SR   # G
	   ###################### G ######################
	   GTT => "V", GTC => "V", GTA => "V", GTG => "V", GTN => "V", # T
	   GCT => "A", GCC => "A", GCA => "A", GCG => "A", GCN => "A", # C
	   GAT => "D", GAC => "D", GAA => "E", GAG => "E",      # DE   # A
	   GGT => "G", GGC => "G", GGA => "G", GGG => "G", GGN => "G"  # G
	   ### T ######### C ######### A ######### G #####
	   );

%CODONS_CANDIDA = %CODONS;
$CODONS_CANDIDA{'CTG'} = 'S';
delete $CODONS_CANDIDA{'CTN'};

%DNA = (
	A => 1,
	T => 1,
	G => 1,
	C => 1,
	N => 1
	);

$TRIPLET = 3;
$START_CODON='ATG';
$STOP_CODON=qr/^TAG|TGA|TAA$/;

$INFINITY = 1e100;
$NONZERO = 1e-20;

$eutils = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=gb&id=";

#########################################
#########################################

# used by Genome objects 

%FEATURES = (
        'INTRON' => {
        	#           +5 *
            # DONOR:  AGgGAa|GTAT.GTTca
            # BRANCH: attaTTTACTAACGAc
            # ACCEPT: ACATtGAAG|AaAtCaC
            #                   * -6
            'LENGTH' => 36,
            'START' => 5,
            'STOP' => -6
            },

        'TELOMERE' => {
            'LENGTH' => 523,
            'START' => 0,
            'STOP' => -0
            },

        'CENTROMERE' => {
        	#           +5 *
            # DONOR:  AGgGAa|GTAT.GTTca
            # BRANCH: attaTTTACTAACGAc
            # ACCEPT: ACATtGAAG|AaAtCaC
            #                   * -6
            'LENGTH' => 0,
            'START' => 0,
            'STOP' => -0
            },

        'ARS' => {
        	#           +5 *
            # DONOR:  AGgGAa|GTAT.GTTca
            # BRANCH: attaTTTACTAACGAc
            # ACCEPT: ACATtGAAG|AaAtCaC
            #                   * -6
            'LENGTH' => 0,
            'START' => 0,
            'STOP' => -0
            }            
        );

#########################################
# ORF objects .. 
#########################################

# if ( $string ~~ @REGEX ) { ... } 
# is currently broken in perl. does not perform regex on each element. 
# must be done with map:  grep {$var ~~ $_} @REGEX
@REGEX = ( 
    qr/([A-Z])([LR])(\d{3})([CW]\-[A-Z])?/, 
    qr/([A-Z]{4})(\d+)([A-Z])(\d+)([grst])?/, 
    qr/([A-Z]{1}[a-z]{2,3})_(\d+)\.(\d+)([a-z])?/ 
    );

%HOMOLOGY = (
    'SGD' => qr/Y[A-P][LR]\d{3}[CW](\-[A-Z])?/, 
#    'SCER' => qr/Y[A-P][LR]\d{3}[CW](\-[A-Z])?/,
    'EGOS' => qr/A[A-G][LR]\d{3}[CW](\-[A-Z])?/,
    # Standard 
    'KWAL' => qr/Kwal_\d+\.\d+[a-z]?/,
    'KPOL' => qr/Kpol_\d+\.\d+[a-z]?/,
    'SCAS' => qr/Scas_\d+\.\d+[a-z]?/,
    'YGOB' => qr/Anc_\d+\.\d+[a-z]?/,
    # Genolevures 
    'KLAC' => qr/KLLA\d+[A-P]\d+[rstg]?/,
    'CGLA' => qr/CAGL\d+[A-P]\d+[rstg]?/,
    'SKLU' => qr/SAKL\d+[A-P]\d+[rstg]?/,
    'ZROU' => qr/ZYRO\d+[A-P]\d+[rstg]?/,
    'KTHE' => qr/KLTH\d+[A-P]\d+[rstg]?/,
    'YLIP' => qr/YALI\d+[A-P]\d+[rstg]?/,
    # messed up 
    'SPOM' => qr/SP[A|B|C|MIT]/,
    'SCER_MITO' => qr/Q\d{4}/,
    );

%SPECIES = (
    'Y' => 'SGD',
    'A' => 'EGOS',
    'ANC' => 'YGOB',
    'CAGL' => 'CGLA',
    'KLLA' => 'KLAC',
    'SAKL' => 'SKLU',
    'ZYRO' => 'ZROU',
    'KLTH' => 'KTHE',
    'YALI' => 'YLIP'
    );

# would be nicer to have this in Annotation::Orf but
# since the settign procedure is run by Annotation::Genome
# it makes more sense to have it here

%EVIDENCE = (    

    MANUAL => {
	INFER => 'FEATURE',
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => 0,
	RANK => -$INFINITY,
	ORDER => 14,
	SOFA => 'sequence_feature',
	SOFA_PARTOF => undef,
	GENBANK_ASN1 => 'misc_feature', 
	GENBANK_ASN1_PARTOF => undef,
	STATUS => undef	
    }, 
    
    NNNN => {
	INFER => 'GAP',
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => 0,
	RANK => -$INFINITY,
	ORDER => 13,
	SOFA => 'assembly_component',
	SOFA_PARTOF => undef,
	GENBANK_ASN1 => 'assembly_gap', # 'gap'
	GENBANK_ASN1_PARTOF => undef, 
	STATUS => undef	
    }, 
    
    # noncoding. centromers etc. 
    HMM => {
	INFER => 'TELOMERE', # should be CHRM_FEAT
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => 0,
	RANK => -$INFINITY,
	ORDER => 12,
	SOFA => 'telomere',
	SOFA_PARTOF => undef,
	GENBANK_ASN1 => undef, # 'centromere' or 'telomere' but both  
	GENBANK_ASN1_PARTOF => undef, # require experimental validation 	
	STATUS => undef	
    }, 
    
    # rna genes 
    RNA => {
	INFER => 'TRNA',     # shoudl be NCRNA 
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => 0,
	RANK => -2,
	ORDER => 11,
	SOFA => 'tRNA',
	SOFA_PARTOF => 'noncoding_exon',
	GENBANK_ASN1 => 'tRNA',
	GENBANK_ASN1_PARTOF => 'exon',
	STATUS => undef	
    }, 	
    
    # STOP/codon ratio > x 
    STOP => {
	INFER => 'PSEUDO', 
	DEF => 0,
	TEST => '>',
	EVAL => 4,
	RANK => -1,
	ORDER => 10,
	SOFA => 'pseudogene',
	#SOFA_PARTOF => 'decayed_exon'	
	SOFA_PARTOF => undef,
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',
	STATUS => undef	
    },
    
    # ka/ks ~ 1
    KAKS2 => {
	INFER => 'FREE',
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => 1,
	RANK => 0,
	ORDER => 9,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'dubious'	
    },
    
    # hits at AA level to ty ORF	
    TY => {
	INFER => 'REPEAT',
	DEF => $INFINITY,
	TEST => '<',
	EVAL => 1e-5,
	RANK => 0,
	ORDER => 8,
	SOFA => 'repeat_region',
	SOFA_PARTOF => undef,
	GENBANK_ASN1 => 'mobile_element', # 'repeat_region'
	GENBANK_ASN1_PARTOF => undef,
	STATUS => undef
    },
    
    # hits at DNA level to ty LTR 
    LTR => {
	INFER => 'REPEAT',
	DEF => $INFINITY,
	TEST => '<',
	EVAL => 1e-5,
	RANK => 0,
	ORDER => 7.5,
	SOFA => 'repeat_region',
	SOFA_PARTOF => undef,
	GENBANK_ASN1 => 'LTR', # 'repeat_region',
	GENBANK_ASN1_PARTOF => undef,
	STATUS => undef
    },

    # hits at AA level to large multi gene family 	
    HCNF => {
	INFER => 'REPEAT',
	DEF => 0,
	TEST => '>',
	EVAL => 0.25, # homology to more than 1/4 of genes in cluster 
	RANK => 0,
	ORDER => 7,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'homology'
    },
    
    # ka/ks << 1
    KAKS => {
	INFER => 'REAL',
	DEF => $INFINITY,
	TEST => '<',
	EVAL => 0.8,
	RANK => 1,
	ORDER => 6,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'conservation'	
    },

    # partial alignment to homolog  
    HSP => {
	#INFER => 'FRAG',
	INFER => 'HSP',
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => .5,
	RANK => 1,
	ORDER => 5,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS', ## Correct? 
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'homology'			
    },

    # HMM hit to YGOB DB 
    YGOB => {
	INFER => 'REAL',
	DEF => $INFINITY,
	TEST => '<',
	EVAL => 1e-5,
	RANK => 1,
	ORDER => 4,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'homology'	
    },
    
    # hits db of real AA 
    AA => {
	INFER => 'REAL',
	DEF => $INFINITY,
	TEST => '<',
	EVAL => 1e-5,
	RANK => 1,
	ORDER => 3,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'homology'			
    },

    # weak hit but synteny support 
    SYNT => {
	INFER => 'REAL',
	DEF => 0,
	TEST => '>',
	EVAL => 2,
	RANK => 1,
	ORDER => 2,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'synteny'	
    }, 
    
    # hits db of real AA 
    NCBI => {
	INFER => 'REAL',
	DEF => $INFINITY,
	TEST => '<',
	EVAL => 1e-10,
	RANK => 1,
	ORDER => 1,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'homology'		
    },
    
    
    # hsp but no valid structure 
    #HSP => {
#	INFER => 'HSP',
#	DEF => $INFINITY,
#	TEST => '<',
#	EVAL => 1e-5,
#	RANK => 1.5,
#	ORDER => 4,
#	SOFA => 'gene'
 #   },
    
    # lda evidence 
    LDA => {
	INFER => 'NOVEL',
	DEF => -$INFINITY,
	TEST => '>',
	EVAL => .8,
	RANK => 2,
	ORDER => 0,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'hypothetical'		
    },  
    
    # no homology but far to long for chance 
    LENGTH => {
	INFER => 'NOVEL', 
	
	DEF => 300, # BEWARE!!!
	# need somewhere to store value of _orf_min and here
	# is ok for the moment-- $orf->data('LENGTH') is never
	# queried, only $orf->length. # BEWARE!!!
	
	TEST => '>',
	EVAL => 600,
	RANK => 2,
	ORDER => -1,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'hypothetical'		
    },              
    
    # count introns as evidence of real 
    INTRONS => {
	INFER => 'NOVEL',
	DEF => 0,
	TEST => '>',
	EVAL => 1,
	RANK => 2.5,
	ORDER => -2,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'dubious'	
    },
    
    # structure 
    STRUCT => {
	INFER => 'HYPO',
	DEF => 0,
	TEST => '=',
	EVAL => 2,
	RANK => 3,
	ORDER => -3,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => 'CDS',
	GENBANK_ASN1_PARTOF => 'exon',	
	STATUS => 'dubious'	
    },
    
    # intergenics
    NONE => {
	INFER => 'INTER',
	RANK => $INFINITY,
	ORDER => -4,
	SOFA => 'gene',
	SOFA_PARTOF => 'exon',
	GENBANK_ASN1 => undef, # **FIX
	GENBANK_ASN1_PARTOF => undef,  # **FIX
	STATUS => 'dubious'
    }
    ); 

@EVIDENCE = grep {!/NONE/} sort { $EVIDENCE{$a}->{ORDER} <=> $EVIDENCE{$b}->{ORDER} } keys %EVIDENCE;

sub _evidence {
    my $self = shift;
    my $args = {@_};
    $args->{'-query'} = uc($args->{'-query'});
    $self->throw("Invalid evidence: ".$args->{'-evidence'})
	unless exists $EVIDENCE{$args->{'-evidence'}};
    $self->throw("Invalid query: ".$args->{'-query'})
	unless exists $EVIDENCE{$args->{'-evidence'}}->{$args->{'-query'}};		
    return $EVIDENCE{$args->{'-evidence'}}->{$args->{'-query'}};
}

sub _set_eval {
    my $self = shift;
    my $args = {@_};		
    $self->throw("Invalid evidence: ".$args->{'-evidence'})
	unless exists $EVIDENCE{$args->{'-evidence'}};
    $self->throw("Cannot set 'EVAL' to ".$args->{'-set'})
	unless $args->{'-set'} =~ /\d/;
    $EVIDENCE{$args->{'-evidence'}}->{'EVAL'} = $args->{'-set'};
    return $self;
}

sub _set_default {
    my $self = shift;
    my $args = {@_};		
    $self->throw("Invalid evidence: ".$args->{'-evidence'})
	unless exists $EVIDENCE{$args->{'-evidence'}};
    $self->throw("Cannot set 'DEFAULT' to ".$args->{'-set'})
	unless $args->{'-set'} =~ /\d/;
    $EVIDENCE{$args->{'-evidence'}}->{'DEF'} = $args->{'-set'};
    return $self;
}

1;
