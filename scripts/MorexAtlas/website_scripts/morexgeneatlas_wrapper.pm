use strict;
use DBI;
use Time::localtime;

package morexgeneatlas_wrapper;
require Exporter;

#------------------------------------------------------------------------------------------------------
# Constructor
#------------------------------------------------------------------------------------------------------

sub new {
  my ($proto, $database) = @_;	
  my $class = ref($proto) || $proto;
  my $self = {};				
  bless ($self, $class);
  $self->__initialisation($database);
  return $self;
}

sub __initialisation {
  my($self, $database) = (@_);
  $self->{'dbh'} = $self->__databaseConnection($database);
}

#-----------------------------------------------------------------------
# Database connection methods
#-----------------------------------------------------------------------

sub __databaseConnection {
  my ($self, $database) = (@_);


  my $dbName = 'morexgeneatlas';
  my $dbDriver = 'mysql';

  # Connect to the database
  $self->{'dbh'} = DBI->connect("DBI:mysql:$database;host=<host IP>", "<database name>", "<password>");



}


#------------------------------------------------------------------------------------
# doSQLStatement(SQL STATEMENT)
# 	Runs the SQL statement that is passed in and returns statements handle.
#------------------------------------------------------------------------------------
sub doSQLStatement {
  my ($self, $stmt) = (@_);
  my($sth) = $self->{'dbh'}->prepare($stmt) || die "Problems preparing the SQL statement:". DBI->errstr;;

  $sth->execute() || die "Problems running the SQL statement:". DBI->errstr;
  return $sth;
}


#////////////////////////////////////////////////////////////////////////////////////
#
# Data retrieval methods
#
#////////////////////////////////////////////////////////////////////////////////////




sub getRelatedGenes { #takes a reference_gene_id and finds all the transcripts related in transcriptome_relationships
  my($self, $gene_id) = @_;

  my $stmt = "SELECT r.query_tran_id, r.query_gene_id, r.transcript_dataset, r.reference_tran_id, r.reference_gene_id, r.reference_dataset, r.reference_match_length, r.query_tran_exon_number, r.query_tran_length, r.major_isoform, r.gffcompare_class_code 
	      FROM transcriptome_relationships as r, gmap_locations as g
	      where r.reference_gene_id = '$gene_id' and r.query_tran_id = g.query_tran_id
	      ";
	
  #print "$stmt;<br>\n";
  
  my $sth = $self->doSQLStatement($stmt);
  
  my %dataStructure;
  my $dataRef = \%dataStructure;
  

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $query_tran_id           = $hashRef ->{query_tran_id};
    my $query_gene_id           = $hashRef ->{query_gene_id};
    my $transcript_dataset      = $hashRef ->{transcript_dataset};
    my $reference_tran_id       = $hashRef ->{reference_tran_id};
    my $reference_gene_id       = $hashRef ->{reference_gene_id};
    my $reference_dataset       = $hashRef ->{reference_dataset};
    my $reference_match_length  = $hashRef ->{reference_match_length};
    my $query_tran_exon_number  = $hashRef ->{query_tran_exon_number};
    my $query_tran_length       = $hashRef ->{query_tran_length};
    my $major_isoform           = $hashRef ->{major_isoform};
    my $gffcompare_class_code   = $hashRef ->{gffcompare_class_code};

    my $joined_details = join("\t", $query_gene_id, $reference_tran_id, $reference_gene_id, $reference_dataset, $reference_match_length, $query_tran_exon_number, $query_tran_length, $major_isoform, $gffcompare_class_code);
    
    $dataStructure{$query_tran_id}{$transcript_dataset} = $joined_details;

  }
  
  return $dataRef;
  
}


sub getTranscriptSequencesByGene {
  my($self, $gene_id, $dataset) = (@_);


  my $stmt = "SELECT transcript_id, seq_length, transcript_sequence
	      FROM transcript_sequences
	      where gene_id = '$gene_id'
	      and dataset_name = '$dataset'
	      ";

  my $sth = $self->doSQLStatement($stmt);

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id      = $hashRef->{transcript_id};
    my $seq_length         = $hashRef->{seq_length};
    my $sequence           = $hashRef->{transcript_sequence};

    my $number = $transcript_id;
    $number =~ s/$gene_id\.//;

    my $joined_details = join("\t", $transcript_id, $seq_length, $sequence);
    
    
    $dataStructure{$number} = $joined_details;
  }
  return $dataRef;
}


sub getProteinSequencesBySeq {
  my($self, $gene_id, $dataset) = (@_);


  my $stmt = "SELECT protein_id, seq_length, protein_sequence
	      FROM protein_sequences
	      where gene_id = '$gene_id'
	      and dataset_name = '$dataset'
	      ";

  my $sth = $self->doSQLStatement($stmt);

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $protein_id         = $hashRef->{protein_id};
    my $sequence           = $hashRef->{protein_sequence};

    my $number = $protein_id;
    $number =~ s/$gene_id\.//;
    
    $dataStructure{$sequence} .= "$number,";
  }
  return $dataRef;
}


sub getPaddedTPMsJS_new{
  
  my ($self, $gene_name) = @_;

  my $stmt = "select v.transcript_name, v.experiment_id, v.tpm_value, s.study_accession, s.run_accession_id, s.cultivar, s.tissue_type, s.bio_replicate, s.expt_conditions, s.dev_stage
              from tpm_values as v, sra_run_accessions as s
              where v.gene_name = '$gene_name' 
              and s.experiment_id = v.experiment_id 
              ";

  #print "$stmt<br>\n";
	      
  my $sth = $self->doSQLStatement($stmt);
  
  my %dataStructure;
  my $dataRef = \%dataStructure;
  
  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id       = $hashRef->{transcript_name};
    my $expt_id             = $hashRef->{experiment_id};
    my $tpm_value           = $hashRef->{tpm_value};
    my $study_accession     = $hashRef->{study_accession};
    my $run_accession       = $hashRef->{run_accession_id};
    my $cultivar            = $hashRef->{cultivar};
    my $tissue              = $hashRef->{tissue_type};
    my $bio_replicate       = $hashRef->{bio_replicate};
    my $expt_condition      = $hashRef->{expt_conditions};
    my $dev_stage           = $hashRef->{dev_stage};

    my $transcript_number = $transcript_id;
    
    $transcript_number =~ s/$gene_name\.//;

    my $joined_details = join("\t", $tpm_value, $study_accession, $run_accession, $cultivar, $tissue, $bio_replicate, $expt_condition, $dev_stage);

    $dataStructure{$transcript_number}{$expt_id} = $joined_details;
  }
  return $dataRef;
}

sub getPaddedTPMsJSSubset{
  
  my ($self, $gene_id, $study_filter) = @_;

  my $stmt = "select v.transcript_id, v.experiment_id, v.tpm_value, s.study_accession, s.sample_accession, s.cultivar, s.tissue_type, s.bio_replicate, s.expt_conditions, s.dev_stage
              from tpm_values as v, genes as g, sra_run_accessions as s
              where g.gene_id = '$gene_id' 
              and v.id = g.id 
              and v.experiment_id = s.experiment_id
              $study_filter
              ";

  #print "$stmt<br>\n";
	      
  my $sth = $self->doSQLStatement($stmt);
  
  my %dataStructure;
  my $dataRef = \%dataStructure;
  
  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id       = $hashRef->{transcript_id};
    my $expt_id             = $hashRef->{experiment_id};
    my $tpm_value           = $hashRef->{tpm_value};
    my $study_accession     = $hashRef->{study_accession};
    my $sample_accession    = $hashRef->{sample_accession};
    my $cultivar            = $hashRef->{cultivar};
    my $tissue              = $hashRef->{tissue_type};
    my $bio_replicate       = $hashRef->{bio_replicate};
    my $expt_conditions     = $hashRef->{expt_conditions};
    my $dev_stage           = $hashRef->{dev_stage};

    my $transcript_number = $transcript_id;
    
    $transcript_number =~ s/$gene_id\.//;

    my $joined_details = join("\t", $tpm_value, $study_accession, $sample_accession, $cultivar, $tissue, $bio_replicate, $expt_conditions, $dev_stage);

    $dataStructure{$transcript_number}{$expt_id} = $joined_details;
  }
  return $dataRef;
}


sub getExperimentMetadata{
    my ($self) = @_;
  
  my $stmt = "SELECT experiment_id, study_accession, sample_accession, cultivar, tissue_type, dev_stage, bio_replicate, expt_condition
	      FROM sra_run_accessions";
  
  my $sth = $self->doSQLStatement($stmt);
  
  my %dataStructure;
  my $dataRef = \%dataStructure;
  
  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $experiment_id   = $hashRef->{experiment_id};
    my $study_accession = $hashRef->{sra_study_accession};
    my $sample_accession   = $hashRef->{sample_accession};
    my $cultivar        = $hashRef->{cultivar};
    my $tissue          = $hashRef->{tissue_type};
    my $dev_stage       = $hashRef->{dev_stage};
    my $bio_replicate   = $hashRef->{bio_replicate};
    my $expt_condition  = $hashRef->{expt_condition};

    my $joined_details = join("\t", $study_accession, $sample_accession, $cultivar, $tissue, $dev_stage, $bio_replicate, $expt_condition);
    
    $dataStructure{$experiment_id} = $joined_details;
  }
  return $dataRef;
}


sub getTranscriptListByGene {
  my($self, $gene_id) = @_;
  
  my $stmt = "SELECT transcript_id 
	      FROM transcript_sequences
	      where gene_id = '$gene_id'
	      ";
	
  
  my $sth = $self->doSQLStatement($stmt);
  
  my @dataStructure;
  my $dataRef = \@dataStructure;
  
  my $count;
  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id = $hashRef->{transcript_id}; 
    
    $dataStructure[$count] = $transcript_id;
    $count++;
  }

  
  return $dataRef;
  
}


sub searchGOAnnotation {

  my($self, $phrase, $dataset) = @_;

  my($stmt) = "SELECT transcript_id, gene_id, dataset_name, pannzer_annotation, go_terms
               FROM gene_annotation where $phrase";

  #print "$stmt;<br>\n";

  my ($sth) = $self->doSQLStatement($stmt);	

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id      = $hashRef->{transcript_id};
    my $gene_id            = $hashRef->{gene_id};
    my $dataset_name       = $hashRef->{dataset_name};
    my $pannzer_annotation = $hashRef->{pannzer_annotation};
    my $go_terms           = $hashRef->{go_terms};

    my $transcript_number = $transcript_id;
    $transcript_number =~ s/$gene_id\.//;

    my($joined_details) = join ("\t", $dataset_name, $pannzer_annotation, $go_terms);

    $dataStructure{$gene_id}{$transcript_number} = $joined_details;
  }
	
  return $dataRef

}

sub searchAllAnnotation {

  my($self, $go_phrase, $blast_phrase) = @_;

  my($stmt) = "SELECT transcript_id, gene_id, dataset_name, pannzer_annotation, go_terms
               FROM gene_annotation where $go_phrase";

  #print "$stmt;<br>\n";

  my ($sth) = $self->doSQLStatement($stmt);	

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id      = $hashRef->{transcript_id};
    my $gene_id            = $hashRef->{gene_id};
    my $dataset_name       = $hashRef->{dataset_name};
    my $pannzer_annotation = $hashRef->{pannzer_annotation};
    my $go_terms           = $hashRef->{go_terms};

    my $transcript_number = $transcript_id;
    $transcript_number =~ s/$gene_id\.//;

    $pannzer_annotation =~ s/\d+\;//;

    my($joined_details) = join ("\t", $dataset_name, $pannzer_annotation, $go_terms);

    $dataStructure{$gene_id}{$transcript_number}{"PANNZER"} = $joined_details;
  }

  my($stmt) = "SELECT query_name, gene_id, dataset_name, blast_db, hit_name, score, percent_id, evalue, description
               FROM representative_blast_hits where $blast_phrase and hit_rank = 1 order by gene_id";

  #print "$stmt;<br>\n";

  my ($sth) = $self->doSQLStatement($stmt);	

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $query_name  = $hashRef->{query_name};
    my $gene_id     = $hashRef->{gene_id};
    my $dataset     = $hashRef->{dataset_name};
    my $blast_db    = $hashRef->{blast_db};
    my $match_name  = $hashRef->{hit_name};
    my $score       = $hashRef->{score};
    my $percent_id  = $hashRef->{percent_id};
    my $evalue      = $hashRef->{evalue}; 
    my $description = $hashRef->{description};

    my $transcript_number = $query_name;
    $transcript_number =~ s/$gene_id\.//;

    my($joined_details) = join ("\t", $dataset, $match_name, "bitscore=$score, \%id=$percent_id, evalue=$evalue, <b>$description</b>");

    $dataStructure{$gene_id}{$transcript_number}{$blast_db} = $joined_details;
  }
	
  return $dataRef

}


sub getGOAnnotation {

  my($self, $gene, $dataset) = @_;

  my($stmt) = "SELECT transcript_id, pannzer_annotation, go_ids, go_terms FROM gene_annotation 
                where gene_id = '$gene'
                and dataset_name = '$dataset'";

  #print "$stmt;<br>";

  my ($sth) = $self->doSQLStatement($stmt);	

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $transcript_id      = $hashRef->{transcript_id};
    my $pannzer_annotation = $hashRef->{pannzer_annotation};
    my $go_ids             = $hashRef->{go_ids};
    my $go_terms           = $hashRef->{go_terms};

    my($joined_details) = join ("\t", $pannzer_annotation, $go_ids, $go_terms);

    my $transcript_number = $transcript_id;
    $transcript_number =~ s/$gene\.//;

    $dataStructure{$transcript_number} = $joined_details;
  }
	
  return $dataRef

}

sub getHomologies {

  my($self, $transcript, $dataset) = @_;

  my($stmt) = "SELECT blast_db, hit_name, score, percent_id, evalue, description
                FROM representative_blast_hits where query_name = '$transcript' and dataset_name = '$dataset' and hit_rank = 1";


  my ($sth) = $self->doSQLStatement($stmt);	

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $blast_db    = $hashRef->{blast_db};
    my $match_name  = $hashRef->{hit_name};
    my $percent_id  = $hashRef->{percent_id};
    my $evalue      = $hashRef->{evalue}; 
    my $score       = $hashRef->{score};
    my $description = $hashRef->{description};

    my($joined_details) = join ("\t", $score, $percent_id, $evalue, $description);

    $dataStructure{$blast_db}{$match_name} = $joined_details;
  }
    
  return $dataRef

}


sub searchContigDescription {

  my($self, $phrase, $dataset) = @_;

  my($stmt) = "SELECT query_name, gene_id, dataset_name, blast_db, hit_name, score, percent_id, evalue, description
               FROM representative_blast_hits where $phrase and hit_rank = 1 order by gene_id";

  #print "$stmt;<br>\n";

  my ($sth) = $self->doSQLStatement($stmt);	

  my %dataStructure;
  my $dataRef = \%dataStructure;

  while(my $hashRef = $sth->fetchrow_hashref('NAME_lc')){

    my $query_name  = $hashRef->{query_name};
    my $gene_id     = $hashRef->{gene_id};
    my $dataset     = $hashRef->{dataset_name};
    my $blast_db    = $hashRef->{blast_db};
    my $match_name  = $hashRef->{hit_name};
    my $score       = $hashRef->{score};
    my $percent_id  = $hashRef->{percent_id};
    my $evalue      = $hashRef->{evalue}; 
    my $description = $hashRef->{description};

    my($joined_details) = join ("\t", $gene_id, $dataset, $match_name, $score, $percent_id, $evalue, $description);

    $dataStructure{$query_name}{$blast_db} = $joined_details;
  }
	
  return $dataRef

}


