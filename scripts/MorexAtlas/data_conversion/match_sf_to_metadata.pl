#!/usr/bin/perl


use strict;

open(META, "MorexAtlas_allProjects_paper_18Aug23_edited_samples.txt") || die;

my $log_file = "match_sf_to_metadata.log";

open (LOG, ">$log_file");

open(METASQL, ">MorexAtlas_allProjects_paper_18Aug23_edited_samples.sql");

print METASQL "use morexgeneatlas;\n";

my %meta;
my %sample_titles;
my $cultivar_count = 0;
my %cultivar_list;
my $tissue_count = 0;
my %tissue_list;
my $devstage_count = 0;
my %devstage_list;
my $exptcond_count = 0;
my %exptcond_list;

while(<META>){

    if($_ =~ /^expt_id/){
        next;
    }

    $_ =~ s/\s+$//;

    $_ =~ s/\"//g;

    my($expt_id, $sample_type, $study_accession, $study_keywords, $run_accession_id, $run_accessions_included, $sample_accession, $tax_id, $sample_description, $sample_title, $cultivar, $tissue_type, $dev_stage, $bio_replicate, $expt_conditions) = split(/\t/, $_, 15);

    if (defined $cultivar_list{$cultivar}){
    } else {
        $cultivar_count++;
        $cultivar_list{$cultivar} = $cultivar_count;
    }
    if (defined $tissue_list{$tissue_type}){
    } else {
        $tissue_count++;
        $tissue_list{$tissue_type} = $tissue_count;
    }
        if (defined $devstage_list{$dev_stage}){
    } else {
        $devstage_count++;
        $devstage_list{$dev_stage} = $devstage_count;
    }
        if (defined $exptcond_list{$expt_conditions}){
    } else {
        $exptcond_count++;
        $exptcond_list{$expt_conditions} = $exptcond_count;
    }

    print METASQL "insert into sra_run_accessions(experiment_id, sample_type, study_accession, study_keywords, run_accession_id, run_accessions_included, sample_accession, taxon_id, sample_description, sample_title, cultivar, tissue_type, dev_stage, bio_replicate, expt_conditions, cultivar_id, tissue_id, dev_stage_id, expt_cond_id) values (\"$expt_id\", \"$sample_type\", \"$study_accession\", \"$study_keywords\", \"$run_accession_id\", \"$run_accessions_included\", \"$sample_accession\", \"$tax_id\", \"$sample_description\", \"$sample_title\", \"$cultivar\", \"$tissue_type\", \"$dev_stage\", \"$bio_replicate\", \"$expt_conditions\", \"$cultivar_list{$cultivar}\", \"$tissue_list{$tissue_type}\", \"$devstage_list{$dev_stage}\", \"$exptcond_list{$expt_conditions}\");\n";

    # Note some samples are named after the run_accession_id, some by the sample_title
    $meta{$run_accession_id} = $_;
    $sample_titles{$sample_title} = $_;
}

# make a file for loading the tab-delimited text into the db
open(SQL, ">load_tpm_data.sql");
print SQL "use morexgeneatlas;\n";

#find the quant.sf salmon output files
my $files = `ls <path to files>/morex_tpm_files/PR*/*/*quant.sf`;

my @files = split(/\s+/, $files);

foreach my $file(@files){

    my @path = split(/\//, $file);

    my $sf_name = $path[11];
    $sf_name =~ s/\.quant\.sf//;
    $sf_name =~ s/\_quant\.sf//;

    my $metadata;

    if(defined $meta{$sf_name}){

        $metadata = $meta{$sf_name};

    } elsif (defined $sample_titles{$sf_name}){

        $metadata = $sample_titles{$sf_name};

    } else {
        print STDERR "$file: $sf_name: cannot find meta $sf_name or sample_title $sf_name\n";
    }

    my($expt_id, $sample_type, $study_accession, $study_keywords, $run_accession_id, $run_accessions_included, $sample_accession, $tax_id, $sample_description, $sample_title, $cultivar, $tissue_type, $dev_stage, $bio_replicate, $expt_conditions) = split(/\t/, $metadata, 15);

    # Want to make a file for the sql bulk load
    my $quant_sql_file = $study_accession . "_" . $run_accession_id . "_load.txt";

    print SQL "load data local infile \"$quant_sql_file\" into table tpm_values;\n";

    open(TXT, ">$quant_sql_file");

    open(SF , "$file");

    my $gene_count = 0;
    my %gene_list;

    while(<SF>){
        if($_ =~ /^Name/){
            next;
        }

        $_ =~ s/\s+$//;

        my ($transcript_name, $length, $effective_length, $tpm, $numreads) = split(/\t/, $_, 5);

        my ($gene_name, $transcript_number) = split(/\.\d+$/, $transcript_name, 2);

        # if the gene has been already assigned a id number - don't do anything
        if (defined $gene_list{$gene_name}){

        # otherwise add to the count and assign to the hash for the gene
        } else {
            $gene_count++;
            $gene_list{$gene_name} = $gene_count;

        }

        print TXT "$expt_id\t$transcript_name\t$gene_list{$gene_name}\t$gene_name\t$tpm\t$numreads\n";

    }

    close TXT;
}

