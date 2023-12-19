#!/usr/bin/perl

use strict;
use CGI;
use GD;

use lib '/var/www/html/morexgeneatlas';
use morexgeneatlas;
use morexgeneatlas_wrapper;



my $query    = CGI->new();
my $sequence    = $query->param("sequence");
my $e_cutoff    = $query->param("e_cutoff");
my $search_type = $query->param("type");
my $database    = $query->param("database");


my $db_dir = '/var/www/html/morexgeneatlas/blast_dbs';
my $temp_dir = '/var/www/html/morexgeneatlas/temp';
my $blast_exe_location = '/var/www/html/morexgeneatlas/ncbi-blast-2.13.0+/bin';


my $db_query = morexgeneatlas_wrapper->new('morexgeneatlas');

my %blast_results;
my %m_align;
my %m_lengths;
my $total_hits;
my $no_aln;
my $db_count;
my $query_length;
my $query_name;

my %dbs = (
  "hvmxrtd" => "iso_rna_hvmx_merged.fasta",
  "hvmxrtdpp" => "iso_rna_hvmx_merged_transuite_v1_transfeat_pep.fasta",
  "morexgenome" => "210316_Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.fasta",
  "BaRT2v18" => "BaRT2v18.fa",
  "180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1" => "180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta", 
  "150831_barley_pseudomolecules" => "150831_barley_pseudomolecules.fasta",
	"Hv_IBSC_PGSB_r1_transcripts" => "Hv_IBSC_PGSB_r1_transcripts_all.fa", 
	"BaRT.1.0_Unsplit_exons" => "BaRT.1.0_Unsplit_exons_renamed.fasta",
    "BaRT.1.0_Peptides" => "BaRT.1.0_.fa.transdecoder.pep_renamed.fasta"
);


my %db_title = (
  "hvmxrtd" => "HvMxRTD", 
  "hvmxrtdpp" => "HvMxRTD Predicted Proteins",
  "morexgenome" => "Morex Genome Reference",
  "BaRT2v18" => "BaRT2v18 Barley Reference Transcripts",
  "180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1" => "The Barke genome assembly", 
  "150831_barley_pseudomolecules" => "Barley MTP Morex Assembly",
	"Hv_IBSC_PGSB_r1_transcripts" => "Barley MTP predicted Transcripts", 
  "BaRT.1.0_Unsplit_exons" => "BaRT.1.0u Barley Reference Transcripts",
  "BaRT.1.0_Peptides" => "BaRT.1.0u Barley Reference Proteins"
                );


my %db_type = (
  "hvmxrtd" => "n",
  "hvmxrtdpp" => "p",
  "morexgenome" => "n",
  "BaRT2v18" => "n",
  "180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1" => "n",
  "150831_barley_pseudomolecules" => "n",
	"Hv_IBSC_PGSB_r1_transcripts" => "n", 
	"BaRT.1.0_Unsplit_exons" => "n",
  "BaRT.1.0_Peptides" => "p"
               );


# Figure out what flavour of blast is needed
my $blast_exe;
if (($db_type{$database} eq "n") && ($search_type eq "nuc")){
  $blast_exe = "blastn";
}
if (($db_type{$database} eq "p") && ($search_type eq "nuc")){
  $blast_exe = "blastx";
}
if (($db_type{$database} eq "n") && ($search_type eq "pro")){
  $blast_exe = "tblastn";
}
if (($db_type{$database} eq "p") && ($search_type eq "pro")){
  $blast_exe = "blastp";
}



print "Content-type: text/html\n\n";

morexgeneatlas::printHeader("BLAST Search Results");


# Check the user input for FASTA format
if (($sequence eq "") || ($sequence =~ /^\s+$/)){
  # if nothing is in the sequence box
  print "Nothing in the sequence<BR>\n";
  morexgeneatlas::printFooter();
  die("No sequence given");

}


print "<div class='container-fluid'>
      <div class='row'>
        <nav class='col-md-2 d-none d-md-block grey sidebar'>
          <div class='sidebar-sticky'>
            <ul class='nav flex-column'>

                <li class='nav-item nav-link'  id='nav-title'><img class='img-fluid' src='images/morexgeneatlas-logo.png' width=150'></li>

                <li class='nav-item nav-link'  id='nav-header'><b>Utilities Menu</b></li>
                <div class='list-group'>
                  <li class='nav-item'><a class='nav-link' href='index.html'>Home</a></li>
                  <li class='nav-item'><a class='nav-link active' href='blast.html'>Homology Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='keyword.html'>Annotation Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='download.html'>Bulk Data Download</a></li>
                  <li class='nav-item'><a class='nav-link active' href='relator.html'>Transcript Relationships</a></li>
                </div>
                <li class='nav-item nav-link'  id='nav-title'></li>
                <li class='nav-item nav-link'  id='nav-header'><b>Links</b></li>
                <li class='nav-item nav-link'  id='nav-title'></li>
                <li class='nav-item'><a class='nav-link ext-links' href='contact.html'>Citation / About Authors</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://ics.hutton.ac.uk/barleyrtd'>BarleyRTD Website</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://ics.hutton.ac.uk'>Information and Computing Sciences \@Hutton</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='http://www.hutton.ac.uk'>The James Hutton Institute</a></li>
                <li class='nav-item'><a class='nav-link ext-links' href='https://www.barleyhub.org/'>International Barley Hub</a></li>
                <li class='nav-item nav-link'  id='nav-title'></li>

              </ul>
          </div>
        </nav>
        <main role='main' class='col-md-9 ml-sm-auto col-lg-10 pt-3 px-4'>\n\n
        
        <div class='data-header'>BLAST Results</div>\n<div class='data-body'>";


my @seqs = split(/\n>/, $sequence);

my $number_of_seqs = @seqs;

if ($number_of_seqs > 30){
  
  print "<h3>Too many sequences submitted ($number_of_seqs). Please submit no more than 30 sequences.</h3><br><a href='javascript:history.go(-1)'> Go Back</a><br><br></div>\n";
  
  morexgeneatlas::printFooter();
  
  die ("morexgeneatlas/blast.cgi: Too many sequences submitted\n");
  
}


# Make a temporary file name for the sequence based on the time/date of submission
my $datestamp = `date '+%d%m%y_%H%M%S'`;
chomp $datestamp;

my $fasta_file = "temp_" . $datestamp . ".fas";

open(TEMP, ">$temp_dir/$fasta_file") || die ("Cannot open temp file '$temp_dir/$fasta_file'\n");
print TEMP "$sequence";
close TEMP;


`chmod 666 $temp_dir/$fasta_file`;


# Run the BLAST search and capture the output


my $blast_out = `$blast_exe_location/$blast_exe -num_threads 4 -db $db_dir/$dbs{$database} -query $temp_dir/$fasta_file -html -evalue $e_cutoff -num_descriptions 50 -num_alignments 50`|| die("Cannot execute blast search : $blast_exe -num_threads 4 -db $db_dir/$dbs{$database} -query $temp_dir/$fasta_file -html -evalue $e_cutoff -num_descriptions 50 -num_alignments 50\n");


# Make a filename for the blast output with added hyperlinks
my $blast_file = "blast_" . $datestamp . ".html";

# Make a html file from the blast output with links to the database
my $linked_html = &create_linked_html($blast_out, $database);

open(TEMP, ">>$temp_dir/$blast_file") || die("Cannot open temp file '$temp_dir/$blast_file'\n");
print TEMP "$linked_html\n";
close TEMP;

`chmod 664 $temp_dir/$blast_file`;


my $seq_type = "nucleotide";
if ($search_type eq "pro"){
  $seq_type = "protein";
}

print "'$number_of_seqs' $seq_type sequences have been submitted<br><br>
<a href='temp/$blast_file'>View full BLAST results</a><br><br>\n";


my @blast_out = split(/<b>Query=<\/b> /, $blast_out);

my $number_results = @blast_out - 1;


for (my $i = 1; $i < @blast_out; $i++){
  
  #print "<pre>Submitting for processing $i : $blast_out[$i]</pre><br>\n\n";
  

  # Calculate how many hits there are to each database and add them up
  my ($query_name, $query_length, $no_hits, $no_aln) = &parse_blast_alignments($blast_out[$i], $database);

  # Don't know if there will be any hits to this sequence, so check that there are hits before drawing a BLAST diagram
  if ($no_hits !=0){

    my $imagefile = $query_name ."_" . $dbs{$database} . "_" . $datestamp . ".png";

    &drawBlastDiagram($imagefile, $query_name, $query_length, $no_hits);
  
  } else {
    
    my $units = "bp";
    if ($search_type eq "pro"){
      $units = "aa";
    }
    print "<b>$query_name ($query_length $units): No hits found</b><br><br>\n";
  }

}

print "</div>\n";


morexgeneatlas::printFooter();


sub parse_blast_alignments {

  my ($result, $db) = @_;
  

  %blast_results = ();
  %m_align = ();
  %m_lengths = ();
  $total_hits = 0;
  $no_aln = 0;
  $db_count = 0;


  my ($search_header, $results) = split(/Sequences producing/, $result);


  my @header = split(/\n/, $search_header);

  ($query_name, my $rest) = split(/\s+/, $header[0]);

  
  #print "header '$search_header'<br>\n";
  foreach my $line(@header) {
    if ($line =~ /Length=/) {
      $query_length = $line;
      $query_length =~ s/\D+//g;

    }

    if ($line =~ /Query=/) {
      (my $poo, $query_name, $rest)  = split(/\s+/, $line, 3);

    }
   }

    my @matches = split(/<a name/, $results);
    my $m_count;
    my $no_hsps;

    if ($result =~ /No hits found/) {

    } else {


    for(my $i = 1; $i < @matches; $i++) {

      $m_count++;

      my ($m_title, @alignments) = split(/Score/, $matches[$i]);

      (my $poo, $m_title) = split(/\<\/a\> /, $m_title);

      my ($m_name, $description) = split(/\s+/, $m_title, 2);
    
      ($description, my $m_length) = split(/Length=/, $description);

      $description =~ s/\'//g;
      $description =~ s/\n//g;
      $description =~ s/\s+/ /g;
    
      my $units = "bp";
      if ($db_type{$database} eq "p"){
        $units = "aa";
      }
    
    
      $m_length =~ s/\s+//g;
    
      $description = "($m_length $units) " . $description;
    
      $m_lengths{$m_name} = $m_length;

      my $aln_count = 0;
      my $highest_score = 0;
      my @m_locs;
      my @query_locs;
      my $dir;

      foreach my $aln(@alignments) {

        $aln_count++;

        my @aln = split(/\n/, $aln);

        my @score_line = split(/\s+/, $aln[0]);

        my $m_score = $score_line[2];
        my $m_evalue = $score_line[7];


        my ($ident, $pos) = split(/,/, $aln[1]);

        ($poo, my $identities) = split(/ = /, $ident);

        # Figure out direction of hit
        # In BLASTN get Strand
        if ($aln[2] =~ /Strand/){

          my @strands = split(/\s+/, $aln[2]);
      
          $dir = $strands[5];

        }
      
        # In TBLASTN (protein seq vs translated nucleotide db) get Frame
        if ($aln[2] =~ /Frame/){
        
          my @frames = split(/ = /, $aln[2]);
      
          if ($frames[1] =~ /-/){
          
            $dir = "Minus";
          
          } else {
          
            $dir = "Plus";
          }
        }


        my $query_aln = "";
        my $subject_aln = "";
      
        foreach my $line(@aln) {

	        if ($line=~ /Query/) {
	          $query_aln = $query_aln . $line;
	        }

	        if ($line=~ /Sbjct/) {
	          $subject_aln = $subject_aln . $line;
	        }

        }

        my @words = split(/\s+/, $query_aln);

        my $query_start = $words[1];
        my $query_stop = $words[@words-1];
      

        @words = split(/\s+/, $subject_aln);
        my $subject_start = $words[1];
        my $subject_stop = $words[@words-1];
      
      
        # Find highest scoring segment and use as a rank for displaying later
        if ($aln_count == 1) {
        
	        $highest_score = $m_score;

        }
      
        #Find the highest and lowest positions in the hsps
        push (@m_locs, $subject_start);
        push (@m_locs, $subject_stop);
      
        push (@query_locs, $query_start);
        push (@query_locs, $query_stop);

        $identities =~ s/, Gaps//;
        $m_evalue =~ s/,$//;

        my $results_details = join("\t", $m_evalue, $query_start, $query_stop);

        $blast_results{$query_name}{$highest_score}{$m_name}{$aln_count} = $results_details;


      }
    
      @m_locs = sort{$a <=>$b}(@m_locs);
      @query_locs = sort{$a <=>$b}(@query_locs);

      $no_hsps =+ $aln_count;
    
      # store the ends of the query and hit for drawing later
    
      $m_align{$query_name}{$highest_score}{$m_name} = "$query_locs[0],$query_locs[@query_locs-1],$m_locs[0],$m_locs[@m_locs-1],$dir,$description";


    }

  }

  return ($query_name, $query_length, $m_count, $no_hsps);

}



sub create_linked_html{
  
  my($blast_out) = @_;
    
  my @linked_html = split(/\n/, $blast_out);
  
  my $linked_html;
  my $query_name;
  
  foreach my $line(@linked_html){
    
    if ($line =~ /Query=<\/b> /){
      
      (my $poo1, $query_name, my $poo2) = split(/\s+/, $line, 3);
      
    }
    
    if ($line =~ /><a name=/){
      
      my ($poo, $match) = split(/<\/a> /, $line);
      
      my ($m_name, $desc) = split(/\s+/, $match, 2);
      
      my $target_name = $query_name . "_" . $m_name;

      if (($m_name =~ /Hv_Mx/) || ($m_name =~ /BaRT/) || ($m_name =~ /BART/) || ($m_name =~ /HORVU/)){

        $line = "$poo</a><a name = $target_name></a><b><a href='https://ics.hutton.ac.uk/morexgeneatlas/search_ids.cgi?seq_name=$m_name'>$m_name</a></b> <a href='javascript:history.go(-1)'>Go Back</a>\n";

      } else  {

        $line = "$poo</a><a name = $target_name></a><b>$m_name</b> <a href='javascript:history.go(-1)'>Go Back</a>\n";

      }
      
      
      
    }
    
    $linked_html .= "$line\n";
    
  }
  
  return $linked_html;
}


sub drawBlastDiagram {
  
  my ($imagefile, $query_name, $query_length, $total_hits) = @_;
  


############################# Make image of the BLAST hits ############################################

  my $diagram_width = 830;
  
  my $scale = $diagram_width / $query_length; # want display to be 600 pixels wide
  
  #Need to be able to change the detail of the scale with different query sequence lengths
  my $main_tick;
  my $minor_tick;
  
  if ($query_length < 1000) {
    
    $main_tick = 100;
    $minor_tick = 20;

  } elsif (($query_length >= 1000) && ($query_length < 2000)) {
    
    $main_tick = 200;
    $minor_tick = 50;

  } elsif (($query_length >= 2000) && ($query_length < 3000)) {
    
    $main_tick = 200;
    $minor_tick = 100;

  } elsif ($query_length >= 3000) {
    
    $main_tick = 500;
    $minor_tick = 100;

  }
  
  my $y_gap = 6;         # How much of a gap is between each match
  my $block_height = 12; # How much height each match takes up on diagram
  my $text_offset = 240; # Leave space for the names of the contigs on the left hand side


  my $y = 30;     # Start placing matches at 30 pixels down from top of image

  my $image_height = ($total_hits * ($block_height + $y_gap)) + 80; # each match consists of the height of the block and the gap under it
  
  my $desc_width = 0;

  my $image_width = ($query_length * $scale) + $text_offset + $desc_width + 10;

  my $image = GD::Image->newPalette($image_width, $image_height);

  #Make colours for use in the image
  my $black     = $image->colorAllocate(0, 0, 0);
  my $lightgray = $image->colorAllocate(236, 242, 244);
  my $darkgray  = $image->colorAllocate(215, 215, 215);
  my $offwhite  = $image->colorAllocate(230, 230, 230);
  my $white     = $image->colorAllocate(255, 255, 255);
  my $red       = $image->colorAllocate(255, 0, 0);
  my $link_col  = $image->colorAllocate(0, 134, 179);
  my $logo_blue       = $image->colorAllocate(0, 122, 194);
  my $logo_green      = $image->colorAllocate(67, 149, 57);
  my $logo_lightblue  = $image->colorAllocate(129,178,154);
  my $logo_lightgreen = $image->colorAllocate(130, 206, 193);
  my $yellow          = $image->colorAllocate(242,204,143); 

  $image->filledRectangle(0, 0, $image_width, $image_height, $lightgray); # Background colour square
  
  $image->string(gdMediumBoldFont, 20 , 10, "$query_name", $logo_blue);
  
  my $label_x = $text_offset + $diagram_width + 10;
  
  my $units = "bp";
  if ($search_type eq "pro"){
        $units = "aa";
  }


  ############################# make the scale along the top ##########################################
  my $scale_width = ($query_length * $scale) + $text_offset;
  $image->filledRectangle($text_offset, 10, $scale_width ,12, $logo_green);
  for(my $n = 0; $n < $query_length; $n++ ){
    if ($n%$main_tick == 0) {
      my $scale_tick = ($n * $scale) + $text_offset;
      ($scale_tick, my $float) =split(/\./, $scale_tick);
      $image->line($scale_tick, 10, $scale_tick ,20, $logo_green);
      my $scale_label = $scale_tick - 5;
      $image->string(gdSmallFont, $scale_label ,20, "$n", $logo_green);
    } elsif ($n%$minor_tick == 0) {
      my $scale_tick = ($n * $scale) + $text_offset;
      ($scale_tick, my $float) =split(/\./, $scale_tick);
      $image->line($scale_tick, 10, $scale_tick ,15, $logo_green);
    }
  }
  ###################################################################################################
  #Print the hits as rectangles

  my $imap_name = "blast_map_" . $query_name;
    
  my $imap = "<map name='$imap_name'>\n";
  
  my $count = 0;
  
  foreach my $q_name(sort keys %blast_results){
    
    #Set the colour for the contig type
    my $assembly_col;
    
    $assembly_col = $image->colorAllocate(217, 102, 255);

    $y = $y + 20;
    
    my $title_height = $y - 5;
    

    
    foreach my $overall_e(sort {$b <=> $a}keys %{$blast_results{$q_name}}){
      foreach my $member(keys %{$blast_results{$q_name}{$overall_e}}){

	      $count++;

	      # Draw background white rectangle and the name of the hit sequence
	      $y = $y + ($block_height + $y_gap);
	      my $y_end = $y + $block_height;
  
        my $end_db_area = $image_width  - 5;

	      my $box_end = $image_width - $desc_width - 10;
	      $image->filledRectangle($text_offset, $y, $box_end ,$y_end, $white);
        
        # This is the rectangle that would be the imap for the link to the sequence
	      $image->filledRectangle(5, $y, $text_offset ,$y_end, $offwhite);

	      my $text_height = $y - 1;

        if (($member =~ /Hv_Mx/) || ($member =~ /BaRT/) || ($member =~ /BART/) || ($member =~ /HORVU/)){ # do not make a hyperlink to the pseudomolecules for Morex or Barke yet

          $imap .= "<area href='search_ids.cgi?seq_name=$member' coords='5, $y, $text_offset, $y_end'>";

        }

        $image->string(gdSmallFont, 5 ,$text_height, "$member", $link_col);


        my($lowest_query_start,$highest_query_stop,$lowest_match_start,$highest_match_stop,$strand,$m_desc) = split(/,/, $m_align{$q_name}{$overall_e}{$member});
        

        my $dir = "->";
        # Protein vs Protein has no direction in output - therefore default is "->"
        if ($strand eq "Minus") {
            $dir = "<-";
        }
        
        $image->string(gdSmallFont, 210 ,$text_height, "$dir", $black);


        # if the hit is a peptide, have to multiply hit length, hit start, and hit stop by 3
        # to map it back onto the query nucleotide sequence
        if (($search_type eq "nuc") && ($db_type{$database} eq "p")) {

            $m_lengths{$member} = $m_lengths{$member} * 3;
            $lowest_match_start = ($lowest_match_start * 3) - 3;
            $highest_match_stop = $highest_match_stop * 3;

        }
        
        # if the hit is a translated nucleotide, the hit length, start and stop need to be divided by 3
        # to map onto a protein query
        if (($search_type eq "pro") && ($db_type{$database} eq "n")) {


            #print "Searching protein vs translated nucleotide - divide hit by three<br>\n";
            
            $m_lengths{$member} = int($m_lengths{$member} / 3);
            $lowest_match_start = int(($lowest_match_start / 3));
            $highest_match_stop = int($highest_match_stop / 3);

        }
        
        
        
        # Relate the ends of the hit back to the query sequence by adding/subtracting the
        # ends of the hit to the query locations
        my $m_start = $lowest_query_start - $lowest_match_start + 1;
        my $m_stop = $highest_query_stop + ($m_lengths{$member} - $highest_match_stop + 1) - 1;


	      # Cope with cases where the hit is to the rev-comp of the query sequence
        if ($dir eq "REV") {
        
          $m_start = $lowest_query_start - ($m_lengths{$member}- $highest_match_stop + 1) + 1;
          $m_stop = $highest_query_stop + $lowest_match_start - 1;
                
	      }
        
        my $diff = $m_stop - $m_start + 1;

        # If the start and stop of the hit sequence is beyond the ends of the query
        # need to set the ends to the ends of the query and later add "+" to show that the match sequence
        # extends beyond this limit
        my $draw_match_start = "";
            
        if ($m_start <= 0) {
          $m_start = 1;
          $draw_match_start = "true";
        }
            
        my $draw_match_stop = "";
        if ($m_stop > $query_length) {
          $m_stop = $query_length;
          $draw_match_stop = "true";
        }
	
        # Calculate where the rectangle should start/stop
	      my $m_end_point = (($m_stop) * $scale) + $text_offset;
	      my $m_start_point = (($m_start) * $scale) + $text_offset;


        $image->filledRectangle($m_start_point, $y, $m_end_point ,$y_end, $yellow);

	      foreach my $aln_count(sort {$b <=> $a} keys %{$blast_results{$q_name}{$overall_e}{$member}}){

	        my ($evalue, $query_start, $query_stop) = split(/\t/, $blast_results{$q_name}{$overall_e}{$member}{$aln_count});

	        my $dir = "for";
	        my $hsp_length = $query_stop - $query_start;
	        my $start = $query_start;
	        my $stop = $query_stop;

	        # Cope with cases where the hit is to the rev-comp of the query sequence
	        if ($query_start > $query_stop) {
	         $dir = "rev";
	         $hsp_length = $query_start - $query_stop;
	          $start = $query_stop;
	          $stop = $query_start;
	        }
	
      	  # Calculate where the rectangle should start/stop
	        my $end_point = (($stop) * $scale) + $text_offset;
	        my $start_point = (($start) * $scale) + $text_offset;

	        # Draw the rectangle for the hit
	        $image->filledRectangle($start_point, $y, $end_point ,$y_end, $logo_lightblue);
          

          
          # Only draw HSP e-value if the HSP is wide enough to take the print
          my $length_hsp = $end_point - $start_point;
          
          if ($length_hsp >= 40) {

            #Draw the total length of the section of hit in the centre of the rectangle
            my $dist_xpos = ((($hsp_length / 2) + $start) * $scale) + $text_offset - 10;
            
            
            if ($evalue =~ /^e/){
              $evalue = "1" . $evalue;
            }
            
            $image->string(gdSmallFont, $dist_xpos ,$text_height, "$evalue", $black);
          }
          
          # If the ends of the hit go beyond the ends of the query, show this with a "+" sign
          my $end_width = $image_width - $desc_width - 20;

          if ($draw_match_start eq "true") {

        	  $image->string(gdSmallFont, $text_offset , $text_height, "+", $black);

          }
            
          if ($draw_match_stop eq "true") {

            $image->string(gdSmallFont, $end_width ,$text_height, "+", $black);

          }
	      }
      }
    }
  }

  $imap .= "</map>\n";

  
  # Open the PNG file and print the GD output to it

  open(PNG, ">/var/www/html/morexgeneatlas/blast_net_images/$imagefile");
  binmode PNG;
  print PNG $image->png();
  close PNG;

  print "<img  src='blast_net_images/$imagefile' usemap='#$imap_name' border=0><br><br>\n";
  print "$imap";

}


