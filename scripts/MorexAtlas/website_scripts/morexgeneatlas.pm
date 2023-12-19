use strict;

package morexgeneatlas;
use morexgeneatlas_wrapper;
use GD;
require Exporter;

my $db_query = morexgeneatlas_wrapper->new('morexgeneatlas');

my $image_dir = "/var/www/html/morexgeneatlas/blast_net_images";


sub printHeader{

  my ($page_title) = @_;

  print "<html lang='en'><head>
<meta http-equiv='content-type' content='text/html; charset=UTF-8'>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>

  <link href='https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.14.0/css/all.min.css' rel='stylesheet'/>
  <script src='https://cdnjs.cloudflare.com/ajax/libs/jquery/3.3.1/jquery.min.js'></script>

  <script src='https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.min.js'></script>
  <link href='https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css' rel='stylesheet' />

  <script src='//unpkg.com/react\@16/umd/react.development.js' ></script>
  <script src='//unpkg.com/react-dom\@16/umd/react-dom.development.js' ></script>
  <script src='//unpkg.com/\@jbrowse/react-linear-genome-view/dist/react-linear-genome-view.umd.development.js' ></script>

  <link rel='stylesheet' href='stylesheet/eorna.css'>

  <script src='test_files/plotly-latest.js'></script><style id='plotly.js-style-global'></style>
  <script> \$(function () { \$('.example-popover').popover({ container: 'button'}) })</script>
  <script>function clearOptions(id){document.getElementById(id).options.length = 0;}</script>


  <title>Morex Gene Atlas - $page_title</title>
  </head>
  
  <body class='body-pos' data-spy='scroll' data-target='#list-data' data-offset='35' id='content'>\n\n";
}

sub printFooter {

  print "
                </main>
            </div>
        </div>
    </body>
</html>";

}

sub transcriptsAccordion{

  my($gene, $dataset) = @_;

  my %transcripts = %{$db_query->getTranscriptSequencesByGene($gene, $dataset)};

  print "
            <div class='anchor data-header' id='trans_seq'>Transcript Sequences from $gene</div>
            <div class='data-body'>\n";

  print "  
            <div class='accordion' id='transcript_accordion'>\n";

  foreach my $number(sort {$a <=> $b} keys %transcripts){
    
    my($transcript_id, $seq_length, $sequence) = split(/\t/, $transcripts{$number}, 3);

    my $ref = "cds" . $number;

    print "  
            <div class='card'>
              <div class='card-header'>
                <button class='btn btn-link btn-nav-accordion collapsed' type='button' data-toggle='collapse' data-target='#$ref'>
                  <span>$transcript_id ($seq_length bp)</span>
                  <i class='fas fa-chevron-up'></i>
                </button>
              </div>
              <div id='$ref' class='collapse'>
                <div class='card-body'><div style='word-wrap:break-word;width: 100%'>
                  $sequence
                </div>
              </div>
            </div>
          </div>\n";
    
  }
  print "
        </div>
      </div>\n";
}



sub aggregateProteinsAccordion{

  my($gene, $dataset) = @_;

  my %proteins = %{$db_query->getProteinSequencesBySeq($gene, $dataset)};

  print "
            <div class='anchor data-header' id='prot_seq'>Protein Sequences from $gene</div>
            <div class='data-body'>\n";

  print "  
            <div class='accordion' id='prot_accordion'>\n";

  my $count_unique_seq;

  foreach my $prot_seq(sort {$a <=> $b} keys %proteins){

    $count_unique_seq++;

    my $seq_length = length($prot_seq);

    my $ref = "prot" . $count_unique_seq;

    $proteins{$prot_seq} =~ s/\,$//;

    my @order = split(/\,/, $proteins{$prot_seq});

    @order = sort {$a <=> $b} (@order);

    my $transcript_list = join(",", @order);

    print "  
            <div class='card'>
              <div class='card-header'>
                <button class='btn btn-link btn-nav-accordion collapsed' type='button' data-toggle='collapse' data-target='#$ref'>
                  <span>$transcript_list ($seq_length aa)</span>
                  <i class='fas fa-chevron-up'></i>
                </button>
              </div>
              <div id='$ref' class='collapse'>
                <div class='card-body'><div style='word-wrap:break-word;width: 100%'>
                  $prot_seq
                </div>
              </div>
            </div>
          </div>\n";
    
  }
  print "
          </div>
        </div>\n";
}



sub printRelatedGenes {
  my ($gene_name) = @_;

  my %related = %{$db_query->getRelatedGenes($gene_name)};

  my @trackDetails = keys(%related); 

    print "
        <div class='anchor data-header' id='related_genes'>Related Genes in Region of $gene_name</div>
        <div class='data-body'>\n";

  if(scalar(@trackDetails)<1){

    print "<h5><font color='#cc3300'>There are no matching genes for this gene.</font></h5></div>\n";

  } else {

    my %gene_relation;
    my %related_datasets;

    print "
          <p>These comparisons between the various barley transcriptomes were produced by first using <a href='https://pubmed.ncbi.nlm.nih.gov/27008021/'>GMAP</a> to map the transcriptomes against the <a download href='downloads/210316_Morex_V3_pseudomolecules_and_unplaced_scaffolds_ENA.zip'>Morex Reference Genome</a> assembly and generate locations. The GMAP mapping was run with parameters restricting the results to the single best match with 90% sequence identity and 80% transcript coverage. Then <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7222033/'>gffcompare</a> was used to compare the GMAP locations of the transcriptomes to the locations of the reference transcript dataset (HvMxRTD annotation). </p>";

    my $full_details = "
        <div class='accordion'>
              <button class='btn btn-link btn-nav-accordion collapsed' type='button' data-toggle='collapse' data-target='#table'>
                  <span>Show Full Details of All Transcripts in Region of $gene_name <i class='fas fa-chevron-up'></i></span>
              </button>
            
            <div id='table' class='collapse'>
              <div style='word-wrap:break-word;width: 100%'>
                <table class='table-bordered' width='100%'><tbody><tr><th>Related Transcript ID</th><th>Dataset</th><th>Reference Transcript ID</th><th>Gffcompare Class Code</th></tr>\n";

    foreach my $related_transcript(sort keys %related){

      foreach my $transcript_dataset(sort {$b <=> $a}keys %{$related{$related_transcript}}){

        my ($query_gene_id, $reference_tran_id, $reference_gene_id, $reference_dataset, $reference_match_length, $query_tran_exon_number, $query_tran_length, $major_isoform, $gffcompare_class_code) = split(/\t/, $related{$related_transcript}{$transcript_dataset}, 9);

        my $transcript_number = $related_transcript;
        $transcript_number =~ s/$query_gene_id\.//;


        $gene_relation{$query_gene_id}{$transcript_dataset}{$transcript_number} = $gffcompare_class_code;

        my $code_link = &get_code_link($query_gene_id, $transcript_number, $gffcompare_class_code);

        my $tran_db_link = &make_db_link($related_transcript);

        $related_datasets{$query_gene_id} = $transcript_dataset;

        $full_details .= "
                  <tr><td><a href='$tran_db_link'>$related_transcript</a></td><td>$transcript_dataset</td><td>$reference_tran_id</td><td>$code_link</td></tr>";
      }
      
    }

    $full_details .= "
                    </table>
                  </div>
                </div>
              </div>
            </div>\n";

    print "
          <table class='table-bordered' width='100%'>
            <tbody>
              <tr><th>Related gene</th><th>Dataset</th><th>Related transcripts (gffcompare overlap code) - hover to see overlap definition, click to view transcript</th></tr>";

    foreach my $query_transcript(sort keys %gene_relation) {

      my $trans_list;

      foreach my $transcript_dataset(sort {$a <=> $b} keys %{$gene_relation{$query_transcript}}){

        foreach my $transcript_number(sort {$a <=> $b} keys %{$gene_relation{$query_transcript}{$transcript_dataset}}){

          my $code_link = &get_code_link($query_transcript, $transcript_number, $gene_relation{$query_transcript}{$transcript_dataset}{$transcript_number});
        
          $trans_list .="$code_link ,";
        }

        $trans_list =~ s/,$//;

        print "
              <tr><td><a href='search_ids.cgi?seq_name=$query_transcript'>$query_transcript</a></td><td>$related_datasets{$query_transcript}</td><td>$trans_list</td></tr>";

      }

    }

    print "
            </tbody>
          </table>\n";
       

    print $full_details;

  }

}


sub make_db_link {

  my($seq_name) = @_;

  my $db_link;

  if($seq_name =~ /^Hv_Mx_/) {
  
    if($seq_name =~ /\./){
      $db_link = "transcript.cgi?seq_name=$seq_name&dataset=HvMxRTD";
    } else {
      $db_link = "gene.cgi?seq_name=$seq_name&dataset=HvMxRTD";
    }
  }

  if($seq_name =~ /^HORVU1Hr/) {
  
    if($seq_name =~ /\./){
      $db_link = "https://ics.hutton.ac.uk/eorna/transcript.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";
    } else {
      $db_link = "https://ics.hutton.ac.uk/eorna/gene.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";
    }
  }

  if($seq_name =~ /^BART/) {
    
    $seq_name =~ s/p/u/;
  
    if($seq_name =~ /\./){
      $db_link = "https://ics.hutton.ac.uk/eorna/transcript.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";
    } else {
      $db_link = "https://ics.hutton.ac.uk/eorna/gene.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";
    }
  }

  if($seq_name =~ /^A/) {
    $db_link = "https://ics.hutton.ac.uk/morexGenes/view_flcdna_sequence2.cgi?seq_name=$seq_name&dataset=assembly3_WGSMorex_rbca.fasta";
  }

    if($seq_name =~ /^MLOC/) {

    if($seq_name =~ /\./){
      $db_link = "https://ics.hutton.ac.uk/morexGenes/view_cds_sequence2.cgi?seq_name=$seq_name&dataset=assembly3_WGSMorex_rbca.fasta";
    } else {
      $db_link = "https://ics.hutton.ac.uk/morexGenes/view_gene2.cgi?seq_name=$seq_name&dataset=assembly3_WGSMorex_rbca.fasta";
    }
  }

  if($seq_name =~ /^JLOC/) {
    
    if($seq_name =~ /\./){
      $db_link = "https://ics.hutton.ac.uk/barleyGenes/view_cds_sequence2.cgi?seq_name=$seq_name&dataset=Cuffmerge_version4";
    } else {
      $db_link = "https://ics.hutton.ac.uk/barleyGenes/view_gene2.cgi?seq_name=$seq_name&dataset=Cuffmerge_version4";
    }
  }

  if($seq_name =~ /^HORVU\.MOREX\.r3/) {

    if($seq_name =~ /\./){
      $db_link = "http://plants.ensembl.org/Hordeum_vulgare/Search/Results?species=Hordeum_vulgare;idx=;q=$seq_name;site=ensemblthis";
    } else {
      $db_link = "http://plants.ensembl.org/Hordeum_vulgare/Search/Results?species=Hordeum_vulgare;idx=;q=$seq_name;site=ensemblthis";
    }
  }

  return $db_link;
}



sub get_code_link {
  my ($gene_name, $trans_number, $code) = @_;

  my $transcript_name = $gene_name . "." . $trans_number;

  if($gene_name =~ /^A/){

    $transcript_name = $gene_name;
    $trans_number = $gene_name;

  }

  my $code_definition_link;

  if ($code eq '=') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='complete, exact intron chain match'>$trans_number (=)</a>"; }
  if ($code eq 'c') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='contained in reference transcript (intron compatible)'>$trans_number (c)</a>";}
  if ($code eq 'k') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='contains reference transcript (reverse containment)'>$trans_number (k)</a>";}
  if ($code eq 'm') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='retained intron(s) compared to reference, full intron chain match everywhere else'>$trans_number (m)</a>";}
  if ($code eq 'n') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='completely overlaps intron from reference transcript, partial or no intron chain match everywhere else'>$trans_number (n)</a>";}
  if ($code eq 'j') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='multi-exon with at least one junction match'>$trans_number (j)</a>";}
  if ($code eq 'e') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='single exon that partially covers an intron from reference'>$trans_number (e)</a>";}
  if ($code eq 'o') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='other same strand overlap with reference exons'>$trans_number (o)</a>";}
  if ($code eq 's') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='intron match on the opposite strand (likely a mapping error)'>$trans_number (s)</a>";}
  if ($code eq 'x') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='exonic overlap on the opposite strand'>$trans_number (x)</a>";}
  if ($code eq 'i') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='fully contained within a reference intron'>$trans_number (i)</a>";}
  if ($code eq 'y') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='contains a reference within its intron(s)'>$trans_number (y)</a>";}
  if ($code eq 'p') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='possible polymerase run-on (close to reference but no overlap)'>$trans_number (p)</a>";}
  if ($code eq 'r') { $code_definition_link = "<a href='search_ids.cgi?seq_name=$transcript_name' data-bs-toggle='tooltip' title='repeat (at least 50% bases are soft-masked)'>$trans_number (r)</a>";}

   return $code_definition_link;
}



sub createJSPlotlyGraph {
  my($gene_id, $study_list) = @_;

  # Get a list of the transcript_ids for this gene
  my @transcripts = @{$db_query->getTranscriptListByGene($gene_id)};

  my $number_of_transcripts = @transcripts;

  print "<div class='anchor data-header' id='tpmgraph'>TPM Values of Transcripts</div><div class='data-body'>\n";

  #convert the unpadded gene id to a padded gene_id to access the tpm_values table
  $gene_id =~ s/-u/-p/;

  # Make a hash from the study list to be able to skip studies that are not in the list
  my %tpms;
  my %study_filter;
  if($study_list eq ""){
    %tpms = %{$db_query->getPaddedTPMsJS_new($gene_id)};
  } else {

    my $study_filter = "";
    my @study_list = split(/,/, $study_list);
    if (@study_list > 0){
      foreach my $item_id(@study_list){
        $study_filter .= "(s.experiment_id = '$item_id') || ";
      }
      $study_filter =~ s/ \|\| $//;
      $study_filter = "and ($study_filter)";

      %tpms = %{$db_query->getPaddedTPMsJSSubset($gene_id,$study_filter)};

    }
  }

  my @trackDetails = keys(%tpms);

  if(scalar(@trackDetails)<1){
    print "<span class='error'>No TPM values available for this gene's transcript models. Try any <a href='#matchingGenes'>matching BART genes below </a>.</span></div>\n";
  } else {

    # Make filename for the plotly page
    my $plotly_filename = $gene_id . "_plotly.html";
  
    print "<a href='plotly/$plotly_filename'>Link to larger plotly expression graph</a>\n";

    # Gather all the plotly javascript into a string to dump out into a separate html document so we can have a link out to it
    my $plotly_text_doc = "
  <div id='plotly_tpm' style='width:100%; height:400px'> 
  </div>
  <script>
  TESTER = document.getElementById('plotly_tpm');\n";

  
    my %tpms_pivot; #going to want to pivot the data matrix to be able to print out as tab-delimited text later
  
    # Make a list of the names of the graph data lists that will be drawn on the plot
    my $trace_list;
    my $traceLabel_text;
    my $report_expt_count;
    my %expt_details;
    my $graph_width;
  

    foreach my $transcript_number(sort {$a <=> $b} keys %tpms){

      my $tracevar = "trace" . $transcript_number;
    
      $trace_list .= "$tracevar,";

      #make the lists of data for the traces
      my $expt_name_list;
      my $tpm_value_list;
      my $label_value_list;
      my $graph_label_text;

      my $expt_count;
      foreach my $expt_id(sort {$a <=> $b} keys %{$tpms{$transcript_number}}){

        my ($tpm_value, $study_accession, $sample_accession, $cultivar, $tissue, $bio_replicate, $expt_condition, $dev_stage) = split(/\t/, $tpms{$transcript_number}{$expt_id}, 8);

        $expt_count++;

        $label_value_list .= "\"0\",";
        $expt_name_list .= "\"$expt_count\","; # Note: this is the column "number" for plotly, not the expt_id in mysql
        $tpm_value_list .= "$tpm_value,";

        $graph_label_text .="\"$expt_id<br>SRA Study: $study_accession<br>Sample ID: $sample_accession<br>$cultivar<br>Tissue: $tissue<br>Dev Stage: $dev_stage<br>Replicate: $bio_replicate &nbsp; &nbsp; Condition: $expt_condition\",\n";

        # Need to pivot the data matrix for printing out in a tab-delimited text file later
        $tpms_pivot{$expt_id}{$transcript_number} = $tpm_value;

        $expt_details{$expt_id} = "$study_accession\t$sample_accession\t$cultivar\t$tissue\t$dev_stage\t$bio_replicate\t$expt_condition";

     }


      if ($transcript_number == 1){
        $report_expt_count = $expt_count;
      }

      # Base the width of the graph on the number of samples shown?
      $graph_width = 100 + ($report_expt_count * 5);
    
      $expt_name_list =~ s/\,$//;
      $tpm_value_list =~ s/\,$//;
      $label_value_list =~ s/\,$//;
      $graph_label_text =~ s/\,$//;

      $traceLabel_text = "var traceLabels = {
    type: 'bar',
    name: '',
    hoverinfo: 'text',
    hoverlabel:{ bgcolor: 'rgb(124,124,124)', font: { color: 'rgb(255,255,255)'}},
    showlegend: false,
    text: [$graph_label_text],
    x: [$expt_name_list],
    y: [$label_value_list]
    };\n";


    $plotly_text_doc .= "var $tracevar = {
    type: 'bar',
    name: '$transcript_number',
    hoverinfo: 'none',
      x: [$expt_name_list],
      y: [$tpm_value_list]
    };\n";

    
    }

    $plotly_text_doc .= $traceLabel_text;
    
    $trace_list =~ s/\,$//;

    $plotly_text_doc .= "
  var data = [traceLabels,$trace_list];
  var layout = {
    legend: {orientation: 'h', side: ''},
    autosize: false,
    width: $graph_width,
    height: 420,
    margin: {
      l: 50,
      r: 50,
      b: 10,
      t: 30,
      pad: 0
    },
    yaxis: {title: 'TPM'},
    xaxis: {
      title: 'Sample',
      showticklabels: false,
    },
    barmode: 'stack',
   colorway: ['#17becf', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
  };

  var config = {
    responsive: true,
    toImageButtonOptions: {scale: 4 }
  }

  Plotly.newPlot('plotly_tpm', data, layout, config);


</script>
\n";

print "<iframe src='https://ics.hutton.ac.uk/morexgeneatlas/plotly/$plotly_filename' height=460 width=100% style='border:none;'>Plotting graph, please wait......</iframe> ";

  #print "There are $report_expt_count samples shown in the graph ($graph_width px)<br>\n";
  
    # Open the stand-alone html file and print the javascript for the plot in it
    open(PLOTLY, ">/var/www/html/morexgeneatlas/plotly/$plotly_filename") || die ("Cannot create file /var/www/html/morexgeneatlas/plotly/$plotly_filename");

    print PLOTLY "<!DOCTYPE html>
  <html lang='en'>
  <head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>

  <link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css'>
  <link rel='stylesheet' href='/morexgeneatlas/stylesheet/eorna.css'>

  <title>Morex Gene Atlas</title>

  <script src='https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js'></script>
  <script src='https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js'></script>
  <script src='https://cdn.plot.ly/plotly-latest.min.js'></script>

  <style>
    .js-plotly-plot .plotly .modebar {
      left: 0%;
      
    }

  </style>
  </head>
  <body>
  <font size=small><b>$gene_id TPM values</b></font>
  $plotly_text_doc\n
  </body>
</html>\n";
  
    close PLOTLY;

# Make the text file download a separate call........?

    # Make a tab-delimited text file of the TPM values for downloading

    my $tpms_filename = $gene_id . "_tpms.txt";


    #print "Creating file '$tpms_filename': ";

    open(TPMS, ">/var/www/html/morexgeneatlas/plotly/$tpms_filename") || die ("Cannot create file /var/www/html/morexgeneatlas/plotly/$tpms_filename");

    # Print headers for the text file
    print TPMS "Experiment ID\tSRA Study\tSRA Run\tCultivar\tTissue\tDev Stage\tBio Replicate\tExpt Condition";
    for (my $i = 1; $i <= $number_of_transcripts; $i++){
      print TPMS"\t$i";
    }
    
    print TPMS "\n";
  
    foreach my $expt_id(sort {$a <=> $b} keys %tpms_pivot){

      print TPMS "$expt_details{$expt_id}";

      foreach my $transcript_number(sort {$a <=> $b} keys %{$tpms_pivot{$expt_id}}) {

        print TPMS "\t$tpms_pivot{$expt_id}{$transcript_number}";

      }
      
      print TPMS "\n";
    }
    close TPMS;

    print "<br><a download href='plotly/$tpms_filename'>Download tab-delimited text file of TPM values</a><br><br>\n";
  
  }
}


sub printGOAnnotation{

  my ($gene_id) = @_;

  my %go_annotation = %{$db_query->getGOAnnotation($gene_id, "HvMxRTD")};

  print "
      <div class='anchor data-header' id='go_annot'>GO Annotation for $gene_id</div>
        <div class='data-body'>";

  my @trackDetails = keys(%go_annotation);

  if(scalar(@trackDetails)<1){

    print "     There is no GO annotation for this gene<br><br></div>\n";
    
  } else {

    print "   
          <h6>Annotation by PANNZER</h6><br>\n";

    print "
          <table class='table-bordered' width='100%'>
          <tbody><tr><th>Transcript ID</th><th>PANNZER annotation</th><th>GO ID</th><th>GO terms</th></tr>";

    foreach my $transcript_number(sort {$a <=> $b} keys %go_annotation){

      my($pannzer_annotation, $go_ids, $go_terms) = split(/\t/, $go_annotation{$transcript_number}, 3);

      my $transcript_id = $gene_id . "." . $transcript_number;

      my @go_ids = split(/\;/, $go_ids);

      my $go_ids_text;
      for (my $i = 0; $i < @go_ids; $i++){

        $go_ids_text .= "<a href='http://amigo.geneontology.org/amigo/term/$go_ids[$i]' data-bs-toggle='tooltip' title='View in AmiGO'>$go_ids[$i]</a> ";

      }

      print "
            <tr><td>$transcript_id</td><td>$pannzer_annotation</td><td>$go_ids_text</td><td>$go_terms</td></tr>";

    }
    print "
        </tbody>
      </table>
    </div>\n";

  }
  
}


sub getHomologies{

  my($gene, $dataset) = @_;

  # Find Blast matches for the longest transcript of this gene

  my ($longest_transcript, $length) = $db_query->{'dbh'}->
         selectrow_array("select transcript_id, seq_length from transcript_sequences where gene_id = '$gene' and dataset_name = '$dataset' order by seq_length desc limit 1");

  my %homologies = %{$db_query->getHomologies($longest_transcript, $dataset)};

  print "
        <div class='anchor data-header' id='homology'>Homologies with Model Species for $gene</div>
          <div class='data-body'>";

  my @trackDetails = keys(%homologies);

  if(scalar(@trackDetails)<1){

    print "     There are no homologies for this gene<br><br></div>\n";
    
  } else {

    print "
          <h6>Top BLASTX hit (E-value cutoff 1e-10) for longest transcript $longest_transcript ($length bp)</h6><br>\n";

    print "
          <table  class='table-bordered' width='100%'><tbody><tr><th>Blast Database</th><th>Match Name</th><th>E-value</th><th>Bit Score</th><th>Percentage Identity</th><th>Description</th></tr>";

    foreach my $blast_db(sort {$a <=> $b} keys %homologies){

      foreach my $match_name(sort {$a <=> $b} keys %{$homologies{$blast_db}}){

        my($score, $percent_id, $evalue, $description) = split(/\t/, $homologies{$blast_db}{$match_name}, 6);

        my $link;
        if($blast_db eq "TAIRPP10"){
          $link = "<a href='http://www.araport.org/locus/$match_name' target='_blank'>$match_name</a>";
        } elsif ($blast_db eq "RICEPP7"){
          $link = "<a href='http://rice.uga.edu/cgi-bin/gbrowse/rice/?name=$match_name' target='_blank'>$match_name</a>";
        } else {
          $link = $match_name;
        }
  

        print "
            <tr><td>$blast_db</td><td>$link</td><td>$evalue</td><td>$score<td>$percent_id</td><td>$description</td></tr>";

      }

    }
    print "
          </tbody>
        </table>
      </div>\n";
  }

}



sub get_JBrowse{
  
  my($seq_name, $dataset) = @_;

  my ($gene_start, $gene_stop, $chr_id) = $db_query->{'dbh'}->selectrow_array("select min(f_start), max(f_end), chr_id from transcript_structures where gene_id = '$seq_name' and dataset_name=\'$dataset\'  group by chr_id");
  
  print "<div class='anchor data-header' id='jbrowse'>MorexGS JBrowse</div><div class='data_body'><br><br>\nLocation on Morex genome - $chr_id:$gene_start..$gene_stop<br><br>\n";

  print "    <div id='jbrowse_linear_genome_view'></div>
    <script type='module'>
      import assembly from './morexgs_jbrowse/assembly.js'
      import tracks from './morexgs_jbrowse/tracks.js'
      const { createViewState, JBrowseLinearGenomeView } = JBrowseReactLinearGenomeView
      const { createElement } = React
      const { render } = ReactDOM

      const state = new createViewState({
        assembly,
        tracks,
        location: '$chr_id:$gene_start..$gene_stop',
        

        defaultSession: {
          name: 'my session',
          view: {
            id: 'linearGenomeView',
            type: 'LinearGenomeView',
            tracks: [
              {
                'type': 'FeatureTrack',
                'configuration': 'iso_rna_hvmx_merged',
                'displays': [
                  {
                    'type': 'LinearBasicDisplay',
                    'height': 300,
                    'trackShowDescriptions': true,
                    'trackShowLabels': true,
                    'trackDisplayMode': 'normal',
                    'renderer': {
                      'type': 'SvgFeatureRenderer'
                    }
                  }
                ]
              },
            ],
            'hideHeader': false,
            'hideHeaderOverview': false,
            'hideNoTracksActive': false,
            'trackSelectorType': 'hierarchical',
            'trackLabels': 'offset',
            'showCenterLine': false,
            'showCytobandsSetting': true,
            'showGridlines': true
          },
        },
      })
      function navTo(event) {
        state.session.view.navToLocString(event.target.dataset.location)
      }
      const buttons = document.getElementsByTagName('button')
      for (const button of buttons) {
        if (button.dataset.type === 'gene_button') {
          button.addEventListener('click', navTo)
        }
      }
      render(
        createElement(JBrowseLinearGenomeView, { viewState: state }),
        document.getElementById('jbrowse_linear_genome_view'),
      )
    </script>";

  print "</div><br><br>\n";
  
}


