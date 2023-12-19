#!/usr/bin/perl


use strict;
use CGI;
use lib '/var/www/html/morexgeneatlas';
use morexgeneatlas;
use morexgeneatlas_wrapper;


my $cgi_query = CGI->new();
my $seq_name  = $cgi_query->param("seq_name");


print "Content-type: text/html\n\n";


my $query = morexgeneatlas_wrapper->new('morexgeneatlas');


if (($seq_name =~ /^Hv_Mx_/) || ($seq_name =~ /^HORVU1Hr/) || ($seq_name =~ /^BaRT/) || ($seq_name =~ /BART/) || ($seq_name =~ /MLOC/) || ($seq_name =~ /^A/) || ($seq_name =~ /JLOC/) || ($seq_name =~ /HORVU\.MOREX\.r3/)){
  
  my $redirect;

  if($seq_name =~ /^Hv_Mx_/) {

    $seq_name =~ s/\.\d+$//;
 
    $redirect = "gene.cgi?seq_name=$seq_name&dataset=HvMxRTD";


  }

  if($seq_name =~ /^HORVU1Hr/) {
  
    if($seq_name =~ /\./){

      $redirect = "https://ics.hutton.ac.uk/eorna/transcript.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";
  
    } else {
 
      $redirect = "https://ics.hutton.ac.uk/eorna/gene.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";

    }
  }

  if($seq_name =~ /^BART/) {
    
    $seq_name =~ s/p/u/;
  
    if($seq_name =~ /\./){

      $redirect = "https://ics.hutton.ac.uk/eorna/transcript.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";
  
    } else {
 
      $redirect = "https://ics.hutton.ac.uk/eorna/gene.cgi?seq_name=$seq_name&dataset=150831_barley_pseudomolecules";

    }
  }

  if($seq_name =~ /^A/) {

    $redirect = "https://ics.hutton.ac.uk/morexGenes/view_flcdna_sequence2.cgi?seq_name=$seq_name&dataset=assembly3_WGSMorex_rbca.fasta";
  

  }

    if($seq_name =~ /^MLOC/) {
    
  
    if($seq_name =~ /\./){

      $redirect = "https://ics.hutton.ac.uk/morexGenes/view_cds_sequence2.cgi?seq_name=$seq_name&dataset=assembly3_WGSMorex_rbca.fasta";
  
    } else {
 
      $redirect = "https://ics.hutton.ac.uk/morexGenes/view_gene2.cgi?seq_name=$seq_name&dataset=assembly3_WGSMorex_rbca.fasta";

    }
  }

  if($seq_name =~ /^JLOC/) {
    
  
    if($seq_name =~ /\./){

      $redirect = "https://ics.hutton.ac.uk/barleyGenes/view_cds_sequence2.cgi?seq_name=$seq_name&dataset=Cuffmerge_version4";
  
    } else {
 
      $redirect = "https://ics.hutton.ac.uk/barleyGenes/view_gene2.cgi?seq_name=$seq_name&dataset=Cuffmerge_version4";

    }
  }

  if($seq_name =~ /^HORVU\.MOREX\.r3/) {
    
  
    if($seq_name =~ /\./){

      $redirect = "http://plants.ensembl.org/Hordeum_vulgare/Search/Results?species=Hordeum_vulgare;idx=;q=$seq_name;site=ensemblthis";

  
    } else {
 
      $redirect = "http://plants.ensembl.org/Hordeum_vulgare/Search/Results?species=Hordeum_vulgare;idx=;q=$seq_name;site=ensemblthis";

    }
  }


  print "<html>
<head>
<meta http-equiv=\"Refresh\" content=\"0;url=$redirect\" />
</head>

<body>

</body>
</html>";

} else {
  
  morexgeneatlas::printHeader("Error");


print "<div class='container-fluid'>
      <div class='row'>
        <nav class='col-md-2 d-none d-md-block grey sidebar'>
          <div class='sidebar-sticky'>
            <ul class='nav flex-column'>

                <li class='nav-item nav-link'  id='nav-title'><img class='img-fluid' src='images/morexgeneatlas-logo.png' width=150'></li>

                <li class='nav-item nav-link'  id='nav-header'><b>Utilities Menu</b></li>
                <div class='list-group'>
                  <li class='nav-item'><a class='nav-link' href='index.html'>Home</a></li>
                  <li class='nav-item'><a class='nav-link' href='blast.html'>Homology Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='keyword.html'>Annotation Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='download.html'>Bulk Data Download</a></li>
                  <li class='nav-item'><a class='nav-link' href='relator.html'>Transcriptome Relationships</a></li>
                </div>
                <li class='nav-item nav-link'  id='nav-title'></li>

            </ul>
          </div>
        </nav>
        <main role='main' class='col-md-9 ml-sm-auto col-lg-10 pt-3 px-4'>\n\n";
  
  print "<div class='data-header'>Error Message</div>
                    <br>
                    <div class='data-body'>
                        <span class='error'>The search term \"$seq_name\" does not match any sequence ID in the database. <br>Please enter a valid sequence identifier for the search.</span>

                    <div class='card-deck  m-3'>
                            <div class='card border-dark'>
                                <div class='card-header option-header'>Search for a Sequence by ID</div>
                                <div class='card-body'>

                                <p>Currently <font color='#5c5f58'><b>Morex Gene Atlas</b></font> contains the HvMxRTD transcripts (<a href=''>publication here</a>).</p>

                                <p>The predicted genes have nomenclature like 'Hv_Mx_chr1HG02979'. Transcript models are named after the the gene loci with a number appended, like 'Hv_Mx_chr1HG02979.10'. 

                                <p>Searching for a gene ID will return all the transcript models from that gene region. You can go directly to a gene or transcript by using this search box:</p>

                                <form action='search_ids.cgi' method='get'>
                                    <div class='input-group mb-3'>
                                        <div class='input-group-prepend'>
                                            <span class='input-group-text'>Search for Sequence ID</span>
                                        </div>
                                        <input type='text' class='form-control' name='seq_name' value='Hv_Mx_chr1HG02979' onfocus=\"this.value=''\">
                                        <div class='input-group-append'><input class='btn btn-secondary red-btn' type='submit' value='Search'></div>
                                    </div>
                                    
                                    </form>
                                </div>
                            </div>
                        </div>
                    </div>\n";
  
  morexgeneatlas::printFooter();
    
}









