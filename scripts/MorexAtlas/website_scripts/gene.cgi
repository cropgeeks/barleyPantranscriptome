#!/usr/bin/perl

use strict;
use CGI;
use GD;
use lib '/var/www/html/morexgeneatlas';
use morexgeneatlas_wrapper;
use morexgeneatlas;

my $cgi_query  = CGI->new();
my $gene_name  = $cgi_query->param("seq_name");
my $dataset    = $cgi_query->param("dataset");
my $study_list = $cgi_query->param("studies");

my $image_dir = "/var/www/html/morexgeneatlas/blast_net_images";
my $db_query = morexgeneatlas_wrapper->new('morexgeneatlas');

print "Content-type: text/html\n\n";

#Get the contig that this gene is from
my ($chr_id) = $db_query->{'dbh'}->
  selectrow_array("select chr_id from transcript_sequences where gene_id = '$gene_name' and dataset_name='$dataset'");
  
  
if ($chr_id eq ""){

  morexgeneatlas::printHeader("Error");

    print "<div class='container-fluid'>
            <div class='row'>


                <nav class='col-md-2 d-none d-md-block grey sidebar'>
                    <div class='sidebar-sticky'>
                        <ul class='nav flex-column'>
                            <li class='nav-item nav-link'  id='nav-title'><img class='img-fluid' src='images/morexgeneatlas-logo.png' width=150'></li>
                            <li class='nav-item nav-link'  id='nav-header'><b>Utilities Menu</b></li>
                            <div class='list-group'>
                            <li class='nav-item'><a class='nav-link active' href='index.html'>Home</a></li>
                            <li class='nav-item'><a class='nav-link' href='blast.html'>Homology Search</a></li>
                            <li class='nav-item'><a class='nav-link' href='keyword.html'>Annotation Search</a></li>
                            <li class='nav-item'><a class='nav-link' href='download.html'>Bulk Data Download</a></li>
                            <li class='nav-item'><a class='nav-link' href='relator.html'>Transcriptome Relationships</a></li>
                            </div>
                        </ul>
                    </div>
                </nav>

                <main role='main' class='col-md-9 ml-sm-auto col-lg-10 pt-3 px-4'>

                    <div class='data-header'>Error Message</div>
                    <br>
                    <div class='data-body'>\n
                      <span class='error'>There are no genes in the database with ID \"$gene_name\". <br>Please enter a valid sequence identifier for the search.</span><br><br>
                            <div class='card-deck  m-3'>
                            <div class='card border-dark'>
                                <div class='card-header option-header'>Search for a Sequence by ID</div>
                                <div class='card-body'>

                                <p>Currently <font color='#5c5f58'><b>Morex Gene Atlas</b></font> contains the HvMxRTD transcripts (<a href=''>publication here</a>).</p>

                                <p>The predicted genes have nomenclature like \"Hv_Mx_chr1HG02979\". Transcript models are named after the the gene loci with a number appended, like \"Hv_Mx_chr1HG02979.10\". 

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
                    </div>
                </main>
            </div>
        </div>
    </body>
</html>";

  
  die("morexgeneatlas/gene.cgi: Lookup contig source for $gene_name, $dataset - No contig in morexgene atlas for this gene id\n");
}

my ($gene_start, $gene_end, $chromosome) = $db_query->{'dbh'}->
  selectrow_array("select min(f_start), max(f_end), chr_id from transcript_structures where gene_id = '$gene_name' and dataset_name=\'$dataset\'  group by chr_id");


morexgeneatlas::printHeader($gene_name);

    print "
    <div class='container-fluid'>
      <div class='row'>
        <nav class='col-md-2 d-none d-md-block grey sidebar'>
          <div class='sidebar-sticky'>
            <ul class='nav flex-column'>

                <li class='nav-item nav-link'  id='nav-title'><img class='img-fluid' src='images/morexgeneatlas-logo.png' width=150'></li>
                <li class='nav-item'></li>
                <li class='nav-item nav-link' id='nav-header'><b>Gene IDs</b></li>
                <li class='nav-item' id='nav-info'>$gene_name</li>
                <li class='nav-item nav-link' id='nav-header'><b>Location</b></li>
                <li class='nav-item' id='nav-info'>$chromosome: $gene_start - $gene_end</li>
                <li class='nav-item nav-link' id='nav-header'><b>Data</b></li>
                <div class='list-group' id='list-data'>
                  <li class='nav-item'><a class='nav-link' href='#tpmgraph'>TPM Values</a></li>
                  <li class='nav-item'><a class='nav-link' href='#jbrowse'>Morex JBrowse</a></li>
                  <li class='nav-item'><a class='nav-link' href='#trans_seq'>Transcripts</a></li>
                  <li class='nav-item'><a class='nav-link' href='#prot_seq'>Proteins</a></li>
                  <li class='nav-item'><a class='nav-link' href='#related_genes'>Related Genes</a></li>
                  <li class='nav-item'><a class='nav-link' href='#homology'>Homology</a></li>
                  <li class='nav-item'><a class='nav-link' href='#go_annot'>GO Annotation</a></li>
                </div>
                <li class='nav-item nav-link'  id='nav-header'><b>Utilities Menu</b></li>
                <div class='list-group'>
                  <li class='nav-item'><a class='nav-link' href='index.html'>Home</a></li>
                  <li class='nav-item'><a class='nav-link' href='blast.html'>Homology Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='keyword.html'>Annotation Search</a></li>
                  <li class='nav-item'><a class='nav-link' href='download.html'>Bulk Data Download</a></li>
                  <li class='nav-item'><a class='nav-link' href='relator.html'>Transcriptome Relationships</a></li>
                </div>
            </ul>
          </div>
        </nav>
        <main role='main' class='col-md-9 ml-sm-auto col-lg-10 pt-3 px-4'>\n\n";


############### Print table of padded TPM values for the transcript of this gene ############

morexgeneatlas::createJSPlotlyGraph($gene_name, $study_list);

############### Get the JBrowse location of the gene ##########################################

morexgeneatlas::get_JBrowse($gene_name, $dataset);

morexgeneatlas::transcriptsAccordion($gene_name, $dataset);

morexgeneatlas::aggregateProteinsAccordion($gene_name, $dataset);


################## Print related genes  ####################################################


morexgeneatlas::printRelatedGenes($gene_name);


############### Get the Rice and TAIR top hits for the longest transcript of the gene #######

morexgeneatlas::getHomologies($gene_name, $dataset);

############### Print the GO annotation for a gene #####################################

morexgeneatlas::printGOAnnotation($gene_name);

morexgeneatlas::printFooter();



