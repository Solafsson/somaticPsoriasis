#!/usr/local/bin/perl
# Function to pull out genomic context for $ARGV[0] bases either side of genomic positions supplied in file $ARGV[1]
# $ARGV[1] file is a tab-delimited text file where each line gives a chromosome and base position
# $ARGV[2] is the name of the output file with <info from $ARGV[0]> \t base \t context (as a string)

use strict;
use warnings;
use lib "/nfs/cangen/lib";
use lib qw(/software/CGP/external-perl/lib/site_perl/5.8.8 
           /software/CGP/external-perl/bioperl
	   /software/pubseq/PerlModules/Ensembl/www_70_1/ensembl/modules/
   	   /software/pubseq/PerlModules/Ensembl/www_70_1/modules
	   /software/pubseq/PerlModules/Ensembl/www_70_1/ensembl-variation/modules  );
				     
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Mapper::Coordinate;
use Bio::Seq;
Bio::EnsEMBL::Registry->load_registry_from_db(
                -host => 'ensembldb.ensembl.org',
		-user => 'anonymous',
		-port => 3306
		);
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor("human", "core", "gene");
my $tr_adaptor    = Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core', 'Transcript' );


open(INPUTLIST, $ARGV[1]);
open(POSLIST, ">", $ARGV[2]);
my $old_fh = select(POSLIST);
$| = 1;
select($old_fh);

my $header = <INPUTLIST>;
chomp $header;
print POSLIST "$header\tWT\tCONTEXT\tGENE\tSTRAND\n";

while (my $posID = <INPUTLIST>) {
	chomp $posID; 
	my @gen_pos = split /\t/, $posID;
	if ($gen_pos[0] eq 23) {$gen_pos[0] = "X"}
	if ($gen_pos[0] eq 24) {$gen_pos[0] = "Y"}
	if ($gen_pos[0] eq 25) {$gen_pos[0] = "MT"}
	my $tr_string = "No_gene\tNo_gene\n";
	my $slice = eval {$slice_adaptor->fetch_by_region( 'chromosome', $gen_pos[0], $gen_pos[1] - $ARGV[0], $gen_pos[1] + $ARGV[0] );};
	if (defined ($slice)) {
		my $seq = $slice->seq();
		my $home_base = substr($seq, $ARGV[0], 1);	
		print POSLIST "$posID\t$home_base\t$seq\t";
		my $transcripts = eval {$tr_adaptor->fetch_all_by_Slice($slice);};
		if (defined ($transcripts)) {
		  my $tr_ct = 0;
		  while ( my $tr = shift @{$transcripts} ) {
		    if ($tr_ct == 0) {
    			$tr_ct++;
			my $dbID      = $tr->dbID();
    			my $start     = $tr->start();
    			my $end       = $tr->end();
    			my $strand    = $tr->strand();
    			my $stable_id = $tr->stable_id();
    			$tr_string = "$stable_id\t$strand\n";
  		    }
		  }

		}	
	
			
	  	print POSLIST "$tr_string";
	}
}

close POSLIST;
close INPUTLIST;
