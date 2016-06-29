#!/usr/bin/perl

#use strict;
#use warnings;
use diagnostics;
use POSIX qw(strftime);

$inputdir = $ENV{'INPUTDIR'};
$obo1file = $ENV{'OBO1FILE'};
$obo2file = $ENV{'OBO2FILE'};

# Program: PRO_ancestors.pl
# Original Author: Mary Dolan
#
# Purpose: Given protein ontology in OBO format
#    generate mapping from PRO ids to MGI ids
#
# Usage: perl map_MGI2PRO.pl
#

#open (PRO_ANCS, ">MGI2PRO.txt");
open (VOCFILE, "> $inputdir/provoc.txt");
open (ANNOTFILE, "> $inputdir/proannot.txt");

%mgi_id = ();
%generic_id = ();
%line_printed = ();

my $date = strftime "%m/%d/%Y", localtime;
#print $date;

######### read in PRO terms
	&parse_obo(%nodes);

######### nodes (PRO ids) to ancestors
	foreach $query_node (keys %nodes) {
	@ancestor_nodes = ();
	&find_ancestors($query_node);
	}
		
######### map MGI ids to PRO ancestors
#print header
        #print PRO_ANCS "mgi_id\tpro_id\tpro_short_name\tpro_name\tpro_synonyms\n";

	foreach $query_node (sort keys %nodes) {
		
		@ancestor_nodes = @ { $nodes{$query_node}->{"ancestors"} };
		foreach $ancestor (@ancestor_nodes) {
			if (exists $mgi_id{$ancestor}) {
				$node_name = $nodes{$query_node}->{"name"};
				if ($node_name =~ /(mouse)/) {
					$short_name = $nodes{$query_node}->{"short_name"};
					if ($short_name eq "") {
						$short_name = $node_name;
						}
					@synonyms = @ { $nodes{$query_node}->{"synonyms"} };
					$synonym_list = join ("|", @synonyms);
					$print_line = "$mgi_id{$ancestor}\t$query_node\t$short_name\t$node_name\t$synonym_list";
					$vocfile_line = "${short_name}\t$query_node\tcurrent\t\t$node_name\t\t$synonym_list\t";
					$annot_line = "$query_node\t$mgi_id{$ancestor}\tJ:232325\tIEA\t\t\tproload\t$date\t\t";

					if (! exists $line_printed{$print_line}) {
						#print PRO_ANCS "$print_line\n";
						print VOCFILE "$vocfile_line\n";
						print ANNOTFILE "$annot_line\n";
						$line_printed{$print_line}++;
						}
					}
				}
			elsif ((exists $generic_id{$ancestor}) and (exists $mgi_id{$generic_id{$ancestor}})){
				$node_name = $nodes{$query_node}->{"name"};
				if ($node_name =~ /(mouse)/) {
					$short_name = $nodes{$query_node}->{"short_name"};
					if ($short_name eq "") {
						$short_name = $node_name;
						}
					@synonyms = @ { $nodes{$query_node}->{"synonyms"} };
					$synonym_list = join ("|", @synonyms);
					$print_line = "$mgi_id{$generic_id{$ancestor}}\t$query_node\t$short_name\t$node_name\t$synonym_list";
					$vocfile_line = "${short_name}\t$query_node\tcurrent\t\t$node_name\t\t$synonym_list\t";
					$annot_line = "$query_node\t$mgi_id{$generic_id{$ancestor}}\tJ:232325\tIEA\t\t\tproload\t$date\t\t";

					if (! exists $line_printed{$print_line}) {
						#print PRO_ANCS "$print_line\n";
						print VOCFILE "$vocfile_line\n";
						print ANNOTFILE "$annot_line\n";
						$line_printed{$print_line}++;
						}
					}
				}
			
			}
		}

#close (PRO_ANCS) || die "PRO_ANCS file was not closed\n";	
close (VOCFILE) || die "VOCFILE file was not closed\n";	
close (ANNOTFILE) || die "ANNOTFILE file was not closed\n";	
exit (0);

############
sub parse_obo()
{
    %nodes = ();

    open (OBO, $obo1file) ||print "ontology open error";
    while (<OBO>) {

    if (/\[Term\]/) {
    #new term
	&get_node_info($node_id, $node_name, $short_name, @parent_nodes, $is_obsolete, $nonmouse);
	if ((! $is_obsolete) and (! $nonmouse)) {
	$term = {
		"name" => $node_name,
		"short_name" => $short_name,
		"parents" => [ @parent_nodes ],
		"ancestors" => [  ],
		};
	
	if ($node_id =~ /PR:/) {
		$nodes{$node_id} = $term;
		}
	}
	} #end term test
     } #end read loop

     close (OBO);

    open (OBO, $obo2file) ||print "ontology open error";
    while (<OBO>) {

    if (/\[Term\]/) {
    #new term
	&get_node_info($node_id, $node_name, $short_name, @synonyms, @parent_nodes, $is_obsolete, $nonmouse);
	if ((! $is_obsolete) and (! $nonmouse)) {
	$term = {
		"name" => $node_name,
		"short_name" => $short_name,
		"synonyms" => [ @synonyms ],
		"parents" => [ @parent_nodes ],
		"ancestors" => [  ],
		};
	
	if ($node_id =~ /PR:/) {
		$nodes{$node_id} = $term;
		}
	}
	} #end term test
     } #end read loop

     close (OBO);
}

############
sub get_node_info()
{
#clear info
    $node_id = "";
    $node_name = "";
    $short_name = "";
    @synonyms = ();
    @parent_nodes = ();
    $organism_gene = 0;
    $is_obsolete = 0;
    $nonmouse = 0;
    $taxon_id = "";

    while (<OBO>) {
    chomp;
    $term_line = $_;

    #node_id
	if ($term_line =~ /^id: (PR:\S+)/) {
	    $node_id = $1;
	    #force uniprot isoform with suffix '-d' to be child of uniprot without suffix
	    if ($node_id =~ /(PR:\S+)-\d+/) {
	    	$parent_id = $1;
	    	$parent_type = "is_a";
	    	$parent = {
			      "id" => $parent_id,
			      "type" => $parent_type,
			      };
	    	push @parent_nodes, $parent;
	    	}
	    }

    #node name
	if ($term_line =~ /^name: (.*)/) {
	    $node_name = $1;
	    }

    #comment: Category=organism-gene.
	if (($term_line =~ /^comment: Category=organism-gene./) and ($node_name =~ /(mouse)/)) {
	    $organism_gene = 1;
	    }
    
    #is_a
    #is_a: PR:000000001 ! protein
	if ($term_line =~ /^is_a: (PR:\S+)/) {
	    $parent_id = $1;
	    $parent_type = "is_a";
	    $parent = {
		      "id" => $parent_id,
		      "type" => $parent_type,
		      };
	#skip "PR:000000001" ! protein and "PR:000029032" ! Mus musculus protein and "PR:000018263" ! amino acid chain
	#PR:000018264 ! proteolytic cleavage product <=== to add?	
	if (($parent_id ne "PR:000029032") and ($parent_id ne "PR:000000001")  and ($parent_id ne "PR:000018263")) {
		push @parent_nodes, $parent;
		}
	    }

    #intersection_of
    if ($term_line =~ /^intersection_of: (PR:\S+)/) {
	$parent_id = $1;
	$parent_type = "intersection_of";
	$parent = {
	      "id" => $parent_id,
	      "type" => $parent_type,
	      };
	#skip "PR:000000001" ! protein and "PR:000029032" ! Mus musculus protein and "PR:000018263" ! amino acid chain
	if (($parent_id ne "PR:000029032") and ($parent_id ne "PR:000000001")  and ($parent_id ne "PR:000018263")) {
		push @parent_nodes, $parent;
		if ($organism_gene) {
		$generic_id{$parent_id} = $node_id;
			}
		}
	}
    
    #derives_from
    #relationship: derives_from PR:000002258 ! BH3-interacting domain death agonist isoform 1
    #intersection_of: derives_from PR:000029649 ! pro-neuregulin-1, membrane-bound isoform
	if ($term_line =~ /: derives_from (PR:\S+)/) {
	    $parent_id = $1;
	    $parent_type = "derives_from";
	    $parent = {
		      "id" => $parent_id,
		      "type" => $parent_type,
		      };
	#skip "PR:000000001" ! protein and "PR:000029032" ! Mus musculus protein and "PR:000018263" ! amino acid chain
	if (($parent_id ne "PR:000029032") and ($parent_id ne "PR:000000001")  and ($parent_id ne "PR:000018263")) {
		push @parent_nodes, $parent;
		}
	    }
  
    #link to gene
    #relationship: has_gene_template MGI:1197518 ! Aatk (mouse)
    #intersection_of: has_gene_template MGI:1274784 ! Adgrv1 (mouse) 
	if ($term_line =~ /: has_gene_template (MGI:\d+)/) {
	    $parent_id = $1;
	    $parent_type = "has_gene_template";
	    $parent = {
		      "id" => $parent_id,
		      "type" => $parent_type,
		      };
	    $mgi_id{$node_id} = $parent_id;
	    }
	    
    #synonym: "TGM3" EXACT PRO-short-label [PRO:DNx]
    #short name
	if ($term_line =~ /^synonym: "(.*)" EXACT PRO-short-label/) {
	    $short_name = $1;
	    }
	    
    #synonym: "TGM3" EXACT PRO-short-label [PRO:DNx]
    #other synonyms
    #synonym: "PNR" RELATED []
    #synonym: "Gm227" RELATED []
    #synonym: "TaR-5" EXACT []
    #synonym: "trace amine receptor 5" EXACT []
	if ($term_line =~ /^synonym: "(.*)" EXACT \[\]/) {
	    $synonym = $1;
	    push @synonyms, $synonym;
	    }	   
	    
    #non-mouse?
	if ($term_line =~ /: only_in_taxon NCBITaxon:(\d+)/) {
	    $taxon_id = $1;
	    if ($taxon_id != 10090) {
	    	$nonmouse = 1;
	    	}
	    }
	    
    #obsolete?
	if ($term_line =~ /^is_obsolete: true/) {
	    $is_obsolete = 1;
	    }

    #end of Term
	if ($term_line eq "") {
	    return();
	    }

    } #end read node
} #end get_node_info routine


#############
sub find_ancestors()
{
	#find parents/ancestors
	&find_parents($query_node);

		# select unique objects
		%seen = ();
    		$seen{$query_node}++;
		foreach $item (@ancestor_nodes) {
    			$seen{$item}++;
			}
		@ancestor_nodes = keys %seen;
	$nodes{$query_node}->{"ancestors"} = [ @ancestor_nodes ];
}

#############
sub find_parents()
{
    local $node_id = shift;
    local @parent_nodes = @ { $nodes{$node_id}->{"parents"} };

    #no more parents (stopping condition)
    if (!@parent_nodes) {
	return();
	}
    #
    foreach $parent (@parent_nodes) {
	$parent_id = $parent->{"id"};
	push @ancestor_nodes, $parent_id;
	&find_parents($parent_id);
	}
} #end find_parents routine

__END__

