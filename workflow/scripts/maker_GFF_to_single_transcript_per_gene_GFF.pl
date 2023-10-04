#!/usr/bin/env perl

=head1 Name

maker_GFF_to_single_transcript_per_gene_GFF.pl -- for each gene, keep the transcript with the lowest
AED score, the largest number of exons and the largest size.

=head1 Description

This program reads GFF maker file and write to stdout a new GFF only with maker records (col2 = "maker") 
and with only one transcript per gene. The selected transcript is the one with the lowest AED score, 
the largest number of exons and the largest size.

=head1 Usage

 maker_GFF_to_single_transcript_per_gene_GFF.pl all.maker.gff
 
=cut

use strict;

my $usage = "\nUsage: maker_GFF_to_single_transcript_per_gene_GFF.pl maker.all.gff > maker.no_alternatives.gff\n\n";

die $usage unless (-f $ARGV[0]);
my $gff = $ARGV[0];
my %genes;   # $genes{'gene'}->['trans1', 'trans2', ...]
my %trans;   # $trans{'transcript'}->{'AED'} # $trans{'transcript'}->{'nb_exons'} # $trans{'transcript'}->{'size'} # $trans{'transcript'}->{'start'} # $trans{'transcript'}->{'end'}
my %contig_length; # $contig_length{'contig'}->contig_length

open (my $fh, $gff) or die "Error, cannot open file $gff";
while (my $line = <$fh>) {
	last if ($line eq "##FASTA\n");
	chomp($line);
	my @cols = split(/\t/, $line);
	$contig_length{$cols[0]} = $cols[4] if ($cols[2] eq 'contig');
	next unless $cols[1] eq 'maker';
	if (($cols[2] eq 'mRNA' || $cols[2] eq 'tRNA') && $line =~ /ID=([^;]+);Parent=([^;]+);.*_AED=(\d\.\d+?)\;/) {
		push(@{$genes{$2}->{mRNA}}, $1);
		$trans{$1}->{AED} = $3;
	}
	if ($cols[2] eq 'exon' && $line =~ /Parent=([^;]+)/) {
		my $list = $1;
		chomp($list);
		foreach my $transcript (split(/,/, $list)) {
			$trans{$transcript}->{nb_exons}++;
			$trans{$transcript}->{size} += $cols[4]-$cols[3]+1;
			$trans{$transcript}->{start} = $cols[3] unless ($trans{$transcript}->{start} && $trans{$transcript}->{start} < $cols[3]);
			$trans{$transcript}->{stop} = $cols[4] unless ($trans{$transcript}->{stop} && $trans{$transcript}->{stop} > $cols[4]);
		}
	}
}
my %selected;
foreach my $gene (keys %genes) {
	my @sorted = sort { $trans{$a}->{AED} <=> $trans{$b}->{AED} || $trans{$b}->{nb_exons} <=> $trans{$a}->{nb_exons} || $trans{$b}->{size} <=> $trans{$a}->{size} } @{$genes{$gene}->{mRNA}};
	$selected{$sorted[0]} = 1;
	$genes{$gene}->{start} = $trans{$sorted[0]}->{start};
	$genes{$gene}->{stop} = $trans{$sorted[0]}->{stop};
}
seek($fh, 0, 0);
print "##gff-version 3\n";
while (my $line = <$fh>) {
	last if ($line eq "##FASTA\n");
	my @cols = split(/\t/, $line);
	next unless $cols[1] eq 'maker';
	# check that element coordinates are not larger than scaffold length
	if (exists($contig_length{$cols[0]}) && $cols[4] > $contig_length{$cols[0]} && $cols[8] =~ /^ID=([^;]+)/) {
		printf STDERR ("Stop moved from %d to %d (%d) for %s %s to avoid exceeding scaffold length\n", $cols[4], $contig_length{$cols[0]}, $contig_length{$cols[0]}-$cols[4], $cols[2], $1);
		$cols[4] = $contig_length{$cols[0]};
	}
	if ($cols[2] eq 'gene' && $cols[8] =~ /^ID=([^;]+)/) {
		# check that gene coordinates cover mRNA coordinates
		if ($cols[3] > $genes{$1}->{start}) { 
			printf STDERR ("Start moved from %d to %d (%d) for gene %s to respect mRNA coordinates\n", $cols[3], $genes{$1}->{start}, $genes{$1}->{start}-$cols[3], $1);
			$cols[3] = $genes{$1}->{start};
		}
		if ($cols[4] < $genes{$1}->{stop}) {
			if (exists($contig_length{$cols[0]}) && $genes{$1}->{stop} <= $contig_length{$cols[0]}) {
				printf STDERR ("Stop moved from %d to %d (%d) for gene %s to respect mRNA coordinates\n", $cols[4], $genes{$1}->{stop}, $genes{$1}->{stop}-$cols[4], $1);
				$cols[4] = $genes{$1}->{stop};
			}
		}
		print join("\t", @cols);
	} elsif (($cols[2] eq 'mRNA' || $cols[2] eq 'tRNA') && $cols[8] =~ /^ID=([^;]+)/) {
		my $id = $1;
		if (exists $selected{$id}) {
			# replace the mRNA name which is not unique for mRNA from Augustus prediction
			# remove maker metrics attributes
			$cols[8] =~ s/Name=.+/Name=$id/;
			print join("\t", @cols);
		}
	} elsif ($cols[8] =~ /Parent=([^;]+)/) {
		my $parentList = $1;
		chomp($parentList);
		foreach my $transcript (split(/,/, $parentList)) {
			if (exists $selected{$transcript}) {
				# replace the list of parents by the single transcript for each exon
				$cols[8] =~ s/$parentList/$transcript/;
				if ($cols[2] eq 'exon') {
					# replace the ID by the single transcript with incremented number of exon for each exon
					$cols[8] =~ s/ID=[^;]+/ID=$transcript:$selected{$transcript}/;
					$selected{$transcript}++;
				}
				print join("\t", @cols);
				last;
			}
		}
	}
}
close $fh;
exit 0;
