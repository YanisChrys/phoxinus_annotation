#!/usr/bin/perl

use Pod::Usage;
use Getopt::Long;

my($opt_help, $opt_man);

GetOptions(
	'help' => \$opt_help,
	'man'  => \$opt_man,
)
or pod2usage( "Try '$0 --help' for more information.");

pod2usage( -verbose => 1 ) if $opt_help;
pod2usage( -verbose => 2 ) if $opt_man;

#pod2usage( -msg => 'Wrong number of params') unless ($ARGV[0] && $ARGV[1]);
#pod2usage( -msg => 'Only param 1 OR 2 could be -') if ($ARGV[0] eq '-' && $ARGV[1] eq '-');

$scope = @ARGV||'all';
for (@ARGV) {$_--};
%c = map { $_, 1 } @ARGV;
while (<STDIN>) {
	chomp;
	@t = split(/\t/,$_);
	for (my $i = 0 ; $i <= $#t ; $i++) {
		next if ($scope ne 'all' && !exists $c{$i});
		next if ($t[$i] =~ /,/);
		$a=reverse($t[$i]);
		$a=~s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
		$t[$i]=scalar reverse $a;
	}
	print join("\t",@t)."\n";
}

=pod

=head1 NAME

 thousandify.pl

=head1 SYNOPSIS

 thousandify.pl <column to thousandify> < file.tsv

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

 Reads a tsv file as STDIN. 
 Print to STDOUT the same tsv file with columns given as params thousandified
 i.e. using comma as thousands separator (9123456 will give 9,123,456)

=head1 AUTHORS

Cedric Cabau

=head1 VERSION

1

=head1 DATE

2013

=head1 KEYWORDS

tsv column thousand

=head1 EXAMPLE

 cat file.tsv | thousandify.pl 2 5 # thousandify columns 2 and 5

=cut
