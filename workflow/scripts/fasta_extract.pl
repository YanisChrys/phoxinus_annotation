#!/usr/bin/perl -w

=head1 Description

  read a FASTA file as STDIN
  if parameter 1 is a file, consider it as a list of accession numbers
  otherwise consider all paramters as accession numbers of sequences
  to be extrated from the FASTA file

=cut

#------------------------------------------------------------

sub usage( $ )
  {
    print STDERR "$_[0]\n";
    system("pod2text $0");
    exit(-1);
  }

#------------------------------------------------------------

sub goodAC( $ \% )
  {
    local ($s, $AC) = @_;
    return 1 if (exists($$AC{$s}));
    return 0;
  }

#------------------------------------------------------------
  
($#ARGV >= 0) || &usage("wrong number of parameters");

if (! -d $ARGV[0] && open(FIN,$ARGV[0])) # file or command providing a list of AC
  {
    while (<FIN>)
      {
	chomp;
	$AC{$_}='';
      }
    close(FIN);
  }
else # AC to be extracted are given as parameters
  {
    map{$AC{$_}='';} @ARGV;
  }

$extract=0;
while (<STDIN>)
  {
    chomp;
    if (($s) = (/^>(\S+)/))
      {
	if ($extract == 1) { $extract = 0;}
	$extract = &goodAC($s,\%AC);
      }
    if ($extract == 1)
      {
	print "$_\n";
      }
  }

