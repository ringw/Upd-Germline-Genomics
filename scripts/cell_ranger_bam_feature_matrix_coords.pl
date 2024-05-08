#!/usr/bin/perl -w
# Extract known Cell Ranger matrix coords (feature and cell barcode) from
# Cell Ranger SAM output. By piping to sort, uniq -c, and awk (change order of
# fields), we will create a Matrix Market sparse matrix of counts.
#
# Arguments:
# STDIN - SAM format
# Cells txt - Cell barcodes, newline separated. Unlike 10X barcodes.tsv, there
# is no header row, and no multiple fields per row!
# 

my $list_of_cells_file = $ARGV[0];
open my $list_of_cells_fh, "<", $list_of_cells_file or die $list_of_cells_file;
my @cells_arr = <$list_of_cells_fh>;
close $list_of_cells_fh;
chomp @cells_arr;
my %cells;
# ONE-BASED Matrix Market output. The %cells hash is not a reverse lookup into a
# Perl array, per se, but has an increment of one added to it.
@cells{@cells_arr} = 1..@cells_arr;

my $list_of_transcripts_file = $ARGV[1];
open my $list_of_transcripts_fh, "<", $list_of_transcripts_file or die $list_of_transcripts_file;
my @transcripts_arr = <$list_of_transcripts_fh>;
chomp @transcripts_arr;
close $list_of_transcripts_fh;
@transcripts{@transcripts_arr} = 1..@transcripts_arr;

my $transcripts_tag = $ARGV[2];

while (<STDIN>) {
  next unless /CB:Z:([ATCG]+-[0-9])/;
  next unless exists($cells{$1});
  my $cell_barcode = $1;
  next unless /xf:i:25/;
  next unless /RE:A:E/;
  for my $transcript_id (/(?:$transcripts_tag:Z:|\G;)([A-z]+\d+)[^; \t]*/go) {
    next unless exists($transcripts{$transcript_id});
    print "$transcripts{$transcript_id}\t$cells{$cell_barcode}\n";
  }
}