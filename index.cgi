#!/usr/bin/perl -w
use strict;
use File::Spec;
use Data::Dumper;
use IPC::Open2;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Bio::SeqIO;
use CGI;
use CGI::Carp qw(fatalsToBrowser warningsToBrowser);
use Time::Piece;
use Time::Seconds;
use URI::Escape;
use Sys::Hostname;
use List::Util qw(min max);

sub nocase { lc($a) cmp lc($b) }

my @tool = qw(blastp blastp blastx tblastn tblastx);

my $q = CGI->new();

my $HOMELINK = $q->p($q->a({-href=>$q->url},'Start a new search'));

print $q->header;

my $cfg = read_config();
#print $q->pre(Dumper($cfg));

my $title = sprintf("%s BLAST+", $cfg->{name});
print $q->start_html(
  -title=>$title,
  -style=>{ -src=>'blast.css' },
);

my $ver = blast_version();
if (my $img = $cfg->{logo}) {
  print $q->img({-src=>$img,-align=>'right'});
}
print $q->h1("$title $ver");
print $HOMELINK;

if ($q->url_param('db') and $q->url_param('id')) {
#  print Dump;
  my $id = $q->url_param('id');
  my $db = $q->url_param('db');
  print $q->h2($id);
  my $exe = File::Spec->catfile($cfg->{exe_path}, 'blastdbcmd');
#  $db = File::Spec->catfile($cfg->{db_path}, $db);
  my $cmd = "$exe -db $db -entry $id";
  print $q->p("Running: <TT>$cmd</TT>");
  print "<PRE>";
  system("$exe -db '$db' -entry '$id' 2>&1");
  print "</PRE>";
}
elsif ($q->param) {
  print $q->h2("Input");
  print "<UL>\n";
#  print $q->pre(Dumper($q->param));
  my $seq = $q->param('seq') or die "need sequence";
  $seq = clean_sequence($seq);
  my $seq_dna = is_dna($seq);
  print $q->li("Treating query sequence as ", $seq_dna ? 'DNA' : 'protein');
  print $q->li("Query sequence has ~", length($seq), " letters");
  my $db_dna = 0;
  my $db = 'NO_DATABASE';
  my $tool = 'NO_TOOL';
  if ($q->param('db') and $q->param('db') =~ m/^(Protein|Nucleotide) (.*)/) {
    $db_dna = $1 eq 'Nucleotide';
    $db = $2;
    if ($db_dna) {
      $tool = $seq_dna ? 'blastn' : 'tblastn';
    } else {
      $tool = $seq_dna ? 'blastx' : 'blastp';
    }
  }
  else {
    die "Must specify a database";
  }
  print $q->li("Database <TT>$db</TT> is ", $db_dna ? 'DNA' : 'protein');
  print $q->li("Using BLAST tool <TT>$tool</TT>");
  my $have_cores = num_cpu();
  my $cores = min($have_cores, 8);
  print $q->li("Identified $have_cores CPU cores, using $cores.");
  my $evalue = $q->param('evalue') || 1E-3;
  print $q->li("Using e-value threshold $evalue");
  print "</UL>\n";

  my $na = $q->param('num_aln') || 50;
  $na =~ m/^(\d+)$/ or die "Invalid num_aln";
  
  print $q->h2("Execute");
  print $q->p("System status: ", $q->tt( qx(uptime) ) );
  my $exe = File::Spec->catfile($cfg->{exe_path}, $tool);
  my $cmd = "nice $exe -db $db -evalue $evalue -num_threads $cores";
  $cmd .= " -num_descriptions $na -num_alignments $na";
  my $style = $q->param('style') || 'HTML';
  $cmd .= " -html" if $style eq 'HTML';

  print $q->p("Running: <TT>$cmd</TT>");
  print $q->p("Please wait...");
  my $t1 = localtime;
  print $q->h2("Result");

  my $pid = open2(\*OUT, \*IN, $cmd);
  print IN $seq;
  close IN;
  
  if ($style eq 'text') {
    print "<PRE>\n";
    print while (<OUT>);
    print "</PRE>\n";
  }
  elsif ($style eq 'HTML') {
    print while (<OUT>);
  }
  else {
    my $in = Bio::SearchIO->new(-format=>'blast', -fh=>\*OUT);
    my $writer = Bio::SearchIO::Writer::HTMLResultWriter->new;
    $writer->start_report( sub { '' } );
    $writer->title( sub { '' } );
    $writer->introduction( sub { '' } );
    $writer->end_report( sub { '' } );
    sub mylink {
      my($self,$hit,$result) = @_;
      #return "hit.name=".$hit->name." db=".$result->database_name;
#      $q->a({-href=>$q->url."?id=".uri_escape($hit->name).
#                    ";db=".uri_escape($result->database_name) }, $hit->name);
      $q->a({-href=>$q->url."?id=".uri_escape($hit->name).
                    ";db=".uri_escape($db) }, $hit->name);
    }
    $writer->hit_link_desc(\&mylink);
    $writer->hit_link_align(\&mylink);
    my $out = Bio::SearchIO->new(-writer=>$writer);
    print "<DIV ID='blastreport'>\n";
    $out->write_result($in->next_result);
    print "</DIV>\n";
 }

  print $q->h2("Timing");
  my $t2 = localtime;
  my $secs = $t2 - $t1;
  print $q->p("Search took ", $secs, " seconds to run on $cores CPU cores.");
  
#  print $q->p($q->a({-href=>$q->url},'Start a new search'));
}
else 
{
  print $q->start_form();
  
  print $q->h2("Query sequence");
  print $q->textarea(-class=>'seqtextarea', -name=>'seq',-rows=>15,-columns=>80,
#    -value=>join('',<DATA>),
  );
  print q{<br><input type="button" value="Clear" onclick="this.form.elements['seq'].value=''">};
  print $q->submit(-value=>'BLAST!');

  print $q->h2("Database");
  my $groups;
  for my $mol ('Protein', 'Nucleotide') {
    my $dbs = get_dbs($cfg->{db_path}, $mol);
    my(@values,%labels);
    for my $db_loc (sort { nocase } keys %{$dbs}) {
      my $n = "$mol $db_loc";
      push(@values, $n);
      $labels{$n} = $dbs->{$db_loc};
    }
    my $group = $q->optgroup(-name=>$mol, -values => \@values, -labels=>\%labels);
    #my $group = $q->optgroup(-name=>$mol, -values => [map {$_ . "_$mol" } @db], -labels=>$db->{$mol});
    $groups .= $group;
  }

  print $q->popup_menu(
    -name=>"db",
    -values =>['(Choose a database)', $groups]
  );

  print $q->h2("E-value");
#  print $q->textfield(-name=>'evalue',-size=>5,-value=>'0.01');
  print $q->popup_menu(-name=>'evalue', -value=>[ map { 10**-$_ } -1..10], -default=>0.01);
  #print $q->pre(Dumper($db));

  print $q->h2("Report");
  print $q->radio_group(-name=>'style', -values=>['fancy','HTML', 'text'], -default=>'fancy');

#  print $q->h2("Sensitivity/Speed trade-off");
#  print $q->radio_group(-name=>'tradeoff', -values=>['default','fast','faster','fastest']);

#  print $q->h2("Tool");
#  print $q->popup_menu(
#   -name => 'tool',
#   -values => \@tool,
#  );

  print $q->h2("Number of alignments");
  print $q->popup_menu(-name=>'num_aln', -value=>[1,10,50,100,250,500], -default=>50);
  
  print $q->h2("Submit");
  print $q->p, $q->submit(-value=>'BLAST!');

#  my @opt = blast_options();
#  print $q->pre(map { "$_\n" } @opt);

  print $q->end_form();
  
  my @core = <core.*>;
#  print $q->p("For admin use only: <TT>@core</TT>");
  unlink @core if @core;

}

print $HOMELINK;

my $admin = $cfg->{admin_name} || 'Administrator';
my $email = $cfg->{admin_email} || 'root@'.hostname;

print $q->hr;
print $q->p("Please contact <A HREF='mailto:$email'>$admin</A> if you have problems.");
print $q->end_html;

#-----------------------------------------------------------------

sub blast_options {
  my $bin = File::Spec->catfile($cfg->{exe_path}, $tool[0]);
  my @opt = grep { m/^\s+-/ } qx($bin -help);
  chomp @opt;
  return @opt;
}

sub blast_version {
  my $bin = File::Spec->catfile($cfg->{exe_path}, $tool[0]);
#  print $q->p("BIN=$bin");
  my($line) = qx($bin -version);
#  print $q->p("LINE=$line");
  chomp $line;
  $line =~ m/(\d.*\d)/;
  return $1 || '0';  
}


#-----------------------------------------------------------------

sub get_dbs {
  my($dirs, $moltype) = @_;
  $moltype = (!$moltype or $moltype =~ m/p/i) ? 'p' : 'n';
  my $EXT = "\\.$moltype(in|al)\$"; # .pin/.pal or .nin/.nal
  my %db;
  for my $dir (split ' ', $dirs) {
    opendir DIR, $dir;
    my @idx = grep { m/$EXT$/ } readdir(DIR);
    closedir DIR;
    for my $idx (@idx) {
      $idx =~ s/(\.\d\d\d?)?$EXT$//; # multi-part blast indices
      my(undef,undef,$fname) = File::Spec->splitpath($idx);
      #$db{$fname} = File::Spec->catpath(undef,$dir,$idx);
      $db{ File::Spec->catpath(undef,$dir,$idx) } = $fname;
    }
  }
  return \%db;
}


#-----------------------------------------------------------------

sub read_config {
  my %cfg;
  open my $fh, '<', 'blast.conf';
  while (<$fh>) {
    chomp;
    next unless m/^(\w+)\s+(.*)$/; # should skip comments as only \w
    $cfg{$1} ||= $2;
  }
  return \%cfg;
}

#-----------------------------------------------------------------

sub num_cpu {
  if ($^O =~ m/linux/i) {
    my($num) = qx(grep -c ^processor /proc/cpuinfo);
    chomp $num;
    return $num if $num =~ m/^\d+/;
  }
  elsif ($^O =~ m/darwin/i) {
    my ($num) = qx(system_profiler SPHardwareDataType | grep Cores);
    $num =~ /.*Cores: (\d+)/;
    $num = $1;
    return $num if $num =~ m/^\d+/;
  }
  return 1; # safe default
}


#-----------------------------------------------------------------

sub is_dna {
  my($s) = @_;
  # http://effectiveperl.blogspot.com/2006/01/idiomatic-perl-counting-number-of.html
  my $nts = uc($s) =~ tr/ATGCN//;  # count nucleotides
  return 1 if $nts > 0.95*length($s);
}


#-----------------------------------------------------------------

sub clean_sequence {
  my($s) = @_;
  $s = ">sequence\n$s\n" if $s !~ m/^\s*>/xms;
  return $s;
}


#-----------------------------------------------------------------

__DATA__
>fadD FAD-dependent dehydrogenase
MIQELQLRVVPEVAGEYELLKTFVSKTLKIGIQEIWHIEILSRSIDARQKTIYFNLKVLV
FIGENYVEKQIALPDFSDVSNAKEVIVVGAGPAGLFACLQLILSGFKPILLERGKDVMKR
PFDLKEVNIHHNVNEDSNYCFGEGGAGTYSDGKLYTRSKKRGNVRQILEWLVGFGANKDI
LVEAHPHIGTNKLPKIVKNIREKIIETGGEIHFEKRVTDLLLNGNQIQGVVTKDGDKVYA
KNIILATGHSARDIFELLHQKGIELELKPIAVGVRVEHQQSLIDSIQYKCEVKSTYLPPS
PYNIVKQVDGRGVYSFCMCPGGVIAACATKPEEIVTNGWSSSKRARPSANSGIVVELKSE
DFKPFAKFGPLAAMEFQKEIERKAWVAGGRTQKAPAQKLMDFVQGKLSTDLPKTSYTPGI
TSVVLGEVLPRFVYQSLQKGFQEFDRSMKGYLTNEAVVHAPETRTSSPVCIPRDPNSLQH
VRIQGLYPCGEGAGYAGGIVSAAMDGIRSAHACVSSMI

