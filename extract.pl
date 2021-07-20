##!/usr/bin/perl -w
### Extract the target site #

my $l=0;
my $Q1=20;
my $Q2=0.2;

parse_command_line();
open FIG,    "<$input";
open OUT,    ">$output";

$cnt=0;

foreach $word(<FIG>){
         chomp($word=$word);
         @file=split(/\s+/,$word);
         $cnt=$cnt+1;
         if($cnt%4==1){$name=substr($word,1);}
         elsif($cnt%4==2){$read=$word;}
         elsif($cnt%4==0){$qual=$word;
              $len=length $read;
              $qual_sum=0;
             # Quality filter
              for($i=0;$i<$len;$i++){
                            $tt=substr($qual,0,$i);
                            $tt=ord($tt)-33;
                            if($tt<$Q1){$qual_sum++;}               
              }
			 # Extract the target site
              if($qual_sum>=$len*$Q2){;}
              else{ 
                   $tag=0;
                   if($read=~/$base/){
                       while($read=~/$base/g ){ 
                       $end=pos($read);
                       $ss=substr($read,$end-$l,$l);
                       $tag++;
                       }
                   print OUT ">$name\n$ss\n";
                   }                               
              }                     
         }
}
           
sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	    if       ($_ =~ /^-i$/)  { $input   = shift  @ARGV; }
        elsif    ($_ =~ /^-o$/)  { $output  = shift  @ARGV; }
        elsif    ($_ =~ /^-b$/)  { $base    = shift  @ARGV; }
        elsif    ($_ =~ /^-l$/)  { $l       = shift  @ARGV; }
        elsif    ($_ =~ /^-Q1$/) { $Q1      = shift  @ARGV; }
        elsif    ($_ =~ /^-Q2$/) { $Q2      = shift  @ARGV; }
        else {
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	    }  
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl reads_filter.pl -i input  -o output -b -Q1  -Q2 [-h]
    i    :input file 
    o    :outputfile
    b    :target site : ATCATC[ATGC]{58}TGTTGC
    l    :length of read
    Q1   :threshold for low quality score [20].
    Q2   :maximum percent of low-quality bases[0.2].  
    h    :display the help information.
EOQ
exit(0);
}
