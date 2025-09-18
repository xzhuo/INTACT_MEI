open IN, "zcat $ARGV[0] |" or die "Cannot open $ARGV[0]\n";
open OUT, ">$ARGV[1]" or die "Cannot write to $ARGV[1]\n";

while (<IN>) {
    chomp;
    next if /^#/;
    @a = split /\t/;
    %info=split /=|;/, $a[7];
    next unless exists $info{"SVLEN"} && abs($info{"SVLEN"})>=50;
    if ($info{"SVTYPE"} eq "INS") {
        print OUT ">$a[2]\n$a[4]\n";
    }elsif ($info{"SVTYPE"} eq "DEL") {
	print OUT ">$a[2]\n$a[3]\n";
    }
}

close IN;
close OUT;

