#!/usr/bin/perl -w

if($#ARGV+1 != 1) {die "requires path of input directory";};

my $dirname = $ARGV[0];

print "input_dir=".$dirname."\n";

my $outdir = "fgs_output";

unless(-e $outdir or mkdir $outdir) {
    die "Unable to create $outdir\n";
}




opendir($indir, "$dirname");
@files = readdir($indir);
closedir $indir;
foreach $folder (@files)
{
#    print "$folder\n";
    my $folder_path = $dirname."/".$folder;
    if(-d  $folder_path)
    {
        
        if (-e "$folder_path"."/contigs" && -e "$folder_path"."/Features/peg/tbl") {
            print "processing $folder\n";            
            $prefix = "$outdir"."/fgs_train_"."$folder";
            $status = `./rast_train_prepare.py -d $folder_path -p $prefix`;
            print "processing $folder done. $status\n";
        }
    }
}
