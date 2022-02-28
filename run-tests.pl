#!/usr/bin/env perl
#
#   Copyright (C) 2022 Genome Research Ltd.
#
#   Author: Petr Danecek <pd3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin";
use File::Temp qw/ tempfile tempdir /;
use Cwd qw/ abs_path /;

my $opts = parse_params();
parse_tests($opts);
find_bams($opts);
if ( $$opts{list_orphans} )
{
    list_orphans($opts);
}
else
{
    run_tests($opts);
}

exit 0;

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print
        "About: test bcftools calling using tests in the tests subdirectory. The test files are formatted as\n",
        "           cmd: command line\n",
        "           fmt: bcftools query format, the default is \"%CHROM:%POS %REF %ALT\", and this part is mandatory\n",
        "           exp: expected result\n",
        "\n",
        "       Unless the command contains the pipe symbol, it will be piped through\n",
        "           bcftools call -m | bcftools view -i'N_ALT>0' | bcftools norm\n",
        "\n",
        "       Small BAM tests can be produced with\n",
        "         ./misc/create-bam-test -b file.bam -f file.fa -r chr:pos -o prefix\n",
        "\n",
        "Usage: run-tests.pl [OPTIONS]\n",
        "Options:\n",
        "   -b, --bcftools PATH     BCFtools executable.\n",
        "   -d, --dry-run           Print commands but do not run anything.\n",
        "   -t, --temp-dir PATH     When given, temporary files will not be removed.\n",
        "   -h, -?, --help          This help message.\n",
        "\n";
    exit -1;
}

sub cygpath {
    my ($path) = @_;
    $path = `cygpath -m $path`;
    $path =~ s/\r?\n//;
    return $path
}
sub safe_tempdir
{
    my $dir = tempdir(CLEANUP=>1);
    if ($^O =~ /^msys/) {
        $dir = cygpath($dir);
    }
    return $dir;
}

sub parse_params
{
    my $opts = { nok=>0, nfailed=>0, bcftools=>'bcftools' };
    while (defined(my $arg=shift(@ARGV)))
    {
        if (                 $arg eq '--list-orphans' ) { $$opts{list_orphans}=1; next }
        if ( $arg eq '-d' or $arg eq '--dry-run' ) { $$opts{dry_run}=1; next }
        if ( $arg eq '-b' or $arg eq '--bcftools' ) { $$opts{bcftools}=shift(@ARGV); next }
        if ( $arg eq '-t' or $arg eq '--temp-dir' ) { $$opts{tmp}=shift(@ARGV); next }
        if ( $arg eq '-?' or $arg eq '-h' or $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{tmp}) ) { $$opts{tmp} = safe_tempdir; }
    $$opts{tests} = $FindBin::RealBin."/tests";
    $$opts{dat}   = $FindBin::RealBin."/dat";
    $$opts{bin}   = $FindBin::RealBin;
    $$opts{bin}   =~ s{run-tests.pl$}{};
    if ( $^O =~ /^msys/ )
    {
        $$opts{tests} = cygpath($$opts{tests});
        $$opts{dat} = cygpath($$opts{dat});
        $$opts{bin} = cygpath($$opts{bin});
    }
    return $opts;
}

sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid)
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    }
    else
    {
        # child
        exec('bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [bash -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub _cmd3
{
    my ($cmd) = @_;

    my $tmp = "$$opts{tmp}/tmp";
    my $pid = fork();
    if ( !$pid )
    {
        exec('bash', '-o','pipefail','-c', "($cmd) 2>$tmp.e >$tmp.o");
    }
    waitpid($pid,0);

    my $status  = $? >> 8;
    my $signal  = $? & 127;

    my (@out,@err);
    if ( open(my $fh,'<',"$tmp.o") )
    {
        @out = <$fh>;
        close($fh) or error("Failed to close $tmp.o");
    }
    if ( open(my $fh,'<',"$tmp.e") )
    {
        @err = <$fh>;
        close($fh) or error("Failed to close $tmp.e");
    }
    unlink("$tmp.o");
    unlink("$tmp.e");

    return ($status,\@out,join('',@err));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed: $cmd\n", $out); }
    return $out;
}

#--------------------------

# Parse tests.txt
sub parse_tests
{
    my ($opts) = @_;
    my @keys = (qw(cmd fmt exp));
    my %jobs = ();
    my @tests = glob("$$opts{tests}/*.txt");
    for my $test (@tests)
    {
        open(my $fh,'<',$test) or error("$test: $!");
        my $num  = 0;
        my $job  = {};
        while (my $line=<$fh>)
        {
            $num++;
            if ( $line=~/^#/ ) { next; }
            for my $key (@keys)
            {
                if ( !($line=~/^$key:/) ) { next; }
                $line =~ s/^$key:\s*//;
                $line =~ s/\s*$//;
                if ( scalar %$job == 0 )
                {
                    $$job{num}  = $num;
                    $$job{file} = $test;
                }
                $$job{$key} = $line;
            }
            my $complete = 1;
            for my $key (@keys)
            {
                if ( !exists($$job{$key}) && $key ne 'fmt' ) { $complete = 0; last; }
            }
            if ( $complete )
            {
                if ( !exists($$job{fmt}) ) { $$job{fmt} = '%CHROM:%POS %REF %ALT'; }
                my $id = "$$job{file}:$$job{num}";
                $id =~ s/\s+/_/g;
                $jobs{$id} = { %$job };
                $job = {};
            }
        }
        close($fh) or error("close failed: $test");
    }
    $$opts{tests} = \%jobs;
}

# Find bams in their subdirectories
sub find_bams
{
    my ($opts,@path) = @_;
    if ( !@path ) { push @path,$$opts{dat}; }
    my $dir = join('/',@path);
    opendir(my $dh,$dir) or error("opendir $dir: $!");
    while (my $item = readdir($dh))
    {
        if ( $item =~ /^\./ ) { next; }
        if ( -d "$dir/$item" ) { find_bams($opts,@path,$item); next; }
        if ( $item =~ /\.fa$/i or $item =~ /\.bam$/i or $item =~ /\.cram$/i )
        {
            my $path = join('/',@path[1..$#path],$item);
            if ( exists($$opts{paths}{$item}) )
            {
                error("Unique relative file names are required, found a clash:\n\t$$opts{paths}{$item}\n\t$path\n\n");
            }
            $$opts{paths}{$item} = $path;
        }
    }
    closedir($dh);
    $$opts{paths}{bcftools} = $$opts{bcftools};
}

sub list_orphans
{
    my ($opts) = @_;
    my %bnames = ();
    for my $file (keys %{$$opts{paths}})
    {
        if ( !($file =~ /\.fa$/i or $file =~ /\.bam$/i or $file =~ /\.cram$/i) ) { next; }
        my $bname = $`;
        $bnames{$bname} = $$opts{paths}{$file};
    }
    my %missing = ();
    for my $key (keys %{$$opts{tests}})
    {
        my $job = $$opts{tests}{$key};
        for my $type (qw(fa bam cram))
        {
            my $tmp = $$job{cmd};
            while ( $tmp=~/(\S+)\.$type/ )
            {
                my $bname = $1;
                $tmp = $';
                if ( exists($bnames{$bname}) ) { $bnames{$bname} = undef; }
            }
        }
    }
    for my $bname (keys %bnames)
    {
        if ( !defined($bnames{$bname}) ) { next; }
        if ( !($bnames{$bname}=~/\.[^.]+$/) ) { error("Could not parse: $bnames{$bname}\n"); }
        my $path = $`;
        print "$path.*\n";
    }
}

sub expand_paths
{
    my ($opts,$test,$args) = @_;
    my $out = $args;
    for my $type (qw(bcftools fa bam cram))
    {
        my $tmp = $out;
        $out = '';
        while ( $tmp=~/\S+\.$type/ )
        {
            my $prev = $`;
            my $file = $&;
            $tmp = $';
            if ( exists($$opts{paths}{$file}) ) { $file = $$opts{paths}{$file}; }
            if ( !-e "$$opts{dat}/$file" ) { error("\n\nNo such file referenced from $test\n\t$$opts{dat}/$file\n\n"); }
            $out .= $prev . $file;
        }
        $out .= $tmp;
    }
    return $out;
}

sub run_test
{
    my ($opts,$exp,$exe) = @_;

    my $stats = $$opts{stats};
    my ($ret,$out,$err) = _cmd3(qq[cd $$opts{dat} && $exe]);

    # Compare the output, the first field must be always "chr:pos ref alt"
    my $has_variant = 0;
    my $has_details = 0;
    my @exp = split(/\s+/,$exp);
    for my $line (@$out)
    {
        chomp($line);
        my @got = split(/\s+/,$line);
        $has_variant = 1;
        for (my $i=0; $i<3; $i++)
        {
            if ( $exp[$i] eq $got[$i] ) { next; }
            $has_variant = 0;
            last;
        }
        if ( !$has_variant ) { next; }
        $has_details = 1;
        for (my $i=3; $i<@exp; $i++)
        {
            if ( $exp[$i] eq $got[$i] ) { next; }
            $has_details = 0;
            last;
        }
        last;
    }

    # Report the status
    if ( $has_variant )
    {
        print ".. site ok\n";
        $$stats{nsites_match}++;
    }
    else
    {
        print ".. site missed\n";
        print "\texpected:$exp\n";
        print "\tfound:   ",join('',@$out),"\n";
        $$stats{nsites_missed}++;
    }
    if ( $has_variant && @exp>3 )
    {
        if ( $has_details )
        {
            print "\t.. details ok\n";
            $$stats{ndetails_match}++;
        }
        else
        {
            print "\t.. details missed\n";
            print "\texpected:$exp\n";
            print "\tfound:   ",join('',@$out),"\n";
            $$stats{ndetails_missed}++;
        }
    }
}

sub run_tests
{
    my ($opts) = @_;
    $$opts{stats} =
    {
        nsites_match    => 0,
        nsites_missed   => 0,
        ndetails_match  => 0,
        ndetails_missed => 0,

    };
    for my $test (sort {$$opts{tests}{$a}{num}<=>$$opts{tests}{$b}{num}} keys %{$$opts{tests}})
    {
        # Create the command line
        my $job = $$opts{tests}{$test};
        my ($chrpos,$ref,$alt,$gt) = split(/\s+/,$$job{exp});
        my @exe = $$job{cmd};
        if ( !($$job{cmd}=~/\|/) )
        {
            if ( !($$job{cmd}=~/\S+\.fa/) ) { error("Could not identify the fasta reference in \"$$job{cmd}\"\n"); }
            my $fa_ref = $&;
            push @exe,qq[ $$opts{bcftools} call -m ];
            push @exe,qq[ $$opts{bcftools} view -i 'N_ALT>0' ];
            push @exe,qq[ $$opts{bcftools} norm -f $fa_ref ];
        }
        push @exe,qq[ $$opts{bcftools} query -f '$$job{fmt}\\n' ];
        my $exe = expand_paths($opts, $test, join('|',@exe));
        print "\n# $test\n$exe\n";

        # Run the test
        if ( $$opts{dry_run} ) { next; }
        run_test($opts,$$job{exp},$exe);
    }

    # Print the summary status
    my $stats = $$opts{stats};
    my $nsites   = $$stats{nsites_match} + $$stats{nsites_missed};
    my $ndetails = $$stats{ndetails_match} + $$stats{ndetails_missed};
    print "\nNumber of sites tested:\n";
    printf "    total  .. %d\n", $nsites;
    printf "    passed .. %d  (%.1f%%)\n", $$stats{nsites_match}, $nsites ? 100.*$$stats{nsites_match}/$nsites : 0;
    printf "    missed .. %d  (%.1f%%)\n", $$stats{nsites_missed}, $nsites ? 100.*$$stats{nsites_missed}/$nsites : 0;
    print "\nNumber of details tested:\n";
    printf "    total  .. %d\n", $ndetails;
    printf "    passed .. %d  (%.1f%%)\n", $$stats{ndetails_match}, $ndetails ? 100.*$$stats{ndetails_match}/$ndetails : $ndetails;
    printf "    missed .. %d  (%.1f%%)\n", $$stats{ndetails_missed}, $ndetails ? 100.*$$stats{ndetails_missed}/$ndetails : $ndetails;
    print "\n";
}


