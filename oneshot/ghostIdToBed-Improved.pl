#!/usr/bin/perl
# need ghost cluster ids from stdin
# scrapes the jgi-id for given est-cluster id from ghost website

use WWW::Mechanize;
use MIME::Base64;
$| = 1;

my @args = (
    Authorization => "Basic " .
    MIME::Base64::encode( "guest" . ':' . "welcome" )
);
my $i = 0;
while (<>) {
    $acc = $_;
    chomp($acc);
    print STDERR $i++."  ".$acc."\n";
    $url = "http://hoya.zool.kyoto-u.ac.jp/cgi-bin/gbrowse/ci?name=in_situ:In_situ-".$acc;
    $m = WWW::Mechanize->new();
    $m->credentials($url, "guest", "welcome");
    $m->get($url, @args);
    $title = $m->title();
    @arr = split(/ /,$title);
    $addr = $arr[4];
    $addr =~ s/s/S/;
    $addr =~ s/:/\t/;
    $addr =~ s/[.][.]/\t/;
    $jgiacc=$addr;
    print $acc."\t".$jgiacc."\n";
}
