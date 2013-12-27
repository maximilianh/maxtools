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
    $url = "http://ghost.zool.kyoto-u.ac.jp/cgi-bin3/txtgetr2.cgi?".$acc;
    $m = WWW::Mechanize->new();
    $m->credentials($url, "guest", "welcome");
    $m->get($url);
    $link=$m->find_link(text => "cDNA browser");
    $gbrowselink=$m->find_link(text => "genome browser");
    if ($link=="") { 
	    $link=$gbrowselink;
	    $url = $link->URI();
	    $m->get($url, @args);
	    $title = $m->title();
	    @arr = split(/ /,$title);
	    $addr = $arr[4];
	    $addr =~ s/s/S/;
	    $addr =~ s/:/\t/;
	    $addr =~ s/[.][.]/\t/;
            $jgiacc=$addr;
	    } 
    else {
            $url = $link->URI();
            print STDERR "cdna browser\n";
	    $m->credentials($url, "guest", "welcome");
	    $m->get($url, @args);
	    $link=$m->find_link(text => "The best-hit gene model (JGI)");
	    if ($link=="") {
		$link=$m->find_link(text => "The best-hit gene model (Kyoto)");
		if ($link=="") {
			$link=$gbrowselink;
			if ($link=="") {
				$jgiacc="FuckingGhost";
			}
			else {
			$url=$link->URI();
			    $m->get($url, @args);
			    $title = $m->title();
			    @arr = split(/ /,$title);
			    $addr = $arr[4];
			    $addr =~ s/s/S/;
			    $addr =~ s/:/\t/;
			    $addr =~ s/[.][.]/\t/;
			    $jgiacc=$addr;
			}

		}
                else {
		$url=$link->URI();
		    $m->get($url, @args);
		    $title = $m->title();
		    @arr = split(/ /,$title);
		    $addr = $arr[4];
		    $addr =~ s/s/S/;
		    $addr =~ s/:/\t/;
		    $addr =~ s/[.][.]/\t/;
		    $jgiacc=$addr;
		}
	}
	    else {
		    $url=$link->URI();
		    @arr = split(/=/,$url);
		    $jgiacc = $arr[2];
		    $jgiacc = "ci01".$jgiacc;

		}
	}
	    print $acc."\t".$jgiacc."\n";
}
