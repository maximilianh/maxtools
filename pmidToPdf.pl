#!/usr/bin/perl
# pubmed crawler
# Maximilian Haeussler, 11/2006

# input: pmid ids on stdin
# output: writes pdfs and other files to current dir
#         creates crawler.$LOG with list of failed pmids

# basic idea of this script is:
# - goto pubmed, get "link out" field from pubmed
# - follow this link to publishers website
# - grab everything that remotely resembles articles or suppl data so
#   everything that has .doc/.xls/.pdf in its name or in the text
#   of links on the page or in pages that contain "pdf" or "Supplement"
#   with special functions for blackwell(javascript)/oxford(frames) and
#   elsevier(scienceDirect-chooser)

use WWW::Mechanize;
use MIME::Base64;
use FileHandle;
$| = 1;

use vars qw($totalSize $MAXSIZE $LOG);
$MAXSIZE=100000000;
$totalSize=0;

# -------------- FUNCTIONS --------------------

sub saveFile() {
    my $agent = shift;
    my $link = shift;
    my $filename = shift;

    $agent->get($link->url());
    open PDF, ">$filename";
    my $content = $agent->response()->content;
    print PDF $content;
    &updateTotalSize($agent);
    close PDF;
    $agent->back();
}

sub getSupplData() {
    my $agent = shift;
    my $pmid = shift;
    my $pdf = shift;

    print STDERR "   SuppData: Searching for urls that points/contain to doc/xls...\n";
    @exts = ("xls", "doc");
    if ($pdf) {
        push(@exts, "pdf");
    }

    my $i = 1;
    foreach $ext (@exts) {
        print STDERR "   SuppData: Searching for $ext...\n";
        @links = $agent->find_all_links(url_regex => qr/\.$ext/i, n=>0 ) ;
        @links2 = $agent->find_all_links(text_regex => qr/\.$ext/i, n=>0 ) ;
        push (@links, @links2);
        foreach $link (@links) {
            my $fname = "$pmid.S$i.$ext";
            print STDERR "   SuppData: Found link with ext $ext, saving as $fname...\n";
            print $LOG " SuppData:$fname ";
            &saveFile($agent, $link, $fname); 
            $i++;
        }
    }
    if ($i==1) {
            print STDERR "   SuppData: No doc/xls found here...\n";
    }
    return ($i-1);
}


# Elsevier is giving you the choice between ScienceDirect or native Journal access
sub elsevierChooser() {
    $agent = shift;

    if ($agent->title() =~ /Locator/) {
        print $LOG " ElsevierChooser(LocatorInTitle) ";
        $agent->follow_link(text_regex => qr/Article via ScienceDirect/i);
    }
}

# JNLP needs another click
sub jnlpFollower() {
    $agent = shift;

    if ($agent->title() eq "resolver") {
        print $LOG " JNLPFollowLink ";
        print STDERR " doing a stupid click for JNLP...\n";
        $agent->follow_link(text_regex => qr/link page/i);
    }
}

# Blackwell tries to code it's url in javascript. We'll still get them.
sub handleBlackwell() {
    $agent = shift;
    $link = 0;
    if ($agent->title() =~ /Blackwell/) {
            print STDERR "blackwell journal, doing strange things now...";
            print $LOG " BlackwellInTitle ";
            $baseurl = "http://www.blackwell-synergy.com".$suffix; 
            $link = $agent->find_link(url_regex => qr/doi\/pdf/i);
            $url= $link->url();
            $url =~ /\/doi\/pdf.*/;
            $url = $&;
            $url =~ s/\'//;
            $url =~ s/\)//;
            $link = WWW::Mechanize::Link->new($url, "", "", "", "", "");
    }
    return $link;
}


sub handleOxford() {
    $agent = shift;

    $link = $agent->find_link(text_regex => qr/here/i);
    if ($link) {
        print $LOG " Oxford(HereLinkFound) ";
        print STDERR "Oxford Journal: Handle frame-websites: here / manual download...\n";
        $agent->get($link->url());
        &updateTotalSize($agent);
        $link = $agent->find_link(url_regex => qr/\.pdf/i ) ;
    }
    return $link;
}

sub renameOldLogs() {
    my $n = shift;
    rename $n.".10",$n.".11";
    rename $n.".9",$n.".10";
    rename $n.".8",$n.".9";
    rename $n.".7",$n.".8";
    rename $n.".6",$n.".7";
    rename $n.".5",$n.".4";
    rename $n.".4",$n.".5";
    rename $n.".3",$n.".4";
    rename $n.".1",$n.".2";
    rename $n."",$n.".1";
}

sub updateTotalSize(){
    my $agent = shift;
    my $content = $agent->response()->content;
    $totalSize+=length($content);
    printf "transferred size=%d, totalSize=%d\n", length($content), $totalSize;
    if ($totalSize > $MAXSIZE) {
        print STDERR "Maximum total size reached, waiting for 3600 seconds...\n";
        sleep 3600;
        $totalSize=0;
    }


}

# -------------- MAIN -------------------------

&renameOldLogs("crawler.log");

$LOG = new FileHandle;
$LOG->open(">crawler.log");
$LOG->autoflush(1);

# $LOG into inist, get cookie
print STDERR "Logging into INIST.fr ...\n";
$agent = WWW::Mechanize->new();
$agent->cookie_jar(HTTP::Cookies->new);
$url = 'http://gate1.inist.fr/login';
$agent->get($url);
$agent->form_number(1);
$agent->field('pass', 'GERS9V29');
$agent->field('user', 'SDVFRC2118');
$agent->click();
$suffix = ".gate1.inist.fr"; #see below

while (<>) {
        chomp;
        $pmid = $_;
        my $i = 0;
        print $LOG "$pmid\t";
        if (-e $pmid.".pdf") {
            print STDERR "$pmid.pdf already exists, skipping this pmid...\n";
            print $LOG " FileExists\n";
            next;
        }
        print STDERR "Getting publisher outlink site from NCBI for pmid $pmid...\n";
	$agent->get("http://www.ncbi.nlm.nih.gov.gate1.inist.fr/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=$pmid&retmode=ref&cmd=prlinks");
        # keeping track of transferred filesize
        &updateTotalSize($agent);

        if ($agent->title eq "Entrez PubMed") {
                print STDERR "No outlink exists for this pmid, giving up...\n";
                print $LOG "NoOutlink\n";
                next;
        }
        if ($agent->title =~ /Restriction de/) {
                print STDERR "Oh damn, now we're blocked again by inist...\n";
                exit 1;
        }

        print STDERR "(check for elsevier and jnlp)\n";
        &elsevierChooser($agent);
        &jnlpFollower($agent);

        $link = &handleBlackwell($agent);

        if (! $link) {
            $link = &handleOxford($agent);
        }

        if (! $link) {

            # try to find url that points to pdf
            print STDERR "Searching for url that points to pdf...\n";
            $link = $agent->find_link(url_regex => qr/pdf$/i ) ;

            if ($link) {
                print STDERR "Found direct link to .pdf ...\n";
            }
            
            # try to find link that contains pdf in description
            else {
                print STDERR "Searching for link-desc that contains pdf (+ test blackwell)...\n";
                $link = $agent->find_link(text_regex => qr/pdf/i ) ;
                if (!$link) {
                    print STDERR "Searching for link-url that contains pdf...\n";
                    $link = $agent->find_link(url_regex => qr/pdf/i ) ;
                }
                if (!$link) {
                    print STDERR "Searching for link-url that contains 'Full-text'...\n";
                    $link = $agent->find_link(url_regex => qr/Full-text/i ) ;
                }
                #print $agent->response()->content;
                if ($link) {
                        print STDERR "Following link, checking if answer is html...\n";
                        $agent->get($link->url());
                        &updateTotalSize($agent);
                        #print STDERR $agent->is_html;
                        if (! $agent->is_html) {
                            print STDERR "Link is not html, we'll treat it as a pdf.\n";
                            my $url = $agent->uri->as_string();
                            $link = WWW::Mechanize::Link->new($url, "", "", "", "", "");
                            print $LOG " PdfTextLinkFound(NoHtml):$url ";
                        }
                        else {
                            print STDERR "Searching again for link that points to pdf...\n";
                            $link = $agent->find_link(url_regex => qr/\.pdf/i ) ;

                            # special case for Oxford Journals' damn framed websites
                            if (! $link) {
                                $link = &handleOxford($agent);
                            } # end special case
                        }
                }
            }
        }

	if (! $link) {
            print STDERR "Link not found. Giving up. \n";
            print $LOG " noPdfLinkFound\n";
            next;
        }
        else {
            print STDERR "Saving pdf to $pmid.pdf.\n";
            print $LOG " PdfSuccess ";
            &saveFile($agent, $link, "$pmid.pdf");
        }
        
        # goback to mainpage and try to collect doc/xls
        print STDERR "Searchin for doc/xls files on main page...\n";
        #$agent->back();
        $nSuppl = &getSupplData($agent, $pmid, 0);

        print STDERR "Searching for link that contains Suppl...\n";
        $link = $agent->find_link(text_regex => qr/Supplementary/i);
        if ($link) {
            print STDERR "link found. Searching for Supplementary data...\n";
            $agent->get($link->url());
            &updateTotalSize($agent);
            $nSuppl += &getSupplData($agent, $pmid, 0);
        }
        else {
            print STDERR "no link with Suppl found\n";
        }

        if ($nSuppl==0) {
            print $LOG "noSuppl ";
        }
        print $LOG "\n";

        print STDERR "----------------------------------\n";
    }

    close $LOG;
