#!/usr/bin/perl
# found on blog.patrick-morgan.net/2005/09/perl-script-to-convert-csv-to-tab.html
#  in the comments

while (<>) {
my @fields = ();
push(@fields, $+) while $_ =~ m{
("[^"]*(?:""[^"]*)*"),?
| ([^,\n\r]+),?
| ,
}gx;
push(@fields, undef) if substr($_, -1,1) eq ',';
print join("\t", @fields);
print "\n";
}
