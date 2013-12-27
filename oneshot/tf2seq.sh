#!/bin/bash
# prints protein sequences for all factors as fasta with matrix id as name
cat /usr/lib/cgi-bin/biobase/transfac/10.2/data/factor.dat | gawk '
/^MX/ {gsub("\\.",""); 
    id = $3; 
    if (seq!="") 
        { print ">"id; print seq; } 
    else 
        { print id" no sequence" > "/dev/stderr"; } 
        }
/^SQ/ { seq = seq $2; inseq=1; } 
/^[/][/]/ {seq=""; inseq=0}' >mat2seq.fa 2> noseq.lst

