#!/usr/bin/gawk -f
# remove non-dna characters from fasta seqs, remove seqs shorter than 9bp (one sequence per line)
BEGIN { RS="\n"; FS=" ";} 

/^>/ { id=$0; next;}  # save fasta id 

/^[^>]/          {  
		  line=$0; 
		  t = gsub (/[^ACTGactg]/,"",line); 
		  if (length(line)>9)  {
		       print id;
		       print line; 
		       }
		}  
