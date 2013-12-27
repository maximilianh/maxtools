#!/usr/bin/gawk -f
# replace lowercase sequence letters with uppercase N
BEGIN { RS="\n"; FS=" ";} 

/^>/ { print; next; } 

/^[^>]/          {  
		  line=$0; 
		  gsub (/[nactg]/,"N",line); 
          print line; 
		}  
