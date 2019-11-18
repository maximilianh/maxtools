# pull out only the CGI error messages from our Apache logs, they're hopefully always above the Stack dump: line.

# example:
#[Tue May 07 14:55:42.876120 2019] [cgi:error] [pid 239917] [client 142.254.109.176:62980] AH01215: [Tue May  7 14:55:42 201] [hgTables] [client 142.254.109.176] [hgsid=725042543_hHaA3SUWUXZZLB] [/cgi-bin/hgTables] mySQL error 1044: Access denied for user 'hguser'@'localhost' to database 'hgTemp' (profile=<noProfile>, host=localhost, db=hg38), referer: http://hgw1.soe.ucsc.edu/cgi-bin/hgTables?hgsid=725042543_hHaA3SUWUXZZLBZ9yCASbGDHHWeq
#[Tue May 07 14:55:42.876133 2019] [cgi:error] [pid 239917] [client 142.254.109.176:62980] AH01215: , referer: http://hgw1.soe.ucsc.edu/cgi-bin/hgTables?hgsid=725042543_hHaA3SUWUXZZLBZ9yCASbGDHHWeq
#[Tue May 07 14:55:42.876170 2019] [cgi:error] [pid 239917] [client 142.254.109.176:62980] AH01215: Stack dump:, referer: http://hgw1.soe.ucsc.edu/cgi-bin/hgTables?hgsid=725042543_hHaA3SUWUXZZLBZ9yCASbGDHHWeq

import fileinput, re
regEx = re.compile(r'\[/cgi-bin/[^\s]+')

lastLine = None
for line in fileinput.input():
    line = line.rstrip("\n")
    if "[/cgi-bin/" in line:
        lastLine = line
    if "Stack dump:" in line:
        print(regEx.split(lastLine)[1].strip())
