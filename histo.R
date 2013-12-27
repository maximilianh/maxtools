args <-  commandArgs(trailingOnly=TRUE)

datafile=args[1]
limit=as.numeric(args[2])
dist=as.numeric(args[3])
label=args[4]
graphfile=args[5]

q <- read.table(datafile, header=F)
attach(q)
# calc relative cg content
#cg = ncg / (ncg+nat)
#png(paste("/tmp/",name,".png", sep = "", collapse = NULL))
png(graphfile, width=1024, height=700)
hist(freq=FALSE, col='blue', V1, breaks=seq(0,max(V1)+dist+1, dist), xlim=c(0,limit), main=label)
#title(factor, "-500 , +500")
detach(q)
#locator()
dev.off()

