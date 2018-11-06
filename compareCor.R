#sanity check
library(synapser)
synLogin()
##1 compare correlation between my counts and Thanneer's counts (logCPM)
#current counts, kelsey
cc <- read.table(synGet("syn17015639")$path, sep = '\t', header=TRUE)
rownames(cc) <- cc[,1]
cc <- cc[,-1]
#former counts, thanneer
fc <- read.table(synGet("syn16847979")$path, sep = '\t', header=TRUE)
rownames(fc) <- fc[,1]
fc <- fc[,-1]

intsct <- intersect(colnames(cc), colnames(fc))
genes <- intersect(rownames(cc), rownames(fc))

c <- cc[genes,intsct]
f <- fc[genes,intsct]

compare <- cor(c,f, use = "complete.obs")
viz <- image(compare)

viz
