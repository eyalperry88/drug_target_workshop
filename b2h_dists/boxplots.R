WORKING_DIR <- "D:/Projects/drug_target_workshop/b2h_dists"
setwd(WORKING_DIR)

# FLT patients
pim <- as.numeric(read.table('FLT3/PIM.txt')[,1])
pi3k <- as.numeric(read.table('FLT3/PI3K.txt')[,1])
mek <- as.numeric(read.table('FLT3/MEK.txt')[,1])

# remove outliers
N <- ceiling(length(pim) * 0.95)

pim <- sort(pim)[1:N]
pi3k <- sort(pi3k)[1:N]
mek <- sort(mek)[1:N]

df <- data.frame(PIM=pim, MEK=mek, PI3K=pi3k)
png(filename="FLT-boxplot.png", units="in", width=8, height=8, res=300)
boxplot(df, main="B2H distribution per known drug target", 
        xlab="Drug Targets", ylab="B2H score",
        boxwex=0.6)
dev.off()



# non FLT patients
pim <- as.numeric(read.table('nonFLT3/PIM.txt')[,1])
pi3k <- as.numeric(read.table('nonFLT3/PI3K.txt')[,1])
mek <- as.numeric(read.table('nonFLT3/MEK.txt')[,1])

# remove outliers
N <- ceiling(length(pim) * 0.95)

pim <- sort(pim)[1:N]
pi3k <- sort(pi3k)[1:N]
mek <- sort(mek)[1:N]

df <- data.frame(PIM=pim, MEK=mek, PI3K=pi3k)
png(filename="nonFLT-boxplot.png", units="in", width=8, height=8, res=300)
boxplot(df, main="B2H distribution per known drug target", 
        xlab="Drug Targets", ylab="B2H score",
        boxwex=0.6)
dev.off()
