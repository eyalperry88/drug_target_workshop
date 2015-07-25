WORKING_DIR <- "D:\\Projects\\drug_target_workshop\\output_validations_no_FLT\\per_knockout"
setwd(WORKING_DIR)

pim <- as.numeric(read.table('PIM.txt')[,1])
akt <- as.numeric(read.table('AKT.txt')[,1])
pim.akt <- as.numeric(read.table('PIM+AKT.txt')[,1])
pi3k <- as.numeric(read.table('PI3K.txt')[,1])
mek <- as.numeric(read.table('MEK.txt')[,1])

pim.akt.expected = pim + akt

N <- 117

pim <- sort(pim)[1:N]
akt <- sort(akt)[1:N]
pim.akt <- sort(pim.akt)[1:N]
pi3k <- sort(pi3k)[1:N]
mek <- sort(mek)[1:N]
pim.akt.expected <- sort(pim.akt.expected)[1:N]

df1 <- data.frame(PIM=pim, AKT=akt, "exp. comb."=pim.akt.expected, "obs. comb."=pim.akt)
boxplot(df1, main="B2H distribution for different knockouts")

df2 <- data.frame(PIM=pim, MEK=mek, PI3K=pi3k)
png(filename="boxplot.png", units="in", width=8, height=8, res=300)
boxplot(df2, main="B2H distribution per known drug target", 
        xlab="Drug Targets", ylab="B2H score",
        boxwex=0.6)
dev.off()
