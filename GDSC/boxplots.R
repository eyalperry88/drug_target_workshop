WORKING_DIR <- "D:/Projects/drug_target_workshop/GDSC"
setwd(WORKING_DIR)

# FLT patients



ap <- as.numeric(read.table('results_per_drug/AP-24534.txt')[,1])
az <- as.numeric(read.table('results_per_drug/AZD7762.txt')[,1])
bo <- as.numeric(read.table('results_per_drug/Bosutinib.txt')[,1])
ce <- as.numeric(read.table('results_per_drug/CEP-701.txt')[,1])
fh <- as.numeric(read.table('results_per_drug/FH535.txt')[,1])
gw <- as.numeric(read.table('results_per_drug/GW441756.txt')[,1])
so <- as.numeric(read.table('results_per_drug/Sorafenib.txt')[,1])
su <- as.numeric(read.table('results_per_drug/Sunitinib.txt')[,1])

# remove outliers
N <- 44

ap <- sort(ap)[1:N]
az <- sort(az)[1:N]
bo <- sort(bo)[1:N]
ce <- sort(ce)[1:N]
fh <- sort(fh)[1:N]
gw <- sort(gw)[1:N]
so <- sort(so)[1:N]
su <- sort(su)[1:N]

drugs <- c('Sorafenib', 'AZD7762', 'AP-24534', 'Sunitinib', 
           'GW441756', 'CEP-701', 'Bosutinib', 'FH535')

png(filename="noFLT-boxplot.png", units="in", width=8, height=8, res=300)
par(las=2)
boxplot(so, az, ap, su, gw, ce, bo, fh,
        main="B2H distribution per drug", 
        names=drugs, xlab="", ylab="B2H score",
        boxwex=0.6)
dev.off()

