##########################
# CB3475

library("karyoploteR")
library("RColorBrewer")

####################
# data coverage

data=read.table("~/Documents/cedar/Tatiana/Strains/strains/CB3475/coverage.tsv", sep="\t", header=1)
mydata=toGRanges(data)


####################
# data SNP

snp=read.table("~/Documents/cedar/Tatiana/Strains/strains/CB3475/snp.tsv", sep="\t", header=1)
mysnp=toGRanges(snp)


####################
# breakpoints

start=read.table("start.tsv", sep="\t")
end=read.table("end.tsv", sep="\t")

mystart=toGRanges(start)
myend=toGRanges(end)

####################
# plot
png("CB3475.png", width=900, height=900)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[1])

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.4)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.8, pch=1, col=brewer.pal(n=8, name="Set2")[1])

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

dev.off()


tiff("CB3475_2.tiff", width=400, height=200)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=1, chromosomes=c("chrI", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.6)
kpLines(kp, data = mydata, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[1])

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.6)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.6, col=brewer.pal(n=8, name="Set2")[1])

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

dev.off()

##########################
# BC1217

library("karyoploteR")


####################
# data coverage

data=read.table("~/Documents/cedar/Tatiana/Strains/strains/BC1217/coverage.tsv", sep="\t", header=1)
mydata=toGRanges(data)


####################
# data SNP

snp=read.table("~/Documents/cedar/Tatiana/Strains/strains/BC1217/SNP.tsv", sep="\t", header=1)
mysnp=toGRanges(snp)


####################
# breakpoints

start=read.table("start.tsv", sep="\t")
end=read.table("end.tsv", sep="\t")

mystart=toGRanges(start)
myend=toGRanges(end)

####################
# plot
png("BC1217.png", width=900, height=900)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "#AAFFAA")
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata, ymax=2, r0=at$r0, r1=at$r1, col="blue")

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "#AAAAFF")
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.4)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.8, pch=1, col="blue")

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

dev.off()





##########################
# KR2839

library("karyoploteR")


####################
# data coverage

data=read.table("coverage.tsv", sep="\t", header=1)
mydata=toGRanges(data)



####################
# data SNP

snp=read.table("~/Documents/cedar/Tatiana/Strains/pilot/visualization/circos/BC986_rufus_het.csv")

temp=paste("chr",snp$V1, sep="")
head(temp)

snp2=data.frame(temp, snp$V2, snp$V3, "value"=snp$V4)
mysnp=toGRanges(snp2)


####################
# breakpoints

start=read.table("start.tsv", sep="\t")
end=read.table("end.tsv", sep="\t")

mystart=toGRanges(start)
myend=toGRanges(end)

####################
# plot

#------------------------------------
png("KR2839.png", width=900, height=900)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "#AAFFAA")
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata, ymax=2, r0=at$r0, r1=at$r1, col="blue")

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "#AAAAFF")
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.4)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.8, pch=1, col="blue")

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

dev.off()



##########################
##########################
# CZ1072

library("karyoploteR")
library("RColorBrewer")

####################
# data coverage

data=read.table("coverage.tsv", sep="\t", header=1)
mydata=toGRanges(data)



####################
# data SNP

snp=read.table("snp.tsv", header=1)
mysnp=toGRanges(snp)


####################
# breakpoints

start=read.table("start.tsv", sep="\t")
end=read.table("end.tsv", sep="\t")

mystart=toGRanges(start)
myend=toGRanges(end)

####################
# plot

#------------------------------------
png("CZ1072.png", width=800, height=800)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=8, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata, ymax=8, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[1])

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.4)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.6, col=brewer.pal(n=8, name="Set2")[1])

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

kpDataBackground(kp, data.panel = 2, col=brewer.pal(n=9, name="Set3")[9], r0=0.1, r1=0.3)
#all
kpPlotRegions(kp, data=c("chrI:10565980-10577022"), col="orange", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrII:15271369-15273155"), col="red", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrIII:3891209-3898488"), col="orange", data.panel=2, r0=0.1, r1=0.3)

dev.off()



##########################
##########################
# RW6002

library("karyoploteR")
library("RColorBrewer")

####################
# data coverage

data=read.table("coverage.tsv", sep="\t", header=1)
mydata=toGRanges(data)

####################
# data SNP

snp=read.table("snp.tsv", header=1)
mysnp=toGRanges(snp)


####################
# breakpoints

start=read.table("start.tsv", sep="\t")
end=read.table("end.tsv", sep="\t")

mystart=toGRanges(start)
myend=toGRanges(end)

####################
# plot

#------------------------------------
png("RW6002.png", width=800, height=800)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=4, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata, ymax=4, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[1])

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.4)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.6, col=brewer.pal(n=8, name="Set2")[1])

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

kpDataBackground(kp, data.panel = 2, col=brewer.pal(n=9, name="Set3")[9], r0=0.1, r1=0.3)
#all
kpPlotRegions(kp, data=c("chrII:3910889-3910983"), col="green", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrIII:11830810-11834895"), col="orange", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrIV:5370371-5371720"), col="blue", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrIV:13847160-13922385"), col="blue", data.panel=2, r0=0.1, r1=0.3)

dev.off()

#------------------------------------
tiff("RW6002_2.tiff", width=400, height=200)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=1, chromosomes=c("chrII","chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.6)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=4, r0=at$r0, r1=at$r1, side=2, cex=0.6)
kpLines(kp, data = mydata, ymax=4, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[1])

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel2")[5])
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.6)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.6, col=brewer.pal(n=8, name="Set2")[1])

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

#kpDataBackground(kp, data.panel = 2, col=brewer.pal(n=9, name="Set3")[9], r0=0.1, r1=0.3)

dev.off()

##########################
# UV5

library("karyoploteR")


####################
# data coverage

data=read.table("coverage.tsv", sep="\t", header=1)
mydata=toGRanges(data)



####################
# data SNP

snp=read.table("~/Documents/cedar/Tatiana/Strains/pilot/visualization/circos/BC986_rufus_het.csv")

temp=paste("chr",snp$V1, sep="")
head(temp)

snp2=data.frame(temp, snp$V2, snp$V3, "value"=snp$V4)
mysnp=toGRanges(snp2)


####################
# breakpoints

start=read.table("start.tsv", sep="\t")
end=read.table("end.tsv", sep="\t")

mystart=toGRanges(start)
myend=toGRanges(end)

####################
# plot

#------------------------------------
png("UV5.png", width=900, height=900)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "#AAFFAA")
kpAddLabels(kp, labels = "Coverage", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata, ymax=2, r0=at$r0, r1=at$r1, col="blue")

at <- autotrack(current.track = 2, total.tracks = 2)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = "#AAAAFF")
kpAddLabels(kp, labels = "SNP", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, numticks = 2, tick.pos = c(0.5, 1), labels = c("Het", "Hom"), r0=at$r0, r1=at$r1, side=2, cex=0.4)
kpPoints(kp, data = mysnp, ymax=2, r0=at$r0, r1=at$r1, cex=0.8, pch=1, col="blue")

kpPlotRegions(kp, mystart, r0=0, r1=0.5, col="#ff8d92")
kpPlotRegions(kp, myend, r0=0, r1=0.5, col="#8d9aff")
kpPlotLinks(kp, data=mystart, data2=myend, col="#fac7ffaa", r0=0.5)

dev.off()



##########################
# FX strains

library("karyoploteR")
library("RColorBrewer")

####################
# data coverage FX30144

data30133=read.table("FX30133/coverage.tsv", sep="\t", header=1)
mydata30133=toGRanges(data30133)

data30144=read.table("FX30144/coverage.tsv", sep="\t", header=1)
mydata30144=toGRanges(data30144)

data30235=read.table("FX30235/coverage.tsv", sep="\t", header=1)
mydata30235=toGRanges(data30235)

data30237=read.table("FX30237/coverage.tsv", sep="\t", header=1)
mydata30237=toGRanges(data30237)

data30257=read.table("FX30257/coverage.tsv", sep="\t", header=1)
mydata30257=toGRanges(data30257)

####################
# breakpoints

start30133=read.table("FX30133/start.tsv", sep="\t")
end30133=read.table("FX30133/end.tsv", sep="\t")
mystart30133=toGRanges(start30133)
myend30133=toGRanges(end30133)

start30235=read.table("FX30235/start.tsv", sep="\t")
end30235=read.table("FX30235/end.tsv", sep="\t")
mystart30235=toGRanges(start30235)
myend30235=toGRanges(end30235)

start30144=read.table("FX30144/start.tsv", sep="\t")
end30144=read.table("FX30144/end.tsv", sep="\t")
mystart30144=toGRanges(start30144)
myend30144=toGRanges(end30144)

start30257=read.table("FX30257/start.tsv", sep="\t")
end30257=read.table("FX30257/end.tsv", sep="\t")
mystart30257=toGRanges(start30257)
myend30257=toGRanges(end30257)

start30237=read.table("FX30237/start.tsv", sep="\t")
end30237=read.table("FX30237/end.tsv", sep="\t")
mystart30237=toGRanges(start30237)
myend30237=toGRanges(end30237)

startFX=read.table("FX_start.tsv", sep="\t")
endFX=read.table("FX_end.tsv", sep="\t")
mystartFX=toGRanges(startFX)
myendFX=toGRanges(endFX)

startFX133_237=read.table("FX30133_FX30237_start.tsv", sep="\t")
endFX133_237=read.table("FX30133_FX30237_end.tsv", sep="\t")
mystartFX133_237=toGRanges(startFX133_237)
myendFX133_237=toGRanges(endFX133_237)

####################
# plot

#------------------------------------
png("FXstrains.png", width=900, height=900)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.5)

at <- autotrack(current.track = 1, total.tracks = 5)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel1")[9])
kpAddLabels(kp, labels = "FX30133", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata30133, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[1])

at <- autotrack(current.track = 2, total.tracks = 5)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel1")[9])
kpAddLabels(kp, labels = "FX30144", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata30144, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[2])

at <- autotrack(current.track = 3, total.tracks = 5)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel1")[9])
kpAddLabels(kp, labels = "FX30235", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata30235, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[3])

at <- autotrack(current.track = 4, total.tracks = 5)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel1")[9])
kpAddLabels(kp, labels = "FX30237", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata30237, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[4])

at <- autotrack(current.track = 5, total.tracks = 5)
kpDataBackground(kp, r0=at$r0, r1=at$r1, color = brewer.pal(n=8, name="Pastel1")[9])
kpAddLabels(kp, labels = "FX30257", r0=at$r0, r1=at$r1, cex=0.8)
kpAxis(kp, ymin=0, ymax=2, r0=at$r0, r1=at$r1, side=2, cex=0.45)
kpLines(kp, data = mydata30257, ymax=2, r0=at$r0, r1=at$r1, col=brewer.pal(n=8, name="Set2")[5])

kpPlotLinks(kp, data=mystart30235, data2=myend30235, col=brewer.pal(n=5, name="Set2")[3], r0=0)
kpPlotLinks(kp, data=mystart30144, data2=myend30144, col=brewer.pal(n=5, name="Set2")[2], r0=0)
kpPlotLinks(kp, data=mystart30257, data2=myend30257, col=brewer.pal(n=5, name="Set2")[5], r0=0)
kpPlotLinks(kp, data=mystart30237, data2=myend30237, col=brewer.pal(n=5, name="Set2")[4], r0=0)
kpPlotLinks(kp, data=mystart30133, data2=myend30133, col=brewer.pal(n=5, name="Set2")[1], r0=0)
#kpPlotLinks(kp, data=mystartFX, data2=myendFX, col="yellow", r0=0)
#kpPlotLinks(kp, data=mystartFX133_237, data2=myendFX133_237, col="black", r0=0)

kpDataBackground(kp, data.panel = 2, col=brewer.pal(n=9, name="Set3")[9], r0=0.1, r1=0.3,)
#all
kpPlotRegions(kp, data=c("chrI:11849276-11851612"), col="blue", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrV:18429890-18465644"), col="blue", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrV:6236672-6239149"), col="orange", data.panel=2, r0=0.1, r1=0.3)
kpPlotRegions(kp, data=c("chrI:12586101-12586223"), col="blue", data.panel=2, r0=0.1, r1=0.3)
#2
kpPlotRegions(kp, data=c("chrX:2932305-2951620"), col="orange", data.panel=2, r0=0.1, r1=0.2)
kpPlotRegions(kp, data=c("chrX:2939368-2940583"), col="orange", data.panel=2, r0=0.2, r1=0.3)
#FX30144
kpPlotRegions(kp, data=c("chrII:5131837-5805620"), col="orange", data.panel=2, r0=0.1, r1=0.2)
kpPlotRegions(kp, data=c("chrII:5136381-5806749"), col="blue", data.panel=2, r0=0.2, r1=0.3)
kpPlotRegions(kp, data=c("chrIV:16464405-16465988"), col="blue", data.panel=2, r0=0.1, r1=0.15)
kpPlotRegions(kp, data=c("chrIV:16464489-16465367"), col="blue", data.panel=2, r0=0.15, r1=0.2)
kpPlotRegions(kp, data=c("chrIV:16469112-16510351"), col="green", data.panel=2, r0=0.2, r1=0.25)
kpPlotRegions(kp, data=c("chrIV:16462855-16478758"), col="green", data.panel=2, r0=0.25, r1=0.3)
#FX30133
kpPlotRegions(kp, data=c("chrV:15826754-15886806"), col="blue", data.panel=2, r0=0.1, r1=0.3)

#start.regs <- toGRanges(data.frame("chrII", 10e6, 10e6))
#end.regs <- toGRanges(data.frame("chrV", 19e6, 19e6))
#kpPlotLinks(kp, data=start.regs, data2=end.regs, col=brewer.pal(n=5, name="Set2")[8], data.panel=2, r0=0.1)

dev.off()



