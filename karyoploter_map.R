##########################
# map

library("karyoploteR")
library(biomaRt)
library(regioneR)
library("yarrr")


#####################
# genes

data=read.table("mart_export_map.tsv", sep="\t", header=T)
head(data)
genes <- toGRanges(data)
genes


####################
# plot
png("map.png", width=900, height=900)

kp <- plotKaryotype(genome = "BSgenome.Celegans.UCSC.ce11", plot.type=2, chromosomes=c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"))
kpAddBaseNumbers(kp, tick.dist=5000000, add.units=TRUE, cex=0.8)

#----------------------------------------------
###chrI

#hT3
kpRect(kp, chr="chrI", x0=1, x1=5245743, y0=0.2, y1=0.8, col="grey", r0=0.2, r1=0)
kpRect(kp, chr="chrX", x0=14773538, x1=17600000, y0=0.2, y1=0.8, col="grey", r0=0.2, r1=0)

#szT1 - translocation and free duplication
kpRect(kp, chr="chrI", x0=1, x1=7631470, y0=0.2, y1=0.8, col=yarrr::transparent("#FF0000FF", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrX", x0=1, x1=2597439, y0=0.2, y1=0.8, col=yarrr::transparent("#FF0000FF", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrI", x0=1, x1=7454116, y0=0.2, y1=0.8, col="#FF0000FF", r0=0.4, r1=0.2)
kpRect(kp, chr="chrX", x0=2144605, x1=17600000, y0=0.2, y1=0.8, col="#FF0000FF", r0=0.4, r1=0.2)

#hIn1 - inversion
kpRect(kp, chr="chrI", x0=11249952, x1=15003739, y0=0.2, y1=0.8, col=yarrr::transparent("#B1FF00FF", trans.val=.6), r0=0.2, r1=0)
kpRect(kp, chr="chrI", x0=11592246, x1=14863482, y0=0.2, y1=0.8, col="#B1FF00FF", r0=0.2, r1=0)

#hDf15 - deletion
kpRect(kp, chr="chrI", x0=12431742, x1=14006044, y0=0.2, y1=0.8, col=yarrr::transparent("#009BFFFF", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrI", x0=13300000, x1=13910000, y0=0.2, y1=0.8, col="#009BFFFF", r0=0.4, r1=0.2)

#hT1 - Translocation
kpRect(kp, chr="chrI", x0=1, x1=8409987, y0=0.2, y1=0.8, col=yarrr::transparent("slateblue", trans.val=.6), r0=0.6, r1=0.4)
kpRect(kp, chr="chrV", x0=1, x1=7207631, y0=0.2, y1=0.8, col=yarrr::transparent("slateblue", trans.val=.6) , r0=0.6, r1=0.4)
#kpRect(kp, chr="chrI", x0=1, x1=8409987, y0=0.2, y1=0.8, col="slateblue", r0=0.6, r1=0.4)
kpRect(kp, chr="chrV", x0=1, x1=6515240, y0=0.2, y1=0.8, col="slateblue", r0=0.6, r1=0.4)

#hT2 - translocation
kpRect(kp, chr="chrI", x0=1, x1=13187133, y0=0.2, y1=0.8, col=yarrr::transparent("violet", trans.val=.6), r0=0.8, r1=0.6)
kpRect(kp, chr="chrIII", x0=4822648, x1=13700000, y0=0.2, y1=0.8, col=yarrr::transparent("violet", trans.val=.6), r0=0.8, r1=0.6)
kpRect(kp, chr="chrI", x0=1, x1=12513420, y0=0.2, y1=0.8, col="violet", r0=0.8, r1=0.6)
kpRect(kp, chr="chrIII", x0=5107330, x1=13500000, y0=0.2, y1=0.8, col="violet", r0=0.8, r1=0.6)

#hDf8
kpRect(kp, chr="chrI", x0=1, x1=11034814, y0=0.2, y1=0.8, col=yarrr::transparent("#009BFFFF", trans.val=.6), r0=1, r1=0.8)
kpRect(kp, chr="chrI", x0=1, x1=6845142, y0=0.2, y1=0.8, col="#009BFFFF", r0=1, r1=0.8)



#----------------------------------------------
###chrII

#mT1 - translocation
kpRect(kp, chr="chrII", x0=6296872, x1=15000000, y0=0.2, y1=0.8, col=yarrr::transparent("violetred", trans.val=.6), r0=0.2, r1=0)
kpRect(kp, chr="chrIII", x0=3635354, x1=13700000, y0=0.2, y1=0.8, col=yarrr::transparent("violetred", trans.val=.6), r0=0.2, r1=0)
kpRect(kp, chr="chrII", x0=6710149, x1=15000000, y0=0.2, y1=0.8, col="violetred", r0=0.2, r1=0)
kpRect(kp, chr="chrIII", x0=3040846, x1=13500000, y0=0.2, y1=0.8, col="violetred", r0=0.2, r1=0)

#mIn1 - inversion
kpRect(kp, chr="chrII", x0=3553628, x1=12704681, y0=0.2, y1=0.8, col=yarrr::transparent("#B1FF00FF", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrII", x0=3982549, x1=12176069, y0=0.2, y1=0.8, col="#B1FF00FF", r0=0.4, r1=0.2)

#mnC1 - complex inversion
kpRect(kp, chr="chrII", x0=4904692, x1=14909258, y0=0.2, y1=0.8, col=yarrr::transparent("#B1FF00FF", trans.val=.6), r0=0.6, r1=0.4)
kpRect(kp, chr="chrII", x0=6710149, x1=14684456, y0=0.2, y1=0.8, col="#B1FF00FF", r0=0.6, r1=0.4)

#----------------------------------------------
###chrIII

#nDf6 - deletion
kpRect(kp, chr="chrIII", x0=3607207, x1=3695411, y0=0.2, y1=0.8, col="deepskyblue2", r0=1, r1=0.8)

#sC1 - complex inversion
kpRect(kp, chr="chrIII", x0=323321, x1=4641137, y0=0.2, y1=0.8, col=yarrr::transparent("#B1FF00FF", trans.val=.6), r0=0.8, r1=0.6)
kpRect(kp, chr="chrIII", x0=491547, x1=3040846, y0=0.2, y1=0.8, col="#B1FF00FF", r0=0.8, r1=0.6)

#qC1 - complex inversion
kpRect(kp, chr="chrIII", x0=1286122, x1=13737952, y0=0.2, y1=0.8, col=yarrr::transparent("#B1FF00FF", trans.val=.6), r0=0.6, r1=0.4)
kpRect(kp, chr="chrIII", x0=2034988, x1=11195201, y0=0.2, y1=0.8, col="#B1FF00FF", r0=0.6, r1=0.4)

#eT1 - reciprocal translocation
kpRect(kp, chr="chrIII", x0=8200764, x1=13700000, y0=0.2, y1=0.8, col=yarrr::transparent("purple", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrV", x0=1, x1=8930675, y0=0.2, y1=0.8, col=yarrr::transparent("purple", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrIII", x0=8193843, x1=13700000, y0=0.2, y1=0.8, col="purple", r0=0.4, r1=0.2)
kpRect(kp, chr="chrV", x0=1, x1=8939684, y0=0.2, y1=0.8, col="purple", r0=0.4, r1=0.2)

#----------------------------------------------
## chrIV

#nT1 - chromoplexy
kpRect(kp, chr="chrIV", x0=1901208, x1=17500000, y0=0.2, y1=0.8, col=yarrr::transparent("plum", trans.val=.6), r0=0.2, r1=0)
kpRect(kp, chr="chrV", x0=1, x1=16832779, y0=0.2, y1=0.8, col=yarrr::transparent("plum", trans.val=.6), r0=0.2, r1=0)
kpRect(kp, chr="chrIV", x0=3618259, x1=17500000, y0=0.2, y1=0.8, col="plum", r0=0.2, r1=0)
kpRect(kp, chr="chrV", x0=1, x1=15080152, y0=0.2, y1=0.8, col="plum", r0=0.2, r1=0)

#sDf22 - deletion
kpRect(kp, chr="chrIV", x0=12526000, x1=13141000, y0=0.2, y1=0.8, col=yarrr::transparent("#009BFFFF", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrIV", x0=12567720, x1=12788344, y0=0.2, y1=0.8, col="#009BFFFF", r0=0.4, r1=0.2)

#mnDp33 - Translocated duplication
kpRect(kp, chr="chrIV", x0=2853000, x1=4028000, y0=0.2, y1=0.8, col=yarrr::transparent("#FF0000FF", trans.val=.6), r0=0.4, r1=0.2)
kpRect(kp, chr="chrIV", x0=2984263, x1=3958552, y0=0.2, y1=0.8, col="#FF0000FF", r0=0.4, r1=0.2)

#----------------------------------------------
##chrV

#sC4
kpRect(kp, chr="chrV", x0=16000000, x1=21000000, y0=0.2, y1=0.8, col=yarrr::transparent("#009BFFFF", trans.val=.6), r0=0.6, r1=0.4)
kpRect(kp, chr="chrV", x0=15071760, x1=21000000, y0=0.2, y1=0.8, col="#009BFFFF", r0=0.6, r1=0.4)

#eDf43 - deletion + CGR
kpRect(kp, chr="chrV", x0=3466642, x1=3714670, y0=0.2, y1=0.8, col="#009BFFFF", r0=0.8, r1=0.6)

#e917 - inversion + CGR
kpRect(kp, chr="chrV", x0=4488422, x1=14577013, y0=0.2, y1=0.8, col="#B1FF00FF", r0=0.8, r1=0.6)

#----------------------------------------------
##chrX

#stDp2
kpRect(kp, chr="chrX", x0=6842000, x1=10309000, y0=0.2, y1=0.8, col=yarrr::transparent("#FF0000FF", trans.val=.6), r0=0.2, r1=0)
kpRect(kp, chr="chrX", x0=6889674, x1=10133182, y0=0.2, y1=0.8, col="#FF0000FF", r0=0.2, r1=0)

#mnDp3 - free duplication
kpRect(kp, chr="chrX", x0=12803000, x1=17600000, y0=0.2, y1=0.8, col=yarrr::transparent("#FF0000FF", trans.val=.6), r0=0.8, r1=0.6)
kpRect(kp, chr="chrX", x0=13201770, x1=17600000, y0=0.2, y1=0.8, col="#FF0000FF", r0=0.8, r1=0.6)

#mnDp1 - translocated duplication
kpRect(kp, chr="chrX", x0=13794000, x1=17600000, y0=0.2, y1=0.8, col="grey", r0=0.6, r1=0.4)


kpPlotMarkers(kp, data=genes, labels=genes$value, data.panel=2, cex=0.8, adjust.label.position = TRUE, text.orientation="vertical", r1=0.5)

dev.off()


