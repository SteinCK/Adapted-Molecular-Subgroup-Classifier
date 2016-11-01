

# install Geoquery, affy, and gcrma

#source("https://bioconductor.org/biocLite.R")
#biocLite( 'GEOquery' )
library( 'GEOquery' )
#biocLite( 'affy' )
library( 'affy' )
#biocLite( 'gcrma' )
library( 'gcrma' )
#install.packages( 'e1071' )
library( e1071 )

# load model 

load( "~/TC-6 model.RData")

# download and prepare data

DAT <- getGEOSuppFiles( 'GSE38627' )  # download CEL files from GEO, large download 735MB

untar( '~/GSE38627/GSE38627_RAW.tar' , exdir="GSE38627/CEL/") # unpack tarball

GZ.NAMES <- dir( "~/GSE38627/CEL/", full.names=TRUE) # paths for compressed CEL files

sapply( 1:length( GZ.NAMES), function(x) gunzip( GZ.NAMES[x] ))  # unzip CEL files

CEL.NAMES <- dir( "~/GSE38627/CEL/", full.names=TRUE)

# read affy and normalize

affy.data <- ReadAffy( filenames = CEL.NAMES )

GCRMA <- gcrma( affy.data , fast=FALSE )

GSE.GCRMA <- exprs( GCRMA )

#### calculate subgroups

TC.RAW <- predict( TC12.MODEL , t(GSE.GCRMA[TC12.PROBES,]) )

TC6 <- sapply( strsplit( as.character( TC.RAW ) , ":" ), "[[" , 1)
TC6[TC6 %in% c("11q13")] <- "CCND1-11q13"
TC6[TC6 %in% c("6p21")] <- "CCND3-6p21"
TC6[TC6 %in% c("16q23 or 8q24","20q11")] <- "MAF"
TC6[TC6 %in% c("4p16")] <- "MMSET"
TC6 <- factor( TC6 , levels=c("D1","D2","CCND1-11q13","CCND3-6p21","MMSET","MAF"))

## scatterplots of subgroups

# CCND2 vs CCND1
plot( GSE.GCRMA["208711_s_at",] ,
      GSE.GCRMA["200953_s_at",], xlim=c(1,16), ylim=c(1,16), 
      bg=c(rainbow( 6 ))[TC6], pch=c(23,22,24,24,25,25)[TC6], xlab="CCND1", ylab="CCND2")
legend( "topright", levels( TC6 ), pch=c(23,22,24,24,25,25), pt.bg=rainbow( 6 ))

# MMSET vs FGFR3
plot( GSE.GCRMA["204379_s_at",] ,
      GSE.GCRMA["222777_s_at",], xlim=c(1,16), ylim=c(1,16), 
      bg=c(rainbow( 6 ))[TC6], pch=c(23,22,24,24,25,25)[TC6], xlab="FGFR3", ylab="MMSET")
legend( "topright", levels( TC6 ), pch=c(23,22,24,24,25,25), pt.bg=rainbow( 6 ))

# MAFB vs MAF
plot( GSE.GCRMA["209347_s_at",] ,
      GSE.GCRMA["218559_s_at",], xlim=c(1,16), ylim=c(1,16), 
      bg=c(rainbow( 6 ))[TC6], pch=c(23,22,24,24,25,25)[TC6], xlab="MAF", ylab="MAFB")
legend( "topright", levels( TC6 ), pch=c(23,22,24,24,25,25), pt.bg=rainbow( 6 ))

# ISL2 vs CCND1
plot( GSE.GCRMA["208711_s_at",] ,
      GSE.GCRMA["232352_at",], xlim=c(1,16), ylim=c(1,16), 
      bg=c(rainbow( 6 ))[TC6], pch=c(23,22,24,24,25,25)[TC6], xlab="MAF", ylab="MAFB")
legend( "topright", levels( TC6 ), pch=c(23,22,24,24,25,25), pt.bg=rainbow( 6 ))


