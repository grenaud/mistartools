#!/usr/bin/env Rscript

library("ggplot2");
library("gplots");
library("reshape2");

args=(commandArgs(TRUE))
#paircoacompute2barplot [dstat file] [sample source] [pdf out prefix]

#print("reading data");
data <- read.table(args[1],header=FALSE);
#print("individuals found");

l<- unlist(  strsplit(  as.character( unlist(strsplit( as.character(data$V1) , c("-") )) ), c("@")) )

ind1<-l[c(TRUE, FALSE,FALSE)]
ind2<-l[c(FALSE,TRUE, FALSE)]
indS<-l[c(FALSE,FALSE,TRUE) ]




m<-as.matrix(data)
m[upper.tri(m)] <- NA

dat<-as.data.frame(cbind(as.matrix(ind1),as.matrix(ind2),as.matrix(indS), m[,seq(2,(length(data)-1))]) ,stringsAsFactors=FALSE)

dat<-dat[dat$V3==args[2],]

cols = seq(3,length(dat)-2);
dat[,cols] = apply(dat[,cols], 2, function(x) as.numeric(x));





val<-acast(dat, V1~V2, value.var="V34")
valmin<-acast(dat, V1~V2, value.var="V36")
valmax<-acast(dat, V1~V2, value.var="V35")
class(val)<-"numeric"
class(valmin)<-"numeric"
class(valmax)<-"numeric"


val<-    round(val*100,2)
valmin<- round(valmin*100,2)
valmax<- round(valmax*100,2)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)



pdf(paste(args[3],"heat.pdf",sep=""),width = 10, height = 10)
#png("filename.png",width=2000,height=800)
#pdf("filename.pdf",width=24,height=8)

p<-heatmap.2(val,
  cellnote = val,  # same data set for cell labels
  main = paste("Pairwise coalescence using source ",args[2],sep=""), # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  key=FALSE, 
  Colv="NA",
          Rowv="NA",
             lwid=c(0.1,4), lhei=c(0.1,0.4)
          )            # turn off column clustering
dev.off();


pdf(paste(args[3],"reorderheat.pdf",sep=""),width = 10, height = 10)
#png("filename.png",width=2000,height=800)
#pdf("filename.pdf",width=24,height=8)

p<-heatmap.2(val,
  cellnote = val,  # same data set for cell labels
  main = paste("Pairwise coalescence using source ",args[2],sep=""), # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
  key=FALSE, 
             lwid=c(0.1,4), lhei=c(0.1,0.4)
          )            # turn off column clustering
dev.off();
