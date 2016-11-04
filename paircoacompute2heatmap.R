#!/usr/bin/env Rscript

library("ggplot2");
library("gplots");
library("reshape2");

args=(commandArgs(TRUE))
#paircoacompute2barplot [mst file] [samples to include, comma delim] [pdf out prefix]

data <- read.table(args[1],header=FALSE);

ind1<-unlist(strsplit( as.character(data$V1) , "-" ))[c(TRUE, FALSE)]
ind2<-unlist(strsplit( as.character(data$V1) , "-" ))[c(FALSE,TRUE) ]


m<-as.matrix(data)

dat<-as.data.frame(cbind(as.matrix(ind1),as.matrix(ind2), m[,seq(2,(length(data)-1))]) ,stringsAsFactors=FALSE)


cols = seq(3,length(dat)-2);
dat[,cols] = apply(dat[,cols], 2, function(x) as.numeric(x));


#lengthindinclude<-length(strsplit( args[2]  , "," )[[1]])

vectorin<-rep(FALSE,length(dat$V1));
#for(indinc in  strsplit( "15302,15304,Andaman,Dinka_A,Mbuti_A,French_A,Papuan_A,Han_A,Yoruba_A,San_A,Mandenka_A,Stuttgart,Loschbour,French_B,Han_B,Mandenka_B,Mbuti_B,Papuan_B,San_B,Yoruba_B,Australian1_B,Dinka_B,Ust"     , "," )[[1]]){
for(indinc in  strsplit( args[2]     , "," )[[1]]){
    print(indinc)
    vectorin <- vectorin | dat$V1==indinc;
}
dat<-dat[vectorin  , ]

vectorin<-rep(FALSE,length(dat$V1));

for(indinc in  strsplit( args[2]     , "," )[[1]]){
    print(indinc)
    vectorin <- vectorin | dat$V2==indinc;
}
dat<-dat[vectorin  , ]

print(dat$V1)
print(dat$V2)


val<-acast(dat, V1~V2, value.var="V65")
valmin<-acast(dat, V1~V2, value.var="V66")
valmax<-acast(dat, V1~V2, value.var="V67")
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
  main = "Pairwise coalescence", # heat map title
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
  main = "Pairwise coalescence", # heat map title
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
