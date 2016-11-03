#!/usr/bin/env Rscript

library("ggplot2");

args=(commandArgs(TRUE))
#paircoacompute2barplot [mst file] [samples to include, comma delim] [ancient samples to exclude, comma delim] [pdf out]

data <- read.table(args[1],header=FALSE);

ind1<-unlist(strsplit( as.character(data$V1) , "-" ))[c(TRUE, FALSE)]
ind2<-unlist(strsplit( as.character(data$V1) , "-" ))[c(FALSE,TRUE) ]


m<-as.matrix(data)

dat<-as.data.frame(cbind(as.matrix(ind1),as.matrix(ind2), m[,seq(2,(length(data)-1))]) ,stringsAsFactors=FALSE)


cols = seq(3,length(dat)-2);
dat[,cols] = apply(dat[,cols], 2, function(x) as.numeric(x));

lengthindinclude<-length(strsplit( args[2]  , "," )[[1]])

vectorin<-rep(FALSE,length(dat$V1));
for(indinc in  strsplit( args[2]  , "," )[[1]]){
#for(indinc in  strsplit( "15302,15304,San_A,San_B" , "," )[[1]]){

    vectorin <- vectorin | dat$V1==indinc;
}
dat<-dat[vectorin  , ]

vectorout<-rep(TRUE,length(dat$V1));
for(indinc in  strsplit( args[3]  , "," )[[1]]){

    vectorout <- vectorout & dat$V2!=indinc;
}


dat<-dat[vectorout , ]

dat$V2 <- factor(dat$V2, levels = unique(as.vector(dat$V2[ order(as.numeric(dat$V1),dat$V65) ] )) ,ordered=TRUE)


ymax<-max(c(dat$V65,dat$V66,dat$V67));
ymin<-min(c(dat$V65,dat$V66,dat$V67));
wiggleroomyaxis<-0.001;

#pdf("divtest.pdf",width=400,height=400)
p<-ggplot(dat,aes(x=V1, y=V65, fill=V2))+ geom_bar(stat = "identity",position="dodge")+coord_cartesian(ylim = c(ymin-wiggleroomyaxis,ymax+wiggleroomyaxis))+ xlab("Ind#2")+ ylab("Div. to anc.")+guides(fill=guide_legend(title="Ind#1"))+geom_errorbar(aes(ymin=V66, ymax=V67), width=.2,position=position_dodge(.9))

ggsave(args[4],p, width=(lengthindinclude*7), height=7, units="in")

#dev.off();
