args <- commandArgs(TRUE)

#open File

directory<-args[1];
print(directory);
filename<-args[2];
print(filename);
start_stats<-args[3];
end_stats<-args[4];
start_param<-args[5];
end_param<-args[6];
numComp<-args[7];
numComp<-as.numeric(numComp);
#print(numComp);
#numComp<-;
print(numComp);


#read file
a<-read.table(paste(directory, filename, sep=""), header=T, nrows=10000, skip=0);
#print (a);
print(names(a));
stats<-a[,c(start_stats:end_stats)]; params<-a[,c(start_param:end_param)]; rm (a);
#stats<-a[,c(12:347)]; params<-a[,c(1:10)]; rm(a);
#stats<-a[,c(25:120)]; params<-a[,c(1:24)]; rm(a);
 
#standardize the params
for(i in 1:length(params)){params[,i]<-(params[,i]-mean(params[,i]))/sd(params[,i]);}

#force stats in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stats)){
	myMax<-c(myMax, max(stats[,i]));
	myMin<-c(myMin, min(stats[,i]));
	stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
}

#transform statistics via boxcox  
library("MASS");	
for(i in 1:length(stats)){		
	d<-cbind(stats[,i], params);
	mylm<-lm(as.formula(d), data=d)			
	myboxcox<-boxcox(mylm, lambda=seq(-50, 80, 1/10), plotit=F, interp=T, eps=1/50);	
	lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);			
	#print(paste(names(stats)[i], myboxcox$x[myboxcox$y==max(myboxcox$y)]));
	myGM<-c(myGM, exp(mean(log(stats[,i]))));			
}

#standardize the BC-stats
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stats)){
	stats[,i]<-(stats[,i]^lambda[i] - 1)/(lambda[i]*myGM[i]^(lambda[i]-1));	
	myBCSDs<-c(myBCSDs, sd(stats[,i]));
	myBCMeans<-c(myBCMeans, mean(stats[,i]));		
	stats[,i]<-(stats[,i]-myBCMeans[i])/myBCSDs[i];
}

#perform pls
library("pls");
#myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp, validation="LOO");
myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp);

#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
for(i in 2:numComp) { myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]); } 
write.table(cbind(names(stats), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs, myPlsrDataFrame), file=paste(directory, "Routput_", filename, sep=""), col.names=F, row.names=F, sep="\t", quote=F);

#make RMSE plot
#pdf(paste(directory, "RMSE_", filename, ".pdf", sep=""),);
pdf(paste(directory, "RMSE_", filename, ".pdf", sep=""), width=12, height=12);
#par(mfrow=c(4,4));
plot(RMSEP(myPlsr));
dev.off();
 


#obsa<-read.table("/mnt/uni/ABC/arvalis/arvalis_both.obs", header=T);
#n<-data.frame(a=1:length(names(obsa)), n=names(obsa));
#pdf(paste(directory, "stats_", filename, ".pdf", sep=""), width=9, height=12);
#par(mfrow=c(5,4), cex=0.5)		
#	for(i in c(1:13,25,26,49:51,63,64,76:80,183:227)){
#	plot(density(stats[,i]), xlim=c(min(stats[,i])-max(stats[,i])+min(stats[,i]),max(stats[,i])+max(stats[,i])-min(stats[,i])), main=names(stats)[i]);	
#	print(paste(n[n[,2]==names(stats)[i],1], obsa[n[n[,2]==names(stats)[i],1]]));
#	lines(c(obsa[,n[n[,2]==names(stats)[i],1]], obsa[,n[n[,2]==names(stats)[i],1]]), c(0,1000), col="red")
#}

#dev.off();

