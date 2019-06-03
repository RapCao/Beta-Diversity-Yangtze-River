library(raster)
library(psych)
library(lme4)
library(ggplot2)

library(ape)

library(tcltk)
library(picante)
library(vegan)
library(fossil)







library(betapart) #BAS法
library(BAT) #POD法
##重命名BAT包中重名的函数
beta.pod.multi<-BAT::beta.multi
beta.pod.pair<-BAT::beta

setwd("C:/_Software/_wd/_WD")

comm <- read.csv(file="temp/comm1.csv")
comm <- as.matrix(comm[,(-1)])

##BAS法
#总体beta多样性指数
BAS_sor_multi <- beta.multi(x=comm,index.family="sorensen")

#成对相异性指数及其组分
BAS_sor_pair <- beta.pair(x=comm,index.family="sorensen")

BAS_sor <- BAS_sor_pair$beta.sor #BAS 法成对相异指数 Sorensen
BAS_sim <- BAS_sor_pair$beta.sim #BAS 法空间周转组分 βsim
BAS_sne <- BAS_sor_pair$beta.sne #BAS 法嵌套组分 βsne

##POD法
#总体beta多样性指数
POD_jac_multi <- beta.pod.multi(comm,func="jaccard")

#成对相异性指数及其组分
POD_jac_pair <- beta.pod.pair(comm,func="jaccard")

POD_Btotal <- POD_jac_pair$Btotal #POD 法成对相异指数 Jaccard
POD_Brepl <- POD_jac_pair$Brepl #POD 法物种替换组分 β-3
POD_Brich <- POD_jac_pair$Brich #POD 法丰富度差异组分 βrich


write.csv(as.matrix(BAS_sor),"beta_out/BAS_sor.csv")
write.csv(as.matrix(BAS_sim),"beta_out/BAS_sim.csv")
write.csv(as.matrix(BAS_sne),"beta_out/BAS_sne.csv")
write.csv(as.matrix(POD_Btotal),"beta_out/POD_Btotal.csv")
write.csv(as.matrix(POD_Brepl),"beta_out/POD_Brepl.csv")
write.csv(as.matrix(POD_Brich),"beta_out/POD_Brich.csv")

png(file="beta_out/BAS_sor.png",bg="transparent",width=1000,height=1000)
plot(BAS_sor)
dev.off()
png(file="beta_out/BAS_sim.png",bg="transparent",width=1000,height=1000)
plot(BAS_sim)
dev.off()
png(file="beta_out/BAS_sne.png",bg="transparent",width=1000,height=1000)
plot(BAS_sne)
dev.off()
png(file="beta_out/POD_Btotal.png",bg="transparent",width=1000,height=1000)
plot(POD_Btotal)
dev.off()
png(file="beta_out/POD_Brepl.png",bg="transparent",width=1000,height=1000)
plot(POD_Brepl)
dev.off()
png(file="beta_out/POD_Brich.png",bg="transparent",width=1000,height=1000)
plot(POD_Brich)
dev.off()

rm(list=ls())
gc()



###beta多样性密度曲线
library(ggplot2)
setwd("C:/_Software/_wd/_WD")

dat_bas_sor<-read.csv("beta_out/BAS_sor.csv")[,-1]
dat_bas_sim<-read.csv("beta_out/BAS_sim.csv")[,-1]
dat_bas_sne<-read.csv("beta_out/BAS_sne.csv")[,-1]
dat_pod_Btotal<-read.csv("beta_out/POD_Btotal.csv")[,-1]
dat_pod_Brepl<-read.csv("beta_out/POD_Brepl.csv")[,-1]
dat_pod_Brich<-read.csv("beta_out/POD_Brich.csv")[,-1]

ft_bas_sor<-c()
ft_bas_sim<-c()
ft_bas_sne<-c()
ft_pod_Btotal<-c()
ft_pod_Brepl<-c()
ft_pod_Brich<-c()
for (i in 1:(length(dat_bas_sor[1,])-1))##不算自身
{
	ft_bas_sor<-c(ft_bas_sor,dat_bas_sor[-1:-i,i])
	ft_bas_sim<-c(ft_bas_sim,dat_bas_sim[-1:-i,i])
	ft_bas_sne<-c(ft_bas_sne,dat_bas_sne[-1:-i,i])
	ft_pod_Btotal<-c(ft_pod_Btotal,dat_pod_Btotal[-1:-i,i])
	ft_pod_Brepl<-c(ft_pod_Brepl,dat_pod_Brepl[-1:-i,i])
	ft_pod_Brich<-c(ft_pod_Brich,dat_pod_Brich[-1:-i,i])
}

beta_all<-data.frame(ft_bas_sor,ft_bas_sim,ft_bas_sne,
	ft_pod_Btotal,ft_pod_Brepl,ft_pod_Brich)
colnames(beta_all)<-c("βsor","βsim","βsne","βcc","βrepl","βrich")

beta_all_1<-beta_all
beta_all_1$βsor[beta_all_1$βsor==1]<-NA
beta_all_1$βsim[beta_all_1$βsim==1]<-NA
beta_all_1$βsne[beta_all_1$βsne==0]<-NA
beta_all_1$βcc[beta_all_1$βcc==1]<-NA

d_bas_sor<-ggplot(beta_all,aes(x=βsor))+geom_density(colour="firebrick3")+geom_rug(aes(x=βsor),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_bas_sor,filename="beta_out/d_bas_sor.png")
d_bas_sim<-ggplot(beta_all,aes(x=βsim))+geom_density(colour="firebrick3")+geom_rug(aes(x=βsim),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_bas_sim,filename="beta_out/d_bas_sim.png")
d_bas_sne<-ggplot(beta_all,aes(x=βsne))+geom_density(colour="firebrick3")+geom_rug(aes(x=βsne),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_bas_sne,filename="beta_out/d_bas_sne.png")


d_pod_Btotal<-ggplot(beta_all,aes(x=βcc))+geom_density(colour="firebrick3")+geom_rug(aes(x=βcc),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_pod_Btotal,filename="beta_out/d_pod_Btotal.png")
d_pod_Brepl<-ggplot(beta_all,aes(x=βrepl))+geom_density(colour="firebrick3")+geom_rug(aes(x=βrepl),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_pod_Brepl,filename="beta_out/d_pod_Brepl.png")
d_pod_Brich<-ggplot(beta_all,aes(x=βrich))+geom_density(colour="firebrick3")+geom_rug(aes(x=βrich),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_pod_Brich,filename="beta_out/d_pod_Brich.png")


##放大部分
d_bas_sor_1<-ggplot(beta_all_1,aes(x=βsor))+geom_density(colour="firebrick3")+geom_rug(aes(x=βsor),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_bas_sor_1,filename="beta_out/d_bas_sor_1.png")
d_bas_sim_1<-ggplot(beta_all_1,aes(x=βsim))+geom_density(colour="firebrick3")+geom_rug(aes(x=βsim),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_bas_sim_1,filename="beta_out/d_bas_sim_1.png")
d_bas_sne_1<-ggplot(beta_all_1,aes(x=βsne))+geom_density(colour="firebrick3")+geom_rug(aes(x=βsne),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_bas_sne_1,filename="beta_out/d_bas_sne_1.png")

d_pod_Btotal_1<-ggplot(beta_all_1,aes(x=βcc))+geom_density(colour="firebrick3")+geom_rug(aes(x=βcc),colour="firebrick3",sides="b")+
geom_line(colour="firebrick3",stat="density")
ggsave(d_pod_Btotal_1,filename="beta_out/d_pod_Btotal_1.png")

rm(list=ls())
gc()



###comm热图
library(pheatmap)
setwd("C:/_Software/_wd/_WD")

comm <- read.csv(file="temp/comm1.csv")
comm <- as.matrix(comm[,(-1)])

png(filename="beta_out/comm.png",height=3500,width=5000,res=500,units="px")
pheatmap(comm,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	legend_breaks=c(0,1),legend_labels=c("无分布","有分布"),
	labels_row="样地编号 ID",
	fontsize=10,fontsize_col=2,
	color=colorRampPalette(c("lightyellow","orange"))(2))
dev.off()

rm(list=ls())
gc()



###beta多样性热图
library(pheatmap)
setwd("C:/_Software/_wd/_WD")

dat_bas_sor<-as.matrix(read.csv("beta_out/BAS_sor.csv"))[,-1]
dat_bas_sim<-as.matrix(read.csv("beta_out/BAS_sim.csv"))[,-1]
dat_bas_sne<-as.matrix(read.csv("beta_out/BAS_sne.csv"))[,-1]
dat_pod_Btotal<-as.matrix(read.csv("beta_out/POD_Btotal.csv"))[,-1]
dat_pod_Brepl<-as.matrix(read.csv("beta_out/POD_Brepl.csv"))[,-1]
dat_pod_Brich<-as.matrix(read.csv("beta_out/POD_Brich.csv"))[,-1]

pheatmap(dat_bas_sor,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	main="βsor",labels_col="",color=colorRampPalette(c("firebrick3","white"))(100),
	filename="beta_out/hm_bas_sor.png")
pheatmap(dat_bas_sim,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	main="βsim",labels_col="",color=colorRampPalette(c("firebrick3","white"))(100),
	filename="beta_out/hm_bas_sim.png")
pheatmap(dat_bas_sne,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	main="βsne",labels_col="",color=colorRampPalette(c("firebrick3","white"))(100),
	filename="beta_out/hm_bas_sne.png")
pheatmap(dat_pod_Btotal,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	main="βcc",labels_col="",color=colorRampPalette(c("firebrick3","white"))(100),
	filename="beta_out/hm_pod_Btotal.png")
pheatmap(dat_pod_Brepl,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	main="β-3",labels_col="",color=colorRampPalette(c("firebrick3","white"))(100),
	filename="beta_out/hm_pod_Brepl.png")
pheatmap(dat_pod_Brich,cluster_row=FALSE,cluster_col=FALSE,border_color=NA,
	main="βrich",labels_col="",color=colorRampPalette(c("firebrick3","white"))(100),
	filename="beta_out/hm_pod_Brich.png")

rm(list=ls())
gc()



###谱系beta多样性
library(picante)
library(ape)

setwd("C:/_Software/_wd/_WD")
comm<-read.csv(file="temp/comm1.csv")
#comm <- as.matrix(comm[,(-1)])
tree<-read.tree(file="tree/tree.txt")

​prune_tree<-prune.sample(comm,tree)
phydist<-cophenetic(prune_tree)

mpd<-ses.mpd(comm,phydist,null.model="taxa.labels",abundance.weighted=T,runs=99)
#计算MPD，得到的mpd.obs即为mpd值，obs.mpd.z为-NRI






####回归
###差异矩阵获取
setwd("C:/_Software/_wd/_WD")
dat_all<-read.csv("temp/alldata.csv")
dat_base<-read.csv("temp/database_output.csv")


##beta
beta1_sor_out<-as.matrix(read.csv("beta_out/BAS_sor.csv")[,-1])
beta1_Btotal_out<-as.matrix(read.csv("beta_out/POD_Btotal.csv")[,-1])

beta1_sim_out<-as.matrix(read.csv("beta_out/BAS_sim.csv")[,-1])
beta1_sne_out<-as.matrix(read.csv("beta_out/BAS_sne.csv")[,-1])

beta1_Brepl_out<-as.matrix(read.csv("beta_out/POD_Brepl.csv")[,-1])
beta1_Brich_out<-as.matrix(read.csv("beta_out/POD_Brich.csv")[,-1])


##生态系统功能
meangpp<-dat_all$meangpp
meannpp<-dat_all$meannpp
sdgpp<-dat_all$sdgpp
sdnpp<-dat_all$sdnpp
ES<-dat_all$ES
rich<-dat_base$richness

meangpp_out<-matrix(0,nrow=length(meangpp),ncol=length(meangpp))
for (i in 1:length(meangpp))
{
	row_out<-c()
	for (j in 1:length(meangpp))
	{
		res<-abs(meangpp[i]-meangpp[j])
		row_out<-c(row_out,res)
	}
	meangpp_out[i,]<-row_out
}
rownames(meangpp_out)<-dat_all$ID
colnames(meangpp_out)<-dat_all$ID
write.csv(meangpp_out,"regr/meangpp.csv")

meannpp_out<-matrix(0,nrow=length(meannpp),ncol=length(meannpp))
for (i in 1:length(meannpp))
{
	row_out<-c()
	for (j in 1:length(meannpp))
	{
		res<-abs(meannpp[i]-meannpp[j])
		row_out<-c(row_out,res)
	}
	meannpp_out[i,]<-row_out
}
rownames(meannpp_out)<-dat_all$ID
colnames(meannpp_out)<-dat_all$ID
write.csv(meannpp_out,"regr/meannpp.csv")

sdgpp_out<-matrix(0,nrow=length(sdgpp),ncol=length(sdgpp))
for (i in 1:length(sdgpp))
{
	row_out<-c()
	for (j in 1:length(sdgpp))
	{
		res<-abs(sdgpp[i]-sdgpp[j])
		row_out<-c(row_out,res)
	}
	sdgpp_out[i,]<-row_out
}
rownames(sdgpp_out)<-dat_all$ID
colnames(sdgpp_out)<-dat_all$ID
write.csv(sdgpp_out,"regr/sdgpp.csv")

sdnpp_out<-matrix(0,nrow=length(sdnpp),ncol=length(sdnpp))
for (i in 1:length(sdnpp))
{
	row_out<-c()
	for (j in 1:length(sdnpp))
	{
		res<-abs(sdnpp[i]-sdnpp[j])
		row_out<-c(row_out,res)
	}
	sdnpp_out[i,]<-row_out
}
rownames(sdnpp_out)<-dat_all$ID
colnames(sdnpp_out)<-dat_all$ID
write.csv(sdnpp_out,"regr/sdnpp.csv")

ES_out<-matrix(0,nrow=length(ES),ncol=length(ES))
for (i in 1:length(ES))
{
	row_out<-c()
	for (j in 1:length(ES))
	{
		res<-abs(ES[i]-ES[j])
		row_out<-c(row_out,res)
	}
	ES_out[i,]<-row_out
}
rownames(ES_out)<-dat_all$ID
colnames(ES_out)<-dat_all$ID
write.csv(ES_out,"regr/ES.csv")

rich_out<-matrix(0,nrow=length(rich),ncol=length(rich))
for (i in 1:length(rich))
{
	row_out<-c()
	for (j in 1:length(rich))
	{
		res<-abs(rich[i]-rich[j])
		row_out<-c(row_out,res)
	}
	rich_out[i,]<-row_out
}
rownames(rich_out)<-dat_all$ID
colnames(rich_out)<-dat_all$ID
write.csv(rich_out,"regr/rich.csv")



##气象
MAP<-dat_all$MAP
MAT<-dat_all$MAT
ATD<-dat_all$ATD

MAP_out<-matrix(0,nrow=length(MAP),ncol=length(MAP))
for (i in 1:length(MAP))
{
	row_out<-c()
	for (j in 1:length(MAP))
	{
		res<-abs(MAP[i]-MAP[j])
		row_out<-c(row_out,res)
	}
	MAP_out[i,]<-row_out
}
rownames(MAP_out)<-dat_all$ID
colnames(MAP_out)<-dat_all$ID
write.csv(MAP_out,"regr/MAP.csv")

MAT_out<-matrix(0,nrow=length(MAT),ncol=length(MAT))
for (i in 1:length(MAT))
{
	row_out<-c()
	for (j in 1:length(MAT))
	{
		res<-abs(MAT[i]-MAT[j])
		row_out<-c(row_out,res)
	}
	MAT_out[i,]<-row_out
}
rownames(MAT_out)<-dat_all$ID
colnames(MAT_out)<-dat_all$ID
write.csv(MAT_out,"regr/MAT.csv")

ATD_out<-matrix(0,nrow=length(ATD),ncol=length(ATD))
for (i in 1:length(ATD))
{
	row_out<-c()
	for (j in 1:length(ATD))
	{
		res<-abs(ATD[i]-ATD[j])
		row_out<-c(row_out,res)
	}
	ATD_out[i,]<-row_out
}
rownames(ATD_out)<-dat_all$ID
colnames(ATD_out)<-dat_all$ID
write.csv(ATD_out,"regr/ATD.csv")


##地理距离
library(geosphere)

lonlat<-cbind(dat_all$LONG,dat_all$LAT)
dists=distm(lonlat, fun=distVincentyEllipsoid)

rownames(dists)<-dat_all$ID
colnames(dists)<-dat_all$ID

write.csv(dists,"regr/dists.csv")



###回归
##差异矩阵转向量
colnames(meangpp_out)<-NULL
colnames(meannpp_out)<-NULL
colnames(sdgpp_out)<-NULL
colnames(sdnpp_out)<-NULL
colnames(ES_out)<-NULL
colnames(rich_out)<-NULL

colnames(MAP_out)<-NULL
colnames(MAT_out)<-NULL
colnames(ATD_out)<-NULL

colnames(dists)<-NULL


ar_meangpp<-c()
ar_meannpp<-c()
ar_sdgpp<-c()
ar_sdnpp<-c()
ar_ES<-c()
ar_rich<-c()

ar_MAP<-c()
ar_MAT<-c()
ar_ATD<-c()

ar_dists<-c()

ar_beta1_sor<-c()
ar_beta1_Btotal<-c()

ar_beta1_sim<-c()
ar_beta1_sne<-c()
ar_beta1_Brepl<-c()
ar_beta1_Brich<-c()

for (i in 1:(nrow(meangpp_out)-1))
{
	n<-length(meangpp_out[i,])

	ar_beta1_sor<-c(ar_beta1_sor,beta1_sor_out[i,(i+1):n])
	ar_beta1_Btotal<-c(ar_beta1_Btotal,beta1_Btotal_out[i,(i+1):n])

	ar_beta1_sim<-c(ar_beta1_sim,beta1_sim_out[i,(i+1):n])
	ar_beta1_sne<-c(ar_beta1_sne,beta1_sne_out[i,(i+1):n])

	ar_beta1_Brepl<-c(ar_beta1_Brepl,beta1_Brepl_out[i,(i+1):n])
	ar_beta1_Brich<-c(ar_beta1_Brich,beta1_Brich_out[i,(i+1):n])


	ar_meangpp<-c(ar_meangpp,meangpp_out[i,(i+1):n])
	ar_meannpp<-c(ar_meannpp,meannpp_out[i,(i+1):n])
	ar_sdgpp<-c(ar_sdgpp,sdgpp_out[i,(i+1):n])
	ar_sdnpp<-c(ar_sdnpp,sdnpp_out[i,(i+1):n])
	ar_ES<-c(ar_ES,ES_out[i,(i+1):n])
	ar_rich<-c(ar_rich,rich_out[i,(i+1):n])

	ar_MAP<-c(ar_MAP,MAP_out[i,(i+1):n])
	ar_MAT<-c(ar_MAT,MAT_out[i,(i+1):n])
	ar_ATD<-c(ar_ATD,ATD_out[i,(i+1):n])

	ar_dists<-c(ar_dists,dists[i,(i+1):n])
}

dat_dif<-data.frame(ar_beta1_sor,ar_beta1_Btotal,ar_beta1_sim,ar_beta1_sne,ar_beta1_Brepl,ar_beta1_Brich,
	ar_meangpp,ar_meannpp,ar_sdgpp,ar_sdnpp,ar_ES,ar_rich,ar_MAP,ar_MAT,ar_ATD,ar_dists)
colnames(dat_dif)<-c("beta1_sor","beta1_Btotal","beta1_sim","beta1_sne","beta1_Brepl","beta1_Brich",
	"meangpp","meannpp","sdgpp","sdnpp","ES","rich","MAP","MAT","ATD","dists")
write.csv(dat_dif,file="regr/dat_dif.csv",row.names=FALSE)

rm(list=ls())
gc()

dat_dif<-read.csv("regr/dat_dif.csv")


##最近八个样地均差
##获取邻近样地的位置和距离
setwd("C:/_Software/_wd/_WD")

dist_source<-read.csv("regr/dists.csv")
ID<-dist_source[,1]
dist_source<-dist_source[,-1]

dist_spot<-data.frame()
for (i in 1:length(dist_source))
{
	dist_mirror<-dist_source
	dist_mirror[i,i]<-NA
	temp_min<-c()
	for (j in 1:8)
	{
		temp_mindist<-min(dist_mirror[i,],na.rm=TRUE)
		temp_minnum<-which.min(dist_mirror[i,])
		temp_min<-c(temp_min,temp_minnum,temp_mindist)
		dist_mirror[i,temp_minnum]<-NA
		j<-j+1
	}
	dist_spot<-rbind(dist_spot,temp_min)
}

dist_spot<-cbind(ID,dist_spot)
colnames(dist_spot)<-c("ID","N1","D1","N2","D2","N3","D3","N4","D4","N5","D5","N6","D6","N7","D7","N8","D8")

write.csv(dist_spot,"spots/dist_spot.csv",row.names=FALSE)

rm(list=ls())
gc()


##临近样地的数据及数据差异
setwd("C:/_Software/_wd/_WD")

dat_all<-cbind(read.csv("temp/alldata.csv"),read.csv("temp/database_output.csv")$richness)
colnames(dat_all)<-c(colnames(dat_all)[-length(colnames(dat_all))],"rich")

dist_spot<-read.csv("spots/dist_spot.csv")
spots<-dist_spot[,seq(from=2,to=16,by=2)]
dist_r<-dist_spot[,seq(from=3,to=17,by=2)]

spots_pick_data <- function(index)
{
	dat_pie<-subset(dat_all,select=index)

	spots_mean<-c()
	spots_dif<-c()
	for (i in 1:nrow(spots))
	{
		num<-spots[i,]
		row_out<-c()
		for (j in num)
		{
			row_out<-c(row_out,dat_pie[j,])
		}
		row_mean<-mean(row_out,na.rm=TRUE)
		row_dif<-dat_pie[i,]-row_mean

		spots_mean<-c(spots_mean,row_mean)
		spots_dif<-c(spots_dif,row_dif)
	}

	res<-data.frame(spots_mean,spots_dif)
	colnames(res)<-c(paste0(index,"_mean"),paste0(index,"_dif"))
	res;
}


spots_pick_beta <- function(index)
{	
	dat_beta<-read.csv(paste0("beta_out/",index,".csv"))[,-1]

	spots_mean<-c()
	spots_dif<-c()
	for (i in 1:nrow(spots))
	{
		num<-spots[i,]
		row_out<-c()
		for (j in num)
		{
			row_out<-c(row_out,dat_beta[i,j])
		}
		row_mean<-mean(row_out,na.rm=TRUE)

		spots_mean<-c(spots_mean,row_mean)
	}

	for (i in 1:nrow(spots))
	{
		num<-spots[i,]
		row_out<-c()
		for (j in num)
		{
			row_out<-c(row_out,spots_mean[j])
		}
		row_dif<-spots_mean[i]-mean(row_out,na.rm=TRUE)

		spots_dif<-c(spots_dif,row_dif)
	}

	res<-data.frame(spots_mean,spots_dif)
	colnames(res)<-c(paste0(index,"_mean"),paste0(index,"_dif"))
	res;
}


dist_mean<-data.frame(rowMeans(dist_r))
colnames(dist_mean)<-c("dist_mean")

dat_spots<-dat_all[,1]
data_var<-c("meangpp","meannpp","sdgpp","sdnpp","ES","rich","MAP","MAT","ATD")
beta_var<-c("BAS_sor","BAS_sim","BAS_sne","POD_Btotal","POD_Brepl","POD_Brich")
for (p_data in data_var)
{
	dat_spots<-cbind(dat_spots,spots_pick_data(p_data))
}
dat_spots<-cbind(dat_spots,dist_mean)
for (p_beta in beta_var)
{
	dat_spots<-cbind(dat_spots,spots_pick_beta(p_beta))
}

write.csv(dat_spots,"spots/dat_spots.csv",row.names=FALSE)




##单项批量内回归
library(ggplot2)
library(GGally)

library(PerformanceAnalytics)

library(glmulti)
library(MuMIn)

setwd("C:/_Software/_wd/_WD")


dat_row1<-read.csv("temp/alldata.csv")
dat_row2<-read.csv("temp/database_output.csv")
dat_pat<-data.frame(dat_row1$meangpp,dat_row1$meannpp,dat_row1$sdgpp,dat_row1$sdnpp,
	dat_row1$ES,dat_row2$richness,dat_row1$MAP,dat_row1$MAT,dat_row1$ATD)
colnames(dat_pat)<-c("meangpp","meannpp","sdgpp","sdnpp","ES","rich","MAP","MAT","ATD")

for(i in 1:length(dat_pat))
{
	dat_pat[i]<-scale(dat_pat[i])
}

dat_dif<-read.csv("regr/dat_dif.csv")


for(i in 1:length(dat_dif))
{
	dat_dif[i]<-scale(dat_dif[i])
}


dat_spots<-read.csv("spots/dat_spots.csv")[,-1]
for(i in 1:length(dat_spots))
{
	dat_spots[i]<-scale(dat_spots[i])
}



GGally::ggpairs(dat_pat,lower = list(continuous = "smooth"))
GGally::ggpairs(dat_pat[,1:6],lower = list(continuous = "smooth"))
GGally::ggpairs(dat_pat[,7:9],lower = list(continuous = "smooth"))

ggpairs(dat_pat[,7:9],lower=list(continuous="smooth",mapping=aes(color="red")),
	diag=list(mapping=aes(color="blue")))



cor.sig <- function(test)
{
    res.cor = cor(test)
    res.sig = res.cor
    res.sig[abs(res.sig) > 0] = NA
    nx = dim(test)[2]
    for (i in 1:nx)
    {
        for (j in 1:nx)
        {
            res.cor1 = as.numeric(cor.test(test[, i], test[, j])$est)
            res.sig1 = as.numeric(cor.test(test[, i], test[, j])$p.value)
            if (res.sig1 <= 0.001){
                sig.mark = "***"
            }
            if (res.sig1 <= 0.01 & res.sig1 > 0.001) {
                sig.mark = "** "
            }
            if (res.sig1 <= 0.05 & res.sig1 > 0.01) {
                sig.mark = "*  "
            }
            if (res.sig1 > 0.05) {
                sig.mark = "   "
            }
            if (res.cor1 > 0) {
                res.sig[i, j] = paste(" ", as.character(round(res.cor1, 3)), sig.mark, sep = "")
            } else {
                res.sig[i, j] = paste(as.character(round(res.cor1, 3)), sig.mark, sep = "")
            }
        }
    }
    as.data.frame(res.sig)
}

cor.sig(dat_pat)

cor.sig(dat_spots[,c(seq(from=2,to=18,by=2),19)])

global.model<-lm(BAS_sor_mean~rich_dif+meannpp_dif+ES_dif+MAT_dif+dist_mean,data=dat_spots)
sor.model<-glmulti(global.model,level=1,crit="aicc")

summary(sor.model)
summary(sor.model)$icvalue



global.model<-lm(BAS_sor_mean~meangpp_dif+meannpp_dif+ES_dif+rich_dif+MAP_dif+MAT_dif+ATD_dif+dist_mean,data=dat_spots)
sor.model<-glmulti(global.model,level=1,crit="aicc")

summary(sor.model)
summary(sor.model)$icvalue

ave_model <- function(index_set,y,data_set)
{
	data_pie<-data_set[,c(index_set,as.character(substitute(y)))]
	npar=length(index_set)
	unit=c(1,0)
	parEst=rep(unit,each=2^(npar-1))
	
	for (i in 2:npar)
	{
		unit=c(i,0)
		parEst.tmp=rep(rep(unit,each=2^(npar-i)),2^(i-1))
		parEst=cbind(parEst,parEst.tmp)
	}

	parMat=cbind(parEst[,1:npar],1)
	dimnames(parMat)=list(1:(2^npar),c(index_set,as.character(substitute(y))))
	
	allModel=list()
	for (i in 1:(dim(parMat)[1]-1))
	{
		data_pie.tmp=data_pie[,parMat[i,]!=0]
		eval(parse(text=paste0("allModel[[i]]=lm(",as.character(substitute(y)),"~.,data=data_pie.tmp)")))
	}

	modelC=lm(y~1,data=data_pie)
	lm.ave<-model.avg(allModel,modelC)
	summary(lm.ave)
	summary(lm.ave)$sw
}



BAS_sor_mean<-dat_spots$BAS_sor_mean
BAS_sim_mean<-dat_spots$BAS_sim_mean
BAS_sne_mean<-dat_spots$BAS_sne_mean
POD_Btotal_mean<-dat_spots$POD_Btotal_mean
POD_Brepl_mean<-dat_spots$POD_Brepl_mean
POD_Brich_mean<-dat_spots$POD_Brich_mean


ave_model(c("meangpp_dif","meannpp_dif","ES_dif","rich_dif","MAP_dif","MAT_dif","ATD_dif","dist_mean"),BAS_sor_mean,dat_spots)




ave_model(colnames(dat_spots)[,c(seq(from=2,to=18,by=2),19)],BAS_sor_mean,dat_spots)


ave_model(c("MAP_dif","ATD_dif","dist_mean"),BAS_sor_mean,dat_spots)


#
'''lm1<-lm(BAS_sor_mean~MAP_dif+ATD_dif+dist_mean,data=dat_spots)
lm2<-lm(BAS_sor_mean~MAP_dif+ATD_dif,data=dat_spots)
lm3<-lm(BAS_sor_mean~MAP_dif+dist_mean,data=dat_spots)
lm4<-lm(BAS_sor_mean~ATD_dif+dist_mean,data=dat_spots)
lm5<-lm(BAS_sor_mean~MAP_dif,data=dat_spots)
lm6<-lm(BAS_sor_mean~ATD_dif,data=dat_spots)
lm7<-lm(BAS_sor_mean~dist_mean,data=dat_spots)
lm8<-lm(BAS_sor_mean~1,data=dat_spots)
lm.ave <- model.avg(lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8)
summary(lm.ave)'''





res.sor.model<-lm(BAS_sor_mean ~ 1 + MAT_dif + dist_mean,data=dat_spots)









global.model<-lm(beta1_sor ~ rich + ES + meangpp + sdnpp, data = dat_dif)
sor.model <- glmulti(global.model, level = 1, crit = "aicc")


#GGally::ggpairs(dat_dif)




##单项差异与beta回归
library(lme4)
library(ggplot2)
library(patchwork)

setwd("C:/_Software/_wd/_WD")
dat_dif<-read.csv("regr/dat_dif.csv")

meangpp<-dat_dif$meangpp
meannpp<-dat_dif$meannpp
sdgpp<-dat_dif$sdgpp
sdnpp<-dat_dif$sdnpp
ES<-dat_dif$ES
rich<-dat_dif$rich
MAP<-dat_dif$MAP
MAT<-dat_dif$MAT
ATD<-dat_dif$ATD
dists<-dat_dif$dists
beta1_sor<-dat_dif$beta1_sor
beta1_sim<-dat_dif$beta1_sim
beta1_sne<-dat_dif$beta1_sne
beta1_Btotal<-dat_dif$beta1_Btotal
beta1_Brepl<-dat_dif$beta1_Brepl
beta1_Brich<-dat_dif$beta1_Brich


lm_eqn <- function(m)
{
	l <- list(#a=format(coef(m)[1],digits=2),
		#b=format(abs(coef(m)[2]),digits=2),
    	r2=format(summary(m)$r.squared,digits=3),
    	p=format(summary(m)$coefficients[2,4],digits=4));

  	'''if (coef(m)[2] >= 0)
  	{
    	#eq <- paste0(as.character(as.expression(substitute(y==a+b%*%x,l))),"\n",as.character(as.expression(substitute(r^2~"="~r2~","~P~"="~p,l))))
    	eq <- as.character(as.expression(substitute(R^2~"="~r2~","~P~"="~p,l)))
 	}
 	else
 	{
 		#eq <- substitute(italic(y)==a-b%*%italic(x)*" "~~italic(r)^2~"="~r2~","~italic(P)~"="~p,l)
 		eq <- as.character(as.expression(substitute(R^2~"="~r2~","~P~"="~p,l)))
 	}'''

 	eq <- as.character(as.expression(substitute(R^2~"="~r2~","~P~"="~p,l)))
 	eq;             
}

p_plot <- function(x,y,y_name=substitute(y),data)
{
	p_value<-summary(lm(y~x,data))$coefficients[2,4]
	if(p_value>=0&&p_value<=0.001)
	{
		x1<-paste(substitute(x),"***",sep="")
	}
	else if(p_value>0.001&&p_value<=0.01)
	{
		x1<-paste(substitute(x),"**",sep="")
	}
	else if(p_value>0.01&&p_value<=0.05)
	{
		x1<-paste(substitute(x),"*",sep="")
	}
	else
	{
		x1<-substitute(x)
	}

	p<-ggplot()+geom_point(aes(x=x,y=y),color="pink",size=0.5)+
	geom_smooth(aes(x=x,y=y),method='lm',data=data)+
	xlab(x1)+
	ylab(y_name)+
	theme(axis.title.x =element_text(size=20), axis.title.y=element_text(size=20))+
	geom_text(aes(x=max(x),y=0.25,label=lm_eqn(lm(y~x,data))),parse=TRUE,size=4,hjust=1);
}

all_plot <- function(y)
{
	y_name<-substitute(y)

	p1<-p_plot(meangpp,y,y_name,dat_dif)
	p2<-p_plot(meannpp,y,y_name,dat_dif)
	p3<-p_plot(sdgpp,y,y_name,dat_dif)
	p4<-p_plot(sdnpp,y,y_name,dat_dif)
	p5<-p_plot(ES,y,y_name,dat_dif)
	p6<-p_plot(rich,y,y_name,dat_dif)
	p7<-p_plot(MAP,y,y_name,dat_dif)
	p8<-p_plot(MAT,y,y_name,dat_dif)
	p9<-p_plot(ATD,y,y_name,dat_dif)
	p10<-p_plot(dists,y,y_name,dat_dif)
	all_p<-p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+plot_layout(ncol=4,heights=c(1,1))

	ggsave(all_p,filename=paste0("regr_out/",substitute(y),"_pic1.png"),width=16,height=9)
}

all_plot(beta1_sor)
all_plot(beta1_sim)
all_plot(beta1_sne)

all_plot(beta1_Btotal)
all_plot(beta1_Brepl)
all_plot(beta1_Brich)



###四川重庆地图
library(maps)
library(mapdata)
library(maptools)
library(ggplot2)
library(ggrepel)

setwd("C:/_Software/_wd/_WD")

map_china<-rgdal::readOGR('map/bou2_4p.shp')
#iconv(map$NAME,from="GBK")##检索省份所在行
map_site<-rbind(map_china[206,],map_china[208,])##四川206，重庆208

#provname<-map_site@data$NAME
provname<-c("四川省","重庆市")
lon<-c(102.71039,107.87505-0.25)
lat<-c(30.61561-0.25,30.05834-0.25)
prov<-data.frame(provname,lon,lat)

map_site<-borders(map_site,colour="gray50",fill="white")

map<-NULL
map<-ggplot()+map_site+ylim(25,35)+geom_text(aes(x=prov$lon,y=prov$lat,label=prov$provname),size=7.5)

COOR <- read.csv("temp/COOR.csv")

map<-map+labs(x="经度 Longitude",y="纬度 Latitude")+
theme(axis.title.x =element_text(size=15),axis.title.y=element_text(size=15))+
geom_point(aes(x=COOR$LONG,y=COOR$LAT),color="darkorange")##坐标轴名+图名居中+样地标点

ggsave(map,filename = "map/Map.png",width = 8,height = 8)



###发育树
##获得发育树
setwd("C:/_Software/_wd/_WD")
datafrag<-read.csv("temp/database_fragment.csv")
species<-data.frame(unique(datafrag$SPECIES))
write.csv(species[,1],"temp/sp_list.csv",row.names=FALSE)


library(plantlist)

sp<-read.csv("temp/sp_list.csv",head=TRUE)
sp<-as.vector(sp$x)

trimsp<-c()
for (i in 1:length(sp))
{
	spl<-strsplit(sp[i],split=" ")

	isvar<-0
	j<-1
	for (j in 1:length(spl[[1]]))
	{
		if (spl[[1]][j]=="var.")
		{
			isvar<-1
			break
		}
	}

	if (is.na(spl[[1]][2])==FALSE)
	{
		if (isvar==1)
		{
			newsp<-paste(spl[[1]][1],spl[[1]][2],spl[[1]][j],spl[[1]][j+1])
		}
		else
		{
			newsp<-paste(spl[[1]][1],spl[[1]][2])
		}
		trimsp<-c(trimsp,newsp)
	}
}

sp2 <- TPL(trimsp)
res <- taxa.table(sp2)
writeLines(res,"temp/taxa_table.txt")

##将taxa_table.txt内的内容复制到http://phylodiversity.net/phylomatic/的taxa中
##网页中storedtree=zanne2014,clean=ture
##得到进化树保存至tree/tree.txt

##发育树绘图
##在http://itol.embl.de/中upload保存的树文件



###统计科属
setwd("C:/_Software/_wd/_WD")

sp<-read.csv("temp/taxa_table1.txt",head=FALSE)

sp_res1<-c()
sp_res2<-c()
for (i in 1:nrow(sp))
{
	res<-strsplit(as.character(sp[i,]),split='/')
	sp_res1<-c(sp_res1,res[[1]][1])
	sp_res2<-c(sp_res2,res[[1]][2])
}
unique(sp_res1)
unique(sp_res2)



###3D分布图
library(rgl)

setwd("C:/_Software/_wd/_WD")

dat_row1<-read.csv("temp/alldata.csv")
dat_row2<-read.csv("temp/database_output.csv")

dat<-data.frame(dat_row1$meangpp,dat_row1$meannpp,dat_row2$richness)
colnames(dat)<-c("meangpp","meannpp","richness")

attach(dat)
plot3d(meangpp,meannpp,richness,xlim=c(0,2500),ylim=c(0,1600),zlim=c(0,20),
	size=5,col="red")



###多元回归模型
##导入数据并标准化
library(glmulti)
library(MuMIn)

setwd("C:/_Software/_wd/_WD")


dat_row1<-read.csv("temp/alldata.csv")
dat_row2<-read.csv("temp/database_output.csv")
dat_pat<-data.frame(dat_row1$meangpp,dat_row1$meannpp,dat_row1$sdgpp,dat_row1$sdnpp,
	dat_row1$ES,dat_row2$richness,dat_row1$MAP,dat_row1$MAT,dat_row1$ATD)
colnames(dat_pat)<-c("meangpp","meannpp","sdgpp","sdnpp","ES","rich","MAP","MAT","ATD")

for(i in 1:length(dat_pat))
{
	dat_pat[i]<-scale(dat_pat[i])
}


dat_dif<-read.csv("regr/dat_dif.csv")
for(i in 1:length(dat_dif))
{
	dat_dif[i]<-scale(dat_dif[i])
}


dat_spots<-read.csv("spots/dat_spots.csv")[,-1]
for(i in 1:length(dat_spots))
{
	dat_spots[i]<-scale(dat_spots[i])
}


##内回归确定因变量
cor.sig <- function(test)
{
    res.cor = cor(test)
    res.sig = res.cor
    res.sig[abs(res.sig) > 0] = NA
    nx = dim(test)[2]
    for (i in 1:nx)
    {
        for (j in 1:nx)
        {
            res.cor1 = as.numeric(cor.test(test[, i], test[, j])$est)
            res.sig1 = as.numeric(cor.test(test[, i], test[, j])$p.value)
            if (res.sig1 <= 0.001){
                sig.mark = "***"
            }
            if (res.sig1 <= 0.01 & res.sig1 > 0.001) {
                sig.mark = "** "
            }
            if (res.sig1 <= 0.05 & res.sig1 > 0.01) {
                sig.mark = "*  "
            }
            if (res.sig1 > 0.05) {
                sig.mark = "   "
            }
            if (res.cor1 > 0) {
                res.sig[i, j] = paste(" ", as.character(round(res.cor1, 3)), sig.mark, sep = "")
            } else {
                res.sig[i, j] = paste(as.character(round(res.cor1, 3)), sig.mark, sep = "")
            }
        }
    }
    as.data.frame(res.sig)
}




cor.sig(dat_pat)

cor.sig(dat_spots[,c(seq(from=2,to=18,by=2),19)])
cor.sig(dat_spots[,c(seq(from=1,to=17,by=2),19)])



##建立总体模型，并用AICc挑选
global.model<-lm(BAS_sor_mean~rich_dif+meannpp_dif+ES_dif+MAT_dif+dist_mean,data=dat_spots)
sor.model<-glmulti(global.model,level=1,crit="aicc")

summary(sor.model)
summary(sor.model)$icvalue



global.model<-lm(BAS_sor_mean~meangpp_dif+meannpp_dif+ES_dif+rich_dif+MAP_dif+MAT_dif+ATD_dif+dist_mean,data=dat_spots)
sor.model<-glmulti(global.model,level=1,crit="aicc")

summary(sor.model)
summary(sor.model)$icvalue


##模型平均
ave_model <- function(index_set,y,data_set)
{
	data_pie<-data_set[,c(index_set,as.character(substitute(y)))]
	npar=length(index_set)
	unit=c(1,0)
	parEst=rep(unit,each=2^(npar-1))
	
	for (i in 2:npar)
	{
		unit=c(i,0)
		parEst.tmp=rep(rep(unit,each=2^(npar-i)),2^(i-1))
		parEst=cbind(parEst,parEst.tmp)
	}

	parMat=cbind(parEst[,1:npar],1)
	dimnames(parMat)=list(1:(2^npar),c(index_set,as.character(substitute(y))))
	
	allModel=list()
	for (i in 1:(dim(parMat)[1]-1))
	{
		data_pie.tmp=data_pie[,parMat[i,]!=0]
		eval(parse(text=paste0("allModel[[i]]=lm(",as.character(substitute(y)),"~.,data=data_pie.tmp)")))
	}

	modelC=lm(y~1,data=data_pie)
	lm.ave<-model.avg(allModel,modelC)
	summary(lm.ave)
	summary(lm.ave)$sw
}



BAS_sor_mean<-dat_spots$BAS_sor_mean
BAS_sim_mean<-dat_spots$BAS_sim_mean
BAS_sne_mean<-dat_spots$BAS_sne_mean
POD_Btotal_mean<-dat_spots$POD_Btotal_mean
POD_Brepl_mean<-dat_spots$POD_Brepl_mean
POD_Brich_mean<-dat_spots$POD_Brich_mean


ave_model(c("meangpp_dif","meannpp_dif","ES_dif","rich_dif","MAP_dif","MAT_dif","ATD_dif","dist_mean"),BAS_sor_mean,dat_spots)




ave_model(colnames(dat_spots)[,c(seq(from=2,to=18,by=2),19)],BAS_sor_mean,dat_spots)


ave_model(c("MAP_dif","ATD_dif","dist_mean"),BAS_sor_mean,dat_spots)


#
'''lm1<-lm(BAS_sor_mean~MAP_dif+ATD_dif+dist_mean,data=dat_spots)
lm2<-lm(BAS_sor_mean~MAP_dif+ATD_dif,data=dat_spots)
lm3<-lm(BAS_sor_mean~MAP_dif+dist_mean,data=dat_spots)
lm4<-lm(BAS_sor_mean~ATD_dif+dist_mean,data=dat_spots)
lm5<-lm(BAS_sor_mean~MAP_dif,data=dat_spots)
lm6<-lm(BAS_sor_mean~ATD_dif,data=dat_spots)
lm7<-lm(BAS_sor_mean~dist_mean,data=dat_spots)
lm8<-lm(BAS_sor_mean~1,data=dat_spots)
lm.ave <- model.avg(lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8)
summary(lm.ave)'''





res.sor.model<-lm(BAS_sor_mean ~ 1 + MAT_dif + dist_mean,data=dat_spots)









global.model<-lm(beta1_sor ~ rich + ES + meangpp + sdnpp, data = dat_dif)
sor.model <- glmulti(global.model, level = 1, crit = "aicc")