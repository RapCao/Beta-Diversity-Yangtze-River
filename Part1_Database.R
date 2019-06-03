setwd("C:\\_Software\\_wd\\_WD")

library(raster)
library(psych)
library(lme4)
library(ggplot2)

library(betapart) #BAS法
#library(BAT) #POD法
library(ape)

library(tcltk)
library(picante)
library(vegan)
library(fossil)



#得到位于四川、重庆的样地的坐标
dat_book2 <- read.csv("src\\Book2.csv")

COOR_sichuan <- subset(dat_book2,name=="四川")
COOR_chongqing <- subset(dat_book2,name=="重庆市")
COOR <- rbind(COOR_sichuan,COOR_chongqing)

dat_cp <- read.csv("src\\Cproject_res2.csv")
cp <- subset(dat_cp,select=ID:y)

COOR <- merge(cp,COOR,all=FALSE)
COOR <- COOR[order(COOR$ID,decreasing=F),]
COOR <- COOR[,c(3,1,2,4)]
colnames(COOR)=c("ID","LONG","LAT","name")

dir.create("temp")
write.csv(COOR,file="temp\\COOR.csv",row.names=FALSE)

rm(list=ls())
gc()


#根据坐标得出样地的生产力、生产力稳定性、降水、温度的数据库


COOR <- read.csv(file="temp\\COOR.csv")
COOR_xy <- COOR[,-4][,-1]

gpp2000<-raster("src\\tif\\MOD17A3_Science_GPP_2000.tif")
gpp2001<-raster("src\\tif\\MOD17A3_Science_GPP_2001.tif")
gpp2002<-raster("src\\tif\\MOD17A3_Science_GPP_2002.tif")
gpp2003<-raster("src\\tif\\MOD17A3_Science_GPP_2003.tif")
gpp2004<-raster("src\\tif\\MOD17A3_Science_GPP_2004.tif")
gpp2005<-raster("src\\tif\\MOD17A3_Science_GPP_2005.tif")
gpp2006<-raster("src\\tif\\MOD17A3_Science_GPP_2006.tif")
gpp2007<-raster("src\\tif\\MOD17A3_Science_GPP_2007.tif")
gpp2008<-raster("src\\tif\\MOD17A3_Science_GPP_2008.tif")
gpp2009<-raster("src\\tif\\MOD17A3_Science_GPP_2009.tif")
gpp2010<-raster("src\\tif\\MOD17A3_Science_GPP_2010.tif")
gpp2011<-raster("src\\tif\\MOD17A3_Science_GPP_2011.tif")
gpp2012<-raster("src\\tif\\MOD17A3_Science_GPP_2012.tif")
gpp2013<-raster("src\\tif\\MOD17A3_Science_GPP_2013.tif")
gpp2014<-raster("src\\tif\\MOD17A3_Science_GPP_2014.tif")
gpp2015<-raster("src\\tif\\MOD17A3_Science_GPP_2015.tif")

npp2000<-raster("src\\tif\\MOD17A3_Science_NPP_2000.tif")
npp2001<-raster("src\\tif\\MOD17A3_Science_NPP_2001.tif")
npp2002<-raster("src\\tif\\MOD17A3_Science_NPP_2002.tif")
npp2003<-raster("src\\tif\\MOD17A3_Science_NPP_2003.tif")
npp2004<-raster("src\\tif\\MOD17A3_Science_NPP_2004.tif")
npp2005<-raster("src\\tif\\MOD17A3_Science_NPP_2005.tif")
npp2006<-raster("src\\tif\\MOD17A3_Science_NPP_2006.tif")
npp2007<-raster("src\\tif\\MOD17A3_Science_NPP_2007.tif")
npp2008<-raster("src\\tif\\MOD17A3_Science_NPP_2008.tif")
npp2009<-raster("src\\tif\\MOD17A3_Science_NPP_2009.tif")
npp2010<-raster("src\\tif\\MOD17A3_Science_NPP_2010.tif")
npp2011<-raster("src\\tif\\MOD17A3_Science_NPP_2011.tif")
npp2012<-raster("src\\tif\\MOD17A3_Science_NPP_2012.tif")
npp2013<-raster("src\\tif\\MOD17A3_Science_NPP_2013.tif")
npp2014<-raster("src\\tif\\MOD17A3_Science_NPP_2014.tif")
npp2015<-raster("src\\tif\\MOD17A3_Science_NPP_2015.tif")

resgpp2000<-0.1*extract(gpp2000, COOR_xy)
resgpp2001<-0.1*extract(gpp2001, COOR_xy)
resgpp2002<-0.1*extract(gpp2002, COOR_xy)
resgpp2003<-0.1*extract(gpp2003, COOR_xy)
resgpp2004<-0.1*extract(gpp2004, COOR_xy)
resgpp2005<-0.1*extract(gpp2005, COOR_xy)
resgpp2006<-0.1*extract(gpp2006, COOR_xy)
resgpp2007<-0.1*extract(gpp2007, COOR_xy)
resgpp2008<-0.1*extract(gpp2008, COOR_xy)
resgpp2009<-0.1*extract(gpp2009, COOR_xy)
resgpp2010<-0.1*extract(gpp2010, COOR_xy)
resgpp2011<-0.1*extract(gpp2011, COOR_xy)
resgpp2012<-0.1*extract(gpp2012, COOR_xy)
resgpp2013<-0.1*extract(gpp2013, COOR_xy)
resgpp2014<-0.1*extract(gpp2014, COOR_xy)
resgpp2015<-0.1*extract(gpp2015, COOR_xy)

resnpp2000<-0.1*extract(npp2000, COOR_xy)
resnpp2001<-0.1*extract(npp2001, COOR_xy)
resnpp2002<-0.1*extract(npp2002, COOR_xy)
resnpp2003<-0.1*extract(npp2003, COOR_xy)
resnpp2004<-0.1*extract(npp2004, COOR_xy)
resnpp2005<-0.1*extract(npp2005, COOR_xy)
resnpp2006<-0.1*extract(npp2006, COOR_xy)
resnpp2007<-0.1*extract(npp2007, COOR_xy)
resnpp2008<-0.1*extract(npp2008, COOR_xy)
resnpp2009<-0.1*extract(npp2009, COOR_xy)
resnpp2010<-0.1*extract(npp2010, COOR_xy)
resnpp2011<-0.1*extract(npp2011, COOR_xy)
resnpp2012<-0.1*extract(npp2012, COOR_xy)
resnpp2013<-0.1*extract(npp2013, COOR_xy)
resnpp2014<-0.1*extract(npp2014, COOR_xy)
resnpp2015<-0.1*extract(npp2015, COOR_xy)

alldata <- data.frame(
	COOR$ID,COOR_xy$LONG,COOR_xy$LAT,

	resgpp2000,resgpp2001,resgpp2002,resgpp2003,
	resgpp2004,resgpp2005,resgpp2006,resgpp2007,
	resgpp2008,resgpp2009,resgpp2010,resgpp2011,
	resgpp2012,resgpp2013,resgpp2014,resgpp2015,

	resnpp2000,resnpp2001,resnpp2002,resnpp2003,
	resnpp2004,resnpp2005,resnpp2006,resnpp2007,
	resnpp2008,resnpp2009,resnpp2010,resnpp2011,
	resnpp2012,resnpp2013,resnpp2014,resnpp2015)

alldata[alldata==6553.5]=NA

alldata$meangpp <- rowMeans(subset(alldata,select=resgpp2000:resgpp2015))

alldata$meannpp <- rowMeans(subset(alldata,select=resnpp2000:resnpp2015))

alldata$sdgpp <- (((alldata$resgpp2000-alldata$meangpp)^2+(alldata$resgpp2001-alldata$meangpp)^2+
	(alldata$resgpp2002-alldata$meangpp)^2+(alldata$resgpp2003-alldata$meangpp)^2+
	(alldata$resgpp2004-alldata$meangpp)^2+(alldata$resgpp2005-alldata$meangpp)^2+
	(alldata$resgpp2006-alldata$meangpp)^2+(alldata$resgpp2007-alldata$meangpp)^2+
	(alldata$resgpp2008-alldata$meangpp)^2+(alldata$resgpp2009-alldata$meangpp)^2+
	(alldata$resgpp2010-alldata$meangpp)^2+(alldata$resgpp2011-alldata$meangpp)^2+
	(alldata$resgpp2012-alldata$meangpp)^2+(alldata$resgpp2013-alldata$meangpp)^2+
	(alldata$resgpp2014-alldata$meangpp)^2+(alldata$resgpp2015-alldata$meangpp)^2)/15)^0.5

alldata$sdnpp <- (((alldata$resnpp2000-alldata$meannpp)^2+(alldata$resnpp2001-alldata$meannpp)^2+
	(alldata$resnpp2002-alldata$meannpp)^2+(alldata$resnpp2003-alldata$meannpp)^2+
	(alldata$resnpp2004-alldata$meannpp)^2+(alldata$resnpp2005-alldata$meannpp)^2+
	(alldata$resnpp2006-alldata$meannpp)^2+(alldata$resnpp2007-alldata$meannpp)^2+
	(alldata$resnpp2008-alldata$meannpp)^2+(alldata$resnpp2009-alldata$meannpp)^2+
	(alldata$resnpp2010-alldata$meannpp)^2+(alldata$resnpp2011-alldata$meannpp)^2+
	(alldata$resnpp2012-alldata$meannpp)^2+(alldata$resnpp2013-alldata$meannpp)^2+
	(alldata$resnpp2014-alldata$meannpp)^2+(alldata$resnpp2015-alldata$meannpp)^2)/15)^0.5

alldata$ES <- (alldata$meannpp/alldata$sdnpp)

g<-ls()
rm(list=g[4:length(g)])
gc()

prec_01 <- raster("src\\tif\\wc2.0_30s_prec_01.tif")
prec_02 <- raster("src\\tif\\wc2.0_30s_prec_02.tif")
prec_03 <- raster("src\\tif\\wc2.0_30s_prec_03.tif")
prec_04 <- raster("src\\tif\\wc2.0_30s_prec_04.tif")
prec_05 <- raster("src\\tif\\wc2.0_30s_prec_05.tif")
prec_06 <- raster("src\\tif\\wc2.0_30s_prec_06.tif")
prec_07 <- raster("src\\tif\\wc2.0_30s_prec_07.tif")
prec_08 <- raster("src\\tif\\wc2.0_30s_prec_08.tif")
prec_09 <- raster("src\\tif\\wc2.0_30s_prec_09.tif")
prec_10 <- raster("src\\tif\\wc2.0_30s_prec_10.tif")
prec_11 <- raster("src\\tif\\wc2.0_30s_prec_11.tif")
prec_12 <- raster("src\\tif\\wc2.0_30s_prec_12.tif")

tavg_01 <- raster("src\\tif\\wc2.0_30s_tavg_01.tif")
tavg_02 <- raster("src\\tif\\wc2.0_30s_tavg_02.tif")
tavg_03 <- raster("src\\tif\\wc2.0_30s_tavg_03.tif")
tavg_04 <- raster("src\\tif\\wc2.0_30s_tavg_04.tif")
tavg_05 <- raster("src\\tif\\wc2.0_30s_tavg_05.tif")
tavg_06 <- raster("src\\tif\\wc2.0_30s_tavg_06.tif")
tavg_07 <- raster("src\\tif\\wc2.0_30s_tavg_07.tif")
tavg_08 <- raster("src\\tif\\wc2.0_30s_tavg_08.tif")
tavg_09 <- raster("src\\tif\\wc2.0_30s_tavg_09.tif")
tavg_10 <- raster("src\\tif\\wc2.0_30s_tavg_10.tif")
tavg_11 <- raster("src\\tif\\wc2.0_30s_tavg_11.tif")
tavg_12 <- raster("src\\tif\\wc2.0_30s_tavg_12.tif")

resprec_01<-extract(prec_01, COOR_xy)
resprec_02<-extract(prec_02, COOR_xy)
resprec_03<-extract(prec_03, COOR_xy)
resprec_04<-extract(prec_04, COOR_xy)
resprec_05<-extract(prec_05, COOR_xy)
resprec_06<-extract(prec_06, COOR_xy)
resprec_07<-extract(prec_07, COOR_xy)
resprec_08<-extract(prec_08, COOR_xy)
resprec_09<-extract(prec_09, COOR_xy)
resprec_10<-extract(prec_10, COOR_xy)
resprec_11<-extract(prec_11, COOR_xy)
resprec_12<-extract(prec_12, COOR_xy)

restavg_01<-extract(tavg_01, COOR_xy)
restavg_02<-extract(tavg_02, COOR_xy)
restavg_03<-extract(tavg_03, COOR_xy)
restavg_04<-extract(tavg_04, COOR_xy)
restavg_05<-extract(tavg_05, COOR_xy)
restavg_06<-extract(tavg_06, COOR_xy)
restavg_07<-extract(tavg_07, COOR_xy)
restavg_08<-extract(tavg_08, COOR_xy)
restavg_09<-extract(tavg_09, COOR_xy)
restavg_10<-extract(tavg_10, COOR_xy)
restavg_11<-extract(tavg_11, COOR_xy)
restavg_12<-extract(tavg_12, COOR_xy)

alldata <- cbind(alldata,

	resprec_01,resprec_02,resprec_03,resprec_04,
	resprec_05,resprec_06,resprec_07,resprec_08,
	resprec_09,resprec_10,resprec_11,resprec_12,

	restavg_01,restavg_02,restavg_03,restavg_04,
	restavg_05,restavg_06,restavg_07,restavg_08,
	restavg_09,restavg_10,restavg_11,restavg_12)

alldata$MAP <- alldata$resprec_01+alldata$resprec_02+alldata$resprec_03+alldata$resprec_04+
	alldata$resprec_05+alldata$resprec_06+alldata$resprec_07+alldata$resprec_08+
	alldata$resprec_09+alldata$resprec_10+alldata$resprec_11+alldata$resprec_12

alldata$MAT <- rowMeans(subset(alldata,select=restavg_01:restavg_12))

ATD <- c()
tavg_temp <- subset(alldata,select=restavg_01:restavg_12)
i=1
while(i<=nrow(alldata))
{
	ATD <- c(ATD,max(tavg_temp[i,])-min(tavg_temp[i,]))
	i=i+1
}
alldata$ATD <- ATD

colnames(alldata)=c("ID","LONG","LAT",

	"gpp2000","gpp2001","gpp2002","gpp2003",
	"gpp2004","gpp2005","gpp2006","gpp2007",
	"gpp2008","gpp2009","gpp2010","gpp2011",
	"gpp2012","gpp2013","gpp2014","gpp2015",

	"npp2000","npp2001","npp2002","npp2003",
	"npp2004","npp2005","npp2006","npp2007",
	"npp2008","npp2009","npp2010","npp2011",
	"npp2012","npp2013","npp2014","npp2015",

	"meangpp","meannpp","sdgpp","sdnpp","ES",

	"prec_01","prec_02","prec_03","prec_04",
	"prec_05","prec_06","prec_07","prec_08",
	"prec_09","prec_10","prec_11","prec_12",

	"tavg_01","tavg_02","tavg_03","tavg_04",
	"tavg_05","tavg_06","tavg_07","tavg_08",
	"tavg_09","tavg_10","tavg_11","tavg_12",

	"MAP","MAT","ATD")

write.csv(alldata,file="temp\\alldata.csv",row.names=FALSE)

rm(list=ls())
gc()



#根据ID得到样地的植物相关统计数据


dat_db <- read.csv("src\\CaoXZ_trees_database1.csv")
dat_all <- read.csv("temp\\alldata.csv")

IDset <- unique(dat_all$ID)

op_ID <- c()
op_numbers <- c()
op_richness <- c()
op_DBH_qua_mean <- c()
op_DBH_arith_mean <- c()
op_H_geo_mean <- c()
op_H_arith_mean <- c()
op_DBH_max <- c()
op_H_max <- c()
op_sabh <- c()
op_H_mean_5th <- c()

for(id_p in IDset)
{
	unit <- subset(dat_db,ID==id_p,select=ID:SPECIES)

	numbers <- nrow(unit)
	richness <- length(unique(unit$SPECIES))

	DBH_qua_mean <- (sum(unit$DBH^2,na.rm=TRUE)/nrow(unit[complete.cases(unit$DBH),]))^0.5
	DBH_arith_mean <- mean(unit$DBH,na.rm=TRUE)
	H_geo_mean <- geometric.mean(unit$H,na.rm=TRUE)
	H_arith_mean <- mean(unit$H,na.rm=TRUE)

	DBH_max <- max(unit$DBH,na.rm=TRUE)
	H_max <- max(unit$H,na.rm=TRUE)

	sabh <- sum(unit$DBH^2*pi/4,na.rm=TRUE)

	H_mean_5th <- mean(unit$H[which(unit$H>=H_max*0.8)])

	op_ID <- c(op_ID,id_p)
	op_numbers <- c(op_numbers,numbers)
	op_richness <- c(op_richness,richness)
	op_DBH_qua_mean <- c(op_DBH_qua_mean,DBH_qua_mean)
	op_DBH_arith_mean <- c(op_DBH_arith_mean,DBH_arith_mean)
	op_H_geo_mean <- c(op_H_geo_mean,H_geo_mean)
	op_H_arith_mean <- c(op_H_arith_mean,H_arith_mean)
	op_DBH_max <- c(op_DBH_max,DBH_max)
	op_H_max <- c(op_H_max,H_max)
	op_sabh <- c(op_sabh,sabh)
	op_H_mean_5th <- c(op_H_mean_5th,H_mean_5th)
}

output <- data.frame(op_ID,op_numbers,op_richness,
	op_DBH_qua_mean,op_DBH_arith_mean,op_H_geo_mean,op_H_arith_mean,
	op_DBH_max,op_H_max,
	op_sabh,
	op_H_mean_5th)
colnames(output)=c("ID","numbers","richness",
	"DBH_qua_mean","DBH_arith_mean","H_geo_mean","H_arith_mean",
	"DBH_max","H_max",
	"sabh",
	"H_mean_5th")

write.csv(output,file="temp\\database_output.csv",row.names=FALSE)

rm(list=ls())
gc()



#数据回归初探


csv_all <- read.csv("temp\\alldata.csv")
csv_op <- read.csv("temp\\database_output.csv")

Regression <- function(csv1=data.frame(default1=0),csv2=data.frame(default2=0),col1="default1",col2="default2",iflog=FALSE)
{
	if(iflog==TRUE)
	{
		dat <- cbind(log(subset(csv1,select=col1)),log(subset(csv2,select=col2)))
	}
	else
	{
		dat <- cbind(subset(csv1,select=col1),subset(csv2,select=col2))
	}

	assign(paste(col1), unlist(dat[1]))
	assign(paste(col2), unlist(dat[2]))

	#print(ES)
	#col1 <- unlist(dat[1])
	#col2 <- unlist(dat[2])
	
	lm_res <- lm(formula=(eval(parse(text=col1)))~(eval(parse(text=col2))))
	print(summary(lm_res))
	points(eval(parse(text=col1)),predict(lm_res,type="response"))

	return(lm_res)
}

#plot(Regression(csv_all,csv_op,"ES","richness",FALSE))
#Regression(csv_all,csv_op,"ES","richness",TRUE)

#Regression(csv_all,csv_op,"ES","richness",FALSE)
Regression(csv_all,csv_op,"meannpp","richness",FALSE)

Regression(csv_all,csv_all,"ES","MAT",FALSE)
Regression(csv_all,csv_all,"meannpp","MAT",FALSE)

Regression(csv_all,csv_all,"ES","MAP",FALSE)
Regression(csv_all,csv_all,"meannpp","MAP",FALSE)

Regression(csv_all,csv_all,"ES","ATD",FALSE)
Regression(csv_all,csv_all,"meannpp","ATD",FALSE)

rm(list=ls())
gc()



#根据ID得到群落物种多样性comm数据
dat_db <- read.csv("src\\CaoXZ_trees_database1.csv")
dat_all <- read.csv("temp\\alldata.csv")

IDset <- unique(dat_all$ID)
IDset_df <- data.frame(IDset)

fragment <- merge(IDset_df,dat_db,by.x="IDset",by.y="ID")
fragment <- subset(fragment,SPECIES!="#N/A")

SP <- unique(fragment$SPECIES)
SP <- SP[order(SP,decreasing=F)]

SP_underlined <- c()
for(sp_word in SP)
{
	word <- strsplit(sp_word," ")
	new_word <- paste(word[[1]],collapse="_")
	SP_underlined <- cbind(SP_underlined,new_word)
}

colnames(fragment)=c("ID","DBH","H","SPECIES")
write.csv(fragment,file="temp\\database_fragment.csv",row.names=FALSE)


comm1 <- matrix(0,nrow=length(IDset), ncol=length(SP))
rownames(comm1) <- IDset
colnames(comm1) <- SP_underlined
comm1 <- as.data.frame(comm1)

#comm2 <- matrix(0,nrow=length(IDset), ncol=length(SP))
#rownames(comm2) <- IDset
#colnames(comm2) <- SP
#comm2 <- as.data.frame(comm2)


i <- 1
for(id_p in IDset)
{
	unit <- subset(fragment,ID==id_p,select=SPECIES)

	unit_SPset <- data.frame(unique(unit$SPECIES))
	unit_SPset <- unit_SPset[order(unit_SPset,decreasing=F),]
	unit_SPset <- head(unit_SPset,length(unit_SPset))

	j <- 1
	for(sp_p in unit_SPset)
	{
		while(j<=length(SP))
		{
			if(sp_p==SP[j])
			{	
				comm1[i,j] <- comm1[i,j]+1
				break
			}
			j <- j+1
		}
	}
	i <- i+1
}

write.csv(comm1,file="temp\\comm1.txt")
#write.csv(comm2,file="temp\\comm2.csv") #丰度表预留

rm(list=ls())
gc()



#BAS法和POD法 总体beta多样性指数、成对相异性指数及其组分的计算


#beta.pod.pair <- BAT::beta
#beta.pod.multi <- BAT::beta.multi

comm <- read.csv(file="temp\\comm1.csv")
comm <- as.matrix(comm[,(-1)])

#BAS法
#总体beta多样性指数
BAS_sor_multi <- beta.multi(x=comm,index.family="sorensen")
BAS_jac_multi <- beta.multi(x=comm,index.family="jaccard")

#成对相异性指数及其组分
BAS_sor_pair <- beta.pair(x=comm,index.family="sorensen")
BAS_jac_pair <- beta.pair(x=comm,index.family="jaccard")

sor_sor <- BAS_sor_pair$beta.sor #BAS 法成对相异指数 Sorensen
sor_sim <- BAS_sor_pair$beta.sim #BAS 法空间周转组分 βsim
sor_sne <- BAS_sor_pair$beta.sne #BAS 法嵌套组分 βsne
jac_jac <- BAS_jac_pair$beta.jac #BAS 法成对相异指数 Jaccard
jac_jtu <- BAS_jac_pair$beta.jtu #BAS 法空间周转组分 βjtu
jac_jne <- BAS_jac_pair$beta.jne #BAS 法嵌套组分 βjne

#POD法
#预留


rm(list=ls())
gc()