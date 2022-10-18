#1.Figure 1
#install.packages("ggpubr")



library(reshape2)
library(ggpubr)
inputFile="input.txt"      
outFile="vioplot.pdf"     
setwd("D:\\biowolf\\bioR\\12.vioplotMulti")    

#读取文件
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
x=colnames(rt)[1]
colnames(rt)[1]="Type"

#把数据转换成ggplot2数据文件
data=melt(rt,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")

#绘制小提琴图
p=ggviolin(data, x="Gene", y="Expression", color = "Type", 
           ylab="Gene expression",
           xlab=x,
           legend.title=x,
           add.params = list(fill="white"),
           palette = c("blue","red"),
           width=1, add = "boxplot")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出
pdf(file=outFile, width=6, height=5)
print(p1)
dev.off()


#2.1 Figure2-5
#2.Survival Analysis
#install.packages("survival")
#install.packages("survminer")

#引用
library(survival)
library(survminer)
inputFile="input.txt"        
outFile="survival.pdf"       
var="CD74"                   #用于生存分析的变量
setwd("D:\\biowolf\\bioR\\35.survivalContinuous")      


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=rt[,c("futime","fustat",var)]

#根据中位值，把样品分为两组
group=ifelse(rt[,3]>median(rt[,3]),"High","Low")
diff=survdiff(Surv(futime, fustat) ~group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = rt) #剩余数目

#绘制
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=TRUE,
                   pval=pValue,
                   pval.size=5,
                   legend.labs=c("High", "Low"),
                   legend.title=var,
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("red", "blue"),
                   risk.table=T,
                   risk.table.height=.25)
pdf(file=outFile,onefile = FALSE,width = 6,height =5)
print(surPlot)
dev.off()

#2.2 Figure2-5
inputFile="input.txt"        
outFile="forest.pdf"         
setwd("D:\\biowolf\\bioR\\38.forest")    


rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"HR.95L")
hrHigh=sprintf("%.3f",rt$"HR.95H")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#出图格式
pdf(file=outFile, width = 6, height =4.5)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#3.pheatmap
#install.packages("pheatmap")

library(pheatmap)          
inputFile="input.txt"      
groupFile="group.txt"  #分组文件   
outFile="heatmap.pdf"      
setwd("D:\\biowolf\\bioR\\17.heatmap")      
rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)     #读取文件
ann=read.table(groupFile,header=T,sep="\t",row.names=1,check.names=F)    #读取样本属性文件


#绘制
pdf(file=outFile,width=6,height=5.5)
pheatmap(rt,
         annotation=ann,
         cluster_cols = T,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_colnames = T,
         scale="row",  #矫正
         #border_color ="NA",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


#4.cor

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#加载包
library(ggplot2)
library(ggpubr)
library(ggExtra)

inputFile="input.txt"      
gene="CD74"             #第一个因素
tumor="BRCA"              #第二个因素
setwd("D:\\biowolf\\bioR\\22.cor")      

#读取输入文件，提取基因表达量
rt=read.table(inputFile,sep="\t",header=T,check.names=F,row.names=1)
x=as.numeric(rt[gene1,])
y=as.numeric(rt[gene2,])

#相关性分析
df1=as.data.frame(cbind(x,y))
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
p1=ggplot(df1, aes(x, y)) + 
  xlab(gene)+ylab(tumor)+
  geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

#出图
pdf(file="cor.pdf",width=5,height=4.8)
print(p1)
dev.off()

#出图2
pdf(file="cor.density.pdf",width=5,height=5)
print(p2)
dev.off()

