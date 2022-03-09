

## SCRIPT


# housekeeping
rm( list = ls())

# data
load ( file.path (PROJHOME, "Output" ,"matrix", "phenotypematrix.rda"))

# libraries
library( nlme )
library(MuMIn)
library(ggplot2)
library(familiarUnfamiliar)


# models

# context independent hierachy
  m = modFU(df, "t.fam" , "t.unfam")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])
ggplot(df, aes ( t.fam , t.unfam , group.num))+
         geom_point () +
         geom_smooth(method = "lm")+
  xlab ( "leader score (familiar flights)")+
  ylab ( "leader score (unfamiliar flights)")

# predictions
# first paragraph.


# speed and cranio-caudal
library(ggplot2)

m = modFU(df, "fam.as" , "t.fam") # switched dependent var and ind var, to account for random intercepts on first two models
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])
m = modFU(df, "fam.re" , "t.fam")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

t.test ( df$fam.re[df$group.num == 1] ,
         df$fam.re[df$group.num == 2])

m = modFU(df, "t.fam" , "fb.fam")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])



m = modFU(df, "t.fam" , "rel.speed")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.fam" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "fam.as" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "fb.fam" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "fb.unfam" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.fam" , "peakfid")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])


# second paragraph
m = modFU(df, "t.unfam" , "fam.rel.re")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])


m = modFU(df, "t.unfam" , "rel.speed")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.unfam" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

# third paragraph


m = modFU(df, "t.fam" , "expl")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.unfam" , "expl")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.fam" , "dave")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.unfam" , "dave")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "dave" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "dave" , "sex")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])



m = modFU(df, "t.unfam" , "rel.speed")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.unfam" , "mean.mass")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

m = modFU(df, "t.unfam" , "fb.fam")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])





# PCA

# data
load ( file.path (PROJHOME, "Output" ,"matrix", "phenotypematrix.rda"))

# libraries
library( nlme )
library(MuMIn)
library(ggplot2)
library(familiarUnfamiliar)

df3 = df[complete.cases(df),]
MORP = df3[, c("ToeLength",  "TarsoMetatarsus", "Wing_length", "Body_width", "mean.mass", "Root_chord.cm.2", "Wing_area" )]
PERS = df3[, c("dave", "bold", "expl", "neo")]
FLIG = df3[,c( "fam.as","fam.re", "peakfid", "fb.fam")]
pcaMORP = princomp(MORP, cor=T, score=T)
pcaPERS = princomp(PERS, cor=T, score=T)
pcaFLIG = princomp(FLIG, cor=T, score=T)
df2 = cbind( df3,
             MORPpc1 = pcaMORP$scores[,1],
             MORPpc2 = pcaMORP$scores[,2],
             MORPpc3 = pcaMORP$scores[,3],
             FLIGpc1 = pcaFLIG$scores[,1],
             FLIGpc2 = pcaFLIG$scores[,2],
             FLIGpc3 = pcaFLIG$scores[,3],
             PERSpc1 = pcaPERS$scores[,1],
             PERSpc2 = pcaPERS$scores[,2],
             PERSpc3 = pcaPERS$scores[,3]
)



# loadings
morp = pcaMORP$loadings
pcaMORP$sdev
pcaMORP$eig

summary(lm(df2$MORPpc1 ~ df2$t.fam))

pcaPERS$loadings
pcaMORP$loadings
view.mat(loadings(pcaFLIG))
capture.output(summary(pcaFLIG), file = file.path (PROJHOME, "Output" , "pcaloadings", "summaryFLIG.csv"))
capture.output(summary(pcaMORP), file = file.path (PROJHOME, "Output" , "pcaloadings", "summaryMORP.csv"))
capture.output(summary(pcaPERS), file = file.path (PROJHOME, "Output" , "pcaloadings", "summaryPERS.csv"))

write.csv(morp, file= file.path ( PROJHOME , "Output" , "pcaloadings" , "morp.csv" ))
flig = pcaFLIG$loadings
summary(pcaFLIG)
write.csv(flig, file= file.path ( PROJHOME , "Output" , "pcaloadings" , "flig.csv" ))
pers = pcaPERS$loadings
pers
write.csv(pers, file= file.path ( PROJHOME , "Output" , "pcaloadings" , "pers.csv" ))


m1 = lme ( t.fam ~
            MORPpc1 +
             MORPpc2 +
             MORPpc3 +
             FLIGpc1 +
             FLIGpc2 +
             FLIGpc3 +
             PERSpc1 +
             PERSpc2 +
             PERSpc3, random = list ( group.num=~1),data = df2)
library(MuMIn)
d1 = dredge(m1)
d1
save( d1 , file = file.path (PROJHOME , "Output" , "dredge", "dredgefam.rda" ))
d1 = d1[1:8,]
write.csv(d1 , file = file.path (PROJHOME , "Output" , "dredge", "dredgefam.csv"))
summary(m1)


m2 = lme ( t.unfam ~
             MORPpc1 +
             MORPpc2 +
             MORPpc3 +
             FLIGpc1 +
             FLIGpc2 +
             FLIGpc3 +
             PERSpc1 +
             PERSpc2 +
             PERSpc3, random = list ( group.num=~1),data = df2)
d2 = dredge(m2)
d2
save( d2 , file = file.path (PROJHOME , "Output" , "dredge", "dredgeunfam.rda" ))
d2 = d2[1:8,]
write.csv(d2 , file = file.path (PROJHOME , "Output" , "dredge", "dredgeunfam.csv"))
summary(m1)



### PLOTSSS

# data
load ( file.path (PROJHOME, "Output" ,"matrix", "phenotypematrix.rda"))

# libraries
library(ggplot2)


{
g1 = ggplot(df , aes ( y = t.unfam , x = fam.re , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g2 = ggplot(df , aes ( y = t.fam , x = fam.re , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g3 = ggplot(df , aes ( y = t.unfam , x = fam.as , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g4 = ggplot(df , aes ( y = t.fam , x = fam.as , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()
}
library(gridExtra)
grid.arrange(g1,g2,g3,g4,ncol = 2)

{
g1 = ggplot(df , aes ( y = t.fam , x = mean.mass , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g2 = ggplot(df , aes ( y = t.fam , x = fb.fam , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g3 = ggplot(df , aes ( y = fb.fam , x = mean.mass , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g4 = ggplot(df , aes ( y = t.unfam , x = mean.mass , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g5 = ggplot(df , aes ( y = t.unfam , x = fb.unfam , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()

g6 = ggplot(df , aes ( y = fb.unfam , x = mean.mass , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point()
}


grid.arrange(g1,g4,g2,g5,g5,g6,ncol = 2)

# Fig. 4

# data
load ( file.path (PROJHOME, "Output" ,"matrix", "phenotypematrix.rda"))

# libraries
library( nlme )
library(MuMIn)
library(ggplot2)
library(familiarUnfamiliar)

# objects
cols1 = c("#8A84E2","#F25757","#F2E863", "green")

df3 = data.frame ( lead = c( df$t.fam, df$t.unfam),
                   expl = rep ( df$expl, 2),
                   famu = c( rep( "Fam", nrow(df)),
                             rep( "Unf", nrow(df))))
g1 = ggplot(df3 , aes ( y = lead , x = expl , col = famu))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=cols1[1:2])+
  geom_point() +
  theme(legend.position = "none")

g2 = ggplot(df , aes ( y = bold , x = neo ))+
  geom_smooth(method = "lm")+
  geom_point()+
  theme(legend.position = "none")

g3 = ggplot(df , aes ( y = dave , x = neo , col = group.num))+
  geom_smooth(method = "lm")+
  geom_point() +
  theme(legend.position = "none")

g4 = ggplot(df , aes ( y = dave , x = mean.mass , col = sex))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=cols1[3:4])+
  geom_point()



boxy.lady( "dave" , "sex", df, jit = 0.1, boxcol = cols1[3:4], colz = "grey50", cexz = 5 )

gridExtra::grid.arrange(g1,g2, g3, g4, ncol = 2)

ggplot(df , aes ( y = t.fam , x = t.unfam , col = group.num, ))+
  geom_smooth(method = "lm")+
  geom_point()


## boldness prediction

# data
load ( file.path (PROJHOME, "Output" ,"matrix", "phenotypematrix.rda"))

# libraries
library(ggplot2)

sum = summary( lm ( df$bold ~ df$t.fam))
sum$coefficients[2,3]
sum = summary( lm ( df$bold ~ df$t.unfam))
sum$coefficients[2,3]
sum = summary( lm ( df$neo ~ df$t.fam))
sum$coefficients[2,3]
sum = summary( lm ( df$neo ~ df$t.unfam))
sum$coefficients[2,3]
sum = summary( lm ( df$expl ~ df$t.fam))
sum$coefficients[2,3]
sum = summary( lm ( df$expl ~ df$t.unfam))
sum$coefficients[2,3]

d1 = df[df$group.num == 1, ]
d2 = df[df$group.num == 2, ]

mod= lm (d2$t.fam~ d2$fam.re )
summary(mod)
ck = cooks.distance(mod)
mean(ck) *3
d2cook = d2[-6,]
mod2= lm (d2cook$t.fam~ d2cook$fam.re )
summary (mod2)
summary ( lm (df$expl ~ df$neo ))
summary ( lm (df$bold ~ df$neo ))
summary ( lm (df$expl ~ df$bold))
summary ( lm (df$dave ~ df$bold))
summary ( lm (df$dave ~ df$neo ))


ggplot( df, aes ( neo , bold))+
  geom_point()+
  geom_smooth(method = "lm")

ggplot( df, aes ( neo , dave))+
  geom_point()+
  geom_smooth(method = "lm")

m = modFU(df2, "t.fam" , "PERSpc2")
summary( m[[1]])
summary( m[[2]])
summary( m[[3]])

