library(vegan)

Birdsnp<-read.csv('satchori_bird_02.csv',header = TRUE)
names(Birdsnp)
env<- Birdsnp[,6:11]
spp<- Birdsnp[,33:59]

ccamodel <- cca(spp~., env)
ccamodelb <- cca(spp~1, env)
add1(ccamodelb, scope = formula(ccamodel), test = "permutation")
ccamodelfinal <- cca(spp~ Canopy_openness + Distance_from_road, env)
add1(ccamodelfinal, scope = formula(ccamodel), test = "permutation")
ano<-anova(ccamodelfinal, by = "terms")
ccamodelfinal
summary(ccamodelfinal)

anova.cca(ccamodelfinal)
anova.cca(ccamodelfinal, by="terms")
anova.cca(ccamodelfinal, by="axis")

plot(ccamodelfinal, xlim=c(-3,3), ylim=c(-3,3), display=c("sp","cn"))

## Common but bad way: use all variables you happen to have in your
## environmental data matrix
plot(ccamodel, xlim=c(-3,3), ylim=c(-3,3), display=c("sp","cn"))

## RDA (ref_Mark)
rdamodel <- rda(spp~., env)
rdamodel2 <- rda(spp~1, env)
add1(rdamodel2, scope = formula(rdamodel), test = "permutation")
rdamodel3 <- rda(spp~ Canopy_openness + Distance_from_road, env)
add1(rdamodel3, scope = formula(ccamodel), test = "permutation")

ano_rda<-anova(rdamodel3, by = "terms")
rdamodel3
summary(rdamodel3)

anova(rdamodel3)
anova(rdamodel3, by="terms")
anova(rdamodel3, by="axis")

plot(rdamodel3, xlim=c(-2,2), ylim=c(-2,2), display=c("sp","cn"))


#Multi-Dimensional Scaling
spp.t = t(Birdsnp)
spp.t2=  spp.t[33:59,]
spp.eu= dist(spp.t2, method = "euclidean")
spp.pco= cmdscale(spp.eu, k=2)
plot(spp.pco)
text(spp.pco,labels = rownames(spp.pco), cex = 1,adj = NULL, pos = 3, offset = 1)


#MDS for species score
spp.bc = vegdist(spp.t2, method = "bray")
spp.pco.bc = cmdscale(spp.bc, eig = TRUE) #Previously we did with Euclidean dist

spp.wa = wascores(spp.pco.bc$points, spp.t2)
plot(spp.pco.bc$points, xlab = "MDS1", ylab = "MDS2")
text(spp.pco.bc$points, labels = rownames(spp.pco.bc$points))
points(spp.wa, pch = "+", col = "darkgreen")


#NMDS
library(MASS)
Birdsnp<-read.csv('satchori_bird_02.csv',header = TRUE)
env<- Birdsnp[,6:11]
bspp<- Birdsnp[,33:59]
bspp.eu = dist(bspp, method = "euclidean") 
bspp.nmds = MASS::isoMDS(bspp.eu)
plot(bspp.nmds$points, xlab = "NMDS1", ylab = "NMDS2", asp = 1)
text(bspp.nmds$points, cex = 1,adj = NULL, pos = 3, offset = 1)

cor(scores(bspp.nmds)[,1], env, method = "spearman")

for (i in 1:5) print (isoMDS(bspp.eu, k=i, trace=FALSE)$stress)


##for species score
library(vegan)
rankindex(env, spp)
spp.nmds2 = metaMDS(spp, distance = "manhattan")
for (i in 1:5) print (metaMDS(spp,distance = "manhattan", k=i, trace=FALSE)$stress*100)


spp.nmds2 = metaMDS(spp, distance = "manhattan", k = 3) #based on  the prevous result
plot(spp.nmds2, type = "t", display = "sites") #produce plot
op = plot(spp.nmds2, type = "p", display = c("site", "species"))
identify(op, what = "species", cex = 0.8, col = "red")

scores(spp.nmds2)
spp.nmds2$data
wascores(spp.nmds2$po, wisconsin(spp), expand = TRUE)
wascores(env, spp)
eigengrad(env, spp)


#PCA with RDA
spp.pca= rda(spp)
scores(spp.pca, choices = 1:2,'sites')
scores(spp.pca, choices = 1:2,'species')
eigenvals(spp.pca)
summary(spp.pca)
screeplot(spp.pca)


spp.nam = abbreviate(colnames(spp), minlength = 5)
plot(spp.pca, display = "sites", type = "text")
text(spp.pca, display = "species", labels = spp.nam, col = "red", cex = 0.6)


op = plot(spp.pca, type= "n")
text(spp.pca, display = "sites")
points(spp.pca, display = "species", pch = "+", col = "red")
sp = identify(op, "species", col = "red", cex = 0.6, labels = spp.nam)
biplot(spp.pca)


#CCA(p390)
library(vegan)
fw.cca = cca(fw.biol ~ TDS + pH + DO2, data = fw.env)