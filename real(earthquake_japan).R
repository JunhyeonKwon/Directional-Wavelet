set.seed(304)

jongmin = read.csv("~/Dropbox/Research/2019_WaveDirec/jongmin/seismological.csv")

time.tmp = strsplit(as.character(jongmin$Origin.Time), " ")

length(time.tmp)

date.tmp = vector()
for(i in 1:length(time.tmp)){
  date.tmp[i] = as.Date(time.tmp[[i]][1], format = "%Y.%m.%d")
}

japan.idx = which(jongmin$Lat >= 20 & jongmin$Lat <= 50 & jongmin$Lon >= 120 & jongmin$Lon <= 150)

c21.idx = which(date.tmp > as.Date("1999-01-01"))
nonzero.idx = which(jongmin$Magn > 0)

japan21 = jongmin[intersect(japan.idx, intersect(c21.idx, nonzero.idx)), ]
length(intersect(japan.idx, intersect(c21.idx, nonzero.idx)))

quilt.plot(japan21$Lon, japan21$Lat, japan21$Magn, nx=128, ny=128, zlim=c(0,9.1))


bgn = list()
bgn$num = round(dim(japan21)[1]*0.4)
bgn$lon = runif(n = bgn$num, min = 120, max = 150)
bgn$lat = runif(n = bgn$num, min = 20, max = 50)
bgn$magn = rnorm(n = bgn$num, mean = 0.1, sd = 0.01)

points(bgn$lon, bgn$lat, cex = 0.5)



#################
japan21.uniq = unique(as.matrix(cbind(japan21$ISC.event, japan21$Lon + rnorm(length(japan21$Lon), sd=0.001), japan21$Lat + rnorm(length(japan21$Lat), sd=0.001), japan21$Magn + rnorm(length(japan21$Magn), sd=0.01))))[,-1]
bgn.dat = cbind(bgn$lon, bgn$lat, bgn$magn)
seis.dat = rbind(japan21.uniq, bgn.dat)


quilt.plot(seis.dat[,1], seis.dat[,2], seis.dat[,3], nx=128, ny=128, zlim=c(0,9.1))


grid.japan.lon = seq(120,150,length.out = 100)
grid.japan.lat = seq(20, 50, length.out = 100)
tpsfit.seis = Tps(seis.dat[,1:2], seis.dat[,3])
out.p = predictSurface(tpsfit.seis, grid.list = list(x=grid.japan.lon, y=grid.japan.lat))
result.tps = t(out.p$z)
quilt.plot(x=rep(grid.japan.lon,each=100), y = rep(grid.japan.lat,100), z = result.tps, nx=100, ny=100, main="TPS with background noise")
points(japan21.uniq[which(japan21.uniq[,3]>5),1], japan21.uniq[which(japan21.uniq[,3]>5),2], pch=16, cex=0.5)


############################
splat.collection = list()
for(i in 1:dim(seis.dat)[1]){
  splat.collection[[i]] = MakeSplat(i, seis.dat)
}

splat.collection = MergeSplatCollection2(splat.collection, splat.num.target = 8, data = seis.dat)
splat.collection.record = splat.collection
splat.collection.new = GetSplitSplatCollectionTotal(splat.collection)
splat.collection.record = c(splat.collection.record, splat.collection.new)
splat.collection.new = GetSplitSplatCollectionTotal(splat.collection.new)

while(!length(splat.collection.new) == 0){
  splat.collection.record = c(splat.collection.record, splat.collection.new)
  # length(splat.collection.record)
  splat.collection.new = GetSplitSplatCollectionTotal(splat.collection.new)
}


# Basis Evaluation
splat.collection.record.2D = MakeSplatCollection2D(splat.collection = splat.collection.record, data = seis.dat, curve.fitting = FALSE)
eval.mat.obs = EvaluateBasis(x = seis.dat[,1:2], splat.collection.2D = splat.collection.record.2D, data = seis.dat, curve.fitting = FALSE)




# Basis Level Division

area.vec = vector()
for(i in 1:length(splat.collection.record.2D)){
  area.vec[i] = pi * sqrt(sum(splat.collection.record.2D[[i]]$major.axis^2)) * sqrt(sum(splat.collection.record.2D[[i]]$minor.axis^2))
}


basis.cand.num = round(length(splat.collection.record) * 2^(0:3) / sum(2^(0:3)))
basis.cand.num[4] = basis.cand.num[4] - (sum(basis.cand.num) - length(splat.collection.record))
basis.cand.num = rev(basis.cand.num)

basis.selected.num = round(basis.cand.num * (dim(seis.dat)[1] /  sum(basis.cand.num)))
basis.selected.num[1] = basis.selected.num[1] - (sum(basis.selected.num) - dim(seis.dat)[1])

basis.cand.cumsum = cumsum(basis.cand.num)


area.rank = rank(log(area.vec))
basis.cand.idx = list()
basis.cand.idx[[1]] = which(area.rank<=basis.cand.cumsum[1])
basis.cand.idx[[2]] = which(area.rank<=basis.cand.cumsum[2] & area.rank > basis.cand.cumsum[1])
basis.cand.idx[[3]] = which(area.rank<=basis.cand.cumsum[3] & area.rank > basis.cand.cumsum[2])
basis.cand.idx[[4]] = which(area.rank<=basis.cand.cumsum[4] & area.rank > basis.cand.cumsum[3]) # coarsest





# Basis selection
resid.dat.list = list()
splat.selected = list()

for(level.now in 4:1){
  cat("Level: ", level.now, "\n")
  
  if(level.now == 4){
    target.now = resid.dat.list[[level.now]] = seis.dat # only for the coarsest level
  }else{
    target.now = resid.dat.list[[level.now]] = cbind(resid.dat.list[[level.now+1]][,1:2], resid.dat.list[[level.now+1]][,3] - est.list[[level.now+1]])  
  }
  
  splat.selected.tmp = SelectSplatForward(basis.cand.idx[[level.now]], eval.mat.obs = eval.mat.obs, basis.max.num = basis.selected.num[level.now], data = target.now, prop.var = 0.8, BIC = FALSE, verbose = TRUE, scaling = FALSE, centering = TRUE, parallelization = TRUE)
  
  splat.selected[[level.now]] = splat.selected.tmp
  selected.cidx.tmp = which(is.element(basis.cand.idx[[level.now]], splat.selected[[level.now]]))
  
  obs.crspd.basis.selected.idx = vector()
  for(splat.idx in splat.selected[[level.now]]){
    # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
    idx.tmp = splat.collection.record[[ splat.idx ]]$idx
    
    trans.tmp = TranslateData(target.now[idx.tmp, 1:2], center = -splat.collection.record.2D[[splat.idx]]$center)
    # dev.tmp = sqrt(sum(trans.tmp^2)) wrong!
    dev.tmp = apply(trans.tmp, 1, function(x) sum(x^2))
    
    obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
  }
  
  obs.idx.tmp = unique(obs.crspd.basis.selected.idx)
  design.matrix.tmp = eval.mat.obs[ obs.idx.tmp, splat.selected[[level.now]] ]
  
  pcr.tmp = pcr.jh(response = target.now[obs.idx.tmp,3], predictor = design.matrix.tmp, prop.var = 0.8, scaling = FALSE, BIC = FALSE, centering = TRUE)
  
  beta.list[[level.now]] = pcr.tmp$beta
  est.list[[level.now]] = beta.list[[level.now]][1] + eval.mat.obs[,splat.selected[[level.now]]] %*% beta.list[[level.now]][-1]
}



################ Level by level bases plotting ###############################

par(mfrow=c(2,2))
for(j in 4:1){
  plot(NULL, xlim=c(120,150), ylim=c(20,50), xlab="", ylab="")
  map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)
  title(j)
  for(i in splat.selected[[j]]){
    rrr = GetEllipseBoundary2D(i, splat.collection.record, data = seis.dat)
    lines(rrr[,1], rrr[,2], col='darkgrey')
  }
}
par(mfrow=c(1,1))



################ Prediction on the grid ###################################
grid.japan.lon = seq(120,150,length.out = 100)
grid.japan.lat = seq(20, 50, length.out = 100)

eval.mat.grid = EvaluateBasis(x = cbind(rep(grid.japan.lon, each=100), rep(grid.japan.lat, 100)), splat.collection.2D = splat.collection.record.2D[unlist(rev(splat.selected))], data = seis.dat, curve.fitting = FALSE)


intercept.list = list()
betac.list = list()
for(level.now in 4:1){
  intercept.list[[level.now]] = beta.list[[level.now]][1]
  betac.list[[level.now]] = beta.list[[level.now]][-1]
}


intercept.vec = unlist(intercept.list)
betac.vec = unlist(rev(betac.list))

splat.cidx4 = 1:length(betac.list[[4]])
pred.grid4 = sum(intercept.vec[4]) + CenteringData(eval.mat.grid[,splat.cidx4]) %*% betac.vec[splat.cidx4]

splat.cidx3 = 1:(length(betac.list[[4]]) + length(betac.list[[3]]))
pred.grid3 = sum(intercept.vec[3:4]) + CenteringData(eval.mat.grid[,splat.cidx3]) %*% betac.vec[splat.cidx3]

splat.cidx2 = 1:(length(betac.list[[4]]) + length(betac.list[[3]]) + length(betac.list[[2]]))
pred.grid2 = sum(intercept.vec[2:4])+ CenteringData(eval.mat.grid[,splat.cidx2]) %*% betac.vec[splat.cidx2]

splat.cidx1 = 1:(length(betac.list[[4]]) + length(betac.list[[3]]) + length(betac.list[[2]]) + length(betac.list[[1]]))
pred.grid1 = sum(intercept.vec[1:4]) + CenteringData(eval.mat.grid[,splat.cidx1]) %*% betac.vec[splat.cidx1]



library("maps")
# quilt.plot(rep(grid.japan.lon, each=100), rep(grid.japan.lat, 100), pred.grid1, nx=100, ny=100, main=4, zlim=c(-2,10))
# map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)
par(mfrow=c(2,2))
quilt.plot(rep(grid.japan.lon, each=100), rep(grid.japan.lat, 100), pred.grid4, nx=100, ny=100, main=4, zlim=c(-1,10))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

quilt.plot(rep(grid.japan.lon, each=100), rep(grid.japan.lat, 100), pred.grid3, nx=100, ny=100, main=3, zlim=c(-1,10))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

quilt.plot(rep(grid.japan.lon, each=100), rep(grid.japan.lat, 100), pred.grid2, nx=100, ny=100, main=2, zlim=c(-1,10))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

quilt.plot(rep(grid.japan.lon, each=100), rep(grid.japan.lat, 100), pred.grid1, nx=100, ny=100, main=1, zlim=c(-1,10))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

par(mfrow=c(1,1))




# used data
quilt.plot(seis.dat[,1], seis.dat[,2], seis.dat[,3], nx=64, ny=64, zlim=c(0,9.1))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

# estimated data
quilt.plot(seis.dat[,1], seis.dat[,2], est.list[[1]]+est.list[[2]]+est.list[[3]]+est.list[[4]], nx=64, ny=64, zlim=c(0,9.1))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

quilt.plot(seis.dat[,1], seis.dat[,2], est.list[[1]]+est.list[[2]]+est.list[[3]]+est.list[[4]], nx=64, ny=64, zlim=c(0,9.1))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)



# true obs
quilt.plot(seis.dat[1:1179,1], seis.dat[1:1179,2], (est.list[[1]]+est.list[[2]]+est.list[[3]]+est.list[[4]])[1:1179], nx=64, ny=64, zlim=c(0,9.1))
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

# error only on true obs
quilt.plot(seis.dat[,1], seis.dat[,2], abs(est.list[[1]]+est.list[[2]]+est.list[[3]]+est.list[[4]]-seis.dat[,3]), nx=64, ny=64)


### Thin-plate spine #####
grid.japan.lon = seq(120,150,length.out = 100)
grid.japan.lat = seq(20, 50, length.out = 100)
tpsfit.seis = Tps(seis.dat[,1:2], seis.dat[,3])
out.p = predictSurface(tpsfit.seis, grid.list = list(x=grid.japan.lon, y=grid.japan.lat))
result.tps = t(out.p$z)
quilt.plot(x=rep(grid.japan.lon,each=100), y = rep(grid.japan.lat,100), z = result.tps, nx=100, ny=100, main="TPS with background noise")
map("world", xlim=c(120, 150), ylim=c(20, 50), add=T)

points(japan21.uniq[which(japan21.uniq[,3]>5),1], japan21.uniq[which(japan21.uniq[,3]>5),2], pch=4, cex=0.8)





# library(rgeos)
# library("ggplot2")
# theme_set(theme_bw())
# library("sf")
# library(rnaturalearth)
# # library("rnaturalearthdata")
# world <- ne_countries(scale = "medium", returnclass = "sf")
# library(ggplot2)
# ggplot(data = world) +
#   geom_sf() +
#   coord_sf(xlim = c(120, 150), ylim = c(20, 50), expand = FALSE)
# 
# 
# seis.df = data.frame(lon = seis.dat[,1], lat = seis.dat[,2], mag = seis.dat[,3])
# pp = ggplot(data = seis.df, aes(x=lon, y=lat, fill=mag))
# pp + geom_tile()

