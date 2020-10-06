source("functions.R")
library(fields)

set.seed(304)
grid.vertex = seq(0, 1, length.out = 33)
grid.coord = numeric(length(grid.vertex)-1)
for(i in 1:length(grid.coord)){
  grid.coord[i] = (grid.vertex[i] + grid.vertex[i+1])/2
}

x.coord = rep(grid.coord, each=length(grid.coord))
y.coord = rep(grid.coord, length(grid.coord))
z.true = GenerateTestImage(x.coord, y.coord) #phi: 0.2766751, sine: 0.09085886
true.image = cbind(x.coord, y.coord, z.true)
quilt.plot(x.coord, y.coord, true.image[,3], nx=32, ny=32, zlim=c(-0.2,1.2))
missing.idx = sample.int(length(true.image[,3]), length(true.image[,3])/2)

sim.image = true.image[-missing.idx,]
quilt.plot(sim.image[,1], sim.image[,2], sim.image[,3], nx=32, ny=32)


# Progressive Splatting ----------------------------
splat.collection = list()
for(i in 1:dim(sim.image)[1]){
  splat.collection[[i]] = MakeSplat(i, sim.image)
}
splat.collection = MergeSplatCollection2(splat.collection, splat.num.target = 10, data = sim.image)

splat.collection.total = splat.collection
splat.collection.new = GetSplitSplatCollectionTotal(splat.collection)
splat.collection.total = c(splat.collection.total, splat.collection.new)
splat.collection.new = GetSplitSplatCollectionTotal(splat.collection.new)

while(!length(splat.collection.new) == 0){
  splat.collection.total = c(splat.collection.total, splat.collection.new)
  # length(splat.collection.total)
  splat.collection.new = GetSplitSplatCollectionTotal(splat.collection.new)
  length(splat.collection.new) == 0
}


# plot(x=NULL, y=NULL, xlim=c(-0.2,1.2), ylim=c(-0.2,1.2))
# abline(v = c(0,1), h = c(0,1), col='grey')
# for(i in 1:length(splat.collection.total)){
#   rrr = GetEllipseBoundary2D(i, splat.collection.total, sim.image, additional.scaler = 0.1)
#   lines(rrr[,1], rrr[,2])
# }



# Basis Evaluation -------------------------------
splat.collection.total.2D = MakeSplatCollection2D(splat.collection = splat.collection.total, data = sim.image, curve.fitting = TRUE)
eval.mat.obs = EvaluateBasis(x = sim.image[,1:2], splat.collection.2D = splat.collection.total.2D, data = sim.image, curve.fitting = TRUE, mm = 1)


# TPS ----------------------------------
tpsfit = Tps(sim.image[,1:2], sim.image[,3])
out.p = predictSurface(tpsfit, grid.list = list(x=grid.coord, y=grid.coord))
result.tps = t(out.p$z)
quilt.plot(x.coord, y.coord, result.tps, nx=32, ny=32)

# Splat Selection -----------------------------
splat.selected.tmp = SelectSplatForward(basis.cand.idx = 1:dim(eval.mat.obs)[2], eval.mat.obs = eval.mat.obs, basis.max.num = dim(sim.image)[1],  data = sim.image, prop.var = 0.8, BIC = FALSE, verbose = TRUE, scaling = FALSE, centering = TRUE, parallelization = TRUE)

# Make Basis Level ------------------------
splat.num.each.level = round(length(splat.selected.tmp) * 2^(3:0) / sum(2^(3:0)))
splat.num.each.level[1] = splat.num.each.level[1] - (sum(splat.num.each.level) - length(splat.selected.tmp))

area.vec = vector()
for(splat.idx in splat.selected.tmp){
  area.vec[length(area.vec)+1] = pi * sqrt(sum(splat.collection.total.2D[[splat.idx]]$major.axis^2)) * sqrt(sum(splat.collection.total.2D[[splat.idx]]$minor.axis^2))
}

area.rank = rank(area.vec)

splat.selected = list()
splat.selected[[1]] = splat.selected.tmp[ which(area.rank <= splat.num.each.level[1]) ]
for(level.idx in 2:4){
  splat.selected[[ level.idx ]] = splat.selected.tmp[ which((area.rank <= cumsum(splat.num.each.level)[level.idx]) & (area.rank > cumsum(splat.num.each.level)[level.idx-1])) ]
}


# Prediction on grid --------------------------
eval.mat.grid = EvaluateBasis(x = true.image[,1:2], splat.collection.2D = splat.collection.total.2D[unlist(rev(splat.selected))], data = sim.image, curve.fitting = TRUE, mm = 1)


pcr.tmp = pcr.jh(response = sim.image[,3], predictor = eval.mat.obs[,unlist(rev(splat.selected))], prop.var = 0.8, scaling = FALSE, BIC = FALSE, centering = TRUE)
est.tmp = cbind(1, eval.mat.grid) %*% pcr.tmp$beta
quilt.plot(x.coord, y.coord, est.tmp, nx=32, ny=32)
proposed.mse = mean((est.tmp - true.image[,3])[-which(is.na(result.tps))]^2)
tps.mse = mean((result.tps - true.image[,3])[-which(is.na(result.tps))]^2)
proposed.mse; tps.mse



