################ Bases Configuration #######################

simulation.list <- readRDS("~/Dropbox/Research/2019_WaveDirec/Jang_rep/sim.result/simulation2020_seed304_14times_list(sine).rds")

form="sine"
grid.tmp = seq(0,1,length.out = 100)
result.true = GenerateTestImage(x.coord=rep(grid.tmp,each=100), y.coord=rep(grid.tmp,100), snr=Inf, form=form, z.flux = FALSE)

sim.idx = 2
foo = simulation.list[[sim.idx]]


area.vec = vector()
curvlen.vec = vector()


for(i in 1:length(foo$splat.collection.total.2D)){
  area.vec[i] = pi * sqrt(sum(foo$splat.collection.total.2D[[i]]$major.axis^2)) * sqrt(sum(foo$splat.collection.total.2D[[i]]$minor.axis^2))
  curvlen.vec[i] = sqrt(sum(foo$splat.collection.total.2D[[i]]$major.axis^2))
}



basis.cand.num = round(length(foo$splat.collection.total) * 2^(0:3) / sum(2^(0:3)))
basis.cand.num[4] = basis.cand.num[4] - (sum(basis.cand.num) - length(foo$splat.collection.total))
basis.cand.num = rev(basis.cand.num)

basis.selected.num = round(basis.cand.num * (dim(foo$data)[1] /  sum(basis.cand.num)))
basis.selected.num[1] = basis.selected.num[1] - (sum(basis.selected.num) - dim(foo$data)[1])

basis.cand.cumsum = cumsum(basis.cand.num)

area.rank = rank(log(area.vec))
basis.cand.idx = list()
basis.cand.idx[[1]] = which(area.rank<=basis.cand.cumsum[1])
basis.cand.idx[[2]] = which(area.rank<=basis.cand.cumsum[2] & area.rank > basis.cand.cumsum[1])
basis.cand.idx[[3]] = which(area.rank<=basis.cand.cumsum[3] & area.rank > basis.cand.cumsum[2])
basis.cand.idx[[4]] = which(area.rank<=basis.cand.cumsum[4] & area.rank > basis.cand.cumsum[3]) # coarsest



################# Wendland's Multiscale ####################################
prop.var = 0.8

splat.selected = list()
beta.list = list()
est.list = list()
resid.dat.list = list()

# 여기서 계산한 eval.mat.obs로 모든 레벨에서 사용!
eval.mat.obs = GetEvaluatedMatrix(foo$data[,1:2], level.end = 1, list(foo$splat.collection.total), data = foo$data, curve.fitting = TRUE)

for(level.now in 4:1){
  cat("Level: ", level.now, "\n")
  
  if(level.now == 4){
    target.now = resid.dat.list[[level.now]] = foo$data # only for the coarsest level
  }else{
    target.now = resid.dat.list[[level.now]] = cbind(resid.dat.list[[level.now+1]][,1:2], resid.dat.list[[level.now+1]][,3] - est.list[[level.now+1]])  
  }
  
  splat.selected.tmp = SelectSplatForward(basis.cand.idx[[level.now]], eval.mat.obs = eval.mat.obs, basis.max.num = basis.selected.num[level.now], data = target.now, prop.var = 0.8, BIC = FALSE, verbose = TRUE, scaling = FALSE)
  
  splat.selected[[level.now]] = splat.selected.tmp
  selected.cidx.tmp = which(is.element(basis.cand.idx[[level.now]], splat.selected[[level.now]]))
  
  obs.crspd.basis.selected.idx = vector()
  for(splat.idx in splat.selected[[level.now]]){
    # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
    idx.tmp = foo$splat.collection.total[[ splat.idx ]]$idx
    
    trans.tmp = GetDeviationFromCurve(target.now[idx.tmp, 1:2], foo$splat.collection.total.2D[[ splat.idx ]]$polygon)
    dev.tmp = sqrt(trans.tmp$lambda^2 + trans.tmp$dev^2)
    
    obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
  }
  
  obs.idx.tmp = unique(obs.crspd.basis.selected.idx)
  design.matrix.tmp = eval.mat.obs[ obs.idx.tmp, splat.selected[[level.now]] ]
  
  pcr.tmp = pcr.jh(response = target.now[obs.idx.tmp,3], predictor = design.matrix.tmp, prop.var = 0.8, scaling = FALSE, BIC = FALSE)
  
  beta.list[[level.now]] = pcr.tmp$beta
  est.list[[level.now]] = eval.mat.obs[,splat.selected[[level.now]]] %*% beta.list[[level.now]]
}




##########################
par(mfrow=c(2,2))
for(j in 4:1){
  plot(NULL, xlim=c(-0.2,1.2), ylim=c(-0.2,1.2), xlab="", ylab="")
  title(j)
  abline(h=c(0,1), v=c(0,1), col='grey')
  for(i in splat.selected[[j]]){
    rrr = GetCurvedSplatBoundary2D(i, foo$splat.collection.total.2D, foo$data, additional.scaler = 1)
    lines(rrr[,1], rrr[,2])
  }
  
}
par(mfrow=c(1,1))


quilt.plot(foo$data[,1:2], est.list[[4]], nx=32, ny=32, zlim=c(-0.35,1.1), main=form)
quilt.plot(foo$data[,1:2], est.list[[4]] + est.list[[3]], nx=32, ny=32, zlim=c(-0.35,1.1), main=form)
quilt.plot(foo$data[,1:2], est.list[[4]] + est.list[[3]] + est.list[[2]], nx=32, ny=32, zlim=c(-0.35,1.1), main=form)
quilt.plot(foo$data[,1:2], est.list[[4]] + est.list[[3]] + est.list[[2]] + est.list[[1]], nx=32, ny=32, zlim=c(-0.35,1.1), main=form)


quilt.plot(foo$data[,1:2], abs(est.list[[4]] - foo$data[,3]), nx=32, ny=32, zlim=c(0,1), main=form)
quilt.plot(foo$data[,1:2], abs(est.list[[4]] + est.list[[3]] - foo$data[,3]), nx=32, ny=32, zlim=c(0,1), main=form)
quilt.plot(foo$data[,1:2], abs(est.list[[4]] + est.list[[3]] + est.list[[2]] - foo$data[,3]), nx=32, ny=32, zlim=c(0,1), main=form)
quilt.plot(foo$data[,1:2], abs(est.list[[4]] + est.list[[3]] + est.list[[2]] + est.list[[1]] - foo$data[,3]), nx=32, ny=32, zlim=c(0,1), main=form)

##########################

eval.mat.grid = GetEvaluatedMatrix(eval.dat = cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), level.end = 1, list(foo$splat.collection.total[unlist(rev(splat.selected))]), data = foo$data, curve.fitting = TRUE)

#eval.mat.reduced.pred = EvaluateBasis(x = cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), splat.collection.2D = foo$splat.collection.total.2D[foo$selected.splat], data = foo$data, curve.fitting = TRUE)


# on grid
beta.vec = unlist(rev(beta.list))

pred.grid4 = eval.mat.grid[,1:length(beta.list[[4]])] %*% beta.vec[1:length(beta.list[[4]])]
pred.grid3 = eval.mat.grid[,1:(length(beta.list[[4]]) + length(beta.list[[3]]))] %*% beta.vec[1:(length(beta.list[[4]]) + length(beta.list[[3]]))]
pred.grid2 = eval.mat.grid[,1:(length(beta.list[[4]]) + length(beta.list[[3]]) + length(beta.list[[2]]))] %*% beta.vec[1:(length(beta.list[[4]]) + length(beta.list[[3]]) + length(beta.list[[2]]))]
pred.grid1 = eval.mat.grid[,1:(length(beta.list[[4]]) + length(beta.list[[3]]) + length(beta.list[[2]]) + length(beta.list[[1]]))] %*% beta.vec[1:(length(beta.list[[4]]) + length(beta.list[[3]]) + length(beta.list[[2]])  + length(beta.list[[1]]))]



quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = result.true, nx=100, ny=100, zlim=c(-0.2,1.2), main="TRUE")  

par(mfrow=c(2,2))
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = pred.grid1, nx=100, ny=100, zlim=c(-0.4,1.4), main="T1, W")
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = pred.grid2, nx=100, ny=100, zlim=c(-0.4,1.4), main="T2, W")  
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = pred.grid3, nx=100, ny=100, zlim=c(-0.4,1.4), main="T3, W")
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = pred.grid4, nx=100, ny=100, zlim=c(-0.4,1.4), main="T4, W")  
par(mfrow=c(1,1))

mean((result.true - pred.grid4)^2)
mean((result.true - pred.grid3)^2)
mean((result.true - pred.grid2)^2)
mean((result.true - pred.grid1)^2)

foo$mse.tps
mean((result.true - foo$result.tps)^2)








