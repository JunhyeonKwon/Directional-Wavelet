nbhd.num = 5
core.num = detectCores()-2

set.seed(304)


form = "line"


cat("Shape:", form, "\n")


############# True Data ####################################
grid.tmp = seq(0,1,length.out = 100)
result.true = GenerateTestImage(x.coord=rep(grid.tmp,each=100), y.coord=rep(grid.tmp,100), snr=Inf, form=form, z.flux = FALSE)


############# Simulation ###################################
# my.cluster = makeCluster(core.num, outfile = "./debug.txt")
my.cluster = makeCluster(core.num)
registerDoParallel(my.cluster)

time_start = Sys.time()

iter.num = 12
simulation.list = foreach(iter = 1:iter.num, .packages = c("fields", "RANN")) %dopar% {
  tryCatch({
    # Data Generation
    x.grid = runif(512, min=0, max=1)
    y.grid = runif(512, min=0, max=1)
    z.val = GenerateTestImage(x.grid, y.grid, snr=10, form=form)
    sim.image = cbind(x.grid, y.grid, z.val)
    
    # Progressive Splatting
    splat.collection = list()
    for(i in 1:dim(sim.image)[1]){
      splat.collection[[i]] = MakeSplat(i, sim.image)
    }
    splat.collection = MergeSplatCollection2(splat.collection, splat.num.target = 8, data = sim.image)
    
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
    
    # Basis Evaluation
    splat.collection.total.2D = MakeSplatCollection2D(splat.collection = splat.collection.total, data = sim.image, curve.fitting = FALSE)
    eval.mat.obs = EvaluateBasis(x = sim.image[,1:2], splat.collection.2D = splat.collection.total.2D, data = sim.image, curve.fitting = FALSE, mm = 1)
    
    
    # TPS
    tpsfit = Tps(sim.image[,1:2], sim.image[,3])
    out.p = predictSurface(tpsfit, grid.list = list(x=grid.tmp, y=grid.tmp))
    result.tps = t(out.p$z)
    
    
    # simulation data and constructed bases from it
    foo = list(data = sim.image, splat.collection.init = splat.collection, splat.collection.total = splat.collection.total, splat.collection.total.2D = splat.collection.total.2D, eval.mat.full = eval.mat.obs, result.tps = result.tps)
    
    
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
    
    # eval.mat.obs를 모든 레벨에서 사용!
    eval.mat.obs = foo$eval.mat.full
    
    for(level.now in 4:1){
      cat("Level: ", level.now, "\n")
      
      if(level.now == 4){
        target.now = resid.dat.list[[level.now]] = foo$data # only for the coarsest level
      }else{
        target.now = resid.dat.list[[level.now]] = cbind(resid.dat.list[[level.now+1]][,1:2], resid.dat.list[[level.now+1]][,3] - est.list[[level.now+1]])  
      }
      
      splat.selected.tmp = SelectSplatForward(basis.cand.idx[[level.now]], eval.mat.obs = eval.mat.obs, basis.max.num = basis.selected.num[level.now], data = target.now, prop.var = prop.var, BIC = FALSE, verbose = FALSE, scaling = FALSE, centering = TRUE)
      
      splat.selected[[level.now]] = splat.selected.tmp
      selected.cidx.tmp = which(is.element(basis.cand.idx[[level.now]], splat.selected[[level.now]]))
      
      obs.crspd.basis.selected.idx = vector()
      
      for(splat.idx in splat.selected[[level.now]]){
        # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
        idx.tmp = foo$splat.collection.total[[ splat.idx ]]$idx
        
        trans.tmp = TranslateData(target.now[idx.tmp, 1:2], center = -foo$splat.collection.total.2D[[splat.idx]]$center)
        # dev.tmp = sqrt(sum(trans.tmp^2))
        dev.tmp = apply(trans.tmp, 1, function(x) sum(x^2))
        
        obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
      }
      
      obs.idx.tmp = unique(obs.crspd.basis.selected.idx)
      
      
      design.matrix.tmp = eval.mat.obs[ obs.idx.tmp, splat.selected[[level.now]] ]
      
      pcr.tmp = pcr.jh(response = target.now[obs.idx.tmp,3], predictor = design.matrix.tmp, prop.var = prop.var, scaling = FALSE, BIC = FALSE, centering = TRUE)
      
      beta.list[[level.now]] = pcr.tmp$beta
      est.list[[level.now]] = CenteringData(eval.mat.obs[,splat.selected[[level.now]]], augmentation = TRUE) %*% beta.list[[level.now]]
    }
    
    
    intercept.list = list()
    betac.list = list()
    for(level.now in 4:1){
      intercept.list[[level.now]] = beta.list[[level.now]][1]
      betac.list[[level.now]] = beta.list[[level.now]][-1]
    }
    
    
    
    ##########################
    eval.mat.grid = EvaluateBasis(x = cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), splat.collection.2D = foo$splat.collection.total.2D[unlist(rev(splat.selected))], data = foo$data, curve.fitting = FALSE, mm = 1)
    
    
    #eval.mat.reduced.pred = EvaluateBasis(x = cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), splat.collection.2D = foo$splat.collection.total.2D[foo$selected.splat], data = foo$data, curve.fitting = TRUE)
    
    
    # on grid
    intercept.vec = unlist(intercept.list)
    betac.vec = unlist(rev(betac.list))
    
    pred.grid.list = list()
    
    splat.cidx4 = 1:length(betac.list[[4]])
    pred.grid.list[[4]] = sum(intercept.vec[4]) + CenteringData(eval.mat.grid[,splat.cidx4]) %*% betac.vec[splat.cidx4]
    
    splat.cidx3 = 1:(length(betac.list[[4]]) + length(betac.list[[3]]))
    pred.grid.list[[3]] = sum(intercept.vec[3:4]) + CenteringData(eval.mat.grid[,splat.cidx3]) %*% betac.vec[splat.cidx3]
    
    splat.cidx2 = 1:(length(betac.list[[4]]) + length(betac.list[[3]]) + length(betac.list[[2]]))
    pred.grid.list[[2]] = sum(intercept.vec[2:4])+ CenteringData(eval.mat.grid[,splat.cidx2]) %*% betac.vec[splat.cidx2]
    
    splat.cidx1 = 1:(length(betac.list[[4]]) + length(betac.list[[3]]) + length(betac.list[[2]]) + length(betac.list[[1]]))
    pred.grid.list[[1]] = sum(intercept.vec[1:4]) + CenteringData(eval.mat.grid[,splat.cidx1]) %*% betac.vec[splat.cidx1]
    
    
    
    
    pcr.mse = mean((result.true - pred.grid.list[[1]])[-which(is.na(foo$result.tps))]^2)
    tps.mse = mean((result.true - foo$result.tps)[-which(is.na(foo$result.tps))]^2)
    
    result = list(pcr.mse = pcr.mse, tps.mse=tps.mse, beta=beta.list, est=est.list, splat.selected = splat.selected, pred.grid = pred.grid.list)
    
    c(foo, result)},
    error = function(e) {
      return(paste0("Iteration ", iter, " caused the error: ", e))
    }, warning = function(w) {
      return(paste0("Iteration ", iter, " caused the warning: ", w))
    })
  
}
time_end = Sys.time()
cat("Elapsed time: ", time_end - time_start,  units(time_end - time_start), "\n")
stopCluster(my.cluster)
############################################################  



### Comparison of proposed method with TPS (Mean squared error) ####
pvec = vector()
tvec = vector()
for(i in 1:14){
  pvec[i] = simulation.list[[i]]$pcr.mse
  tvec[i] = simulation.list[[i]]$tps.mse
}

mean(pvec) # line: 0.01492715 (no iteration), sine: 0.02300716 (no iteration)
mean(tvec) # line: 0.0173149 (no iteration), sine: 0.01403779 (no iteration)





### Plotting ##################################################

# TRUE
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = result.true, nx=100, ny=100, zlim=c(-0.2,1.2), main=form)  

# observed value
quilt.plot(foo$data[,1],foo$data[,2],foo$data[,3], nx=100, ny=100, zlim=c(-0.2,1.2), main=form)

# predicted value
par(mfrow=c(2,2))
for(j in 4:1){
  quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = pred.grid.list[[j]], nx=100, ny=100, zlim=c(-0.2,1.2), main=paste0("predicted value of ", form, ", Level ", j))   
}
par(mfrow=c(1,1))

# abs(error) of the proposed method on the grid
par(mfrow=c(2,2))
for(j in 4:1){
  quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = abs(pred.grid.list[[j]] - result.true), nx=100, ny=100, zlim=c(-0.2,1.2), main=paste0("abs(error) of ", form, ", Level ", j))   
}
par(mfrow=c(1,1))

# TPS result
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = foo$result.tps, nx=100, ny=100, zlim=c(-0.2,1.2), main=form)

# abs(error) of TPS on the grid
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = abs(foo$result.tps-result.true), nx=100, ny=100, main=form, zlim = c(0,1))


# Bases disposition of each level
par(mfrow=c(2,2))
for(j in 4:1){
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", main=paste0(form, ", Level ", j))
  abline(h=c(0,1), v=c(0,1), col='grey')
  for(i in splat.selected[[j]]){
    rrr = GetEllipseBoundary2D(i, foo$splat.collection.total.2D, data = foo$data)
    lines(rrr[,1], rrr[,2])
  }
}
par(mfrow=c(1,1))
