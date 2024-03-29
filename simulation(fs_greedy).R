source("functions.R")

############# Simulation ###################################
nbhd.num = 5
grid.tmp = seq(0,1,length.out = 100)

core.num = detectCores()-1

set.seed(304)

curve.fitting = FALSE
z.flux = FALSE
iter.num = 7
# verbose.sim = FALSE

time_start = Sys.time()
simulation.list.total = list()

# simulated.direction = c("line", "sine")
simulated.direction = c("line", "sine", "circle", "cross", "phi")


for(curve.fitting in c(FALSE, TRUE)){
  
  for(form in simulated.direction){
    if(which(simulated.direction == form) == 1){
      cat("** SIMULATION SETTING **\n")
      cat("Start time:", as.character(Sys.time()), "\n")
      cat("Direction:", simulated.direction, "\n")
      cat("Iteration:", iter.num, "times for each direction\n")
      cat("CPU threads:", core.num, "\n")
      cat("Number of nbhd obs for bases construction:", nbhd.num, "\n")
      cat("z-axis fluctuation:", z.flux, "\n")
      cat("Non-linear basis:", curve.fitting, "\n\n")    
    }
    
    time_tmp1 = Sys.time()
    cat("Simulation for the shape '", form, "' has started (", as.character(time_tmp1), ").\n")
    
    # True Data
    result.true = GenerateTestImage(x.coord=rep(grid.tmp,each=100), y.coord=rep(grid.tmp,100), snr=Inf, form=form, z.flux = z.flux)
    
    # my.cluster = makeCluster(core.num, outfile = "./debug.txt")
    # my.cluster = makeCluster(core.num)
    # registerDoParallel(my.cluster)
    registerDoParallel(cores = core.num)
    
    # simulation.result = foreach(iter = 1:iter.num, .packages = c("fields", "RANN", "princurve")) %dopar% {
    simulation.result = foreach(iter = 1:iter.num) %dopar% {
      tryCatch({
        # Data Generation
        x.grid = runif(512, min=0, max=1)
        y.grid = runif(512, min=0, max=1)
        z.val = GenerateTestImage(x.grid, y.grid, snr=10, form=form, z.flux = z.flux)
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
        
        # if(verbose.sim){
        #   cat("Data generation and progressive splatting are done.\n")  
        # }
        
        # Basis Evaluation
        splat.collection.total.2D = MakeSplatCollection2D(splat.collection = splat.collection.total, data = sim.image, curve.fitting = curve.fitting)
        eval.mat.obs = EvaluateBasis(x = sim.image[,1:2], splat.collection.2D = splat.collection.total.2D, data = sim.image, curve.fitting = curve.fitting, mm = 1)
        
        # if(verbose.sim){
        #   cat("Basis evaluation on each observation is done.\n")  
        # }
        
        # TPS
        tpsfit = Tps(sim.image[,1:2], sim.image[,3])
        out.p = predictSurface(tpsfit, grid.list = list(x=grid.tmp, y=grid.tmp))
        result.tps = t(out.p$z)
        
        # if(verbose.sim){
        #   cat("Thin-plate spline for the data is done.\n")  
        # }
        
        # simulation data and constructed bases from it
        foo = list(data = sim.image, splat.collection.init = splat.collection, splat.collection.total = splat.collection.total, splat.collection.total.2D = splat.collection.total.2D, eval.mat.full = eval.mat.obs, result.tps = result.tps)
        
        
        splat.selected = list()
        splat.selected.tmp = SelectSplatForward(basis.cand.idx = 1:dim(eval.mat.obs)[2], eval.mat.obs = eval.mat.obs, basis.max.num = dim(foo$data)[1],  data = foo$data, prop.var = prop.var, BIC = FALSE, verbose = FALSE, scaling = FALSE, centering = TRUE, parallelization = FALSE)
        
        splat.num.each.level = round(length(splat.selected.tmp) * 2^(0:3) / sum(2^(0:3)))
        splat.num.each.level[4] = splat.num.each.level[4] - (sum(splat.num.each.level) - length(splat.selected.tmp))
        
        area.vec = vector()
        for(splat.idx in splat.selected.tmp){
          area.vec[length(area.vec)+1] = pi * sqrt(sum(foo$splat.collection.total.2D[[splat.idx]]$major.axis^2)) * sqrt(sum(foo$splat.collection.total.2D[[splat.idx]]$minor.axis^2))
        }
        
        area.rank = rank(area.vec)
        
        splat.selected[[1]] = splat.selected.tmp[ which(area.rank <= splat.num.each.level[1]) ]
        for(level.idx in 2:4){
          splat.selected[[ level.idx ]] = splat.selected.tmp[ which((area.rank <= cumsum(splat.num.each.level)[level.idx]) & (area.rank > cumsum(splat.num.each.level)[level.idx-1])) ]
        }
        
        
        beta.list = list()
        est.list = list()
        resid.dat.list = list()
        
        prop.var = 0.8
        
        for(level.now in 4:1){
          
          if(level.now == 4){
            target.now = resid.dat.list[[level.now]] = foo$data # only for the coarsest level
          }else{
            target.now = resid.dat.list[[level.now]] = cbind(resid.dat.list[[level.now+1]][,1:2], resid.dat.list[[level.now+1]][,3] - est.list[[level.now+1]])  
          }
          
          
          obs.crspd.basis.selected.idx = vector()
          
          if(curve.fitting){
            # for basis with non-linear direction
            for(splat.idx in splat.selected[[level.now]]){
              # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
              idx.tmp = foo$splat.collection.total[[ splat.idx ]]$idx
              
              trans.tmp = GetDeviationFromCurve(target.now[idx.tmp, 1:2], foo$splat.collection.total.2D[[ splat.idx ]]$polygon)
              dev.tmp = sqrt(trans.tmp$lambda^2 + trans.tmp$dev^2)
              
              obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
            }
          }else{
            # for basis with linear direction
            for(splat.idx in splat.selected[[level.now]]){
              # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
              idx.tmp = foo$splat.collection.total[[ splat.idx ]]$idx
              
              trans.tmp = TranslateData(target.now[idx.tmp, 1:2], center = -foo$splat.collection.total.2D[[splat.idx]]$center)
              # dev.tmp = sqrt(sum(trans.tmp^2))
              dev.tmp = apply(trans.tmp, 1, function(x) sum(x^2))
              
              obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
            }
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
        eval.mat.grid = EvaluateBasis(x = cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), splat.collection.2D = foo$splat.collection.total.2D[unlist(rev(splat.selected))], data = foo$data, curve.fitting = curve.fitting, mm = 1)
        
        
        #eval.mat.reduced.pred = EvaluateBasis(x = cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), splat.collection.2D = foo$splat.collection.total.2D[foo$selected.splat], data = foo$data, curve.fitting = curve.fitting)
        
        
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
        }
        # , warning = function(w) {
        #   return(paste0("Iteration ", iter, " caused the warning: ", w))
        # }
      )
      
    }
    # stopCluster(my.cluster)
    
    simulation.list.total[[which(simulated.direction == form)]] = simulation.result
    
    time_tmp2 = Sys.time()
    cat("Simulation for the shape: '", form, "' has ended (", as.character(time_tmp2), ").\n")
    cat("(elapsed time: ", time_tmp2 - time_tmp1, units(time_tmp2 - time_tmp1), ")\n\n")
    
    saveRDS(simulation.list.total, file = paste0("sim0915_",curve.fitting, "_", which(simulated.direction == form), ".rds"))
    
    if(which(simulated.direction == form) == length(simulated.direction)){
      time_end = Sys.time()
      cat("Simulation for the every shape has ended (", as.character(time_end), ").\n")
      cat("(total elapsed time: ", time_end - time_start, units(time_end - time_start), ")\n\n")
    }
  }  
  
}




############################################################  





########### Results #####################################
form = "cross"
result.tmp = simulation.list.total[[which(simulated.direction == form)]]

# True Data
result.true = GenerateTestImage(x.coord=rep(grid.tmp,each=100), y.coord=rep(grid.tmp,100), snr=Inf, form=form, z.flux = z.flux)

### Comparison of proposed method with TPS (Mean squared error) ####
pvec = vector()
tvec = vector()
for(i in 1:iter.num){
  if(is.character(result.tmp[[i]])){
    pvec[i] = NA
    tvec[i] = NA
  }else{
    pvec[i] = result.tmp[[i]]$pcr.mse
    tvec[i] = result.tmp[[i]]$tps.mse  
  }
}

mean(pvec, na.rm = TRUE) # line: 0.01492715 (no iteration), sine: 0.02300716 (no iteration)
mean(tvec, na.rm = TRUE) # line: 0.0173149 (no iteration), sine: 0.01403779 (no iteration)



### Plotting ##################################################
sim.idx = 1
foo = result.tmp[[sim.idx]]

# TRUE
quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = result.true, nx=100, ny=100, zlim=c(-0.2,1.2), main=form)  

# observed value
quilt.plot(foo$data[,1],foo$data[,2],foo$data[,3], nx=100, ny=100, zlim=c(-0.2,1.2), main=form)

# predicted value
par(mfrow=c(2,2))
for(j in 4:1){
  quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = foo$pred.grid[[j]], nx=100, ny=100, zlim=c(-0.2,1.2), main=paste0("predicted value of ", form, ", Level ", j))   
}
par(mfrow=c(1,1))

# abs(error) of the proposed method on the grid
par(mfrow=c(2,2))
for(j in 4:1){
  quilt.plot(x=rep(grid.tmp,each=100), y = rep(grid.tmp,100), z = abs(foo$pred.grid[[j]] - result.true), nx=100, ny=100, zlim=c(-0.2,1.2), main=paste0("abs(error) of ", form, ", Level ", j))   
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
