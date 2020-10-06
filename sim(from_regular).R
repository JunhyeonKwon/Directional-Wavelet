source("functions.R")

library(progress)

iter.num = 20

# Generating simulation data ------------------------------------
set.seed(304)
grid.vertex = seq(0, 1, length.out = 33)
grid.coord = numeric(length(grid.vertex)-1)
for(i in 1:length(grid.coord)){
  grid.coord[i] = (grid.vertex[i] + grid.vertex[i+1])/2
}

x.coord = rep(grid.coord, each=length(grid.coord))
y.coord = rep(grid.coord, length(grid.coord))

z.true = GenerateTestImage(x.coord, y.coord, noise.sd = 0, form = "phi")
true.image = cbind(x.coord, y.coord, z.true)


complete.image.list = list()
sim.image.list = list()
missing.idx.list = list()
for(i in 1:iter.num){
  z.value = GenerateTestImage(x.coord, y.coord, form = "phi")
  complete.image = cbind(x.coord, y.coord, z.value)
  missing.idx = sample.int(length(complete.image[,3]), length(complete.image[,3])/2)
  sim.image = complete.image[-missing.idx,]
  
  complete.image.list[[i]] = complete.image
  missing.idx.list[[i]] = missing.idx
  sim.image.list[[i]] = sim.image
}
# phi: 0.2766751 -> 0.1243541(sd=0.1) / 0.06217707(sd=0.05)
# sine: 0.09085886 -> 0.1243541(sd=0.1) / 0.06217707(sd=0.05)


# Iteration starts! --------------------------------------

pb = progress_bar$new(total = iter.num)

simulation.result = list()
for(iter in 1:iter.num){
  pb$tick()
  
  if(iter == 1){
    time_start = Sys.time()
    cat("Start time:", as.character(time_start), "\n")
  }
  
  
  complete.image = complete.image.list[[iter]]
  missing.idx = missing.idx.list[[iter]]
  sim.image = sim.image.list[[iter]]
  
  # Generating Bases through Progressive Splatting --------------------------
  splat.collection.total = GenerateSplat(sim.image)
  splat.collection.total.2D = MakeSplatCollection2D(splat.collection = splat.collection.total, data = sim.image, curve.fitting = TRUE)
  
  # Basis Evaluation -------------------------------
  eval.mat.obs = EvaluateBasis(x = sim.image[,1:2], splat.collection.2D = splat.collection.total.2D, data = sim.image, curve.fitting = TRUE)
  
  # Splat Selection -----------------------------
  splat.selected.tmp = SelectSplatForward(basis.cand.idx = 1:dim(eval.mat.obs)[2], basis.init.cidx = NULL, eval.mat.obs = eval.mat.obs, basis.max.num = dim(sim.image)[1],  data = sim.image, prop.var = 0.8, verbose = FALSE, parallelization = TRUE)
  
  # Split Splats to Each Level ------------------------
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
  
  
  
  # Bases evaluation on grid --------------------------
  eval.mat.grid = EvaluateBasis(x = complete.image[,1:2], splat.collection.2D = splat.collection.total.2D[unlist(rev(splat.selected))], data = sim.image, curve.fitting = TRUE)
  
  # TPS ----------------------------------
  tpsfit = Tps(sim.image[,1:2], sim.image[,3])
  out.p = predictSurface(tpsfit, grid.list = list(x=grid.coord, y=grid.coord))
  est.tps = t(out.p$z)
  mse.tps = mean((est.tps - true.image[,3])[-which(is.na(est.tps))]^2)
  
  # PCR with entire bases -----------------------
  eval.mat.rearranged = eval.mat.obs[, unlist(rev(splat.selected))]
  pcr.tmp = pcr.jh(response = sim.image[,3], predictor = eval.mat.rearranged, prop.var = 0.8, scaling = FALSE, BIC = FALSE, centering = TRUE)
  est.pcr.entire = CenteringData(eval.mat.grid, augmentation = TRUE) %*% pcr.tmp$beta
  mse.pcr.entire = mean((est.pcr.entire - true.image[,3])[-which(is.na(est.tps))]^2)
  
  # PCR with entire bases + orthogonalization + thsd'ing -----------------------
  basis.num = unlist(rev(lapply(splat.selected, length)))
  basis.cidx = list()
  for(i in 1:4){
    basis.cidx[[i]] = ( c(1, (cumsum(basis.num)+1)[-4])[i] ):( cumsum(basis.num)[i] )
  }
  
  eval.mat.rearranged.ortho = CenteringData(eval.mat.rearranged)
  for(i in 1:3){
    proj.dat = eval.mat.rearranged.ortho[,unlist(basis.cidx[1:i])]
    proj.mat = proj.dat %*% solve(t(proj.dat) %*% proj.dat) %*% t(proj.dat)
    resid.mat = diag(1, dim(proj.mat)[1]) - proj.mat
    eval.mat.rearranged.ortho[, basis.cidx[[i+1]] ] = resid.mat %*% eval.mat.rearranged[, basis.cidx[[i+1]] ]
  }
  
  pcr.tmp2 = pcr.jh(response = sim.image[,3], predictor = eval.mat.rearranged.ortho, prop.var = 0.8, scaling = FALSE, BIC = FALSE, centering = TRUE)
  
  pcr.coef = pcr.coef2 = pcr.tmp2$beta
  for(i in 4:1){
    pcr.coef2[-1][basis.cidx[[i]]] = ebayesthresh(pcr.coef[-1][basis.cidx[[i]]])
  }
  
  eval.mat.grid.ortho = CenteringData(eval.mat.grid)
  for(i in 1:3){
    proj.dat = eval.mat.grid.ortho[,unlist(basis.cidx[1:i])]
    proj.mat = proj.dat %*% solve(t(proj.dat) %*% proj.dat) %*% t(proj.dat)
    resid.mat = diag(1, dim(proj.mat)[1]) - proj.mat
    eval.mat.grid.ortho[, basis.cidx[[i+1]]] = resid.mat %*% eval.mat.grid[, basis.cidx[[i+1]]]
  }
  
  est.pcr.entire.ortho = CenteringData(eval.mat.grid.ortho, augmentation = TRUE) %*% pcr.coef
  est.pcr.entire.ortho.thsd = CenteringData(eval.mat.grid.ortho, augmentation = TRUE) %*% pcr.coef2
  
  mse.pcr.entire.ortho = mean((est.pcr.entire.ortho - true.image[,3])[-which(is.na(est.tps))]^2)
  mse.pcr.entire.ortho.thsd = mean((est.pcr.entire.ortho.thsd - true.image[,3])[-which(is.na(est.tps))]^2)
  
  # SVD on design matrix (using entire bases) -----------------------
  svd.tmp = svd(eval.mat.rearranged)
  reduced.dim = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > 0.8)[1]
  desvd.beta = svd.tmp$v[,1:reduced.dim] %*% diag(1/svd.tmp$d[1:reduced.dim])  %*% t(svd.tmp$u[,1:reduced.dim]) %*% sim.image[,3]
  
  est.desvd.entire = eval.mat.grid %*% desvd.beta
  mse.desvd.entire = mean((est.desvd.entire - true.image[,3])[-which(is.na(est.tps))]^2)
  
  desvd.beta.one = svd.tmp$v[,1:reduced.dim] %*% diag(1/svd.tmp$d[1:reduced.dim])  %*% t(svd.tmp$u[,1:reduced.dim]) %*% rep(1, dim(sim.image)[1])
  est.desvd.one = eval.mat.grid %*% desvd.beta.one
  est.desvd.rescaled = est.desvd.entire / est.desvd.one
  mse.desvd.rescaled = mean((est.desvd.rescaled - true.image[,3])[-which(is.na(est.tps))]^2)
  
  # SVD on design matrix (using entire bases) + orthogonalization + thsd'ing -----------------------
  eval.mat.rearranged.ortho = eval.mat.rearranged
  for(i in 1:3){
    proj.dat = eval.mat.rearranged.ortho[,unlist(basis.cidx[1:i])]
    proj.mat = proj.dat %*% solve(t(proj.dat) %*% proj.dat) %*% t(proj.dat)
    resid.mat = diag(1, dim(proj.mat)[1]) - proj.mat
    eval.mat.rearranged.ortho[, basis.cidx[[i+1]] ] = resid.mat %*% eval.mat.rearranged[, basis.cidx[[i+1]] ]
  }
  
  svd.tmp = svd(eval.mat.rearranged.ortho)
  reduced.dim = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > 0.8)[1]
  desvd.beta.ortho = svd.tmp$v[,1:reduced.dim] %*% diag(1/svd.tmp$d[1:reduced.dim])  %*% t(svd.tmp$u[,1:reduced.dim]) %*% sim.image[,3]
  
  desvd.beta.ortho.thsd = desvd.beta.ortho
  for(i in 4:1){
    desvd.beta.ortho.thsd[basis.cidx[[i]]] = ebayesthresh(desvd.beta.ortho[basis.cidx[[i]]])
  }
  
  eval.mat.grid.ortho = eval.mat.grid
  for(i in 1:3){
    proj.dat = eval.mat.grid.ortho[,unlist(basis.cidx[1:i])]
    proj.mat = proj.dat %*% solve(t(proj.dat) %*% proj.dat) %*% t(proj.dat)
    resid.mat = diag(1, dim(proj.mat)[1]) - proj.mat
    eval.mat.grid.ortho[, basis.cidx[[i+1]]] = resid.mat %*% eval.mat.grid[, basis.cidx[[i+1]]]
  }
  
  est.desvd.entire.ortho = eval.mat.grid.ortho %*% desvd.beta.ortho
  est.desvd.entire.ortho.thsd = eval.mat.grid.ortho %*% desvd.beta.ortho.thsd
  
  mse.desvd.entire.ortho = mean((est.desvd.entire.ortho - true.image[,3])[-which(is.na(est.tps))]^2)
  mse.desvd.entire.ortho.thsd = mean((est.desvd.entire.ortho.thsd - true.image[,3])[-which(is.na(est.tps))]^2)
  

  
  # SVD on design matrix + Wendland's multiscale approach + orthogonalization--------------------
  beta.desvd.wendland.list = list()
  est.desvd.wendland.list = list()
  resid.dat.list = list()
  
  prop.var = 0.8
  
  for(level.now in 4:1){
    
    if(level.now == 4){
      target.now = resid.dat.list[[level.now]] = sim.image # only for the coarsest level
    }else{
      target.now = resid.dat.list[[level.now]] = cbind(resid.dat.list[[level.now+1]][,1:2], resid.dat.list[[level.now+1]][,3] - est.desvd.wendland.list[[level.now+1]])  
    }
    
    obs.crspd.basis.selected.idx = vector()
    
    if(curve.fitting){
      # for basis with non-linear direction
      for(splat.idx in splat.selected[[level.now]]){
        # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
        idx.tmp = splat.collection.total[[ splat.idx ]]$idx
        
        trans.tmp = GetDeviationFromCurve(target.now[idx.tmp, 1:2], splat.collection.total.2D[[ splat.idx ]]$polygon)
        dev.tmp = sqrt(trans.tmp$lambda^2 + trans.tmp$dev^2)
        
        obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
      }
    }else{
      # for basis with linear direction
      for(splat.idx in splat.selected[[level.now]]){
        # level.now의 각 basis에 대응되는 obs 중 principal curve의 중심과 거리가 가장 가까운 obs를 찾아냄
        idx.tmp = splat.collection.total[[ splat.idx ]]$idx
        
        trans.tmp = TranslateData(target.now[idx.tmp, 1:2], center = -splat.collection.total.2D[[splat.idx]]$center)
        # dev.tmp = sqrt(sum(trans.tmp^2))
        dev.tmp = apply(trans.tmp, 1, function(x) sum(x^2))
        
        obs.crspd.basis.selected.idx[length(obs.crspd.basis.selected.idx)+1] = idx.tmp[which.min(dev.tmp)]
      }
    }
    
    obs.idx.tmp = unique(obs.crspd.basis.selected.idx)
    
    design.matrix.tmp = eval.mat.obs.ortho[ obs.idx.tmp, splat.selected[[level.now]] ]
    
    
    svd.tmp = svd(design.matrix.tmp)
    reduced.dim = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > 0.8)[1]
    beta.desvd.wendland = svd.tmp$v[,1:reduced.dim] %*% diag(1/svd.tmp$d[1:reduced.dim])  %*% t(svd.tmp$u[,1:reduced.dim]) %*% target.now[,3]
    est.desvd.wendland = design.matrix.tmp %*% beta.desvd.wendland
    
    
    beta.desvd.wendland.list[[level.now]] = beta.desvd.wendland
    est.desvd.wendland.list[[level.now]] = est.desvd.wendland
    
  }
  
  beta.desvd.wendland.list2 = beta.desvd.wendland.list
  for(i in 4:1){
    beta.desvd.wendland.list2[[i]] = ebayesthresh(beta.desvd.wendland.list[[i]])
  }
  
  
  est.desvd.wendland.ortho  = eval.mat.grid.ortho %*% unlist(rev(beta.desvd.wendland.list))
  est.desvd.wendland.ortho.thsd  = eval.mat.grid.ortho %*% unlist(rev(beta.desvd.wendland.list2))
  
  mse.desvd.wendland.ortho = mean((est.desvd.wendland.ortho - true.image[,3])[-which(is.na(est.tps))]^2)
  mse.desvd.wendland.ortho.thsd = mean((est.desvd.wendland.ortho.thsd - true.image[,3])[-which(is.na(est.tps))]^2)
  
  
  # returning list ----------------------
  result.list = list()
  result.list$sim.image = sim.image
  result.list$complete.image = complete.image
  result.list$true.image = true.image
  result.list$missing.idx = missing.idx
  
  result.list$splat.collection.total = splat.collection.total
  result.list$splat.collection.total.2D = splat.collection.total.2D
  
  result.list$eval.mat.obs = eval.mat.obs
  result.list$eval.mat.grid = eval.mat.grid
  
  result.list$splat.selected = splat.selected
  
  result.list$est.tps = est.tps
  result.list$mse.tps = mse.tps
  
  result.list$est.pcr.entire = est.pcr.entire
  result.list$mse.pcr.entire = mse.pcr.entire
  
  result.list$est.pcr.entire.ortho = est.pcr.entire.ortho
  result.list$mse.pcr.entire.ortho = mse.pcr.entire.ortho
  result.list$est.pcr.entire.ortho.thsd = est.pcr.entire.ortho.thsd
  result.list$mse.pcr.entire.ortho.thsd = mse.pcr.entire.ortho.thsd
  
  result.list$est.desvd.entire = est.desvd.entire
  result.list$mse.desvd.entire = mse.desvd.entire
  result.list$est.desvd.rescaled = est.desvd.rescaled
  result.list$mse.desvd.rescaled = mse.desvd.rescaled
  
  result.list$est.desvd.entire.ortho = est.desvd.entire.ortho
  result.list$mse.desvd.entire.ortho = mse.desvd.entire.ortho
  result.list$est.desvd.entire.ortho.thsd = est.desvd.entire.ortho.thsd
  result.list$mse.desvd.entire.ortho.thsd = mse.desvd.entire.ortho.thsd
  
  result.list$est.desvd.wendland.ortho = est.desvd.wendland.ortho
  result.list$mse.desvd.wendland.ortho = mse.desvd.wendland.ortho
  result.list$est.desvd.wendland.ortho.thsd = est.desvd.wendland.ortho.thsd
  result.list$mse.desvd.wendland.ortho.thsd = mse.desvd.wendland.ortho.thsd
  
  simulation.result[[iter]] = result.list
  
  
  saveRDS(simulation.result, file=paste0("sim1007_zflat_greedy_curvedonly_phi", ".rds"))
  
  if(iter == iter.num){
    time_end = Sys.time()
    cat("End time:", as.character(time_end), "\n")
    cat("(total elapsed time: ", time_end - time_start, units(time_end - time_start), ")\n\n")
  }

}