# nbhd 미리 설정 잘하자.

if(!require(fields)){
  install.packages("fields")
}
library(fields)

if(!require(RANN)){
  install.packages("RANN")
}
library(RANN)

if(!require(igraph)){
  install.packages("igraph")
}
library(igraph)

if(!require(princurve)){
  install.packages("princurve")
}
library(princurve)

if(!require(foreach)){
  install.packages("foreach")
}
library(foreach)

if(!require(doParallel)){
  install.packages("doParallel")
}
library(doParallel)

if(!require(Rcpp)){
  install.packages("Rcpp")
}
library(Rcpp)

if(!require(RcppArmadillo)){
  install.packages("RcppArmadillo")
}
library(RcppArmadillo)

if(!require(EbayesThresh)){
  install.packages("EbayesThresh")
}
library(EbayesThresh)


registerDoParallel(cores=detectCores()-1)

nbhd.num = 8

splat.unit = list(idx=NULL, center=NULL, normal.vector=NULL, major.axis=NULL, minor.axis=NULL, splat.left=NULL, splat.right=NULL, idx2=NULL, cov.2d=NULL, area.2d=NULL, ev.ratio.2d = NULL, activated=TRUE)
# cov가 3차원, 2차원 혼용되고 있지 않은지 체크할 필요 있음


GenerateTestShape = function(x, form="sine"){
  if(form=="line"){
    return(0.6*x + 0.2)
  }else if(form=="sine"){
    return(sin(2*pi*x)/4 + 0.5)  
  }else if(form=="sine5"){
    return(sin(5*2*pi*x)/4 + 0.5)  
  }else if(form == "exp"){
    return(exp(3*x)/20-1/20)
  }else if(form == "quadratic" | form == "quad"){
    return(2 * (x - 0.5)^2 + 0.25)
  }else{
    return(0.05 / abs(x - 0.5) + 0.5)
  }
}

GenerateTestImage = function(x.coord, y.coord, direction.width=1/60, noise.sd=0.05, form="sine", z.flux = FALSE){
  x.len = length(x.coord)
  y.len = length(y.coord)
  
  if(x.len != y.len){
    warning("lengths of x.coord and y.coord should be same.")
  }
  
  
  if(form == "inv"){
    return.vec.nf = exp(-(abs(y.coord-0.5) + 0.5 - GenerateTestShape(x.coord, form))^2 / (2*direction.width^2))
  }else if(form == "circle"){
    return.vec.nf = exp(-((x.coord-0.5)^2 + (y.coord-0.5)^2 - 1/16)^2 / (0.5*direction.width^2))
  }else if(form == "cross"){
    vec.nf.tmp1 = exp(-(x.coord - 0.5)^2 / (direction.width^2))
    vec.nf.tmp2 = exp(-(y.coord - 0.5)^2 / (direction.width^2))
    return.vec.nf = apply(cbind(vec.nf.tmp1, vec.nf.tmp2), 1, max)
  }else if(form == "phi"){
    vec.nf.tmp1 = exp(-((x.coord-0.5)^2 + (y.coord-0.5)^2 - 1/9)^2 / (0.5*direction.width^2))
    vec.nf.tmp2 = exp(-(x.coord - GenerateTestShape(y.coord, "line"))^2 / (2*direction.width^2))
    return.vec.nf = apply(cbind(vec.nf.tmp1, vec.nf.tmp2), 1, max)
  }else{
    return.vec.nf = exp(-(y.coord - GenerateTestShape(x.coord, form))^2 / (2*direction.width^2))
  }
  
  if(z.flux){
    return.vec.nf = return.vec.nf  * (0.3 * cos(x.coord * (2 * pi)) + 1)
    # (0.5 * sin(x.coord * (2 * pi)) + 0.6 + ifelse(x.coord<0.5,-0.2,0.4))
    # (sin(x.coord * (2 * pi)) + ifelse(x.coord<0.5,0,1) + 0.5)*2/3
    # (16*ifelse(x.coord < 0.5, -1, 1) * (x.coord - ifelse(x.coord < 0.5, 0.25, 0.75))^2 + ifelse(x.coord < 0.5, 1, 0))    
  }
  
  if(noise.sd == 0){
    noise.vec = 0
  }else{
    noise.vec = rnorm(x.len, sd = noise.sd)  
  }
  return.vec = return.vec.nf + noise.vec
  
  return(return.vec)
}
# https://math.stackexchange.com/questions/1956699/getting-a-transformation-matrix-from-a-normal-vector

RotateNorm = function(data, normal.vector, inverse=TRUE){
  n1 = normal.vector[1]
  n2 = normal.vector[2]
  n3 = normal.vector[3]
  
  u1 = c(n2/sqrt(n1^2 + n2^2), -n1/sqrt(n1^2 + n2^2), 0)
  u2 = c(n1*n3/sqrt(n1^2 + n2^2), n2*n3/sqrt(n1^2 + n2^2), -sqrt(n1^2 + n2^2))
  u3 = c(n1, n2, n3)
  
  rotMat = rbind(u1, u2, u3)
  
  if(inverse){ # rotate data to be on a plane containing the origin whose normal vector is normal.vector
    return(data %*% rotMat)
  }else{ # make normal.vector be perpendicular to the x,y-plane
    return(data %*% t(rotMat))
  }
  
}


# RotateNorm = function(data, normal.vector, inverse=TRUE){
#   n1 = normal.vector[1]
#   n2 = normal.vector[2]
#   n3 = normal.vector[3]
#   
#   u1 = c(n2/sqrt(n1^2 + n2^2), -n1/sqrt(n1^2 + n2^2), 0)
#   u2 = c(n1*n3/sqrt(n1^2 + n2^2), n2*n3/sqrt(n1^2 + n2^2), -sqrt(n1^2 + n2^2))
#   u3 = c(n1, n2, n3)
#   
#   rotMat = rbind(u1, u2, u3)
#   
#   center = apply(data, 2, mean)
#   data.centered = t(apply(data, 1, function(x) x - center))
#   
#   if(inverse){ # rotate data to be on a plane containing the origin whose normal vector is normal.vector
#     return(center + (data.centered %*% rotMat))
#   }else{ # make normal.vector be perpendicular to the x,y-plane
#     return(center + (data.centered %*% t(rotMat)))
#   }
#   
# }
# 


MakeEllipse2D = function(major.axis, minor.axis, resol=100){
  major.len = sqrt(sum(major.axis^2))
  minor.len = sqrt(sum(minor.axis^2))
  
  x.points = seq(-major.len, major.len, length.out = resol)
  y.points = minor.len * sqrt(1 - (x.points/major.len)^2)
  
  return( rbind(cbind(x.points, y.points), cbind(rev(x.points), -rev(y.points))) )
}

# plot(MakeEllipse2D(major.axis, minor.axis), type='l')

TranslateData = function(data, center=c(0,0,0)){
  
  if(is.vector(data)){
    return(matrix(data + center, ncol=length(center)))
  }else{
    return(t(apply(data, 1, function(x) x+center))  )
  }
}

GetAxesScaler = function(data.target, center, major.axis, minor.axis){
  
  major.len = sqrt(sum(major.axis^2))
  minor.len = sqrt(sum(minor.axis^2))
  
  major.axis.unit = major.axis / major.len
  minor.axis.unit = minor.axis / minor.len
  
  data.major = TranslateData(data.target, -center) %*% major.axis.unit
  data.minor = TranslateData(data.target, -center) %*% minor.axis.unit
  
  dev.scaled = sqrt((data.major / major.len)^2 + (data.minor / minor.len)^2)
  return(max(dev.scaled))
}


MakeSplat = function(idx, data){
  splat.tmp = splat.unit
  
  nbhd.num = 8
  
  if(length(idx)==1){
    center = data[idx, ]
    nbhd.idx = nn2(data, matrix(center, ncol=3), k=nbhd.num)$nn.idx
    idx2 = union(idx, nbhd.idx)
    
    svd.tmp = svd(cov(data[idx2,]))
    major.axis = svd.tmp$v[,1]
    minor.axis = svd.tmp$v[,2]
    normal.vector = svd.tmp$v[,3]
    
  }else if(length(idx) < nbhd.num){
    center = apply(data[idx,], 2, mean)
    nbhd.idx = nn2(data, matrix(center, ncol=3), k=nbhd.num)$nn.idx
    idx2 = union(idx, nbhd.idx)
    
    svd.tmp = svd(cov(data[idx2, ]))
    major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1] ### need to check
    minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2] ### need to check
    normal.vector = svd.tmp$v[,3]
    
  }else{
    center = apply(data[idx,], 2, mean)
    idx2 = idx
    
    svd.tmp = svd(cov(data[idx2, ]))
    major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1]
    minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2]
    normal.vector = svd.tmp$v[,3]
    
  }
  
  cov.3d.tmp = cov(data[idx2,1:3])
  cov.2d.tmp = cov.3d.tmp[1:2,1:2]
  
  idx.aux = setdiff(idx2, idx)
  if(length(idx.aux)==0){
    axes.scaler = GetAxesScaler(data[idx,], center, major.axis, minor.axis)
  }else{
    idx.aux.min = idx.aux[ which.min( mahalanobis(data[idx.aux,], center = center, cov = cov.3d.tmp) ) ]
    axes.scaler = GetAxesScaler(data[union(idx, idx.aux.min),], center, major.axis, minor.axis)
  }
  # axes.scaler = axes.scaler * 0.5 # modification
  major.axis = axes.scaler * major.axis
  minor.axis = axes.scaler * minor.axis
  
  
  
  splat.tmp$idx = idx
  splat.tmp$center = center
  splat.tmp$normal.vector = normal.vector * sign(normal.vector[3])
  splat.tmp$major.axis = major.axis * sign(major.axis[3])
  splat.tmp$minor.axis = minor.axis * sign(minor.axis[3])
  splat.tmp$idx2 = idx2
  splat.tmp$cov.2d = cov.2d.tmp
  svd.tmp = svd(cov.2d.tmp)
  splat.tmp$ev.ratio.2d = sqrt(svd.tmp$d[2] / svd.tmp$d[1])
  
  return(splat.tmp)
}

# idx2 대신 splat 의 중심에서 가장 가까운 데이터 하나까지만 포함하도록 장단축 길이 조정하게 바꿈
# MakeSplat = function(idx, data){
#   splat.tmp = splat.unit
#   
#   nbhd.num = 8
#   
#   if(length(idx)==1){
#     
#     center = data[idx, ]
#     nbhd.idx = nn2(data, matrix(center, ncol=3), k=nbhd.num)$nn.idx
#     idx2 = union(idx, nbhd.idx)
#     
#     svd.tmp = svd(cov(data[idx2,]))
#     major.axis = svd.tmp$v[,1]
#     minor.axis = svd.tmp$v[,2]
#     normal.vector = svd.tmp$v[,3]
#     
#     axes.scaler = GetAxesScaler(data[nbhd.idx,], center, major.axis, minor.axis)
#     major.axis = axes.scaler * major.axis
#     minor.axis = axes.scaler * minor.axis
#     
#   }else if(length(idx)==2){
#     center = apply(data[idx,], 2, mean)
#     
#     # nbhd.idx1 = nn2(data, matrix(data[idx[1],], ncol=3), k=nbhd.num)$nn.idx
#     # nbhd.idx2 = nn2(data, matrix(data[idx[2],], ncol=3), k=nbhd.num)$nn.idx
#     # nbhd.idx = union(nbhd.idx1, nbhd.idx2)
#     nbhd.idx = nn2(data, matrix(center, ncol=3), k=nbhd.num)$nn.idx
#     idx2 = union(idx, nbhd.idx)
#     
#     svd.tmp = svd(cov(data[idx2, ]))
#     major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1] ### need to check
#     minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2] ### need to check
#     normal.vector = svd.tmp$v[,3]
#     
#     axes.scaler = GetAxesScaler(data[idx2,], center, major.axis, minor.axis)
#     major.axis = axes.scaler * major.axis
#     minor.axis = axes.scaler * minor.axis
#     
#   }else{
#     center = apply(data[idx,], 2, mean)
#     nbhd.idx = nn2(data, matrix(center, ncol=3), k=nbhd.num)$nn.idx
#     idx2 = union(idx, nbhd.idx)
#     
#     svd.tmp = svd(cov(data[idx2, ]))
#     major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1]
#     minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2]
#     normal.vector = svd.tmp$v[,3]
#     
#     axes.scaler = GetAxesScaler(data[idx2,], center, major.axis, minor.axis)
#     major.axis = axes.scaler * major.axis
#     minor.axis = axes.scaler * minor.axis
#   }
#   
#   splat.tmp$idx = idx
#   splat.tmp$center = center
#   splat.tmp$normal.vector = normal.vector * sign(normal.vector[3])
#   splat.tmp$major.axis = major.axis * sign(major.axis[3])
#   splat.tmp$minor.axis = minor.axis * sign(minor.axis[3])
#   splat.tmp$idx2 = idx2
#   
#   return(splat.tmp)
# }


ProjectWithNormal = function(data, normal.vector){
  x = data
  n = normal.vector
  projection.matrix = diag(length(normal.vector)) -  n %*% t(n) / sum(n^2)
  
  return(x %*% t(projection.matrix))
  
}

MergeSplat = function(splat.left, splat.right, data, error.metric = "L2"){
  splat.merge = splat.unit
  idx.tmp = union(splat.left$idx, splat.right$idx)
  
  nbhd.num = 8
  
  if(error.metric == "L2"){
    splat.merge = MakeSplat(idx.tmp, data)  
  }else if(error.metric == "Cohen-Steiner"){
    major.axis.left = splat.left$major.axis * sign(splat.left$major.axis[3])
    minor.axis.left = splat.left$minor.axis * sign(splat.left$minor.axis[3])
    major.len.left = sqrt(sum(splat.left$major.axis^2))
    minor.len.left = sqrt(sum(splat.left$minor.axis^2))
    splat.area.left = pi * major.len.left * minor.len.left
    splat.center.left = splat.left$center
    normal.vector.left = splat.left$normal.vector * sign(splat.left$normal.vector[3])
    
    major.axis.right = splat.right$major.axis * sign(splat.right$major.axis[3])
    minor.axis.right = splat.right$minor.axis * sign(splat.right$minor.axis[3])
    major.len.right = sqrt(sum(splat.right$major.axis^2))
    minor.len.right = sqrt(sum(splat.right$minor.axis^2))
    splat.area.right = pi * major.len.right * minor.len.right
    splat.center.right = splat.right$center
    normal.vector.right = splat.right$normal.vector * sign(splat.right$normal.vector[3])
    
    normal.vector.merge = (splat.area.left * normal.vector.left + splat.area.right * normal.vector.right) / (splat.area.left + splat.area.right)
    splat.center.merge = (splat.area.left * splat.center.left + splat.area.right * splat.center.right) / (splat.area.left + splat.area.right)
    
    resol.scaler = 50 / min(c(splat.area.left, splat.area.right))
    splat.points.left = TranslateData(RotateNorm(cbind(MakeEllipse2D(major.axis.left, minor.axis.left,resol=resol.scaler*splat.area.left), 0), normal.vector.left), splat.center.left)
    splat.points.right = TranslateData(RotateNorm(cbind(MakeEllipse2D(major.axis.right, minor.axis.right,resol=resol.scaler*splat.area.right), 0), normal.vector.right), splat.center.right)
    
    splat.points.merge = rbind(ProjectWithNormal(splat.points.left, normal.vector.merge),
                               ProjectWithNormal(splat.points.right, normal.vector.merge))
    
    splat.points.new = RotateNorm(splat.points.merge, normal.vector.merge, inverse = FALSE)[,-3]
    svd.tmp = svd(cov(splat.points.new))
    major.axis.tmp = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1] ### need to check
    minor.axis.tmp = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2] ### need to check
    
    axes.scaler = GetAxesScaler(splat.points.new, apply(splat.points.new, 2, mean), major.axis.tmp, minor.axis.tmp)
    major.axis.tmp = axes.scaler * major.axis.tmp
    minor.axis.tmp = axes.scaler * minor.axis.tmp
    
    major.axis.merge = RotateNorm(c(major.axis.tmp,0), normal.vector.merge, inverse = TRUE)
    major.axis.merge = major.axis.merge * sign(major.axis.merge[3])
    minor.axis.merge = RotateNorm(c(minor.axis.tmp,0), normal.vector.merge, inverse = TRUE)
    minor.axis.merge = minor.axis.merge * sign(minor.axis.merge[3])
    
    
    splat.merge$idx = idx.tmp
    splat.merge$center = as.vector(splat.center.merge)
    splat.merge$normal.vector = as.vector(normal.vector.merge)
    splat.merge$major.axis = as.vector(major.axis.merge)
    splat.merge$minor.axis = as.vector(minor.axis.merge)
   
    nbhd.idx = nn2(data, matrix(splat.merge$center, ncol=3), k=nbhd.num)$nn.idx
    idx2 = union(idx.tmp, nbhd.idx)
    splat.merge$idx2 = idx2 
  }
  
  splat.merge$splat.left = splat.left
  splat.merge$splat.right = splat.right
  
  return(splat.merge)
  
}

MergeSplatCollection = function(splat.collection, splat.num.target, data, error.metric = "L2", verbose = FALSE){
  splat.num = length(splat.collection)
  if(verbose) cat("From", splat.num, "...\n")
  while(splat.num > splat.num.target){
    edge.tmp = MakeEdges(splat.collection, data = data, error.metric=error.metric)
    splat.graph = graph(edges = edge.tmp$edge.list, n=length(splat.collection), directed = FALSE)
    edge_attr(splat.graph, "error.metric", index = E(splat.graph)) = edge.tmp$edge.value
    vertex_attr(splat.graph, "splat", index = V(splat.graph)) = splat.collection
    
    eidx.min = which.min(edge_attr(splat.graph, "error.metric"))
    vidx.left = edge.tmp$edge.list[2*eidx.min - 1]
    vidx.right = edge.tmp$edge.list[2*eidx.min]
    splat.merge = MergeSplat(splat.collection[[vidx.left]], splat.collection[[vidx.right]], data, error.metric)
    
    splat.collection = splat.collection[-c(vidx.left, vidx.right)]
    splat.collection[[length(splat.collection)+1]] = splat.merge
    
    splat.num = length(splat.collection)
    if(verbose) cat(splat.num, "\n")
  }
  
  return(splat.collection)
}


MergeSplatCollection2 = function(splat.collection, splat.num.target, data, error.metric = "L2", verbose = FALSE){
  
  nbhd.num = 8
  splat.num = length(splat.collection)
  center.mat = GetSplatCenters(splat.collection)
  # cost.mat doesn't need to be symmetric
  # cost.mat[,j] consists of merging cost for nearest nhbd for the j-th splat
  
  # Make initial cost.mat
  cost.mat = matrix(Inf, nrow = splat.num, ncol = splat.num)
  for(splat.idx in 1:splat.num){
    # cost.mat의 각 칼럼에 nn merging cost를 채워넣는다
    nbhd.idx = nn2(center.mat, matrix(center.mat[splat.idx,], ncol=3), k=nbhd.num+1)$nn.idx[-1]
    for(j in 1:nbhd.num){
      cost.mat[nbhd.idx[j], splat.idx] = CalculateErrorMetric(splat.collection[[splat.idx]], splat.collection[[ nbhd.idx[j] ]], data, error.metric)
    }
  }
  
  # 본격적 merging 시작
  if(verbose) cat("From", splat.num, "...\n")  
  splat.collection.tmp = splat.collection
  while(length(splat.collection.tmp) > splat.num.target){
    merge.result = MergeSplatOneStep(splat.collection.tmp, data, error.metric, cost.mat, center.mat)
    
    splat.collection.tmp = merge.result$splat.collection
    cost.mat = merge.result$cost.mat
    center.mat = merge.result$center.mat
    
    if(verbose) cat(length(splat.collection.tmp), "\n")
  }
  
  return(splat.collection.tmp)
}





MergeSplatOneStep = function(splat.collection, data, error.metric = "L2", cost.mat, center.mat){
  # rough version (nearest nbhd is found only for the initial splat collection)
  
  # Find which pair of splats to merge
  min.arr.idx = which(cost.mat == min(cost.mat), arr.ind = TRUE)
  if(dim(min.arr.idx)[1] > 1){
    chosen.idx = sample(1:dim(min.arr.idx)[1], size  = 1)
    merging.idx = min.arr.idx[chosen.idx, ]
  }else{
    merging.idx = min.arr.idx
  }
  
  
  # Merge the splats
  splat.merge = MergeSplat(splat.collection[[merging.idx[1]]], splat.collection[[merging.idx[2]]], data, error.metric)
  
  # new splat collection
  splat.collection.new = c(splat.collection[-merging.idx], list(splat.merge))
  splat.num.new = length(splat.collection.new)
  
  # new center.mat
  center.mat.new = rbind(center.mat[-merging.idx,], splat.collection.new[[splat.num.new]]$center)
  
  # new cost.mat
  cost.mat.new = cost.mat[-merging.idx, -merging.idx]
  cost.mat.new = rbind(cbind(cost.mat.new, Inf), Inf)
  
  # nbhd.idx of newly merged splat
  nbhd.idx = nn2(center.mat.new, matrix(center.mat.new[splat.num.new,], ncol=3), k=nbhd.num+1)$nn.idx[-1]
  
  # Update nn merging cost for the newly merged splat and update the existing costs
  for(i in 1:length(nbhd.idx)){
    # for the newly merged splat...
    cost.mat.new[nbhd.idx[i], splat.num.new] = CalculateErrorMetric(splat.collection.new[[splat.num.new]], splat.collection.new[[ nbhd.idx[i] ]], data, error.metric)
    
    # to update the existing costs...
    # splat pair with the largest cost will be replaced if the new cost is smaller
    if(sum(cost.mat.new[,nbhd.idx[i]] < cost.mat.new[nbhd.idx[i], splat.num.new]) < nbhd.num){
      rmv.idx = which.max(cost.mat.new[which(!is.infinite(cost.mat.new[,nbhd.idx[i]])), nbhd.idx[i]])
      rmv.idx = which(!is.infinite(cost.mat.new[,nbhd.idx[i]]))[rmv.idx]
      cost.mat.new[rmv.idx, nbhd.idx[i]] = Inf
      cost.mat.new[splat.num.new, nbhd.idx[i]] = cost.mat.new[nbhd.idx[i], splat.num.new]
    }
    
  }
  
  return(list(splat.collection = splat.collection.new, cost.mat = cost.mat.new, center.mat = center.mat.new))
}




SplitSplat = function(splat.merge, data){
  
  splat.left = splat.merge$splat.left
  splat.right = splat.merge$splat.right
  
  return(list(splat.left = splat.left, splat.right = splat.right))
}


CalculatePointPlaneDistance = function(idx, splat, data){
  
  points = matrix(data[idx,], ncol=3)
  
  center = splat$center
  normal.vector = splat$normal.vector
  
  numerator = abs( TranslateData(points, center) %*% normal.vector )
  denom = sqrt(sum(normal.vector^2))
  
  return(as.vector(numerator / denom))
}


# CalculateErrorMetric = function(splat.merge, data){
#   splat.left = splat.merge$splat.left
#   splat.right = splat.merge$splat.right
#   
#   dist.between.centers = sqrt(sum((splat.left$center - splat.right$center)^2))
#   
#   error.tmp = sum(CalculatePointPlaneDistance(splat.merge$idx, splat.merge, data)^2)
#   
#   dist.between.centers * error.tmp
# }


CalculateErrorMetric = function(splat.left, splat.right, data, error.metric="L2"){
  
  dist.between.centers = sqrt(sum((splat.left$center - splat.right$center)^2))
  
  if(error.metric == "L2"){
    major.len.left = sqrt(sum(splat.left$major.axis^2))
    minor.len.left = sqrt(sum(splat.left$minor.axis^2))
    splat.area.left = pi * major.len.left * minor.len.left
    
    major.len.right = sqrt(sum(splat.right$major.axis^2))
    minor.len.right = sqrt(sum(splat.right$minor.axis^2))
    splat.area.right = pi * major.len.right * minor.len.right
    
    splat.merge = MergeSplat(splat.left, splat.right, data)
    error.tmp = sum(CalculatePointPlaneDistance(splat.merge$idx, splat.merge, data)^2)
    
    return(dist.between.centers * error.tmp * splat.area.left * splat.area.right)
    
  }else if(error.metric == "Cohen-Steiner"){
    major.len.left = sqrt(sum(splat.left$major.axis^2))
    minor.len.left = sqrt(sum(splat.left$minor.axis^2))
    splat.area.left = pi * major.len.left * minor.len.left
    
    major.len.right = sqrt(sum(splat.right$major.axis^2))
    minor.len.right = sqrt(sum(splat.right$minor.axis^2))
    splat.area.right = pi * major.len.right * minor.len.right
    
    normal.vector.left = splat.left$normal.vector
    normal.vector.right = splat.right$normal.vector
    
    return(dist.between.centers * splat.area.left * splat.area.right * sum((normal.vector.left - normal.vector.right)^2))
    
  }
  
}

# 아래의 L2 에서 면적 계산 반영하도록 바꿈
# CalculateErrorMetric = function(splat.left, splat.right, data, error.metric="L2"){
#   
#   dist.between.centers = sqrt(sum((splat.left$center - splat.right$center)^2))
#   
#   if(error.metric == "L2"){
#     splat.merge = MergeSplat(splat.left, splat.right, data)
#     error.tmp = sum(CalculatePointPlaneDistance(splat.merge$idx, splat.merge, data)^2)
#     
#     return(dist.between.centers * error.tmp)  
#     
#   }else if(error.metric == "Cohen-Steiner"){
#     major.len.left = sqrt(sum(splat.left$major.axis^2))
#     minor.len.left = sqrt(sum(splat.left$minor.axis^2))
#     splat.area.left = pi * major.len.left * minor.len.left
#     
#     major.len.right = sqrt(sum(splat.right$major.axis^2))
#     minor.len.right = sqrt(sum(splat.right$minor.axis^2))
#     splat.area.right = pi * major.len.right * minor.len.right
#     
#     normal.vector.left = splat.left$normal.vector
#     normal.vector.right = splat.right$normal.vector
#     
#     return(dist.between.centers * splat.area.left * splat.area.right * sum((normal.vector.left - normal.vector.right)^2))
#   
#     }
#   
# }


GetSplatCenters = function(splat.collection){
  splat.num = length(splat.collection)
  center.mat = matrix(0, ncol=length(splat.collection[[1]]$center), nrow=splat.num)
  
  for(i in 1:splat.num){
    center.mat[i,] = splat.collection[[i]]$center
  }
  
  return(center.mat)
}


MakeEdgeList = function(splat.collection, nbhd.num = 8){
  center.mat = GetSplatCenters(splat.collection)
  splat.num = length(splat.collection)
  
  edge.list.mat = matrix(0, ncol = splat.num, nrow = nbhd.num * 2)
  
  for(i in 1:splat.num){
    
    edge.list.mat[(1:nbhd.num)*2 - 1, i] = i
    edge.list.mat[(1:nbhd.num)*2, i] = nn2(center.mat, matrix(center.mat[i,], ncol=3), k=nbhd.num+1)$nn.idx[-1]
  }
  
  return(as.vector(edge.list.mat))
}



MakeEdgeListUnique = function(edge.list.mat){
  edge.list.vec = as.vector(edge.list.mat)
  
  # make 1st column always smaller than the 2nd column
  for(i in 1:(length(edge.list.vec)/2)){
    if(edge.list.vec[2*i-1] > edge.list.vec[2*i]){
      tmp = edge.list.vec[2*i]
      edge.list.vec[2*i] = edge.list.vec[2*i-1]
      edge.list.vec[2*i-1] = tmp
    }
  }
  
  # sort the 1st column in ascending order
  edge.list.mat.tmp = matrix(edge.list.vec, nrow=2)
  edge.list.mat.sorted = t(edge.list.mat.tmp[,order(edge.list.mat.tmp[1,])])
  
  # sort the 2nd column in ascending order
  range.llim = c(1,which(diff(edge.list.mat.sorted[,1])!=0)+1)
  range.rlim = c(which(diff(edge.list.mat.sorted[,1])!=0), dim(edge.list.mat.tmp)[2])
  for(i in 1:length(unique(edge.list.mat.sorted[,1]))){
    range.tmp = range.llim[i]:range.rlim[i]
    if(length(range.tmp)==1){
      next
    }
    tmp = edge.list.mat.sorted[range.tmp,]
    tmp = tmp[order(tmp[,2]),]
    edge.list.mat.sorted[range.tmp,] = tmp
  }
  
  # make the elements in the 2nd column unique
  edge.list.mat.unique = edge.list.mat.sorted[-(which(diff(edge.list.mat.sorted[,2])==0) + 1),]
  
  return(edge.list.mat.unique)
}

MakeEdges = function(splat.collection, nbhd.num = 8, data, error.metric = "L2"){
  center.mat = GetSplatCenters(splat.collection)
  splat.num = length(splat.collection)
  
  edge.list.mat = matrix(0, ncol = splat.num, nrow = nbhd.num * 2)
  edge.val.vec = numeric(nbhd.num * splat.num)
  
  for(i in 1:splat.num){
    edge.list.mat[(1:nbhd.num)*2 - 1, i] = i
    edge.list.mat[(1:nbhd.num)*2, i] = nn2(center.mat, matrix(center.mat[i,], ncol=3), k=nbhd.num+1)$nn.idx[-1]
    
    knn.idx = edge.list.mat[(1:nbhd.num)*2, i]
    error.tmp = numeric(nbhd.num)
    for(j in 1:nbhd.num){
      error.tmp[j] = CalculateErrorMetric(splat.collection[[i]], splat.collection[[ knn.idx[j] ]], data, error.metric)
    }
    
    edge.val.vec[(i-1)*nbhd.num + 1:nbhd.num] = error.tmp
  }
  
  return(list(edge.list = as.vector(edge.list.mat), edge.value = edge.val.vec))
}

MakeEdges = function(splat.collection, nbhd.num = 8, data, error.metric = "L2"){
  center.mat = GetSplatCenters(splat.collection)
  splat.num = length(splat.collection)
  
  edge.list.mat = matrix(0, ncol = splat.num, nrow = nbhd.num * 2)
  edge.val.vec = numeric(nbhd.num * splat.num)
  
  for(i in 1:splat.num){
    edge.list.mat[(1:nbhd.num)*2 - 1, i] = i
    edge.list.mat[(1:nbhd.num)*2, i] = nn2(center.mat, matrix(center.mat[i,], ncol=3), k=nbhd.num+1)$nn.idx[-1]
  }
  
  edge.list.mat.unique = MakeEdgeListUnique(edge.list.mat)
  
  edge.val.vec = vector()
  for(i in 1:dim(edge.list.mat.unique)[1]){
    edge.val.vec[i] = CalculateErrorMetric(splat.collection[[ edge.list.mat.unique[i,1] ]], splat.collection[[  edge.list.mat.unique[i,2] ]], data, error.metric)
  }
  
  return(list(edge.list = as.vector(t(edge.list.mat.unique)), edge.value = edge.val.vec))
}


# MakeSplat2D = function(idx, data){
#   splat.tmp = splat.unit
#   
#   nbhd.num = 8
#   
#   if(length(idx)==1){
#     
#     center = data[idx, 1:2]
#     nbhd.idx = nn2(data[,1:2], matrix(center, ncol=2), k=nbhd.num)$nn.idx
#     idx2 = union(idx, nbhd.idx)
#     
#     cov.tmp = cov(data[idx2,1:2])
#     svd.tmp = svd(cov.tmp)
#     major.axis = svd.tmp$v[,1]
#     minor.axis = svd.tmp$v[,2]
#     
#     axes.scaler = GetAxesScaler(data[idx2,1:2], center, major.axis, minor.axis)
#     major.axis = axes.scaler * major.axis
#     minor.axis = axes.scaler * minor.axis
#     
#   }else if(length(idx)==2){
#     center = apply(data[idx,1:2], 2, mean)
#     
#     # nbhd.idx1 = nn2(data, matrix(data[idx[1],], ncol=3), k=nbhd.num)$nn.idx
#     # nbhd.idx2 = nn2(data, matrix(data[idx[2],], ncol=3), k=nbhd.num)$nn.idx
#     # nbhd.idx = union(nbhd.idx1, nbhd.idx2)
#     nbhd.idx = nn2(data[,1:2], matrix(center, ncol=2), k=nbhd.num)$nn.idx
#     idx2 = union(idx, nbhd.idx)
#     
#     cov.tmp = cov(data[idx2, 1:2])
#     svd.tmp = svd(cov.tmp)
#     major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1] ### need to check
#     minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2] ### need to check
#     
#     axes.scaler = GetAxesScaler(data[idx2,1:2], center, major.axis, minor.axis)
#     major.axis = axes.scaler * major.axis
#     minor.axis = axes.scaler * minor.axis
#     
#   }else{
#     center = apply(data[idx,1:2], 2, mean)
#     nbhd.idx = nn2(data[,1:2], matrix(center, ncol=2), k=nbhd.num)$nn.idx
#     idx2 = union(idx, nbhd.idx)
#     
#     cov.tmp = cov(data[idx2, 1:2])
#     svd.tmp = svd(cov.tmp)
#     major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1]
#     minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2]
#     
#     axes.scaler = GetAxesScaler(data[idx2,1:2], center, major.axis, minor.axis)
#     major.axis = axes.scaler * major.axis
#     minor.axis = axes.scaler * minor.axis
#   }
#   
#   splat.tmp$idx = idx
#   splat.tmp$center = center
#   splat.tmp$major.axis = major.axis * sign(major.axis[2])
#   splat.tmp$minor.axis = minor.axis * sign(minor.axis[2])
#   splat.tmp$idx2 = idx2
#   splat.tmp$cov = cov.tmp
#   
#   return(splat.tmp)
# }

MakeSplat2D = function(splat, data){
  splat.tmp = splat
  idx2 = splat.tmp$idx2
  idx = splat.tmp$idx
  cov.2d = splat.tmp$cov.2d
  center = splat.tmp$center[1:2]
  
  nbhd.num = 8
  
  svd.2d = svd(cov.2d)
  major.axis = sqrt(svd.2d$d[1]) * svd.2d$v[,1]
  minor.axis = sqrt(svd.2d$d[2]) * svd.2d$v[,2]
  
  idx.aux = setdiff(idx2, idx)
  if(length(idx.aux)==0){
    axes.scaler = GetAxesScaler(data[idx,1:2], center, major.axis, minor.axis)
  }else{
    idx.aux.min = idx.aux[ which.min( mahalanobis(data[idx.aux,1:2], center = center, cov = cov.2d) ) ]
    axes.scaler = GetAxesScaler(data[union(idx, idx.aux.min),1:2], center, major.axis, minor.axis)
  }
  # axes.scaler = axes.scaler * 2 # modification
  major.axis = axes.scaler * major.axis
  minor.axis = axes.scaler * minor.axis
  
  # axes.scaler = GetAxesScaler(data[idx2,1:2], center, major.axis, minor.axis)
  # major.axis = axes.scaler * major.axis
  # minor.axis = axes.scaler * minor.axis
  
  splat.tmp$center = center
  splat.tmp$major.axis = major.axis * sign(major.axis[2])
  splat.tmp$minor.axis = minor.axis * sign(minor.axis[2])
  
  return(splat.tmp)
}
# idx2 대신 splat 의 중심에서 가장 가까운 데이터 하나까지만 포함하도록 장단축 길이 조정하게 바꿈
# MakeSplat2D = function(splat, data){
#   splat.tmp = splat
#   idx2 = splat.tmp$idx2
#   
#   nbhd.num = 8
#   
#   center = apply(data[idx2,1:2], 2, mean)
#   cov.tmp = cov(data[idx2, 1:2])
#   svd.tmp = svd(cov.tmp)
#   major.axis = sqrt(svd.tmp$d[1]) * svd.tmp$v[,1]
#   minor.axis = sqrt(svd.tmp$d[2]) * svd.tmp$v[,2]
#   
#   axes.scaler = GetAxesScaler(data[idx2,1:2], center, major.axis, minor.axis)
#   major.axis = axes.scaler * major.axis
#   minor.axis = axes.scaler * minor.axis
#   
#   splat.tmp$center = center
#   splat.tmp$major.axis = major.axis * sign(major.axis[2])
#   splat.tmp$minor.axis = minor.axis * sign(minor.axis[2])
#   splat.tmp$cov = cov.tmp
#   
#   return(splat.tmp)
# }


MakeSplatCollection2D = function(splat.collection, data, curve.fitting = FALSE){
  splat.collection.2D = list()
  splat.num = length(splat.collection)
  if(!curve.fitting){
    for(i in 1:splat.num){
      splat.collection.2D[[i]] = MakeSplat2D(splat.collection[[i]], data)
    }  
  }else{
    for(i in 1:splat.num){
      splat.collection.2D[[i]] = MakeSplatAlongCurve(splat.collection[[i]], data)
    }
  }
  return(splat.collection.2D)
}


EvalWendland = function(r){
  phi.out = r
  
  idx.tmp = which(abs(r) <= 1)
  
  phi.out[idx.tmp] = (1-r[idx.tmp])^6 * (3 + 18*r[idx.tmp] + 35*r[idx.tmp]^2) / 3
  phi.out[which(!abs(r) <= 1)] = 0
  
  return(phi.out)
}

EvalHaar = function(r){
  phi.out = r
  
  phi.out[which(abs(r) <= 1)] = 1
  phi.out[which(!abs(r) <= 1)] = 0
  
  return(phi.out)
}



EvaluateBasis = function(x, splat.collection.2D, data, curve.fitting = FALSE, kernel.scaling = 0.1){
  
  if(is.matrix(x)){
    if(ncol(x) == 3) x = x[,1:2]
  }else if(is.vector(x)){
    if(length(x) == 3) x = x[1:2]
  }
  
  splat.num = length(splat.collection.2D)
  
  r.mat = matrix(0, ncol=splat.num, nrow=length(x)/2)
  # ncol = length(x)/2 는 x가 x, y 좌표를 갖고 있음을 이용한 것
  
  for(splat.idx in 1:splat.num){
    # curve.fitting = TRUE라면 함수에 처음부터 변환된 공간에서의 중심점과 공분산이 들어와야 함
    splat.data.idx = splat.collection.2D[[splat.idx]]$idx
    splat.data.idx2 = splat.collection.2D[[splat.idx]]$idx2
    
    splat.data.idx.num = length(splat.data.idx)
    
    if(curve.fitting){
      x.transformed.tmp = GetDeviationFromCurve(x, splat.collection.2D[[splat.idx]]$polygon)
      data.transformed.tmp = GetDeviationFromCurve(data[,1:2], splat.collection.2D[[splat.idx]]$polygon)
      
      x.transformed = cbind(x.transformed.tmp$lambda, x.transformed.tmp$dev)
      data.transformed = cbind(data.transformed.tmp$lambda, data.transformed.tmp$dev)
      
      splat.center = apply(matrix(data.transformed[splat.data.idx,], ncol=2), 2, mean)
      splat.data.cov = cov(data.transformed[splat.data.idx2,])
      
      idx.aux = setdiff(splat.data.idx2, splat.data.idx)
      if(length(idx.aux) != 0){
        idx.aux.min = idx.aux[ which.min( mahalanobis(data.transformed[idx.aux,1:2], center = splat.center, cov = cov(data.transformed[splat.data.idx2,1:2])) ) ]
        splat.data.idx = union(splat.data.idx, idx.aux.min)
      }
      delta = 1 / max(sqrt(mahalanobis(x = data.transformed[splat.data.idx,], center = splat.center, cov = splat.data.cov)))
      
      r.tmp = delta * mahalanobis(x = x.transformed, center = splat.center, cov = splat.data.cov)
    }else{
      splat.center = splat.collection.2D[[splat.idx]]$center
      splat.data.cov = splat.collection.2D[[splat.idx]]$cov.2d
      
      idx.aux = setdiff(splat.data.idx2, splat.data.idx)
      if(length(idx.aux) != 0){
        idx.aux.min = idx.aux[ which.min( mahalanobis(data[idx.aux,1:2], center = splat.center, cov = splat.data.cov) ) ]
        splat.data.idx = union(splat.data.idx, idx.aux.min)
      }
      
      delta = 1 / max(sqrt(mahalanobis(x = data[splat.data.idx,1:2], center = splat.center, cov = splat.data.cov)))
      
      r.tmp = delta * mahalanobis(x = x, center = splat.center, cov = splat.data.cov)
    }
    
    r.mat[,splat.idx] = r.tmp * kernel.scaling
  }
  
  
  return( apply(r.mat, 2, function(x) EvalWendland(x)) )
}

# idx2 대신 splat 의 중심에서 가장 가까운 데이터 하나까지만 포함하도록 장단축 길이 조정하게 바꿈
# curve fitting 도 옵션으로 추가
# EvaluateBasis = function(x, splat.collection.2D, data){
#   
#   if(is.matrix(x)){
#     if(ncol(x) == 3) x = x[,1:2]
#   }else if(is.vector(x)){
#     if(length(x) == 3) x = x[1:2]
#   }
#   
#   splat.num = length(splat.collection.2D)
#   
#   r.mat = matrix(0, nrow=splat.num, ncol=length(x)/2)
#   # ncol = length(x)/2 는 x가 x, y 좌표를 갖고 있음을 이용한 것
#   
#   for(splat.idx in 1:splat.num){
#     
#     splat.data.idx2 = splat.collection.2D[[splat.idx]]$idx2
#     splat.data.cov = splat.collection.2D[[splat.idx]]$cov
#     
#     splat.center = splat.collection.2D[[splat.idx]]$center
#     
#     delta = 1 / max(sqrt(mahalanobis(x = data[splat.data.idx2, 1:2], center = splat.center, cov = splat.data.cov)))
#     r.tmp = delta * mahalanobis(x = x, center = splat.center, cov = splat.data.cov)
#     
#     r.mat[splat.idx,] = r.tmp
#   }
#   
#   
#   return( apply(r.mat, 2, function(x) EvalWendland(x)) )
# }


GetBasisMemberIdx = function(splat, data){
  splat.tmp = splat
  major.axis.tmp = splat.tmp$major.axis
  minor.axis.tmp = splat.tmp$minor.axis
  normal.vector.tmp = splat.tmp$normal.vector
  splat.center.tmp = splat.tmp$center
  
  splat.points.tmp.proj.centered = RotateNorm(cbind(MakeEllipse2D(major.axis.tmp, minor.axis.tmp,resol=100), 0), normal.vector.tmp)[,-3]
  points.tmp = splat.points.tmp.proj.centered
  point.idx = sample(1:dim(points.tmp)[1], size=3)
  
  mat.tmp = cbind(points.tmp[point.idx,1]^2, points.tmp[point.idx,2]^2, points.tmp[point.idx,1]*points.tmp[point.idx,2])
  
  data.tmp = TranslateData(data, center = -splat.center.tmp)
  data.tmp = cbind(data.tmp[,1]^2, data.tmp[,2]^2, data.tmp[,1]*data.tmp[,2])
  return( which(data.tmp %*% solve(mat.tmp, rep(-1,3)) + 1 >= 0) )
}


ProjectOnPlane = function(data, normal.vector, origin){
  x = data
  k.vec = TranslateData(-x, center=origin) %*% normal.vector / sum(normal.vector^2)
  return(x + k.vec %*% normal.vector)
}

# ConvertCurve2Polygon = function(principal.curve){
#   fit = principal.curve
#   point.seg = fit$s[fit$ord, ]
#   
#   idx.duplicated = which(apply(diff(point.seg), 1, function(x) sum(x^2))==0)
#   if(length(idx.duplicated) != 0){
#     point.seg = point.seg[-idx.duplicated,]  
#   }
#   
#   idx.first = 1
#   tmp = point.seg[idx.first,] - point.seg[idx.first+1,]
#   x.first = point.seg[idx.first,1] + sign(tmp[1])*5
#   slope1 = tmp[2] / tmp[1]
#   y.first = slope1 * (x.first - point.seg[idx.first,1]) + point.seg[idx.first,2]
#   
#   idx.last = dim(point.seg)[1]
#   tmp = point.seg[idx.last-1,] - point.seg[idx.last,]
#   x.last = point.seg[idx.last,1] + sign(-tmp[1])*5
#   slope2 = tmp[2] / tmp[1]
#   y.last = slope2 * (x.last - point.seg[idx.last, 1]) + point.seg[idx.last, 2]
#   
#   return( rbind(c(x.first, y.first), point.seg, c(x.last, y.last)) )
# }

ConvertCurve2Polygon = function(principal.curve){
  fit = principal.curve
  point.seg = fit$s[fit$ord, ]
  
  idx.duplicated = which(apply(diff(point.seg), 1, function(x) sum(x^2))==0)
  if(length(idx.duplicated) != 0){
    point.seg = point.seg[-idx.duplicated,]  
  }
  
  idx.first = 1
  tmp = point.seg[idx.first,] - point.seg[idx.first+1,]
  slope1 = tmp[2] / tmp[1]
  
  if(sign(tmp[1])<0){
    if(slope1>0){
      x.first = min(c((0-point.seg[idx.first,2]) / slope1 + point.seg[idx.first,1]-1e-5, 0))
    }else{
      x.first = min(c((1-point.seg[idx.first,2]) / slope1 + point.seg[idx.first,1]-1e-5, 0))
    }
  }else{
    if(slope1>0){
      x.first = max(c((1-point.seg[idx.first,2]) / slope1 + point.seg[idx.first,1]+1e-5, 1))
    }else{
      x.first = max(c((0-point.seg[idx.first,2]) / slope1 + point.seg[idx.first,1]+1e-5, 1))
    }
  }
  
  y.first = slope1 * (x.first - point.seg[idx.first,1]) + point.seg[idx.first,2]
  
  idx.last = dim(point.seg)[1]
  tmp = point.seg[idx.last-1,] - point.seg[idx.last,]
  slope2 = tmp[2] / tmp[1]
  
  if(sign(tmp[1])<0){
    if(slope2>0){
      x.last = max(c((1-point.seg[idx.last,2]) / slope2 + point.seg[idx.last,1]+1e-5, 1))
    }else{
      x.last = max(c((0-point.seg[idx.last,2]) / slope2 + point.seg[idx.last,1]+1e-5, 1))
    }
  }else{
    if(slope2>0){
      x.last = min(c((0-point.seg[idx.last,2]) / slope2 + point.seg[idx.last,1]-1e-5, 0))
    }else{
      x.last = min(c((1-point.seg[idx.last,2]) / slope2 + point.seg[idx.last,1]-1e-5, 0))
    }
  }
  
  y.last = slope2 * (x.last - point.seg[idx.last, 1]) + point.seg[idx.last, 2]
  
  return( rbind(c(x.first, y.first), point.seg, c(x.last, y.last)) )
}



# getDistToCurve = function(p, point.seg){
#   
#   point0 = p
#   dist.min = Inf
#   sign.min = 1
#   for(p.idx in 2:dim(point.seg)[1]){
#     point1 = point.seg[p.idx-1,]
#     point2 = point.seg[p.idx,]
#     
#     if(sum((point1 - point0)^2)*sum((point2 - point0)^2)==0){
#       dist.min = 0
#       break
#     }
#     
#     r.tmp = sum((point0 - point1) * (point2 - point1)) / sum((point2 - point1)^2)
#     
#     if(r.tmp < 0){
#       dist.tmp = sum((point0 - point1)^2)
#     }else if(r.tmp > 1){
#       dist.tmp = sum((point0 - point2)^2)
#     }else{
#       dist.tmp = sum((point0 - point1)^2) - r.tmp^2 * sum((point1 - point2)^2)
#     }
#     dist.tmp = sqrt(dist.tmp)
#     sign.tmp = sign((point0[2] - point1[2]) - ((point2[2] - point1[2])/(point2[1] - point1[1])*(point0[1]-point1[1])))
#     
#     if(dist.tmp < dist.min){
#       dist.min = dist.tmp; sign.min = sign.tmp
#     }
#     
#   }
#   return(dist.min * sign.min)
# }

GetPointDeviationFromCurve = function(p, point.seg){
  
  point0 = p
  dist.min = Inf
  sign.min = 1
  lambda.min = NULL
  lambda.vec = cumsum(c(0,sqrt(apply(diff(point.seg)^2,1,sum)))) - sqrt(sum(diff(point.seg)[1,]^2))
  
  for(seg.idx in 2:dim(point.seg)[1]){
    point1 = point.seg[seg.idx-1,]
    point2 = point.seg[seg.idx,]
    
    if(sum((point1 - point0)^2)==0){
      dist.min = 0
      sign.min = 0
      lambda.min = lambda.vec[seg.idx-1]
      break
    }
    if(sum((point2 - point0)^2)==0){
      dist.min = 0
      sign.min = 0
      lambda.min = lambda.vec[seg.idx]
      break
    }
    
    r.tmp = sum((point0 - point1) * (point2 - point1)) / sum((point2 - point1)^2)
    # r.tmp = round(r.tmp, 10)
    if(sum((point2 - point1)^2) == 0){ # point1 == point2인 경우
      dist.tmp = sum((point0 - point1)^2)
      lambda.tmp = lambda.vec[seg.idx - 1]
    }else if(r.tmp <= 0){ # point1과 가장 가까운 경우
      dist.tmp = sum((point0 - point1)^2)
      lambda.tmp = lambda.vec[seg.idx - 1]
    }else if(r.tmp >= 1){ # point2와 가장 가까운 경우
      dist.tmp = sum((point0 - point2)^2)
      lambda.tmp = lambda.vec[seg.idx]
    }else{ # point1과 point2 사이의 점이 가장 가까운 경우
      dist.tmp = sum((point0 - point1)^2) - r.tmp^2 * sum((point1 - point2)^2)
      lambda.tmp = (1 - r.tmp) * lambda.vec[seg.idx - 1] + r.tmp * lambda.vec[seg.idx]
    }
    
    if(dist.tmp<0){
      dist.tmp = 0
    }else{
      dist.tmp = sqrt(dist.tmp)  
    }
    sign.tmp = sign((point0[2] - point1[2]) - ((point2[2] - point1[2])/(point2[1] - point1[1])*(point0[1]-point1[1])))
    
    if(dist.tmp < dist.min){
      dist.min = dist.tmp; sign.min = sign.tmp; lambda.min = lambda.tmp;
    }
    
  }
  return(list(lambda = lambda.min, dev = dist.min * sign.min))
}


GetDeviationFromCurve = function(x, point.seg, arma = FALSE){
  # to remove duplicated points
  idx.duplicated = which(apply(diff(point.seg), 1, function(x) sum(x^2))==0)
  if(length(idx.duplicated) != 0){
    point.seg = point.seg[-idx.duplicated,]  
  }
  
  if(length(x)==2){
    x = as.matrix(x, ncol=2)
    if(arma){
      result.tmp = GetPointDeviationFromCurve_arma(x, point.seg)  
    }else{
      result.tmp = GetPointDeviationFromCurve(x, point.seg)  
    }
    
    lambda.vec = result.tmp$lambda
    dev.vec = result.tmp$dev
  }else{
    if(arma){
      result.tmp = apply(x, 1, function(x) GetPointDeviationFromCurve_arma(as.matrix(x, ncol=2), point.seg))  
    }else{
      result.tmp = apply(x, 1, function(x) GetPointDeviationFromCurve(as.matrix(x, ncol=2), point.seg))  
    }
    
    result.len = length(result.tmp)
    
    dev.vec = numeric(result.len)
    lambda.vec = numeric(result.len)
    
    for(result.idx in 1:result.len){
      dev.vec[result.idx] = result.tmp[[result.idx]]$dev
      lambda.vec[result.idx] = result.tmp[[result.idx]]$lambda
    }
  }
  
  return(list(lambda = lambda.vec, dev = dev.vec))
}



GetRotationMatrix.2d = function(theta){
  return(matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol=2))
}

RotateMat = function(normal.vector){
  n1 = normal.vector[1]
  n2 = normal.vector[2]
  n3 = normal.vector[3]
  
  u1 = c(n2/sqrt(n1^2 + n2^2), -n1/sqrt(n1^2 + n2^2), 0)
  u2 = c(n1*n3/sqrt(n1^2 + n2^2), n2*n3/sqrt(n1^2 + n2^2), -sqrt(n1^2 + n2^2))
  u3 = c(n1, n2, n3)
  
  rotMat = rbind(u1, u2, u3)
  
  return(rotMat)
  
}


FindClosestSplat2D = function(x, splat.collection.2D){
  splat.center = GetSplatCenters(splat.collection.2D)
  dist.from.center = apply(TranslateData(splat.center, center = -x), 1, function(x) sum(x^2))
  return(which.min(dist.from.center))
}

MakeEvaluationList = function(data, splat.collection.2D){
  eval.list = list()
  for(splat.idx in 1:length(splat.collection.2D)){
    eval.list[[splat.idx]] = vector()
  }
  
  for(i in 1:dim(data)[1]){
    splat.idx = FindClosestSplat2D(data[i,1:2], splat.collection.2D)
    tmp = eval.list[[splat.idx]] 
    eval.list[[splat.idx]] = c(tmp, i)
  }
  
  return(eval.list)
}





GetEllipseBoundary = function(splat.idx, splat.collection, split=FALSE){
  splat.tmp = splat.collection[[splat.idx]]
  
  if(split & !is.null(splat.tmp$splat.left)){
    splat.left.tmp = splat.tmp$splat.left
    splat.right.tmp = splat.tmp$splat.right
    
    major.axis.left = splat.left.tmp$major.axis
    minor.axis.left = splat.left.tmp$minor.axis
    normal.vector.left = splat.left.tmp$normal.vector
    center.left = splat.left.tmp$center
    
    major.axis.2d.left = as.vector(RotateMat(normal.vector.left) %*% major.axis.left)[-3]
    theta.left.tmp = atan2(major.axis.2d.left[2], major.axis.2d.left[1])
    rotMat.left = GetRotationMatrix.2d(theta.left.tmp)
    boundary.left = TranslateData(RotateNorm(cbind(MakeEllipse2D(major.axis.left, minor.axis.left)%*%t(rotMat.left), 0), normal.vector.left), center.left)
    
    
    major.axis.right = splat.right.tmp$major.axis
    minor.axis.right = splat.right.tmp$minor.axis
    normal.vector.right = splat.right.tmp$normal.vector
    center.right = splat.right.tmp$center
    
    major.axis.2d.right = as.vector(RotateMat(normal.vector.right) %*% major.axis.right)[-3]
    theta.right.tmp = atan2(major.axis.2d.right[2], major.axis.2d.right[1])
    rotMat.right = GetRotationMatrix.2d(theta.right.tmp)
    boundary.right = TranslateData(RotateNorm(cbind(MakeEllipse2D(major.axis.right, minor.axis.right)%*%t(rotMat.right), 0), normal.vector.right), center.right)
    
    boundary = cbind(boundary.left, boundary.right)
  }else{
    major.axis = splat.tmp$major.axis
    minor.axis = splat.tmp$minor.axis
    normal.vector = splat.tmp$normal.vector
    center = splat.tmp$center
    
    major.axis.2d = as.vector(RotateMat(normal.vector) %*% major.axis)[-3]
    theta.tmp = atan2(major.axis.2d[2], major.axis.2d[1])
    rotMat = GetRotationMatrix.2d(theta.tmp)
    boundary = TranslateData(RotateNorm(cbind(MakeEllipse2D(major.axis, minor.axis)%*%t(rotMat), 0), normal.vector), center)      
  }
  
  return(boundary)
}


GetEllipseBoundary2D = function(splat.idx, splat.collection, data, additional.scaler = 1){
  splat.tmp = splat.collection[[splat.idx]]
  
  splat.center = splat.tmp$center[1:2]
  splat.data.idx = splat.tmp$idx
  splat.data.idx2 = splat.tmp$idx2
  splat.cov.2d = splat.tmp$cov.2d
  
  idx.aux = setdiff(splat.data.idx2, splat.data.idx)
  if(length(idx.aux) != 0){
    idx.aux.min = idx.aux[ which.min( mahalanobis(data[idx.aux,1:2], center = splat.center, cov = splat.cov.2d) ) ]
    splat.data.idx = union(splat.data.idx, idx.aux.min)
  }
  splat.data = data[splat.data.idx,1:2]
  
  svd.tmp = svd(splat.tmp$cov.2d)
  major.axis.tmp = svd.tmp$v[,1] * sqrt(svd.tmp$d[1])
  minor.axis.tmp = svd.tmp$v[,2] * sqrt(svd.tmp$d[2])
  # axes.scaler = GetAxesScaler(splat.data, apply(matrix(splat.data, ncol=2), 2, mean), major.axis.tmp, minor.axis.tmp)
  axes.scaler = additional.scaler * GetAxesScaler(splat.data, splat.center, major.axis.tmp, minor.axis.tmp)
  major.axis.tmp = axes.scaler * major.axis.tmp
  minor.axis.tmp = axes.scaler * minor.axis.tmp
  
  sss = MakeEllipse2D(major.axis.tmp, minor.axis.tmp)
  sss = sss %*% GetRotationMatrix.2d(-atan2(major.axis.tmp[2], major.axis.tmp[1]))
  ellipse = TranslateData(sss, center = splat.center)
  
  return(ellipse)
}





PrintSplatIdxNum = function(splat.collection){
  splat.num = length(splat.collection)
  for(i in 1:splat.num){
    cat(paste0(i, " : "))
    print(sort(splat.collection[[i]]$idx))
  }
}


GetSplitSplatCollection = function(splat.collection){
  splat.num = length(splat.collection)
  splat.collection.new = list()
  
  for(splat.idx in 1:splat.num){
    splat.num.new = length(splat.collection.new)
    if(length(splat.collection[[splat.idx]]$idx) > 3){
      splat.tmp = SplitSplat(splat.collection[[splat.idx]])
      splat.collection.new[[splat.num.new + 1]] = splat.tmp$splat.left
      splat.collection.new[[splat.num.new + 2]] = splat.tmp$splat.right
    }
  }
  
  return(splat.collection.new)
}


GetEvaluatedMatrix  = function(eval.dat, level.end, splat.multilevel, data, curve.fitting = FALSE, orthogonalization = FALSE){
  eval.mat = matrix(0, ncol=1, nrow=dim(eval.dat)[1])
  for(i in 1:level.end){
    splat.level.tmp = splat.multilevel[[i]]
    
    if(i == 1){
      splat.level.2D = MakeSplatCollection2D(splat.collection = splat.level.tmp, data, curve.fitting)
      eval.tmp = EvaluateBasis(x = eval.dat, splat.collection.2D = splat.level.2D, data, curve.fitting)
      # curve.fitting = TRUE 이면 각 스플랫에서의 데이터 변환이 MakeSplatCollection2D()와 EvaluateBasis()에서 중복되어 이루어짐
      # 중복되는 결과가 같음은 확인하였다.
      eval.mat = cbind(eval.mat, eval.tmp)
      eval.mat = eval.mat[,-1]
    }else{
      splat.level.2D = MakeSplatCollection2D(splat.collection = splat.level.tmp, data)
      eval.tmp = EvaluateBasis(eval.dat, splat.level.2D, data)
      if(orthogonalization){
        hat.tmp = eval.mat %*% solve(t(eval.mat) %*% eval.mat) %*% t(eval.mat)
        eval.tmp = (diag(dim(eval.tmp)[1]) - hat.tmp) %*% eval.tmp
      }
      eval.mat = cbind(eval.mat, eval.tmp)
    }
    
  }
  
  return(eval.mat)
}



# MakeSplatAlongCurve = function(splat, data){
#   splat.tmp = splat
#   idx2 = splat.tmp$idx2
#   
#   dat.tmp = data[idx2, 1:2]
#   
#   curvefit = principal_curve(dat.tmp, approx_points = length(idx2)-1) # default 설정은 dev~=0으로 만들어서..?
#   curvefit.poly = ConvertCurve2Polygon(curvefit)
#   if(dim(which(is.nan(curvefit.poly), arr.ind = TRUE))[1] != 0){
#     nan.idx = which(is.nan(curvefit.poly), arr.ind = TRUE)[,1]
#     curvefit.poly = curvefit.poly[-nan.idx,]
#   }
#   
#   curvefit.dev = GetDeviationFromCurve(dat.tmp, curvefit.poly)
#   
#   dat.transformed = cbind(curvefit.dev$lambda, curvefit.dev$dev)
#   
#   center = apply(dat.transformed, 2, mean)
#   cov.tmp = cov(dat.transformed)
#   
#   major.axis.tmp = c(sd(curvefit.dev$lambda), 0)
#   minor.axis.tmp = c(0, sd(curvefit.dev$dev))
#   
#   axes.scaler = GetAxesScaler(dat.transformed, center, major.axis = major.axis.tmp, minor.axis = minor.axis.tmp)
#   
#   splat.tmp$center = center
#   splat.tmp$major.axis = axes.scaler * major.axis.tmp
#   splat.tmp$minor.axis = axes.scaler * minor.axis.tmp
#   splat.tmp$cov.2d = cov.tmp
#   splat.tmp$polygon = curvefit.poly
#   
#   return(splat.tmp)
# }

# dat.transformed, dat.transformed2 구분하기
MakeSplatAlongCurve = function(splat, data){
  splat.tmp = splat
  idx2 = splat.tmp$idx2
  idx = splat.tmp$idx

  dat.tmp2 = data[idx2, 1:2]

  curvefit = principal_curve(dat.tmp2, approx_points = length(idx2)-1) # default 설정은 dev~=0으로 만들어서..?
  curvefit.poly = ConvertCurve2Polygon(curvefit)
  if(dim(which(is.nan(curvefit.poly), arr.ind = TRUE))[1] != 0){
    nan.idx = which(is.nan(curvefit.poly), arr.ind = TRUE)[,1]
    curvefit.poly = curvefit.poly[-nan.idx,]
  }

  curvefit.dev2 = GetDeviationFromCurve(dat.tmp2, curvefit.poly)
  dat.transformed2 = cbind(curvefit.dev2$lambda, curvefit.dev2$dev)
  dat.transformed = matrix(dat.transformed2[1:length(idx),], ncol=2)

  center = apply(matrix(dat.transformed, ncol=2), 2, mean)

  cov.tmp = cov(dat.transformed2)
  major.axis.tmp = c(sd(curvefit.dev2$lambda), 0)
  minor.axis.tmp = c(0, sd(curvefit.dev2$dev))

  #
  dat.transformed.xtd = dat.transformed
  if(length(idx2) > length(idx)){
    # idx.aux here means lil bit different from other idx.aux's in the other ftn's
    idx.aux = (length(idx)+1):length(idx2)
    idx.aux.min = idx.aux[ which.min( mahalanobis(dat.transformed2[idx.aux,1:2], center = center, cov = cov.tmp) ) ]
    dat.transformed.xtd = dat.transformed2[c(1:length(idx), idx.aux.min),]
  }
  #
  axes.scaler = GetAxesScaler(dat.transformed.xtd, center, major.axis = major.axis.tmp, minor.axis = minor.axis.tmp)

  splat.tmp$center = center
  splat.tmp$major.axis = axes.scaler * major.axis.tmp
  splat.tmp$minor.axis = axes.scaler * minor.axis.tmp
  splat.tmp$cov.2d = cov.tmp
  splat.tmp$polygon = curvefit.poly

  return(splat.tmp)
}



# FindIdxSplat = function(splat.collection){
#   obs.num = 0
#   splat.num = length(splat.collection)
#   for(i in 1:splat.num){
#     obs.num = obs.num + length(splat.collection[[i]]$idx)
#   }
#   
#   idx.splat.mat = matrix(nrow=obs.num, ncol=2)
#   idx.splat.mat[,1] = 1:obs.num
#   
#   for(i in 1:splat.num){
#     idx.splat.mat[splat.collection[[i]]$idx,2] = i
#   }
#   
#   return(idx.splat.mat)
# }

FindIdxSplat = function(splat.collection){
  obs.num = 0
  obs.idx = vector()
  
  splat.num = length(splat.collection)
  for(i in 1:splat.num){
    obs.num = obs.num + length(splat.collection[[i]]$idx)
    obs.idx = c(obs.idx, splat.collection[[i]]$idx)
  }
  obs.idx = sort(unique(obs.idx))
  
  idx.splat.mat = matrix(nrow=obs.num, ncol=2)
  idx.splat.mat[,1] = 1:obs.num
  
  for(i in 1:splat.num){
    for(j in splat.collection[[i]]$idx){
      idx.splat.mat[which(j==obs.idx), 2] = i
      idx.splat.mat[which(j==obs.idx), 1] = j
    }
  }
  
  return(idx.splat.mat)
}


RegressWithPC = function(splat.collection, n.init.splat, n.level, prop.var = 0.99, prop.resid = 0.8, intercept=FALSE, curve.fitting=FALSE, data){
  
  splat.collection = MergeSplatCollection(splat.collection = splat.collection, splat.num.target = n.init.splat, data = data)
  
  splat.multilevel = list()
  eg.vector.list = list()
  beta.list = list()
  
  ## splat level 1
  splat.multilevel[[1]] = splat.collection
  eval.mat = GetEvaluatedMatrix(eval.dat = data, level.end = 1, splat.multilevel = splat.multilevel[1], data = data, curve.fitting=curve.fitting)
  svd.tmp = svd(cov(eval.mat))
  pc.num = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > prop.var)[1]
  eg.vector = svd.tmp$v[,1:pc.num]
  eg.vector.list[[1]] = eg.vector
  eval.mat.pc = eval.mat %*% eg.vector
  
  # linear regression model without the intercept term
  wave.mat = data.frame(z.val = data[,3], basis = eval.mat.pc)
  x = model.matrix(z.val~., data=wave.mat)[,-1]
  y = wave.mat$z.val
  
  if(intercept){
    lmfit = lm(y~x)
  }else{
    lmfit = lm(y~x+0)  
  }
  
  beta.list[[1]] = coef(lmfit)
  yhat = lmfit$fitted.values
  resid.tmp = data[,3] - yhat
  
  thsd.tmp = quantile(abs(resid.tmp), prob=prop.resid)
  obs.idx.big.resid = which(abs(resid.tmp) > thsd.tmp)
  splat.idx.big.resid = unique( FindIdxSplat(splat.collection)[obs.idx.big.resid, 2] )
  
  ## splat level 2+
  if(n.level > 1){
    for(level.idx in 2:n.level){
      
      if(length(splat.idx.big.resid)==0){
        break
      }
      
      splat.level.tmp = GetSplitSplatCollection(splat.multilevel[[level.idx-1]][splat.idx.big.resid])
      if(length(splat.level.tmp) == 0){
        break
      }else{
        splat.multilevel[[level.idx]] = splat.level.tmp
      }
      eval.mat = GetEvaluatedMatrix(eval.dat = data, level.end = level.idx, splat.multilevel = splat.multilevel, data = data)
      svd.tmp = svd(cov(eval.mat))
      pc.num = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > prop.var)[1]
      eg.vector = svd.tmp$v[,1:pc.num]
      eg.vector.list[[level.idx]] = eg.vector
      eval.mat.pc = eval.mat %*% eg.vector
      
      # linear regression model without the intercept term
      wave.mat = data.frame(z.val = data[,3], basis = eval.mat.pc)
      x = model.matrix(z.val~., data=wave.mat)[,-1]
      y = wave.mat$z.val
      
      if(intercept){
        lmfit = lm(y~x)
      }else{
        lmfit = lm(y~x+0)  
      }
      
      beta.list[[level.idx]] = coef(lmfit)
      yhat = lmfit$fitted.values
      
      resid.tmp = data[,3] - yhat
      thsd.tmp = quantile(abs(resid.tmp), prob=prop.resid)
      obs.idx.big.resid = which(abs(resid.tmp) > thsd.tmp)
      
      ttt = FindIdxSplat(splat.multilevel[[level.idx]])
      obs.idx.big.resid = intersect(obs.idx.big.resid, ttt[,1])
      
      splat.idx.big.resid = vector()
      for(i in obs.idx.big.resid){
        splat.idx.big.resid = c(splat.idx.big.resid, ttt[which(ttt[,1] == i), 2])
      }
      splat.idx.big.resid = sort(unique(splat.idx.big.resid))
    }  
  }
  
  return(list(splat.multilevel = splat.multilevel, eg.vector = eg.vector.list, beta = beta.list))
}

RegressWithPC = function(x, eval.mat, prop.var = 0.99, prop.resid = 0.8, intercept=FALSE, data){

  eg.vector = NULL
  beta = NULL

  svd.tmp = svd(cov(eval.mat))
  pc.num = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > prop.var)[1]
  eg.vector = svd.tmp$v[,1:pc.num]

  eval.mat.pc = eval.mat %*% eg.vector

  # linear regression model without the intercept term
  wave.mat = data.frame(z.val = x[,3], basis = eval.mat.pc)
  x = model.matrix(z.val~., data=wave.mat)[,-1]
  y = wave.mat$z.val

  if(intercept){
    lmfit = lm(y~x)
  }else{
    lmfit = lm(y~x+0)
  }

  beta = coef(lmfit)

  return(list(eg.vector = eg.vector, beta = beta))
}



GetSplatArea2D = function(splat.idx, splat.collection, data){
  splat.tmp = splat.collection[[splat.idx]]
  
  splat.center = splat.tmp$center[1:2]
  splat.data.idx = splat.tmp$idx
  splat.data.idx2 = splat.tmp$idx2
  splat.cov.2d = splat.tmp$cov.2d
  
  idx.aux = setdiff(splat.data.idx2, splat.data.idx)
  if(length(idx.aux) != 0){
    idx.aux.min = idx.aux[ which.min( mahalanobis(data[idx.aux,1:2], center = splat.center, cov = splat.cov.2d) ) ]
    splat.data.idx = union(splat.data.idx, idx.aux.min)
  }
  splat.data = sim.image[splat.data.idx,1:2]
  
  svd.tmp = svd(splat.tmp$cov.2d)
  major.axis.tmp = svd.tmp$v[,1] * sqrt(svd.tmp$d[1])
  minor.axis.tmp = svd.tmp$v[,2] * sqrt(svd.tmp$d[2])
  axes.scaler = GetAxesScaler(splat.data, splat.center, major.axis.tmp, minor.axis.tmp)
  major.axis.tmp = axes.scaler * major.axis.tmp
  minor.axis.tmp = axes.scaler * minor.axis.tmp
  
  area.tmp = sqrt(sum(major.axis.tmp^2) * sum(minor.axis.tmp^2)) * pi
  return(area.tmp)
}



GetSplitSplatCollectionTotal = function(splat.collection){
  splat.num = length(splat.collection)
  splat.collection.new = list()
  
  for(splat.idx in 1:splat.num){
    splat.num.new = length(splat.collection.new)
    splat.tmp = SplitSplat(splat.collection[[splat.idx]])
    if(is.null(splat.tmp$splat.left)){
      next
    }
    splat.collection.new[[splat.num.new + 1]] = splat.tmp$splat.left
    splat.collection.new[[splat.num.new + 2]] = splat.tmp$splat.right
  }
  
  return(splat.collection.new)
}




GetCurvedSplatBoundary2D = function(splat.idx, splat.collection.2D, data, additional.scaler = 1){
  if(is.null(splat.collection.2D[[splat.idx]]$polygon)){
    cat("Use splat.collection.2D with principal curve polgon.\n")
  }
  
  idx.tmp = splat.collection.2D[[splat.idx]]$idx
  idx2.tmp = splat.collection.2D[[splat.idx]]$idx2
  
  center.trans = splat.collection.2D[[splat.idx]]$center
  major.trans = splat.collection.2D[[splat.idx]]$major.axis
  minor.trans = splat.collection.2D[[splat.idx]]$minor.axis
  major.trans.len = sqrt(sum(major.trans^2))# * additional.scaler
  minor.trans.len = sqrt(sum(minor.trans^2)) * additional.scaler
  
  poly.tmp = splat.collection.2D[[splat.idx]]$polygon
  poly.tmp = poly.tmp[-c(1,dim(poly.tmp)[1]),]
  
  lambda.vec = GetDeviationFromCurve(poly.tmp, splat.collection.2D[[splat.idx]]$polygon)$lambda
  lambda.vec = lambda.vec[which(lambda.vec < center.trans[1] + major.trans.len & lambda.vec > center.trans[1] - major.trans.len)]
  lambda.vec = c(center.trans[1] - major.trans.len, lambda.vec, center.trans[1] + major.trans.len)
  
  expsn.ratio = (lambda.vec[2] - lambda.vec[1]) / (lambda.vec[3] - lambda.vec[2])
  poly.start = poly.tmp[1,] + (poly.tmp[1,] - poly.tmp[2,]) * expsn.ratio
  
  lambda.vec.rev = rev(lambda.vec)
  expsn.ratio = (lambda.vec.rev[2] - lambda.vec.rev[1]) / (lambda.vec.rev[3] - lambda.vec.rev[2])
  poly.dim.tmp = dim(poly.tmp)
  poly.end = poly.tmp[poly.dim.tmp[1],] + (poly.tmp[poly.dim.tmp[1],] - poly.tmp[poly.dim.tmp[1]-1,]) * expsn.ratio
  
  poly.supp = rbind(poly.start, poly.tmp, poly.end)
  
  ellipse.x = lambda.vec - center.trans[1]
  ellipse.x[1] = -major.trans.len # weird part
  ellipse.x[length(ellipse.x)] = major.trans.len # weird part
  lambda.vec.diff = diff(lambda.vec)
  
  poly.prop = lambda.vec.diff / (range(lambda.vec)[2] - range(lambda.vec)[1])
  resol.tmp = 200
  lambda.num.vec = as.vector(round(resol.tmp * poly.prop)) # number of lambdas for evaluation on each segment
  
  curved.supp.plus = vector()
  curved.supp.minus = vector()
  
  for(seg.idx in 1:length(lambda.num.vec)){
    
    lambda.num.tmp = lambda.num.vec[seg.idx] + 2
    # this means that we will divide the segment into (lambda.num.tmp - 1) pieces
    stick.len.tmp = lambda.vec.diff[seg.idx] / (lambda.num.tmp - 1)
    
    poly.x = seq(from=poly.supp[seg.idx,1], to=poly.supp[seg.idx+1,1], length.out = lambda.num.tmp)
    poly.y = seq(from=poly.supp[seg.idx,2], to=poly.supp[seg.idx+1,2], length.out = lambda.num.tmp)
    
    poly.x = poly.x[-length(poly.x)]
    poly.y = poly.y[-length(poly.y)]
    
    # poly.supp.new = rbind(poly.supp.new, cbind(poly.x, poly.y))
    
    slope.tmp = poly.supp[seg.idx,] - poly.supp[seg.idx+1,]
    perp.vec.tmp = c(-slope.tmp[2], slope.tmp[1])
    perp.vec.tmp = perp.vec.tmp / sqrt(sum(perp.vec.tmp^2))
    
    xxx = ellipse.x[seg.idx] + stick.len.tmp * (0:lambda.num.vec[seg.idx])
    yyy.plus = center.trans[2] + minor.trans.len * sqrt(1 - xxx^2/major.trans.len^2)
    yyy.minus = center.trans[2] - minor.trans.len * sqrt(1 - xxx^2/major.trans.len^2)
    
    
    curved.supp.plus = rbind(curved.supp.plus, cbind(poly.x, poly.y) + yyy.plus %o% perp.vec.tmp)
    curved.supp.minus = rbind(curved.supp.minus, cbind(poly.x, poly.y) + yyy.minus %o% perp.vec.tmp)
    
  }
  
  curved.supp = cbind(c(curved.supp.plus[,1], rev(curved.supp.minus[,1]), curved.supp.plus[1,1]), c(curved.supp.plus[,2], rev(curved.supp.minus[,2]), curved.supp.plus[1,2]))
  
  return(curved.supp)
}



SelectSplatForward = function(basis.cand.idx, basis.init.cidx = NULL, eval.mat.obs, basis.max.num, data, prop.var = 0.8, verbose = FALSE, parallelization = FALSE){
  
  AIC.old = AIC.new = 0
  
  if(is.null(basis.init.cidx)){
    
    resid.tmp = numeric(length(basis.cand.idx))
    
    for(basis.cidx in 1:length(basis.cand.idx)){
      resid.tmp.tmp = diag(x=1, nrow=dim(eval.mat.obs)[1]) - 1/sum(eval.mat.obs[,basis.cand.idx[basis.cidx]]^2) * eval.mat.obs[,basis.cand.idx[basis.cidx]] %*% t(eval.mat.obs[,basis.cand.idx[basis.cidx]])
      resid.tmp.tmp = resid.tmp.tmp %*% data[,3]
      resid.tmp[basis.cidx] = sum(resid.tmp.tmp^2)
    }
    
    basis.init.cidx = which.min(resid.tmp)
    
  }
  
  basis.selected.idx.old = basis.selected.idx.new = basis.cand.idx[basis.init.cidx]
  
  design.matrix.tmp = eval.mat.obs[, basis.selected.idx.old]
  
  lmfit.tmp = lm(data[,3]~design.matrix.tmp)
  AIC.old = AIC(lmfit.tmp)
  
  n.var = length(basis.selected.idx.old)
  
  while(n.var < basis.max.num){
    if(verbose){
      cat("n.var:", n.var, "/ AIC:", AIC.old,"\n")
    }
    
    remained.cidx  = which(!is.element(basis.cand.idx, basis.selected.idx.old))
    
    if(parallelization){
      AIC.vec = foreach(basis.cidx = remained.cidx, .combine = 'c') %dopar% {
        basis.selected.idx.tmp = c(basis.selected.idx.old, basis.cand.idx[basis.cidx])
        
        design.matrix.tmp = eval.mat.obs[, basis.selected.idx.tmp]
        lmfit = lm(data[,3]~design.matrix.tmp)
        AIC(lmfit)
      }
    }else{
      
      AIC.vec = vector()
      
      for( basis.cidx in remained.cidx ){
        
        basis.selected.idx.tmp = c(basis.selected.idx.old, basis.cand.idx[basis.cidx])
        
        design.matrix.tmp = eval.mat.obs[, basis.selected.idx.tmp]
        
        lmfit = lm(data[,3]~design.matrix.tmp)
        AIC.vec[length(AIC.vec) + 1] = AIC(lmfit)
      }  
    }
    
    
    AIC.new = min(AIC.vec)
    
    remained.cidx.min = remained.cidx[which.min(AIC.vec)]
    basis.selected.idx.new = c(basis.selected.idx.old, basis.cand.idx[remained.cidx.min])
    
    if(AIC.new > AIC.old){
      if(verbose){
        cat("no more variable will be added.\n\n")
      }
      break
    }else{
      AIC.old = AIC.new
      basis.selected.idx.old = basis.selected.idx.new
      if(is.infinite(AIC.old)){
        break
      }
    }
    
    n.var = length(basis.selected.idx.old)
  }
  
  return(basis.selected.idx.old)
}


EstimateSeq = function(splat.collection.2D, selected.splat.idx, level.partition, eval.mat.obs=NULL, eval.mat.pred=NULL, prop.var = 0.99, curve.fitting=FALSE, data){
  
  grid.tmp = seq(from=0, to=1, length.out = 100)
  
  level.partition.vec = vector()
  for(i in length(level.partition):1){
    level.partition.vec = c(level.partition.vec, level.partition[[i]])
  }
  
  if(is.null(eval.mat.pred)){
    eval.mat.pred = GetEvaluatedMatrix(cbind(rep(grid.tmp,each=100), rep(grid.tmp,100)), level.end = 1, splat.multilevel = list(splat.collection.record[selected.splat.idx]), data = data, curve.fitting = curve.fitting)
  }
  if(is.null(eval.mat.obs)){
    eval.mat.obs = GetEvaluatedMatrix(data[,1:2], level.end = 1, list(splat.collection.record[selected.splat.idx]), data = data, curve.fitting = curve.fitting)
  }
  
  eval.mat.pred = eval.mat.pred[,level.partition.vec]
  eval.mat.obs = eval.mat.obs[,level.partition.vec]
  
  svd.tmp = svd(cov(eval.mat.obs))
  pc.num = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > prop.var)[1]
  eg.vector = svd.tmp$v[,1:pc.num]
  eval.mat.pc = eval.mat.obs %*% eg.vector
  wave.mat = data.frame(z.val = data[,3], basis = eval.mat.pc)
  x = model.matrix(z.val~., data=wave.mat)[,-1]
  y = wave.mat$z.val
  lmfit = lm(y~x+0)
  
  partition.size = list()
  for(i in 1:length(level.partition)){
    partition.size[[i]] = length(level.partition[[i]])
  }
  
  level.partition.new = list()
  idx.tmp = 1:partition.size[[length(partition.size)]]
  for(i in length(partition.size):1){
    level.partition.new[[i]] = idx.tmp
    if(i > 1){
      idx.tmp = rev(idx.tmp)[1] + (1:partition.size[[i-1]])
    }
  }
  
  beta.list = list()
  beta.list[[1]] = eg.vector %*% coef(lmfit)
  
  if(length(level.partition.new)>1){
    A.mat.old = t(eval.mat.obs) %*% (eval.mat.obs)
    for(i in 2:length(level.partition.new)){
      B.mat = A.mat.old[-level.partition.new[[i-1]], level.partition.new[[i-1]]]
      A.mat.new = A.mat.old[-level.partition.new[[i-1]], -level.partition.new[[i-1]]]
      E.mat = solve(A.mat.new) %*% B.mat
      beta.list[[i]] = cbind(diag(rev(level.partition.new[[i-1]])[1] - dim(E.mat)[2]), E.mat) %*% beta.list[[i-1]]
      A.mat.old = A.mat.new
    }
  }
  
  pcr.pred.list = list()
  eval.mat.pred.tmp = matrix(nrow=dim(eval.mat.pred))
  for(i in length(beta.list):1){
    eval.mat.pred.tmp = cbind( eval.mat.pred.tmp, eval.mat.pred[,level.partition.new[[i]] ])
    if(i == length(beta.list)){
      eval.mat.pred.tmp = eval.mat.pred.tmp[,-1]
    }
    pcr.pred.list[[i]] = eval.mat.pred.tmp %*% beta.list[[i]]
  }
  
  return(list(pcr.pred.list = pcr.pred.list, beta.list = beta.list))
}

pcr.jh = function(response, predictor, prop.var = 0.8, scaling = FALSE, BIC = FALSE, centering = FALSE){
  
  if(scaling){
    svd.tmp = svd(cor(predictor))
  }else{
    svd.tmp = svd(cov(predictor))
  }
  
  if(centering){
    pc.num = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > prop.var)[1]
    eg.vector = svd.tmp$v[,1:pc.num]
    
    predictor.centered = CenteringData(predictor, augmentation = FALSE)
    
    design.mat.pc = predictor.centered %*% eg.vector
    design.mat.pc.y = data.frame(response = response, basis = design.mat.pc)
    x = model.matrix(response~., data=design.mat.pc.y)[,-1]
    y = design.mat.pc.y$response
    
    lmfit = lm(y~x)
    beta = c(lmfit$coefficients[1], as.matrix(eg.vector) %*% lmfit$coefficients[-1])
    est = lmfit$fitted.values
    
  }else{
    pc.num = which(cumsum(svd.tmp$d)/sum(svd.tmp$d) > prop.var)[1]
    eg.vector = svd.tmp$v[,1:pc.num]
    design.mat.pc = apply(predictor, 2, mean) %*% eg.vector
    design.mat.pc.y = data.frame(response = response, basis = design.mat.pc)
    x = model.matrix(response~., data=design.mat.pc.y)[,-1]
    y = design.mat.pc.y$response
    
    lmfit = lm(y~x+0)
    beta = eg.vector %*% lmfit$coefficients
    est = lmfit$fitted.values
  }
  
  # AIC = length(lmfit$fitted.values) * log(2 * pi * sum(lmfit$residuals^2) / length(lmfit$fitted.values)) + length(lmfit$fitted.values) + 2 * pc.num # 이러면 더하나 마나 한 것도 basis로 넣게될 듯?
  
  if(!BIC){
    AIC = length(lmfit$fitted.values) * log(2 * pi * sum(lmfit$residuals^2) / length(lmfit$fitted.values)) + length(lmfit$fitted.values) + 2 * (length(lmfit$coefficients) + 1)
  }else{
    AIC = length(lmfit$fitted.values) * log(2 * pi * sum(lmfit$residuals^2) / length(lmfit$fitted.values)) + length(lmfit$fitted.values) + log(length(lmfit$fitted.values)) * (length(lmfit$coefficients) + 1)
  }
  
  return(list(beta = beta, est = est, info.crit = AIC))
  
}


CenteringData <- function(data, augmentation = FALSE) {
  data.centered = data - rep(1, nrow(data)) %*% t(apply(data, 2, mean))
  if(augmentation){
    return(cbind(1, data.centered))  
  }else{
    return(data.centered)
  }
}


GetLOOD = function(splat.2D, data, curve.fitting = FALSE){
  # LOOD stands for Leave-One-Out Deviation
  splat.tmp = splat.2D
  idx2 = splat.tmp$idx2
  idx = splat.tmp$idx
  
  dat.tmp2 = data[idx2, 1:2]
  
  dev.vec = vector()
  
  if(curve.fitting){
    for(i in 1:nrow(dat.tmp2)){
      curvefit = principal_curve(dat.tmp2[-i,], approx_points = length(idx2)-2)
      
      curvefit.poly = ConvertCurve2Polygon(curvefit)
      if(dim(which(is.nan(curvefit.poly), arr.ind = TRUE))[1] != 0){
        nan.idx = which(is.nan(curvefit.poly), arr.ind = TRUE)[,1]
        curvefit.poly = curvefit.poly[-nan.idx,]
      }
      
      curvefit.dev2 = GetDeviationFromCurve(dat.tmp2, curvefit.poly)
      dev.vec[i] = mean(curvefit.dev2$dev^2)
    }
  }else{
    for(i in 1:nrow(dat.tmp2)){
      
      dat.tmp2.loo = dat.tmp2[-i,]
      center.loo = apply(dat.tmp2.loo, 2, mean)
      svd.loo = svd(cov(dat.tmp2.loo))
      
      major.axis.loo = svd.loo$v[,1]
      
      dev.vec.loo = vector()
      for(j in 1:nrow(dat.tmp2)){
        dev.vec.loo[j] = sum((dat.tmp2[j,] - center.loo)^2) - sum((dat.tmp2[j,] - center.loo) * major.axis.loo)^2
      }
      
      dev.loo = mean(dev.vec.loo)
      
      dev.vec[i] = dev.loo
      
    }
  }
  
  return(mean(dev.vec))
}


GenerateSplat = function(data, splat.min.num = 10){
  splat.collection = list()
  for(i in 1:dim(data)[1]){
    splat.collection[[i]] = MakeSplat(i, data)
  }
  splat.collection = MergeSplatCollection2(splat.collection, splat.num.target = splat.min.num, data = data)
  splat.collection.total = splat.collection
  
  splat.collection.new = GetSplitSplatCollectionTotal(splat.collection)
  while(length(splat.collection.new) > 0){
    splat.collection.total = c(splat.collection.total, splat.collection.new)
    splat.collection.new = GetSplitSplatCollectionTotal(splat.collection.new)
  }
  
  return(splat.collection.total)
}