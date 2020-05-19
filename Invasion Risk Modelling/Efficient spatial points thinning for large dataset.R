## for occurrence data, remove points with distance <dist to decrease potential data biase


dis_filter_largeData = function(pts, dist, ncore){
  
  n = length(pts)
  ni = floor(n/10000)
  
  crds= pts@coords
  km = kmeans(crds, centers= ni + 1)
  pts$group=km$cluster
  pts=pts[order(pts$group),]
  
   
  if (ni < ncore){ncore = ni}
  
  cl=makeCluster(ncore)
  registerDoParallel(cl)
  pts_thin = foreach(i=1:ni, .combine=rbind,.packages=c("geosphere"),.export=c("pts","dist")) %dopar%{
    
    id1 = (i-1) * 10000 + 1
    id2 = i * 10000
    
    pts_p = pts[id1:id2,]
    disMa = distm(pts_p)
    disMa2 = disMa/1000
    disMa2[upper.tri(disMa2, diag=T)] = 99
    disMa2[disMa2<dist] = NA
    pts_p2 = pts_p[complete.cases(disMa2),]
    
    return(pts_p2) 
  }
  stopCluster(cl)
  
  pts_thin
  
  id1 = ni * 10000 +1
  id2 = length(pts)
  pts_p = pts[id1:id2,]
  
  disMa = distm(pts_p)
  disMa2 = disMa/1000
  disMa2[upper.tri(disMa2, diag=T)] = 99
  disMa2[disMa2<dist] = NA
  pts_p2 = pts_p[complete.cases(disMa2),]
  
  pts_thin2 = rbind(pts_thin, pts_p2)
  
  return(pts_thin2)
}
