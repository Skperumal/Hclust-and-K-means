               # K- Means and hierarchical Clustering implementation


hierarchial_Clustering = function(data_input, method=c("single","complete","average","centroid"),cut_tree)
{
  if(!is.matrix(data_input)) data_input = as.matrix(data_input)
  
  type_of_hclust = switch (match.arg(method),
                     single   = min,
                     complete = max,
                     average  = mean)
  whole_N = nrow(data_input)
  diag(data_input)=Inf
  # HEIGHT OF CLUSTER
  height = rep(0,whole_N-1)                   
  # FOR GROUPING
  group_n = -(1:whole_N)                       
  # MERGING MATRIX
  merge_matrix = matrix(0,nrow=whole_N-1, ncol=2)   
  for(j in seq(1,whole_N-1))
  {
    #COMPUTE SMALLEST OF ALL FOR EACH ROW
    height[j] = min(data_input)
    i = which(data_input - height[j] == 0, arr.ind=TRUE)
    #TAKE ONE PER PAIR
    i = i[1,,drop=FALSE]
    p = group_n[i]
    p = p[order(p)]
    merge_matrix[j,] = p
    # Agglomerate GROUPING into Jth group
    grp = c(i, which(group_n %in% group_n[i[1,group_n[i]>0]]))
    group_n[grp] = j
    # CONCAT REPLACEMENT DISTANCES 
    if (method=="centroid")
    { 
      #r=tapply(merge_matrix[1:2], list(rep(cut_tree, ncol(merge_matrix)), col(merge_matrix)), mean)
      r = apply(data_input[i,],2,mean)
    }
    else{
      r = apply(data_input[i,],2,type_of_hclust)
    }
    # DIST MATRIX
    data_input[min(i),] = data_input[,min(i)] = r
    data_input[min(i),min(i)]        = Inf
    data_input[max(i),] = data_input[,max(i)] = Inf
  }
  
  label_for_dendo <- unlist(strsplit(input_lables, split = "\r\n"))
  
  structure(list(merge = merge_matrix, height = height, order = ordering_for_dendo(merge_matrix),
                 labels = trimws(label_for_dendo), method = method, 
          call = match.call(), dist.method = "euclidean"),class = "hclust")
}

ordering_for_dendo = function(merge_matrix)
{
  whole_N = nrow(merge_matrix) + 1
  order_in = rep(0,whole_N)
  order_in[1] = merge_matrix[whole_N-1,1]
  order_in[2] = merge_matrix[whole_N-1,2]
  point = 2
  for(i in seq(whole_N-2,1))
  {
  for(j in seq(1,point))
  {
    if(order_in[j] == i)
    {
      order_in[j] = merge_matrix[i,1]
        if(j==point)
        {
          point = point + 1
          order_in[point] = merge_matrix[i,2]
        } else
        {
          point = point + 1
          for(k in seq(point, j+2)) order_in[k] = order_in[k-1]
          order_in[j+1] = merge_matrix[i,2]
        }
      }
    }
  }
  -order_in
}


#Compute distance of each row and column
dis_computation = function(x,is_centroid)
{
  x = as.matrix(x)
  eval_dis = apply(x*x,1,sum) %*% matrix(1.0,1,nrow(x))
  if (is_centroid==1){
  eval_dis=eval_dis-mean(eval_dis)}
  absol_of_matrix=abs(eval_dis + t(eval_dis) - 2 * x %*% t(x))
  sqrt(absol_of_matrix)
}



# READ DATA FROM FILE
attributes<-read.table("nci.data.txt", header = FALSE, sep = "", dec = ".")
input_lables<-readChar("label.txt", file.info("label.txt")$size)
features<-t(data.frame(attributes))

#COMPUTE THE DISTANCE OF MATRIX 
distance_of_matrix=dis_computation(scale(features),is_centroid=0)

#IMPLEMENTATION OF HIERARCHIAL CLUSTERING
cut_tree=0
clust_data_single = hierarchial_Clustering(distance_of_matrix, method="single",cut_tree)
plot(clust_data_single,hang = -2,labels = NULL,main = "hierarchial Clustering implementation", cex = 0.7,xlab = "Label(Disease)")

clust_data_compl = hierarchial_Clustering(distance_of_matrix, method="complete",cut_tree)
plot(clust_data_compl,hang = -2,labels = NULL,main = "hierarchial Clustering implementation", cex = 0.7,xlab = "Label(Disease)")

clust_data_avg = hierarchial_Clustering(distance_of_matrix, method="average",cut_tree)
plot(clust_data_avg,hang = -2,labels = NULL,main = "hierarchial Clustering implementation", cex = 0.7,xlab = "Label(Disease)")

# For centroid cut tree
cut_tree=cutree(clust_data_avg, 2)
distance_of_matrix=dis_computation(scale(features),is_centroid=1)
clust_data_centro = hierarchial_Clustering(distance_of_matrix, method="centroid",cut_tree)
plot(clust_data_centro,hang = -2,labels = NULL,main = "hierarchial Clustering implementation", cex = 0.7,xlab = "Label(Disease)")

#TASK 4
km.out <- kmeans(matrix(scale(features)), 2, nstart=100)
km.out$cluster
plot(scale(features), col=(km.out$cluster+1), main="K-Means Clustering Results
     with K=2", xlab="", ylab="")
km.out <- kmeans(matrix(scale(features)), 3, nstart=100)
km.out$cluster
plot(scale(features), col=(km.out$cluster+1), main="K-Means Clustering Results
     with K=3", xlab="", ylab="")
km.out <- kmeans(matrix(scale(features)), 14, nstart=100)
km.out$cluster
plot(scale(features), col=(km.out$cluster+1), main="K-Means Clustering Results
     with K=14", xlab="", ylab="")