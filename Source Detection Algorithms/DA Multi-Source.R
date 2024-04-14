library(igraph)
library(gtools)
library(e1071)
library(RSpectra)


graph_names=c("Jazz","Karate","Lesmis","Facebook")

GN=1
while(GN<=length(graph_names)){
load(file=paste("./Code, Datasets and Results/Graphs/",graph_names[GN],".RData",sep=''))



if(graph_names[GN]=="Facebook"){
	nos=c(10,12,14,16)
}else{
nos=c(2,3,4,5)}
nos_index=1
while(nos_index<=length(nos)){
path=paste("./Code, Datasets and Results/Multiple Sources Infected Nodes/",nos[nos_index],"/",graph_names[GN],"_Hetero_10_",nos[nos_index],"S.RData",sep='')
object=load(file=path)
infected_nodes_object=get(object)



DA_FB_Ht_2=list() 
DA_sum=0
est_sum=0
outer_sum_original=c()
outer_sum_radiusPenalty=c()
det_rate=c()
time_vector=c()

i=1
while(i<=100){
	time=system.time({ 
	sources=V(graph)[infected_nodes_object[[i]][1:nos[nos_index]]]$name
	sg=induced_subgraph(graph, infected_nodes_object[[i]], impl = c("copy_and_delete"))
	A=as.matrix(get.adjacency(sg)) 
	l=list()
	j=1
	ev=eigs(A,1,which="LM")$values
	while(j<=nrow(A)){
		tryCatch({
		l[j]=abs(ev-eigs(A[-j,-j],1,which="LM")$values)/ev},
		error = function(e){

		l[j]=abs(ev-eigen((A[-j,-j]))$values[1])/ev})
		#print(j)
		j=j+1
		#print(j)
	}
	l=unlist(l)
	estimated_DA=V(sg)[as.numeric(which(max(l)==l))][1]$name

	seeds=c()
	for(r in 1:nos[nos_index]){
		seeds=append(seeds,V(sg)[as.numeric(which(rev(sort(l))[r]==l)[1])]$name)
	}
	}) # time ends
time_vector=append(time_vector,as.numeric(time[3]))
	
if(graph_names[GN]!="Facebook"){
flag=0
	outer_sum=10000000000
	if(length(seeds)<=length(sources)){
		flag=1
	perm=gtools::permutations(n=length(sources),r=length(seeds),v=sources)}
	else{
		perm=gtools::permutations(n=length(seeds),r=length(sources),v=seeds)}

	pit=1
	while(pit<=nrow(perm))
	{
		#print(perm[pit,])
		sit=1
		inner_sum=0
		while(sit<=length(perm[pit,])){
			if(flag==1){
			inner_sum=inner_sum+as.numeric(distances(graph, v = as.character(seeds[sit]), to = as.character(perm[pit,sit]), mode = c("all"), weights=NULL, algorithm = c("unweighted")))
			}
			else{
			inner_sum=inner_sum+as.numeric(distances(graph, v = as.character(sources[sit]), to = as.character(perm[pit,sit]), mode = c("all"), weights=NULL, algorithm = c("unweighted")))
			}

			sit=sit+1
		}
		#print(inner_sum)
		inner_sum=inner_sum/length(sources)
		#print(inner_sum)
		if(inner_sum<outer_sum){
			outer_sum=inner_sum
		}
		pit=pit+1

	}
	#print(k)
#print(outer_sum)
outer_sum_original=append(outer_sum_original,outer_sum)
outer_sum_radiusPenalty=append(outer_sum_radiusPenalty,(outer_sum+radius(sg)*abs(length(sources)-length(seeds))))
}

else{
			outer_sum_original=NA
			outer_sum_radiusPenalty=NA
			outer_sum=NA
		}


dr=length(intersect(as.character(sources),as.character(seeds)))/length(sources)
det_rate=append(det_rate,dr)

est_sum=est_sum+outer_sum

	print(paste("i:",i))

	
	i=i+1
}
print(length(time_vector))

df=data.frame(outer_sum_original,outer_sum_radiusPenalty,det_rate,time_vector)

write.csv(df,paste("./Code, Datasets and Results/Results/DA/",object,".csv",sep=''),sep='\t',row.names=F)
print(paste("***********nos_index: ",nos_index))
nos_index=nos_index+1
}
print(paste("++++++++++++++GN: ",GN))
GN=GN+1
}
	