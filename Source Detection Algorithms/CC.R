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



CC_Karate_Ht_10_5S_DR=list()
CC_Karate_Ht_10_5S_AED=list()

outer_sum_original=c()
outer_sum_radiusPenalty=c()
det_rate=c()
time_vector=c()
est_sum=0
out_min_sum=0
f1_sum=0
i=1
while(i<=100){
	time=system.time({ 
	sources=V(graph)[infected_nodes_object[[i]][1:nos[nos_index]]]$name
	sg=induced_subgraph(graph, infected_nodes_object[[i]], impl = c("copy_and_delete"))
	radius=radius(sg)
	nnodes<<-length(V(sg))
	old_max=0
	conv=F
	set.seed(100)
	old_centers=as.character(c((sample(V(sg))[1:nos[nos_index]]$name)))
	convergence=F
	while(!convergence){
		part=list()
		h=1
		while(h<=length(old_centers)){
			part[[h]]=old_centers[h]
			h=h+1
		}
		for(node in as.character(V(sg)$name))
		{
			dist=distances(sg, v = node, to = old_centers, mode = c("all"), weights=NULL, algorithm = c("unweighted"))
			which=which(min(dist)==dist)[1]
			
				part[[which]]=c(part[[which]],node)
			
		}
	new_centers=list()
	for(xd in 1:length(part)){
		if(length(part[[xd]])>1){part[[xd]]=part[[xd]][-1]}
	}
		for(pt in part){
			ssg=induced_subgraph(sg, pt, impl = c("copy_and_delete"))
			if(V(ssg)!=1){
				tr=closeness(ssg, vids = V(ssg), mode = c( "all"), weights = NULL, normalized = FALSE)
			}
			else{tr=0}
			nc=which(max(tr)==tr)[1]
			new_centers[length(new_centers)+1]=pt[nc]}
		new_centers=unlist(new_centers)
		p=(old_centers==new_centers)
		convergence=prod(p)
		old_centers=new_centers
	}

}) # time ends
time_vector=append(time_vector,as.numeric(time[3]))

indices=unlist(new_centers)
seeds=indices
print(paste("number of estimated sources: ",length(indices)))





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

		inner_sum=inner_sum/length(sources)

		if(inner_sum<outer_sum){
			outer_sum=inner_sum
		}
		pit=pit+1

	}

outer_sum_original=append(outer_sum_original,outer_sum)
outer_sum_radiusPenalty=append(outer_sum_radiusPenalty,(outer_sum+radius(sg)*abs(length(sources)-length(seeds))))
dr=length(intersect(as.character(sources),as.character(seeds)))/length(sources)
det_rate=append(det_rate,dr)



est_sum=est_sum+outer_sum


	print(paste("i:",i))
	i=i+1
	
}
print(length(time_vector))
df=data.frame(outer_sum_original,outer_sum_radiusPenalty,det_rate,time_vector)

write.csv(df,paste("./Code, Datasets and Results/Results/CC/",object,".csv",sep=''),sep='\t',row.names=F)
print(paste("***********nos_index: ",nos_index))
nos_index=nos_index+1
}
print(paste("++++++++++++++GN: ",GN))
GN=GN+1
}
	