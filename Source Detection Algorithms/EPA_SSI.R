l2=list()
temp <<- list()
level=1
score<-function(Du1,Di,l,level,temp){
		if(level>(radius(sg))){
			return(0)
		}
	for(item in l){
		l2=append(l2,attributes(which(Ai[item,]==1))$names)

	}
	l2=unlist(l2)
	l3=setdiff(l2,l)
	l3=setdiff(l3,intersect(l3,temp))
	#print(l3)
	l3=as.character(l3)
	scr=0
	if(length(l3)==0){
		return (0)
	}
	for(item in l){
		scr=scr+(Di[item,item]/Du1[item,item])/(1/(1+log(Du1[item,item]))) 
	}
	alpha=0
	temp=unique(unlist(append(temp,l)))
	if(length(l3)==0){
		return (0)
	}
	else{
		return (scr+score(Du1,Di,l3,level+1,temp))
	}
}

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



Au=as.matrix(get.adjacency(graph))
Du=diag(rowSums(Au))
rownames(Du)=as.character(V(graph)$name)
colnames(Du)=as.character(V(graph)$name)


EPA_SSI_Karate_Ht_10_5S_DR=list()
EPA_SSI_Karate_Ht_10_5S_AED=list()
f1_sum=0
out_min_sum=0
sys=0
est_sum=0
outer_sum_original=c()
outer_sum_radiusPenalty=c()
det_rate=c()
time_vector=c()
EPA_SSI_Karate_Ht_10_5S_DR=list()
i=1
while(i<=100){
	time=system.time({ 
	sources=V(graph)[infected_nodes_object[[i]][1:nos[nos_index]]]$name
	sg=induced_subgraph(graph, infected_nodes_object[[i]], impl = c("copy_and_delete"))
	radius=radius(sg)
	nnodes<<-length(V(sg))
	Ai<<-as.matrix(get.adjacency(sg))
	Di=diag(rowSums(Ai))
	colnames(Di)=as.character(V(sg)$name)
	rownames(Di)=as.character(V(sg)$name)
	Du1=Du[rownames(Di),colnames(Di)]
	kk=1
	seeds=list()
	est_sum_in=0
	old=0
	while(kk<=nos[nos_index]){
		scr_list=list()
		k=1
		for(node in as.character(V(sg)$name) ){
			scr=score(Du1,Di,as.character(node),level,temp)
			penalty=as.numeric(eccentricity(sg, vids = node, mode = c("all")))
			scr_list[k]=scr/penalty
			k=k+1
		} 
		scr_list=unlist(scr_list)
		index=which(max(scr_list)==scr_list)[1]
		seeds=append(seeds,as.character(V(sg)[index]$name))
		seeds=unlist(seeds)
		Ai[which(rownames(Ai)==(V(sg)[index]$name)),]=c(rep(0,length(rownames(Ai))))
		Ai[,which(rownames(Ai)==(V(sg)[index]$name))]=c(rep(0,length(rownames(Ai))))
		Di[which(rownames(Di)==(V(sg)[index]$name)),]=c(rep(0,length(rownames(Di))))
		Di[,which(rownames(Di)==(V(sg)[index]$name))]=c(rep(0,length(rownames(Di))))
		Du1=Du[rownames(Di),colnames(Di)]
		if(length(Di)==0){
		#break
		}
		kk=kk+1
		new=max(scr_list)
		if(new==old){
		
		}
		else{
			old=new
		}
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
		#print(inner_sum)
		if(inner_sum<outer_sum){
			outer_sum=inner_sum
		}
		pit=pit+1

	}
	
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

write.csv(df,paste("./Code, Datasets and Results/Results/EPA_SSI/",object,".csv",sep=''),sep='\t',row.names=F)
print(paste("***********nos_index: ",nos_index))
nos_index=nos_index+1
}
print(paste("++++++++++++++GN: ",GN))
GN=GN+1
}
