library(igraph)
library(gtools)
library(e1071)



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

 "%^%" <- function(S, power) 
    with(eigen(S), vectors %*% (values^power * t(vectors))) 

	if(graph_names[GN]!="Facebook"){
	WP=get.adjacency(graph,attr='weight',sparse=F)
	A=as.matrix(get.adjacency(graph))
	A1=A
	D=diag(rowSums(A))}
	WP1=(D%^%-1/2)%*%(WP)%*%(D%^%-1/2)
	
	S=WP1
	print("S computed")
	}
	else{
		load(file=paste("./Code, Datasets and Results/FROST_S/S_",graph_names[GN],"_FROST.RData",sep=''))
		print("S loaded")
		A=as.matrix(get.adjacency(graph))
		A1=A
	}

	
	rownames(S)=rownames(A)
	colnames(S)=colnames(A)
	S1=S
	
	out_min_sum=0
s=list()
alpha=0.5 # constant throughout the network
seeds_no=c()
frost_Karate_Ht_2=list()
frost_sum=0
k=1
f1_sum=0
est_sum=0
f1=c()
rec=c()
prec=c()
mcc=c()
TP=c()
TN=c()
FP=c()
FN=c()
det_rate=c()
outer_sum_original=c()
outer_sum_radiusPenalty=c()
time_vector=c()
	
	candidates=list()
	while(k<=100){
		A=A1
		S=S1
		frontiers=c()
		for(i in as.character(infected_nodes_object[[k]])) {
			nbr=as.numeric(which(A[i,]==1))
			frontiers=append(frontiers,nbr)
		}
		
		candidates[[k]]=unique(c(infected_nodes_object[[k]],frontiers))
		k=k+1}

		k=1
	while(k<=100){
	time=	system.time({ 
	
	
	 A=A1
	 S=S1
	

	Y=c(rep(-1,length(candidates[[k]])))
	inds=which(infected_nodes_object[[k]] %in% candidates[[k]])
	Y[inds]=c(rep(1,length(infected_nodes_object[[k]])))
	df=NULL
	A=A[candidates[[k]],candidates[[k]]]
	S=S[candidates[[k]],candidates[[k]]]

	
	df=as.data.frame(rbind(Y))
	iter=2
	while(iter<=5){ 
		G=c(rep(0,length(candidates[[k]])))
		for(i in sample(1:length(candidates[[k]]))) {
			summation=0

			nbr=as.numeric(which(A[i,]==1))
			

			for(j in nbr){
				summation=summation+unlist(S[i,j])*unlist(df[iter-1,j]) 
			}
			G[[i]]=alpha*summation+(1-alpha)*unlist(Y[i]) 
			
		}
		
		df=as.data.frame(rbind(df,unlist(G)))
		
		iter=iter+1
	}
	

	indices=list()
	n_s=0
	it=1
	Z=as.numeric(df[5,])
	while(it<=length(Z)) {
		nbr=as.numeric(which(A[it,]==1))
		if(Y[it]==1){
		if(length(which(Z[nbr]>=Z[it]))==0)
		{
			n_s=n_s+1
			indices=append(indices,rownames(A)[it])
		}
	}
		it=it+1

	}
	}) # time ends
	time_vector=append(time_vector,as.numeric(time[3]))
	print(paste("number of estimated sources: ",n_s))
	indices=unlist(indices)
	s=append(s,n_s)

	#sources=V(graph)[infected_nodes_object[[k]][1:4]]$name
	sources=V(graph)[infected_nodes_object[[k]][1:nos[nos_index]]]$name
	print(paste("number of actual sources: ",length(sources)))
	seeds=indices
seeds=as.numeric(seeds)
sg=induced_subgraph(graph, infected_nodes_object[[k]], impl = c("copy_and_delete"))
if(length(seeds)==0){
	est_sum=est_sum+radius(sg)
	f1=append(f1,0)
	prec=append(prec,0)
	rec=append(rec,0)
	seeds_no=append(seeds_no,length(seeds))
	det_rate=append(det_rate,0)
	mcc=append(mcc,-1)
	TP=append(TP,0)
	TN=append(TN,0)
	FP=append(FP,0)
	FN=append(FN,0)
	
	outer_sum_original=append(outer_sum_original,radius(sg)*nos[nos_index])
	outer_sum_radiusPenalty=append(outer_sum_radiusPenalty,radius(sg)*nos[nos_index])

	k=k+1
	next
}

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

	est_sum=est_sum+outer_sum

	

	
	tp=length(intersect(sources,seeds))

	pr=tp/length(indices)
	prec=append(prec,pr) 
	rc=tp/length(sources)
	rec=append(rec,rc)
	# code for mcc begins
	fp=length(indices)-tp
	fn=length(sources)-tp
	actual_negatives=setdiff(as.character(infected_nodes_object[[k]]),as.character(sources))
	tn=length(actual_negatives)-length(intersect(as.character(actual_negatives),as.character(indices)))
	MCC=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
	mcc=append(mcc,MCC)
	TP=append(TP,tp)
	TN=append(TN,tn)
	FP=append(FP,fp)
	FN=append(FN,fn)
	# code for mcc ends
	if(pr==0&&rc==0){
	f1=append(f1,0)
	}
	else{
	f1=append(f1,((2*pr*rc)/(pr+rc)))
	}

seeds_no=append(seeds_no,length(seeds))
dr=length(intersect(as.character(sources),as.character(seeds)))/length(sources)
det_rate=append(det_rate,dr)



print(paste("frost, k:",k))


	k=k+1
	}
print(paste('time_vector: ',length(time_vector)))

df=data.frame(outer_sum_original,outer_sum_radiusPenalty,seeds_no,f1,prec,rec,mcc,det_rate,time_vector,TP,TN,FP,FN)

write.csv(df,paste("./Code, Datasets and Results/Results/FROST/",object,".csv",sep=''),sep='\t',row.names=F)
print(paste("***********nos_index: ",nos_index))
nos_index=nos_index+1
}
print(paste("++++++++++++++GN: ",GN))
GN=GN+1
}
