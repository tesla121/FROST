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


DC_Karate_Ht_10_5S_DR=list()
DC_Karate_Ht_10_5S_AED=list()

est_sum=0
out_min_sum=0
f1_sum=0
outer_sum_original=c()
outer_sum_radiusPenalty=c()
det_rate=c()
time_vector=c()
out_min_sum=0
seeds_no=c()
s=list()
i=1
f1_sum=0
est_sum=0
f1=c()
rec=c()
prec=c()
det_rate=c()
mcc=c()
TP=c()
TN=c()
FP=c()
FN=c()
outer_sum_original=c()
outer_sum_radiusPenalty=c()
time_vector=c()

i=1
while(i<=100){
	time=system.time({ 
	if (nos[nos_index]==2 && graph_names[GN]=='Jazz'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 6))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==3 && graph_names[GN]=='Jazz'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 7))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==4 && graph_names[GN]=='Jazz'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 7))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==5 && graph_names[GN]=='Jazz'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 2, 8))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==2 && graph_names[GN]=='Karate'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 4))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==3 && graph_names[GN]=='Karate'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 5))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==4 && graph_names[GN]=='Karate'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 5))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==5 && graph_names[GN]=='Karate'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 5))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==2 && graph_names[GN]=='Lesmis'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 5))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==3 && graph_names[GN]=='Lesmis'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 5))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==4 && graph_names[GN]=='Lesmis'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 5))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==5 && graph_names[GN]=='Lesmis'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109)
		seed=round(runif(1, 100, 109))
		set.seed(seed)
		rand=round(runif(1, 1, 7))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==10 && graph_names[GN]=='Facebook'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120)
		seed=round(runif(1, 100, 120))
		set.seed(seed)
		rand=round(runif(1, 4, 15))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==12 && graph_names[GN]=='Facebook'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120)
		seed=round(runif(1, 100, 120))
		set.seed(seed)
		rand=round(runif(1, 6, 17))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==14 && graph_names[GN]=='Facebook'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120)
		seed=round(runif(1, 100, 120))
		set.seed(seed)
		rand=round(runif(1, 8, 18))
		print(paste("rand is: ",rand))
	}
	else if (nos[nos_index]==16 && graph_names[GN]=='Facebook'){
		seed_set=c(100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120)
		seed=round(runif(1, 100, 120))
		set.seed(seed)
		rand=round(runif(1, 10, 20))
		print(paste("rand is: ",rand))
	}
	sources=V(graph)[infected_nodes_object[[i]][1:nos[nos_index]]]$name
	sg=induced_subgraph(graph, infected_nodes_object[[i]], impl = c("copy_and_delete"))
	radius=radius(sg)
	nnodes<<-length(V(sg))
	old_max=0
	conv=F
	set.seed(100)
	old_centers=as.character(c((sample(V(sg))[1:rand]$name)))
	iterations=1
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
		tr=list()
		for(p in pt)
		{
			tr[length(tr)+1]=sum(distances(sg, v = p, to = pt, mode = c("all"), weights=NULL, algorithm = c("unweighted")))
		}
		tr=unlist(tr)
		nc=which(min(tr)==tr)
		nc=which(min(tr)==tr)[1]
		new_centers[length(new_centers)+1]=pt[nc]
	}
	new_centers=unlist(new_centers)
	p=(old_centers==new_centers)
	convergence=prod(p)
	if(iterations>=10){
		convergence=T
	}
	old_centers=new_centers
	iterations=iterations+1
}
}) # time ends
time_vector=append(time_vector,as.numeric(time[3]))

indices=unlist(new_centers)
seeds=indices
seeds=as.numeric(seeds)
print(paste("number of estimated sources: ",length(indices)))
print(paste("number of actual sources: ",length(sources)))

if(length(seeds)==0){
	est_sum=est_sum+radius(sg)
	f1=append(f1,0)
	prec=append(prec,0)
	rec=append(rec,0)
	mcc=append(mcc,-1)
	TP=append(TP,0)
	TN=append(TN,0)
	FP=append(FP,0)
	FN=append(FN,0)
	seeds_no=append(seeds_no,length(seeds))
	det_rate=append(det_rate,0)
	outer_sum_original=append(outer_sum_original,radius(sg)*nos[nos_index])
	outer_sum_radiusPenalty=append(outer_sum_radiusPenalty,radius(sg)*nos[nos_index])
	i=i+1
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
	actual_negatives=setdiff(as.character(infected_nodes_object[[i]]),as.character(sources))
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


	print(paste("DC, i:",i))
	i=i+1	
}
print(length(time_vector))

df=data.frame(outer_sum_original,outer_sum_radiusPenalty,seeds_no,f1,prec,rec,mcc,det_rate,time_vector,TP,TN,FP,FN)

write.csv(df,paste("./Code, Datasets and Results/Results/DC Randomised/",object,".csv",sep=''),sep='\t',row.names=F)
print(paste("***********nos_index: ",nos_index))
nos_index=nos_index+1
}
print(paste("++++++++++++++GN: ",GN))
GN=GN+1
}

