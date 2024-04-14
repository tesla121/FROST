library(igraph)
library(gtools)
library(e1071)

graph_names=c("Jazz","Karate","Lesmis","Facebook")
graph_names_KC=c("Jazz_KC","Karate_KC","Lesmis_KC","Facebook_KC")


GN=1
while(GN<=length(graph_names)){
load(file=paste("./Code, Datasets and Results/Graphs/",graph_names_KC[GN],".RData",sep=''))

if(graph_names[GN]=="Facebook"){
	nos=c(10,12,14,16)
}else{
nos=c(2,3,4,5)}
nos_index=1
while(nos_index<=length(nos)){
path=paste("./Code, Datasets and Results/Multiple Sources Infected Nodes/",nos[nos_index],"/",graph_names[GN],"_Hetero_10_",nos[nos_index],"S.RData",sep='')
object=load(file=path)
infected_nodes_object=get(object)



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

while(i<=100){
	time=system.time({ 

sources=V(graph)[infected_nodes_object[[i]][1:nos[nos_index]]]$name
	sg=induced_subgraph(graph, infected_nodes_object[[i]], impl = c("copy_and_delete"))
radius=radius(sg)
nnodes<<-length(V(sg))

old_max=0
conv=F



k=1

while(k<=10){



set.seed(100)
old_centers=as.character(c((sample(V(sg))[1:k]$name)))

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
	dist=distances(sg, v = node, to = old_centers, mode = c("all"), weights=NULL, algorithm = c("bellman-ford"))
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
		tr[length(tr)+1]=sum(distances(sg, v = p, to = pt, mode = c("all"), weights=NULL, algorithm = c("bellman-ford")))
	}
	tr=unlist(tr)
nc=which(min(tr)==tr)[1]

new_centers=append(new_centers,pt[nc])
}

new_centers=unlist(new_centers)

p=(old_centers==new_centers)
convergence=prod(p)
old_centers=new_centers

}



 new_max=0
 for(pt in part){
 	ssg=induced_subgraph(graph, pt, impl = c("copy_and_delete"))
 	ecc_max=as.numeric(max(eccentricity(ssg, vids = V(ssg), mode = c("all"))))
 	if(ecc_max>new_max)
 	{
 		new_max=ecc_max
 	}
 }


 	k=k+1
 	if(old_max==new_max){
 		break
 	}
 	old_max=new_max
 	

 }




}) # time ends
time_vector=append(time_vector,as.numeric(time[3]))
indices=unlist(new_centers)
seeds=indices
seeds=as.numeric(seeds)
print(paste("number of estimated sources: ",length(indices)))

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


	print(paste("k-center, i:",i))
	i=i+1	
}
print(length(time_vector))
df=data.frame(outer_sum_original,outer_sum_radiusPenalty,seeds_no,f1,prec,rec,mcc,det_rate,time_vector,TP,TN,FP,FN)

write.csv(df,paste("./Code, Datasets and Results/Results/K-Center/",object,".csv",sep=''),sep='\t',row.names=F)
print(paste("***********nos_index: ",nos_index))
nos_index=nos_index+1
}
print(paste("++++++++++++++GN: ",GN))
GN=GN+1
}

