library(bnlearn)
library(Rgraphviz)
library(visNetwork)
library(igraph)
data("asia")
set.seed(123456789)


# Assignment 1 ------------------------------------------------------------

plot.network <- function(structure, ht = "400px",my_title){
  # https://www.r-bloggers.com/bayesian-network-example-with-the-bnlearn-package/
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  
  return(visNetwork(nodes, edges, height = ht, width = "100%"
                    ,background="#eeefff",main=my_title))
}

#---------------------------------------------------------------
dag = hc(asia,restart = 3,score="aic")
dag
plot.network(dag,my_title="")
arcs(dag)
vstructs(dag)
cpdag(dag)
#---------------------------------------------------------------
dag2 = hc(asia,restart = 3,score="bic")
dag2
plot.network(dag2,my_title = "")
arcs(dag2)
vstructs(dag2)
cpdag(dag2)
#---------------------------------------------------------------
initial_dag<-model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")
rnd <- random.graph(nodes = colnames(asia))
dag3=hc(asia,restart=2,score="bde",max.iter = 20,star=rnd)
plot.network(dag3,my_title = "")
arcs(dag3)
vstructs(dag3)
cpdag(dag3)
score(dag3,asia)
#---------------------------------------------------------------


# graphviz.plot(dag, shape = "ellipse", highlight = list(arcs = wl))


# Assignment 2 ------------------------------------------------------------

library(gRain)

## 80% of the sample size
smp_size <- floor(0.80 * nrow(asia))

## set the seed 
set.seed(12345)


asia_character<-data.frame(lapply(asia, as.character), 
                           stringsAsFactors=FALSE) # create a df only for to use in the setEvidence function

indx <- sample(nrow(asia), size = smp_size)


asia_test<-asia_character[-indx,] # use this only for the states arg in setEvidence

train_data <- asia[indx, ]
test_data <- asia[-indx, ]


dag.fit<-hc(train_data,restart = 5,score="bde",max.iter = 10) # learn structure

modelstring(dag.fit)

plot.network(dag.fit,my_title = "Network with brute force")

# dag.igraph<-igraph.from.graphNEL(as.graphNEL(bn.model)) # convert to igraph object
# plot(dag.igraph, frame = TRUE,
#      main = " BN Network Asia",
#      vertex.color = sample(colors(),1), vertex.label.color = sample(colors(),1),
#      vertex.label.cex = 3, vertex.size = 50,
#      edge.arrow.size = 0.8, edge.color = "black") # plot the igraph object 


bn.model<-bn.fit(dag.fit,data=train_data,method = "mle") # fit the model-learn the parameters 
bn.model

bn.model.grain<-compile(as.grain(bn.model)) # compile as grain object

col.param=length(colnames(asia))/2

par(mfrow=c(2,col.param))
for(i in 1:length(bn.model)){
  obj<-bn.model[[i]]$prob
  if(names(bn.model)[i]%in%c("A","S","T")){
  barplot(obj,col=sample(colors()),density=60,
       main=paste("Conditional for",names(bn.model)[i]))
  }else{
    plot(obj,col=sample(colors()),
         main=paste("Conditional for",names(bn.model)[i]))
    }
  ""
}

pred_mat<-double(dim(test_data)[1])


for(i in 1:dim(test_data)[1]){
  evidence.obj<-setEvidence(bn.model.grain,nodes=colnames(test_data[-2])
                     ,states=asia_test[i,-2])
  q<-querygrain(evidence.obj,nodes=c("S"),type="marginal")
  pred_mat[i]<-ifelse(q$S[1]>q$S[2],"no","yes")
}


head(pred_mat)

table1<-prop.table(table(test_data[,2],pred_mat,
                         dnn=c("true", "prediction")))
table1

require(caret)
fourfoldplot(table1, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_as2<-sum(pred_mat==test_data[,2])/dim(test_data)[1]
cat("The accuracy of the model is:",accuracy_as2*100,"%")


# Now for true network ----------------------------------------------------

dag.true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
plot.network(dag.true,my_title = "True Network")

# Assignment 3 ------------------------------------------------------------
  
#bnlearn::parents((bn.model),"S") # parents of "S"

#bnlearn::children(bn.model,"S")
#sapply(bnlearn::children(dag.fit,"S"),function(y){bnlearn::parents(dag.fit,node=y)}) # parents of children


mv.blanket<-bnlearn::mb(dag.fit,"S") # markov blancket for S

indexes<-c()
for(obj in mv.blanket){
  ind<-which(obj==colnames(test_data))
  print(ind)
  indexes<-c(indexes,ind)
}  

pred_mat1<-double(dim(test_data)[1])

for(i in 1:dim(test_data)[1]){
  evidence.obj<-setEvidence(bn.model.grain,nodes=colnames(test_data[indexes])
                            ,states=asia_test[i,indexes])
  q1<-querygrain(evidence.obj,nodes=c("S"),type="marginal")
  pred_mat1[i]<-ifelse(q1$S[1]>q1$S[2],"no","yes")
}


head(pred_mat1)

table2<-prop.table(table(test_data[,2],pred_mat1,
                         dnn=c("true", "prediction")))
table2

require(caret)
fourfoldplot(table2, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")

accuracy_as3<-sum(pred_mat1==test_data[,2])/dim(test_data)[1]
cat("The accuracy of the model is:",accuracy_as3*100,"%")

# Assignment 4 ------------------------------------------------------------

string<-""
for(i in 1:length(colnames(asia))){
  if(colnames(asia)[i]!="S"){
    str2<-paste("[",colnames(asia)[i],"|S]",sep="")
  }else{
    str2<-"[S]"
  }
  string<-paste(string,str2,sep="")
}

string

naive.dag<- model2network(string)
plot.network(naive.dag,my_title="Naive Bayes Network")

naive.model<-bn.fit(naive.dag,data=train_data,method = "mle") # fit the model-learn the parameters 
naive.model

naive.model.grain<-compile(as.grain(naive.model))

pred_mat3<-double(dim(test_data)[1])

for(i in 1:dim(test_data)[1]){
  naive.obj<-setEvidence(naive.model.grain,nodes=colnames(test_data[-2])
                        ,states=asia_test[i,-2])
 
  q3<-querygrain(naive.obj,nodes="S",type="marginal")
  pred_mat3[i]<-ifelse(q3$S[1]>q3$S[2],"no","yes")
  
}

head(pred_mat3)

table3<-prop.table(table(test_data[,2],pred_mat3,
                         dnn=c("true", "prediction")))
table3

fourfoldplot(table3, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_as4<-sum(pred_mat3==test_data[,2])/dim(test_data)[1]
cat("The accuracy of the model is:",accuracy_as4*100,"%")

library(knitr)
x<-rbind(accuracy_as2,accuracy_as3,accuracy_as4)
x<-as.data.frame(x) ; colnames(x)<-c("Accuracy")

kable(x, digits = 2, caption = "A table produced by printr.")
  

# Code to plot with graphviz ----------------------------------------------

arc.set <- arcs(naive.dag)

highlight.opts<-list(nodes=colnames(asia),col =sample(colors(),1), fill = sample(colors(),1),arcs=arc.set,lty=5,lwd=2)
graphviz.plot(naive.model,highlight = highlight.opts,layout = "neato")
  


# Prediction function -----------------------------------------------------


prediction_func<-function(fit.dag,train_data,test_data,method,index,node){
  fit.model<-bn.fit(fit.dag,data=train_data,method = method) # fit the model-learn the parameters 
  
  fit.model.grain<-compile(as.grain(fit.model))
  
  pred_vector<-double(dim(test_data)[1])
  
  for(i in 1:dim(test_data)[1]){
    evidence.obj<-setEvidence(fit.model.grain,nodes=colnames(test_data[-index])
                           ,states=asia_test[i,-index])
    query.obj<-querygrain(evidence.obj,nodes=node,type="marginal")
    pred_vector[i]<-ifelse(query.obj[[1]][1]>query.obj[[1]][2],"no","yes")
  }
  return(pred_vector)
}

pred_naive<-prediction_func(xs[[best_index]],train_data,test_data,"mle", 2,"S")


tab<-prop.table(table(test_data[,2],pred_naive,
                         dnn=c("true", "prediction")))
tab

fourfoldplot(tab, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_naive<-sum(pred_naive==test_data[,2])/dim(test_data)[1]
cat("The accuracy of the model is:",accuracy_naive*100,"%")




# Find best parameters -----------------------------------------------------------

set.seed(123456789)
restartGrid<-seq(1,10,1) ; scoreGrid<-c("loglik","aic","bic","bdla","bdj","bde","bds","mbde")
max.iterGrid<-seq(1,10,1) ; rnd<-random.graph(nodes = colnames(asia))

gridMatrix<-expand.grid(restartGrid,scoreGrid,max.iterGrid)
colnames(gridMatrix)<-c("restart","score","max.iter")
gridMatrix[,2]<-as.character(gridMatrix[,2])

xs<-apply(gridMatrix,1,function(x){hc(asia,restart=as.numeric(x[1]),score=x[2],
                             max.iter=as.numeric(x[3]),star=rnd)})
xx<-lapply(xs,function(y){score(y,asia)})
best_index=which.max(unlist(xx))

best_combination<-gridMatrix[best_index,]
cat("The best combination is :\n")
best_combination

plot.network(xs[[best_index]],my_title = "Network with Brute Force")








