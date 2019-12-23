# load libraries and data 
library(bnlearn)
library(Rgraphviz)
library(visNetwork)
library(igraph)
library(gRain)
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


restartGrid<-seq(1,10,1) ; scoreGrid<-c("loglik","aic","bic","bdla","bdj","bde","bds","mbde")
max.iterGrid<-seq(1,10,1) ; rnd<-random.graph(nodes = colnames(asia))

gridMatrix<-expand.grid(restartGrid,scoreGrid,max.iterGrid)
colnames(gridMatrix)<-c("restart","score","max.iter")
gridMatrix[,2]<-as.character(gridMatrix[,2])

xs<-apply(gridMatrix,1,function(x){hc(asia,restart=as.numeric(x[1]),score=x[2],
                             max.iter=as.numeric(x[3]),star=rnd)})
xx<-lapply(xs,function(y){bnlearn::score(y,asia)})
best_index=which.max(unlist(xx))

best_combination<-gridMatrix[best_index,]


knitr::kable(best_combination)


cat("The string model is :\n")
modelstring(xs[[best_index]])



# Assignment 2 ------------------------------------------------------------

## 80% of the sample size
smp_size <- floor(0.80 * nrow(asia))

asia_character<-data.frame(lapply(asia, as.character), 
                           stringsAsFactors=FALSE) # create a df only for to use in the setEvidence function

indx <- sample(nrow(asia), size = smp_size)


asia_test<-asia_character[-indx,] # use this only for the states arg in setEvidence

train_data <- asia[indx, ]
test_data <- asia[-indx, ]


dag.fit<-hc(train_data,restart = 5,score="bde",max.iter = 10) # learn structure

bn.model<-bn.fit(dag.fit,data=train_data,method = "mle") # fit the model-learn the parameters 

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

pred_bf<-prediction_func(dag.fit,train_data,test_data,"mle", 2,"S")

tab_bf<-prop.table(table(test_data[,2],pred_bf,
                         dnn=c("true", "prediction")))
cat("The confusion matrix for the brute force model is :\n")
tab_bf


fourfoldplot(tab_bf, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_bf<-sum(pred_bf==test_data[,2])/dim(test_data)[1]
cat("The accuracy of brute force model is:",accuracy_bf*100,"%")

dag.true = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")



pred_true<-prediction_func(dag.true,train_data,test_data,"mle", 2,"S")

tab_true<-prop.table(table(test_data[,2],pred_true,
                         dnn=c("true", "prediction")))
cat("The confusion matrix for the brute force model is :\n")
tab_true


fourfoldplot(tab_true, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_true<-sum(pred_true==test_data[,2])/dim(test_data)[1]
cat("The accuracy of true model is:",accuracy_true*100,"%")


# Assignment 3 ------------------------------------------------------------

mv.blanket<-bnlearn::mb(dag.fit,"S") # markov blancket for S

indexes<-c()
for(obj in mv.blanket){
  ind<-which(obj==colnames(test_data))
  indexes<-c(indexes,ind)
}  
r<-seq(1:8) ; indexes<-r[-indexes]

# ----------------------------------------------------------
pred_mb<-prediction_func(dag.fit,train_data,test_data,"mle", indexes,"S")

tab_mb<-prop.table(table(test_data[,2],pred_mb,
                         dnn=c("true", "prediction")))
cat("The confusion matrix for the markov blancket model is :\n")
tab_mb



fourfoldplot(tab_mb, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_mb<-sum(pred_mb==test_data[,2])/dim(test_data)[1]
cat("The accuracy of markov blancket model is:",accuracy_mb*100,"%")


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

naive.dag<- model2network(string)



pred_naive<-prediction_func(naive.dag,train_data,test_data,"mle", 2,"S")

tab_naive<-prop.table(table(test_data[,2],pred_naive,
                         dnn=c("true", "prediction")))

cat("The confusion matrix for the naive bayes model is :\n")
tab_naive


fourfoldplot(tab_naive, color = c(sample(colours(),1), sample(colours(),1)),
             conf.level = 0, margin = 1, main = "Confusion Matrix")


accuracy_naive<-sum(pred_naive==test_data[,2])/dim(test_data)[1]
cat("The accuracy of the naive bayes model is:",accuracy_naive*100,"%")


library(knitr)
accuracy_df<-rbind(accuracy_bf,accuracy_true,accuracy_mb,accuracy_naive)
accuracy_df<-as.data.frame(accuracy_df) 
rownames(accuracy_df)<-c("Accuracy BF","Accuracy True","Accuracy MB","Accuracy Naive")
colnames(accuracy_df)<-"Summary Table"

kable(accuracy_df, digits = 2, caption = "Accuracy Table.")

# Code to plot with graphviz ----------------------------------------------

arc.set <- arcs(naive.dag)

highlight.opts<-list(nodes=colnames(asia),col =sample(colors(),1), fill = sample(colors(),1),arcs=arc.set,lty=5,lwd=2)
graphviz.plot(naive.model,highlight = highlight.opts,layout = "neato")


