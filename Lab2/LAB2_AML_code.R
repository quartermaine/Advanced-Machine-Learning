
# Libraries ---------------------------------------------------------------

library(HMM)
library(seqHMM) #another library for HMMs
library(diagram)
library(entropy)
library(dplyr)


# -------------------------------------------------------------------------

# Question 1 --------------------------------------------------------------

names=1:10
States = c() # initialize the hidden states z
for(i in 1:10){
  States=c(States,paste("State_",names[i],sep=""))
}

States # print the vector with the States x
Symbols = letters[1:10] # initialie the observed states 
# transmission matrix
transProbs = matrix(0,10,10) # transmisison matrix-hidden
diag(transProbs) = 0.5
diag(transProbs[,-1]) = 0.5 ; transProbs[10,1]=0.5
emissionProbs = toeplitz(c(0.2,0.2,0.2,rep(0,5),0.2,0.2)) # emission matrix-observable
# print transition probabilities
TP=as.data.frame(transProbs) ; colnames(TP)=States ; knitr::kable(TP) 
# print emission probabilities
EP=as.data.frame(emissionProbs) ; colnames(EP)=Symbols ; knitr::kable(EP)
# initialize HMM 
hmm = initHMM(States, Symbols, 
              transProbs = transProbs, emissionProbs = emissionProbs)
# using the diagram package 
plotmat(transProbs,box.size = 0.05,box.col = "lightyellow",
        arr.length = 0.4,main="Transition - Matrix") # plot the transition-matrix 
plotmat(emissionProbs,box.size = 0.05,box.col = "lightyellow",
        arr.length = 0.4,main="Emission - Matrix") # plot the emission-matrix 

# -------------------------------------------------------------------------

# Question 2 --------------------------------------------------------------

nSim=100
sim = simHMM(hmm, nSim) # simulate 100 timesteps

# -------------------------------------------------------------------------

# Question 3 --------------------------------------------------------------

sim.obs=sim$observation # discard the hidden states 
# filter 
f = forward(hmm, sim.obs)
f_tab= prop.table(exp(f),margin=2)  # marginalize over the columns so thats why we have 2 to prop table
# smooth
b = backward(hmm, sim.obs)
#b_tab=prop.table(exp(b),2) 
smooth_= posterior(hmm,sim.obs) # or we can use the prop.table(exp(f)*exp(b),2)
# viterbi algorithm for most probable path
vit = viterbi(hmm, sim.obs) 

# -------------------------------------------------------------------------

# Question 4 --------------------------------------------------------------

# compute accuracy for filter
t.index=apply(f_tab,2,which.max) # marginalize over the columns 
t=sapply(t.index,function(y){States[y]})
# accuracy for filter
ac1=sum(t==sim$states)/length(sim$states)
# compute accuracy for smooth
t1.index=apply(smooth_,2,which.max)
t1=sapply(t1.index,function(y){States[y]})
# accuracy for smooth
ac2=sum(t1==sim$states)/length(sim$states)
# compute accuracy for viterbi
ac3=sum(vit==sim$states)/length(sim$states)

accuracy_tab=cbind(ac1,ac2,ac3) 
colnames(accuracy_tab)=c("Filter","Smooth","Viterbi")
rownames(accuracy_tab)="Accuracy"
accuracy_tab

knitr::kable(as.data.frame(accuracy_tab))

# -------------------------------------------------------------------------

# Question 5 --------------------------------------------------------------

accuracy_func=function(user_hmm,number_sim){
  # --function that returns the filter ,smooth,viterbi accuracy for a hmm , and number of time steps --
  simHMM = simHMM(user_hmm,number_sim)
  sim_hmm_obs=simHMM$observation # discard the hidden states 
  # filter  
  forward_pass = forward(user_hmm, sim_hmm_obs)
  forward_table= prop.table(exp(forward_pass),2)  # marginalize over the cols Margin=2 to prop table
  # smooth
  back_pass = backward(user_hmm, sim_hmm_obs)
  # back_table=prop.table(exp(back_pass),2)
  smooth_= posterior(user_hmm,sim_hmm_obs)  # or prop.table(exp(forward_pass)*exp(back_pass),2)
  # viterbi algorithm for most probable path
  viterbi_res = viterbi(user_hmm, sim_hmm_obs)
  
  # compute accuracy for filter
  t.index=apply(forward_table,2,which.max) # marginalize over the columns 
  t=sapply(t.index,function(y){user_hmm$`States`[y]})
  # accuracy for filter
  filter_accuracy=sum(t==simHMM$states) / length(simHMM$states)
  # compute accuracy for smooth
  t1.index=apply(smooth_,2,which.max)
  t1=sapply(t1.index,function(y){user_hmm$`States`[y]})
  # accuracy for smooth
  smooth_accuracy=sum(t1==simHMM$states) / length(simHMM$states)
  # compute accuracy for viterbi
  viterbi_accuracy=sum(viterbi_res==simHMM$states) / length(simHMM$states)
  # compute entropy for filter
  filter_entropy=apply(forward_table, 2, FUN = entropy.empirical)
  
  return(list(simHMM,accuracies=
                c(filter_accuracy,smooth_accuracy,viterbi_accuracy),
              entropy=list(filt_ent=filter_entropy,smooth_ent=smooth_entropy)))
}

# replicate the results from accuracy function 10 times with 100 timespters every time
nreps=100
rep_res=replicate(nreps,accuracy_func(hmm,100))
# create a matrix with the results from the rep_res --the second row of rep_res
results_mat=matrix(unlist(rep_res[2,]), nrow=3,ncol = nreps, byrow = F)
# colnames-rownames
colnames(results_mat)=paste("replicate n=",1:nreps)
rownames(results_mat)=c("filter_accuracy","smooth_accuracy",
                    "viterbi_accuracy")
results_mat
# make table 
# knitr::kable(as.data.frame(results_mat)) # don't fit page
# make tables 
t1 <- kable(as.data.frame(results_mat[,1:10]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t2 <- kable(as.data.frame(results_mat[,11:20]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t3<-kable(as.data.frame(results_mat[,21:30]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t4<-kable(as.data.frame(results_mat[,31:40]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t5<-kable(as.data.frame(results_mat[,41:50]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t6<-kable(as.data.frame(results_mat[,51:60]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t7<-kable(as.data.frame(results_mat[,61:70]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t8<-kable(as.data.frame(results_mat[,71:80]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t9<-kable(as.data.frame(results_mat[,81:90]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t10<-kable(as.data.frame(results_mat[,91:100]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")

t1 ; t2 ; t3 ; t4 ; t5 ; t6 ; t7 ; t8 ; t9 ; t10

# plot of the accuracies 
cols=c("gray30","mediumorchid","chartreuse" )
plot(1:nreps,results_mat[2,],type="l",ylim=c(0,0.9),col=cols[1],lwd=2,
     main=paste("Smooth~Filter~Viterbi for\n ",nreps," replications for 100 timesteps "),
     panel.first = grid(25,25),ylab = "Accuracy",xlab="replications")
lines(1:nreps,results_mat[1,],col=cols[2],lwd=2)
lines(1:nreps,results_mat[3,],col=cols[3],lwd=2)
legend("bottomleft",legend = c("Smooth","Filter","Viterbi"),
       col = cols,lty=1,lwd=2)

# -------------------------------------------------------------------------

# Question 6 --------------------------------------------------------------

sample_grid=seq(10,300,5) # generate grid of timesteps
# apply the accuracy funciton to the grid 
accuracy_mat=sapply(sample_grid,function(y){accuracy_func(hmm,y)})
# create a matrix with the results from the accuracy_mat --the second row of raccuracy_mat
res_acc_mat=matrix(unlist(accuracy_mat[2,]), nrow=3,ncol = length(sample_grid), byrow = F)
# row-column names
colnames(res_acc_mat)=paste("nSample=",sample_grid)
rownames(res_acc_mat)=c("filter_accuracy","smooth_accuracy","viterbi_accuracy")
#  print the matrix
res_acc_mat
# construct a table 
# knitr::kable(data.frame(res_acc_mat)) # dont't fit page

t11 <- kable(as.data.frame(res_acc_mat[,1:10]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t12 <- kable(as.data.frame(res_acc_mat[,11:20]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t13<-kable(as.data.frame(res_acc_mat[,21:30]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t14<-kable(as.data.frame(res_acc_mat[,31:40]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t15<-kable(as.data.frame(res_acc_mat[,41:50]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")
t16<-kable(as.data.frame(res_acc_mat[,51:59]), format = "latex", booktabs = TRUE)%>%
  kable_styling(latex_options = "scale_down")

t11 ; t12 ; t13 ; t14 ; t15 ; t16

# plot of the accuracies 
cols1=c("grey73","deeppink2","peachpuff4")
plot(sample_grid,res_acc_mat[2,],type="l",ylim=c(0,0.9),col=cols1[1],lwd=2,
     main=paste("Smooth~Filter~Viterbi Accuracy \nfor n=",sample_grid[1],":",rev(sample_grid)[1]),
     panel.first = grid(25,25),ylab = "Accuracy")
lines(sample_grid,res_acc_mat[1,],col=cols1[2],lwd=2)
lines(sample_grid,res_acc_mat[3,],col=cols1[3],lwd=2)
legend(x=10,y=0.18,legend = c("Smooth","Filter","Viterbi"),
       col = cols1,lty=1,lwd=2)


sample_indx<-sample(1:50,9) ; sample_indx

par(mfrow=c(3,3))
for(i in sample_indx){
plot(accuracy_mat[,i]$entropy,type="o",col=sample(colors(),1),pch=20,
     ylab = "Entropy",xlab="nSample",panel.first = grid(25,25),
     main=paste("Entropy for nSampe=", sample_grid[i]," \nfor 10 timesteps"))
}

# take a sample with 100 timesteps from 4
samp1=rep_res[,60]

plot(samp1$entropy$filt_ent,main = "Entropy filtered and smoothed",
     type = "l", col = "mediumspringgreen",
     ylab = "Entropy", ylim = c(0,2),panel.first = grid(25,25))
lines(samp1$entropy$smooth_ent, col = "mediumorchid1")
legend(x = "topright", lty=1,legend = c("Smoothed","Filtered"),
       col = c("mediumspringgreen","mediumorchid1"))

# -------------------------------------------------------------------------

# Question 7 --------------------------------------------------------------

rndSample=rep_res[1,sample(1:dim(rep_res)[2],1)]

post=posterior(hmm,rndSample[[1]]$observation)
prob_step101=transProbs%*%post[,100]

df=data.frame(prob_step101) ; rownames(df)=States

knitr::kable(df)

# -------------------------------------------------------------------------

# Other  ------------------------------------------------------------------

par(mfrow=c(5,5))
for(i in 5:14){
plot.ts(f_tab[,i])
pch1=as.character(1:10)
points(f_tab[,i], pch=pch1, cex=1.25, font=4, col=1:10)
legend("topleft",legend = States,col=1:10,pch=pch1)
}







