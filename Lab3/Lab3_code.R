
# Question1 ---------------------------------------------------------------

set.seed(1234)

transition_func=function(zt,sd_trans){
  # function that samples a component and returns 
  # random sample from transition model
  U=sample(1:3,1)
  #
  if(U==1){
    zt= rnorm(1,zt,sd_trans)
    }else if(U==2){
      zt = rnorm(1,zt+1,sd_trans)
    }else{
      zt= rnorm(1,zt+2,sd_trans)
    }
  return(zt)
}

emisssion_func=function(zt,sd_emiss){
  # function that samples a component and returns 
  # random sample from emission model
  U=sample(1:3,1)
  #
  if(U==1){
      xt= rnorm(1,zt,sd_emiss)
    }else if(U==2){
      xt = rnorm(1,zt+1,sd_emiss)
    }else{
      xt= rnorm(1,zt+2,sd_emiss)
    }
  return(xt)
}
  

make_model=function(nTimes,upper_limit,sd_trans,sd_emiss){
  # function to simulate states and observations 
  zt=double(nTimes)
  #
  xt=double(nTimes)
  #
  for(time in 1:nTimes){
    if(time==1){
      zt[time]=runif(1,0,upper_limit) 
    }else{
    zt[time]=transition_func(zt[(time-1)],sd_trans)
    }
    xt[time]=emisssion_func(zt[time],sd_emiss)
  }
  return(list(states=zt,observations=xt))
}


x.sim=make_model(100,100,1,1)

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(1:100,x.sim$states,col="darkcyan",pch=3,
     panel.first = grid(25,25),ylab="",xlab="timesteps",
     main="Plot of States and Observations")
lines(1:100,x.sim$observations,col="red",lwd=2)
legend("topleft",legend=c("States","Observations"),
       col = c("darkcyan","red"),lty=c(NA,1),pch=c(3,NA))
#
d_states=density(x.sim$states)
plot(y=d_states$x,x=d_states$y,col="cornflowerblue",lwd=3,
     main = "Density States",type="l")
polygon(y=d_states$x,x=d_states$y,col="red",density=15)
d_observations=density(x.sim$observations)
plot(y=d_observations$x,x=d_observations$y,col="cornflowerblue",lwd=3,
     main="Density Observations",type="l")
polygon(y=d_observations$x,x=d_observations$y,col="red",density=15)


# Partical Filter ------------

observations=x.sim$observations # observations of the simulation
states=x.sim$states # states from the simulation

emission_model=function(zt,mu,sd_emiss){
  # function that returns the probability-weight from the emission model
  result=( dnorm(zt,mu,sd_emiss)+dnorm(zt,(mu-1),sd_emiss)+dnorm(zt,(mu+1),sd_emiss) ) / 3
  return(result) 
}

particle_func=function(obs,nTimes,nParts,sd_trans,sd_emiss,I){
  # function that returns the particles and weights given 
  # observations, timesteps, #particles, sd of transition,sd of emission,
  # switch betwwen update and not update weights (I=1 means update) 
  init_particles=runif(nTimes,0,100) 
  particle_mat=matrix(0,nrow=nTimes+1,
                      ncol=length(init_particles))
  particle_mat[1,]=init_particles
  X_bar=double(nTimes)
  wt=matrix(0,nrow=nTimes,ncol=nParts)
  # 
  for(i in 1:(nTimes)){
    #cat("---------iteration: ",i)
    #cat("\n")
    for( j in 1:length(init_particles)){
      X_bar[j]=transition_func(particle_mat[i,j],sd_trans)
      wt[i,j]=emission_model(obs[i], X_bar[j], sd_emiss)
      #X_bar[j]=future_particle
    }
    if(I==1){
    wt[i,]=wt[i,]/sum(wt[i,])
    }else{
      wt[i,]=rep(1,nTimes)
    }
    # 
    particle_mat[(i+1),]=sample(X_bar,100,prob=wt[i,],replace = T)
    
  }
  return(list(particles=particle_mat,weigts=wt))  
  
}

plot_st=function(true_states,predicted_states,sd_title){
  # plot function for the true and predicted staates for given emisison sd 
  plot(true_states, pch="x",xlab="step", 
       ylab="timestep", 
       main=(paste("Particle filter preds \nwith emission sd=",sd_title)),lwd=2,
       panel.first = grid(25,25))
  lines(predicted_states, type="l", col="blue")
  legend("topleft",legend=c("true","predicted"), 
         col=c("black","blue"),pch=c("x",NA), lty=c(NA,1))
}

plot_parts=function(particles,timestep,obs,states){
  # plot for the particles for given timestep
  plot(particles[timestep,],rep(0, length(particles[timestep,])),pch="|",yaxt='n',
     main = paste("Time step", timestep), xlab = "position", ylab = "")
  points(obs[timestep], 0,col="red",pch="S",cex=1.1)
  points(states[timestep], 0,col="blue",pch="O")
  
  }


sd1=1 # set the emission sd
# sample from particle func
res=particle_func(obs=observations,
                  nTime=100,nParts=100,sd_trans =1,sd_emiss=sd1,I=1)
# calculate the row mean for each timestep and discard the first (initial)
pred_states<-rowMeans(res$particles)[-1]

# plot for states and filter particle preds
plot_st(x.sim$states,pred_states,sd1)
# plot of particles for different timesteps
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE),heights=c(3, 3),
       respect=F)
plot_parts(res$particles,1,observations,states)
plot_parts(res$particles,sample(2:99,1),observations,states)
plot_parts(res$particles,sample(2:99,1),observations,states)
plot_parts(res$particles,100,observations,states)
par(mai=c(0,0,0,0))
plot.new()
legend("center",legend = c("particle","state","observation",paste("Emission sd=",sd1)),
       col = c("black","red","blue",NA),pch=c("|","S","O",NA))


# -------------------------------------------------------------------------

# Question 2 --------------------------------------------------------------

# For emission sd = 5 -------------
sd5=5 # emission sd =5 
# 
res_sd5=particle_func(obs=observations,
                      nTime=100,nParts=100,sd_trans =1,sd_emiss=sd5,I=1)
# filter particles predictions with emission sd=5
pred_states_sd5<-rowMeans(res_sd5$particles)[-1]
# 
plot_st(x.sim$states,pred_states_sd5,sd5)
# 
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE),heights=c(3, 3),
       respect=F)
plot_parts(res_sd5$particles,1,observations,states)
plot_parts(res_sd5$particles,sample(2:99,1),observations,states)
plot_parts(res_sd5$particles,sample(2:99,1),observations,states)
plot_parts(res_sd5$particles,100,observations,states)
par(mai=c(0,0,0,0))
plot.new()
legend("center",
       legend = c("particle","state","observation",paste("Emission sd=",sd5)),
       col = c("black","red","blue",NA),pch=c("|","S","O",NA))

# For emission sd = 50 --------------

sd50=50 # emission sd =50
res_sd50=particle_func(obs=observations,
                       nTime=100,nParts=100,sd_trans =1,sd_emiss=sd50,I=1)
# filter particle predictions eith emission sd=50
pred_states_sd50<-rowMeans(res_sd50$particles)[-1]
# 
plot_st(x.sim$states,pred_states_sd50,sd50)
# 
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE),heights=c(3, 3),
       respect=F)
plot_parts(res_sd50$particles,1,observations,states)
plot_parts(res_sd50$particles,sample(2:99,1),observations,states)
plot_parts(res_sd50$particles,sample(2:99,1),observations,states)
plot_parts(res_sd50$particles,100,observations,states)
par(mai=c(0,0,0,0))
plot.new()
legend("center",
       legend = c("particle","state","observation",paste("Emission sd=",sd50)),
       col = c("black","red","blue",NA),pch=c("|","S","O",NA))

# -------------------------------------------------------------------------

# Question 3 --------------------------------------------------------------
# use non updated weights for particle filtering with emission sd=1
res_sd1_w=particle_func(obs=observations,
                        nTime=100,nParts=100,sd_trans =1,sd_emiss=sd1,I=2)
# make predictions for each time step(-->rows)
pred_states_sd1_w<-rowMeans(res_sd1_w$particles)[-1]
# 
plot_st(x.sim$states,pred_states_sd1_w,sd1)
# 
layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE),heights=c(3, 3),
       respect=F)
plot_parts(res_sd1_w$particles,1,observations,states)
plot_parts(res_sd1_w$particles,sample(2:99,1),observations,states)
plot_parts(res_sd1_w$particles,sample(2:99,1),observations,states)
plot_parts(res_sd1_w$particles,100,observations,states)
par(mai=c(0,0,0,0))
plot.new()
legend("center",
       legend = c("particle","state","observation",paste("Emission sd=",sd1)),
       col = c("black","red","blue",NA),pch=c("|","S","O",NA))
  
# -------------------------------------------------------------------------
# END 
# Other -------------------------------------------------------------------
# another implementation of particle filter starting from t=2,...
particle_func_v2=function(obs,nTimes,nParts,sd_trans,sd_emiss,I){
  # function that returns the particles and weights given 
  # observations, timesteps, #particles, sd of transition,sd of emission,
  # switch betwwen update and not update weights (I=1 means update) 
  init_particles=runif(nTimes,0,100) 
  particle_mat=matrix(0,nrow=nTimes+1,
                      ncol=length(init_particles))
  particle_mat[1,]=init_particles
  X_bar=double(nTimes)
  wt=matrix(0,nrow=nTimes+1,ncol=nParts)
  # 
  for(i in 2:(nTimes+1)){
    #cat("---------iteration: ",i)
    #cat("\n")
    for( j in 1:length(init_particles)){
      X_bar[j]=transition_func(particle_mat[(i-1),j],sd_trans)
      wt[i,j]=emission_model(obs[(i-1)], X_bar[j], sd_emiss)
      #X_bar[j]=future_particle
    }
    if(I==1){
      wt[i,]=wt[i,]/sum(wt[i,])
    }else{
      wt[i,]=rep(1,nTimes)
    }
    # 
    particle_mat[i,]=sample(X_bar,100,prob=wt[i,],replace = T)
    
  }
  return(list(particles=particle_mat,weigts=wt))  
  
}


# sd=1 -----------------------------------

sd1=1
res21=particle_func_v2(obs=observations,
                       nTime=100,nParts=100,sd_trans =1,sd_emiss=sd1,I=1)

pred_states_21<-rowMeans(res21$particles)[-1]

# plot for states and filter particle preds
plot_st(x.sim$states,pred_states_21,sd1)

# sd=5 -----------------------------------

sd1=5
res22=particle_func_v2(obs=observations,
                       nTime=100,nParts=100,sd_trans =1,sd_emiss=sd1,I=1)

pred_states_22<-rowMeans(res22$particles)[-1]

# plot for states and filter particle preds
plot_st(x.sim$states,pred_states_22,sd1)

# sd=50 -----------------------------------

sd1=50
res23=particle_func_v2(obs=observations,
                       nTime=100,nParts=100,sd_trans =1,sd_emiss=sd1,I=1)

pred_states_23<-rowMeans(res23$particles)[-1]

# plot for states and filter particle preds
plot_st(x.sim$states,pred_states_23,sd1)

# -------------------------------------------------------------------------





















  
