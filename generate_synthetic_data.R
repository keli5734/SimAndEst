#todo: set to detect on specivic days: 1,4,8,11...
sim.hh.func.fixed <- function(N,
                              hh.size = sample(3:7,1),
                              tests.per.week = 2,
                              
                              p.comm.base.infant.fix = 0.001,
                              p.comm.multiplier.sibling = 1,
                              p.comm.multiplier.parent = 1,
                              p.comm.multiplier.elder = 1,
                              
                              p.hh.base.infant = 0.1,
                              p.hh.multiplier.sibling = 1,
                              p.hh.multiplier.parent = 1,
                              p.hh.multiplier.elder = 1,
                              
                              p.imm.base.sibling = 0.5,
                              p.imm.base.parent = 0.9, 
                              p.imm.base.elder = 0.6,
                              
                              partial.immunity.infant = 0.8,
                              partial.immunity.sibling = 0.5,
                              partial.immunity.parent = 0.2,
                              partial.immunity.elder = 0.3,
                              
                              delay = T,
                              duration.latent = 4,
                              duration.infect.inf = 10,
                              multiplier.dur.sibpar = 0.5,
                              #What is probability that an infected person is detected on each day 
                              p.detect = 0.999, #what is overall probability that an infection is detected?
                              
                              amplitude = 2.6581,
                              phase = -0.408
                              
){
  
  #p.detect.day <- 1 - exp(log(1-p.detect)/duration_infectious) 
  
  
  time.steps = -as.integer(as.Date("2024-09-21") - as.Date("2025-04-17"))  #272 # ~8 month study from October-April 212 days (30.3*7)
  
  test.days1 <- seq(from=1, by=7, to=time.steps) 
  test.days2 <- seq(from=5, by=7, to=time.steps) 
  test.days3 <- seq(from=3, by=7, to=time.steps) 
  
  if(tests.per.week==1){
    test.days = test.days1
  }else if(tests.per.week==2){
    test.days = c(test.days1,test.days2)
  }else if(tests.per.week==3){
    test.days = c(test.days1,test.days2,test.days3)
  }
  test.days = sort(test.days)
  baseline.test.days <- test.days
  
  # household composition -------------------------------------------------
  # one infant, at least two adults, optional older siblings, up to two elders
  n.adult.base <- 2
  remaining <- hh.size - (1 + n.adult.base)
  n.elder <- ifelse(remaining > 0, sample(0:min(2, remaining), 1), 0)
  remaining <- remaining - n.elder
  n.sib <- max(remaining, 0)
  
  hh.roles <- c("infant",
                rep("adult", n.adult.base),
                rep("sibling", n.sib),
                rep("elder", n.elder))
  
  
  latent <- array(0, dim=c(4,time.steps,hh.size))
  
  infectious <- array(0, dim=c(4,time.steps,hh.size))
  
  immune <- matrix(NA, nrow=time.steps, ncol=hh.size)
  #vax <- matrix(NA, nrow=time.steps, ncol=hh.size)
  
  
  immune[1, hh.roles=="infant"] <- 0 # all infants are susceptible
  immune[1, hh.roles %in% c("adult")] <- rbinom(sum(hh.roles %in% c("adult")), 1, p.imm.base.parent)
  immune[1, hh.roles %in% c("elder")] <- rbinom(sum(hh.roles %in% c("elder")), 1, p.imm.base.elder)
  immune[1, hh.roles=="sibling"] <- rbinom(sum(hh.roles=="sibling"), 1, p.imm.base.sibling)
  
  
  
  p.comm.init <- p.comm.base.infant.fix * exp(amplitude * cos((2*pi*(1+40*7)/365.25) + phase))
  p.comm.init.vec <- ifelse(hh.roles=="infant", p.comm.init,
                            ifelse(hh.roles %in% c("adult"), p.comm.init * p.comm.multiplier.parent, 
                                   ifelse(hh.roles == "elder", p.comm.init * p.comm.multiplier.elder,
                                   p.comm.init * p.comm.multiplier.sibling)))
  latent[1,1,] <- ifelse(immune[1,] >0, 0, rbinom(hh.size, 1, p.comm.init.vec))   #1st of the 4 latent classes
  
  infectious[1,1,] <- 0 #first of the 4 infectious classes
  
  
  
  
  #vax[1,] <- 0
  
  for(i in 2:time.steps){
    
    p.comm.base.infant = p.comm.base.infant.fix * exp(amplitude * cos((2*pi*(i+40*7)/365.25) + phase))
    
    for(j in 1:hh.size){
      
      #Contribution of community infection
      if(hh.roles[j]=="infant"){
        p.comm = p.comm.base.infant
        partial.immunity = partial.immunity.infant
        duration.infect = duration.infect.inf
        
      } else if(hh.roles[j] %in% c("adult")){
        p.comm = p.comm.base.infant * p.comm.multiplier.parent
        partial.immunity = partial.immunity.parent
        duration.infect = duration.infect.inf*multiplier.dur.sibpar
      }else if(hh.roles[j] %in% c("elder")){
        p.comm = p.comm.base.infant * p.comm.multiplier.elder
        partial.immunity = partial.immunity.elder
        duration.infect = duration.infect.inf*multiplier.dur.sibpar
      } else{
        p.comm = p.comm.base.infant * p.comm.multiplier.sibling
        partial.immunity = partial.immunity.sibling
        duration.infect = duration.infect.inf*multiplier.dur.sibpar
        
      }
      
      
      # contribution of household infection by type of infectious contact.
      N.infectious.infant  <- sum(infectious[,(i-1),hh.roles=="infant"])
      N.infectious.sibling <- sum(infectious[,(i-1),hh.roles=="sibling"])
      N.infectious.parent <- sum(infectious[,(i-1),hh.roles %in% c("adult")])
      N.infectious.elder <- sum(infectious[,(i-1),hh.roles %in% c("elder")])
      
      
      # log.p.comm <- log(p.comm.base.infant) + vax.status*log((1-VE_susceptibility/100))
      # p.comm <- exp(log.p.comm)
      
      #Contribution of HH infection from unvaccinated contact, per contact
      # log.p.hh.unvax <- log(p.hh.base.infant) + vax.status*log((1-VE_susceptibility/100))
      # p.hh.unvax <- exp(log.p.hh.unvax) 
      
      #Contribution of HH infection from vaccinated contact, per contact
      # log.p.hh.vax <- log(p.hh.base.infant) + vax.status*log((1-VE_susceptibility/100)) +
      #   log(1-VE_txn/100)
      # p.hh.vax <- exp(log.p.hh.vax)
      
      #need to exponentiate these by the number of infectious HH members in each category
      
      log.qi <- log(1 - p.comm) + N.infectious.infant*log(1-p.hh.base.infant) + 
        N.infectious.sibling*log(1-p.hh.base.infant* p.hh.multiplier.sibling) + 
        N.infectious.parent*log(1-p.hh.base.infant * p.hh.multiplier.parent) + 
        N.infectious.elder*log(1-p.hh.base.infant * p.hh.multiplier.elder)
      
      qi <- exp(log.qi)
      
      pri <- (1- qi)  #probability of being infected at the time point
      
      susceptible <- (1-sum(infectious[,(i-1),j])) * (1-immune[(i-1),j]) * (1-sum(latent[,(i-1),j])) #is person susceptible at previous time step
      
      new.inf <-  rbinom(1,1,pri)*susceptible + rbinom(1,1,pri*partial.immunity)* immune[(i-1),j]* ((((i-1) - max(c(0,which(immune[1:(i-1),j]==0))))>25)+((i-1)<26))
      
      #is a person who was previously susceptible now latent? Is a person previously immune now latent?
      
      
      
      #If a person has been latent for <4 time points, stay latent in current time point
      # stay.latent <-   (sum(latent[1:(i-1),j])<4) 
      # 
      # #If a person has been infectious for <5 time points, stay infectious
      # stay.infectious <- (sum(infectious[1:(i-1),j])<5) 
      # 
      
      #Can use same value regardless of whether person was in latent 1,2 3 etc
      stay.latent <-   max(latent[,(i-1),j])*(rbinom(1,1, (1-1/duration.latent)^1/4))  #Transition to infectious if latent previously
      
      stay.infectious <- max(infectious[,(i-1),j])*(rbinom(1,1, (1-1/duration.infect)^1/4))
      
      #transition between latent states
      latent[1,i,j] <- (new.inf  + #if susceptible at previous time step, do they become infected (latent)
                          latent[1,(i-1),j]*stay.latent ) #If latent previously, do they stay latent?
      
      latent[2,i,j] <- latent[1,(i-1),j]*(1-stay.latent ) + #If latent previously, do they stay latent?
        latent[2,(i-1),j]*stay.latent 
      
      latent[3,i,j] <- latent[2,(i-1),j]*(1-stay.latent ) + #If latent previously, do they stay latent?
        latent[3,(i-1),j]*stay.latent 
      
      latent[4,i,j] <- latent[3,(i-1),j]*(1-stay.latent ) + #If latent previously, do they stay latent?
        latent[4,(i-1),j]*stay.latent 
      
      infectious[1,i,j]   <-     latent[4,(i-1),j]*(1-stay.latent) +  #if latent at previous time step, do they become infectious?
        infectious[1,(i-1),j]*(stay.infectious)     #if infectious at previous step, do they become immune
      
      infectious[2,i,j]   <-    infectious[1,(i-1),j]*(1-stay.infectious) +  #if latent at previous time step, do they become infectious?
        infectious[2,(i-1),j]*(stay.infectious)     #if infectious at previous step, do they become immune
      
      infectious[3,i,j]   <-    infectious[2,(i-1),j]*(1-stay.infectious) +  #if latent at previous time step, do they become infectious?
        infectious[3,(i-1),j]*(stay.infectious)
      
      infectious[4,i,j]   <-    infectious[3,(i-1),j]*(1-stay.infectious) +  #if latent at previous time step, do they become infectious?
        infectious[4,(i-1),j]*(stay.infectious)
      
      
      immune[i,j]   <-     infectious[4,(i-1),j] *(1-stay.infectious) + #if previously infectious, do they become immune?
        immune[i-1,j] - #if previously immune, stay immune
        immune[i-1,j]*new.inf #unless newly infected from the immune state.
      #vax[i,j]<- ifelse(vax.status[j]==1,1,0)
    }
  }
  
  #Observation probability with dynamic testing ---------------------------------
  latent2 <- apply(latent,c(2,3),sum) # sum across the 4 subclasses
  infectious2 <- apply(infectious,c(2,3),sum) # sum across the 4 subclasses
  
  test.schedule <- rep(FALSE, time.steps)
  test.schedule[baseline.test.days] <- TRUE
  detect.inf_full <- matrix(0, nrow=time.steps, ncol=hh.size)
  daily.testing <- FALSE
  neg.streak <- rep(0, hh.size)
  
  for(d in 1:time.steps){
    if(!test.schedule[d]) next
    detect.inf_full[d,] <- infectious2[d,]
    if(!daily.testing && any(detect.inf_full[d,] > 0)){
      daily.testing <- TRUE
      neg.streak <- ifelse(detect.inf_full[d,]==0,1,0)
      if(d < time.steps) test.schedule[(d+1):time.steps] <- TRUE
    } else if(daily.testing){
      neg.streak <- neg.streak + (detect.inf_full[d,]==0)
      neg.streak[detect.inf_full[d,]>0] <- 0
      if(all(neg.streak >= 2)){
        daily.testing <- FALSE
        if(d < time.steps) test.schedule[(d+1):time.steps] <- FALSE
        remaining <- baseline.test.days[baseline.test.days > d]
        test.schedule[remaining] <- TRUE
      }
    }
  }
  
  test.days <- which(test.schedule)
  detect.inf <- detect.inf_full[test.days,,drop=FALSE]
  test <- apply(infectious,c(3),sum) # sum across the 4 subclasses
  
  n.true.infection <- apply(infectious2, 2, function(x) {
    b<-rle(x)
    length(b$values[b$values == 1])
  })
  n.detected.infection <- apply(detect.inf, 2, function(x){
    b<-rle(x)
    length(b$values[b$values == 1])
  })
  
  
  first.infection.detected.start <- apply(detect.inf, 2, function(x){
    if(sum(x)>0){ #if that person had an infection
      z <- test.days[which(x==1)[1] ] #first detected infection day is the first of test = 1
    }else{
      z=-9999
    }
    return(z)
  })
  
  second.infection.detected.start <- apply(detect.inf, 2, function(x){
    b<-  rle(x)
    if(length(b$values[b$values == 1])>1){ #if that person had 2 infections
      z <- test.days[cumsum(b$lengths)[which(b$values == 1)[2]-1]+1]
    }else{
      z=-9999
    }
    return(z)
  })
  
  first.infection.detected.end <- apply(detect.inf, 2, function(x){
    if(sum(x)>0){
      b<- rle(x)
      z <- test.days[cumsum(b$lengths)[which(b$values == 1)[1]]]
    }else{
      z=-9999
    }
    return(z)
  })
  
  second.infection.detected.end <- apply(detect.inf, 2, function(x){
    b<-  rle(x)
    if(length(b$values[b$values == 1])>1){
      z <- test.days[cumsum(b$lengths)[which(b$values == 1)[2]]]
    }else{
      z=-9999
    }
    return(z)
  })
  
  
  #when is person actually INFECTED?
  first.infection.true.date <- apply(latent2, 2, function(x){
    if(sum(x)>0){
      z <- which(x==1)[1] 
    }else{
      z=-9999
    }
    return(z)
  })
  second.infection.true.date <- apply(latent2, 2, function(x){
    b<-  rle(x)
    if(length(b$values[b$values == 1])>1){
      z <- cumsum(b$lengths)[which(b$values == 1)[2]-1]+1
    }else{
      z=-9999
    }
    return(z)
  })
  
  
  
  first.infection.true.duration <- apply(infectious2, 2, function(x){
    if(sum(x)>0){
      b<-  rle(x)
      z <- b$lengths[which(b$values == 1)][1]
    }else{
      z=-9999
    }
    return(z)
  })
  
  second.infection.true.duration <- apply(infectious2, 2, function(x){
    b<-  rle(x)
    if(length(b$values[b$values == 1])>1){
      z <- b$lengths[which(b$values == 1)][2] 
    }else{
      z=-9999
    }
    return(z)
  })
  
  
  first.infection.infectious.day <- apply(infectious2, 2, function(x){
    if(sum(x)>0){
      z2 <- rle(x)
      z <- paste(which(x==1)[1]:(which(x==1)[1]+z2$lengths[which(z2$values==1)[1]]-1) , collapse=',')
    }else{
      z='9999'
    }
    return(z)
  })
  second.infection.infectious.day <- apply(infectious2, 2, function(x){
    z2<-  rle(x)
    if(length(z2$values[z2$values == 1])>1){
      z <- paste((cumsum(z2$lengths)[which(z2$values == 1)[2]-1]+1):(cumsum(z2$lengths)[which(z2$values == 1)[2]] ), collapse=',')
    }else{
      z='9999'
    }
    return(as.character(z))
  })
  
  
  
  
  
  out.df <- cbind.data.frame('indiv.index'=1:hh.size,
                             role=hh.roles,
                             n.true.infection,n.detected.infection,
                             first.infection.detected.start,first.infection.detected.end,first.infection.true.date, first.infection.true.duration,first.infection.infectious.day,
                             second.infection.detected.start,second.infection.detected.end,second.infection.true.date,second.infection.true.duration,second.infection.infectious.day,
                             'HH'=N)
  out.df$infected <- 1*( out.df$n.true.infection >0)
  out.df$detected.infected <- 1*( out.df$n.detected.infection >0)
  
  return(out.df)
  
  #out.ls <- list('latent'=latent,'infectious'=infectious,'immune'=immune ,'first.infection.date'=first.infection.date,'vax.status'=vax.status)
  
}