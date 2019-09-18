
## SCRIPT FOR CREATING TABLE OF CUMULATIVE SURVIVAL, BASED ON PATTERSON ET AL. 2002, AND MODIFIED BY G. CHAPRON. 
## RUN THE FOLLOWING SCRIPT, THEN ATTACH THE TABLE OF SURVIVAL TO THE CURRENT DATA FRAME. THE FIELDS MUST HAVE THE NAMES
## "ent", "exit", "event", and "cause", ALL IN LOWER CASE LETTERS. ONCE THE SCRIPT IS RUN, THERE WILL BE A FUNCTION IN YOUR 
## CURRENT CHAPTER CALLED "cause.survival". UAGE OF THIS FUNCTION IS cause.survival(table, causenumber), WHERE table IS THE 
## TABLE WITH THE DATA, AND CAUSENUMBER AN INTEGER USED TO IDENTIFY THE CAUSE OF MORTALITY. 

library(survival)


"cause.survival" = function(table, p)

{
	assign("p", p)
	
	# create tables to hold the results of GKM survival estimates for all events and for cause specific event
	temp.all <- summary(survfit(Surv(ent, exit, event) ~ 1, data = table))										## NB: NOTE THAT "CAUSE" IS CHANGED TO "cause3; CHANGE TO cause TO MAKE ORIGINAL FORMULATIONS!!!!!
	temp.s <- summary(survfit(Surv(ent, exit, cause3 == p) ~ 1, data = table))
	
	#combine the two tables so survival of all events can be combined with those of the cause specific events
	s.df <- data.frame(time = temp.s$time, n.event = temp.s$n.event, n.risk = temp.s$n.risk, survival = temp.s$surv)
	all.df <- data.frame(time = temp.all$time, n.event = temp.all$n.event, n.risk = temp.all$n.risk, survival = temp.all$surv)
	all.s.df <- merge(all.df, s.df, by.x = "time", by.y = "time", all.x = F, suffixes = c(".all", ".s"))					## Creates NA where s.df holds no data
	assign("n", all.s.df)
	x <- length(n[,1])
	
	#create temporary placeholders for the calculation of the mortality rate and the cause-specific cumulative incidence func.
	tmp.string <- numeric(x)					
	tmp.string2 <- numeric(x)					
	t <- 1

	#cycle through the records of the table, including all events to calculate mortality rate and CIF
	while(t <= x) {
		tmp.string[1] <- n$n.event.s[1]/n$n.risk.s[1]
		if (t == 1) tmp.string[t] <- NA else tmp.string[t] <- (n$survival.all[t-1] * n$n.event.s[t])/n$n.risk.s[t]			## The key to the estimation
		if(is.na(tmp.string[t])) tmp.string2[t] <- NA else tmp.string2[t] <- sum(tmp.string[1:t], na.rm = T)
		t = t + 1
	}

	
	MORT <- data.frame(mort.rate = tmp.string)
	CIF2 <- data.frame(CIF = tmp.string2)
	CIF.s.all <- cbind(all.s.df, MORT, CIF2)

# Calculate the variance, standard error and the Confidence Intervals around CIF

SE <- numeric(x)					
totvar.t <- numeric(x)				

#Reset all temporary variables
t <- 1
j <- 1
Ij <- 0
cumvar.p1 <- 0
cumvar.p2 <- 0
cumvar.p3 <- 0

#loop for the total number of records
while (t <= x) 
{
   It <- CIF.s.all$CIF[t]
   if(is.na(It)) {		
		CIF.s.all$cumvar[t] <- "NA"
		CIF.s.all$StdErr[t] <- "NA"
		CIF.s.all$CI.u[t] <- "NA" 
		CIF.s.all$CI.l[t] <- "NA"
		t = t + 1
	}
   else 
	{
		while (j < t) 
		{
		if(is.na(CIF.s.all$CIF[j]))
			Ij <- Ij			
	     else	
			Ij <- CIF.s.all$CIF[j]		
       cumvar.p1 <- cumvar.p1 + (It - Ij)^2 * (CIF.s.all$n.event.all[j]/(CIF.s.all$n.risk.all[j] * (CIF.s.all$n.risk.all[j] - CIF.s.all$n.event.all[j])))

		if(!is.na(CIF.s.all$CIF[j]))
		{
			if(j == 1)
				Sj3 <- 1
			else
			Sj3 <- CIF.s.all$survival.all[j-1]		
			Ijc <- CIF.s.all$CIF[j]
			cumvar.p3 <- cumvar.p3 + (It - Ijc)*(Sj3)*(CIF.s.all$n.event.all[j] / (CIF.s.all$n.risk.all[j])^2)
		}
		j <- j + 1
   		}

		if (t == 1) 
			Sj2 <- 1  
		else
			Sj2 <- CIF.s.all$survival.all[t-1]  
					
  		cumvar.p2 <- (Sj2)^2 * (((CIF.s.all$n.event.all[t])*(CIF.s.all$n.risk.all[t] - CIF.s.all$n.event.all[t]))/(CIF.s.all$n.risk.all[t])^3) + cumvar.p2

		#total all three components of the variance equation to get the final variance,  generate std. err and confidence intervals
		#Assign all results to the output table
		totvar.t[t] <- cumvar.p1 + cumvar.p2 - (2 * cumvar.p3)
		CIF.s.all$cumvar[t] <- totvar.t[t] 
		SE[t] <- sqrt(totvar.t[t])
		CIF.s.all$StdErr[t] <- SE[t]
		CIF.s.all$CI.u[t] <- CIF.s.all$CIF[t] + (1.645 * SE[t]) 
		CIF.s.all$CI.l[t] <- CIF.s.all$CIF[t] - (1.645 * SE[t]) 
		t = t + 1
		j <- 1
   }

cumvar.p1 <- 0
cumvar.p3 <- 0
Ij <- 0
It <- 0
}

#Variance calculations end here ----------
	
	return(CIF.s.all)
}