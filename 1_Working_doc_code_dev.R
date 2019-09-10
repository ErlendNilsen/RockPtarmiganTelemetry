

#### Descriptive statistics; 

library(tidyverse)
library(jsonlite)
library(lubridate)

d <- as_tibble(read_delim("data/Occurrence_data.txt", delim="\t", quote = ""))

d <- d %>% 
  mutate(dynamicProperties = purrr::map(dynamicProperties, ~ fromJSON(.) %>% as.data.frame())) %>% 
  unnest(dynamicProperties) %>%
  filter(organismName!=4272265) %>%
  mutate(Year=year(eventDate))

##############################
##### Descriptive stats; 

Deaths <- d %>% filter(state!="alive") %>%
          group_by(organismName, state) %>%
          count() %>%
          mutate(n=1) %>%
          group_by(state) %>%
          count()

###############################
##### NUMBER OF RELOCATIONS; 

temp1 <- d %>% group_by(organismName) %>%
        count()


############################
#### Setting up capture history for each bird: 

N_Ind <- d %>% count(organismName) %>%
          filter(n>1)

Death_cause <-  d %>% select(organismName, state, eventDate) %>%
                 filter(state=="dead (non-harvest mortality)" | state=="dead (harvest)") %>%
                 group_by(organismName) %>%
                 summarise(death=min(eventDate), cause=state[1]) %>%
                select(organismName, cause)

CH_first <- d %>% group_by(organismName, Year) %>%
            filter(event=="capture") %>%
            summarise(Capture=min(eventDate))
      
CH_lastAlive <- d %>% group_by(organismName, Year) %>%
  filter(state=="alive") %>%
  summarise(LastAlive=max(eventDate))

CH_firstAlive <- d %>% group_by(organismName, Year) %>%
  filter(state=="alive") %>%
  summarise(firstAlive=min(eventDate))

CH_dead <- d %>% group_by(organismName, Year) %>%
  filter(event!="capture" & state!="alive") %>%
  summarise(Dead=min(eventDate))


CH <- full_join(CH_first, CH_lastAlive) %>%
      full_join(., CH_firstAlive) %>%
      full_join(., CH_dead) %>%
      full_join(., Death_cause) %>%
      right_join(., N_Ind)

######

CH <- CH %>% mutate(cap_year=year(Capture), year_lastAlive=year(LastAlive), year_dead=year(Dead)) %>%
             mutate(we_cap=week(Capture), we_lastAlive=week(LastAlive), we_firstAlive=week(firstAlive), week_death=week(Dead)+1) %>%
             mutate(Entry_temp=we_firstAlive, Exit_temp=ifelse(is.na(week_death), we_lastAlive, round((week_death+we_lastAlive)/2)))  %>%
             mutate(event_temp=if_else(is.na(week_death), 0, 1)) %>%
             mutate(cause_temp=if_else(event_temp==0, "censored", cause)) 


## THIS SECTION WILL PREPARE DATA FOR THE SPECIFIC ANALYSIS CONDUCTED HERE, ASSUMING THTAT 
## THE CAPTURE HISTORY IS CONSTRUCTED FOR FEBRUARY 1ST - JULY 31. 
Censor_week <- week(dmy("01-08-2012"))
Begin_week <- week(dmy("01-02-2012"))

CH <-        CH %>% filter(!is.na(Entry_temp) & Year<2014) %>%
             mutate(entry0=Entry_temp-Begin_week, exit0=if_else(Exit_temp>Censor_week, Censor_week-Begin_week, Exit_temp-Begin_week)) %>%
             mutate(entry=if_else(entry0<1,1,entry0), exit=exit0+1) %>%
             mutate(event=if_else(Exit_temp>Censor_week, 0, event_temp), cause_temp2=if_else(Exit_temp>Censor_week, "censored", cause_temp)) %>%
             mutate(cause_spec=recode(cause_temp2, "censored"=0, "dead (harvest)"=1, "dead (non-harvest mortality)"=2))

## Adding age and sex 

AgeSex <- d %>% mutate(Year=year(eventDate)) %>%
          mutate(age=recode(lifeStage, "juvenile"=1, "adult"=2)) %>%
          group_by(organismName, Year, sex) %>%
          summarise(age=min(age))

CH <- left_join(CH, AgeSex)
## Data set ready for analysis

################################################################
## ESTIMATING SURVIVAL POBABILITY USING K-M MODELS, 
## TESTING FOR DIFFERNCE IN SURVIVAL BETWEEN YEARS, AGE AND SEX

library(survival)
library(ggplot2)
library(survminer)

M1 <- survfit(Surv(entry, exit, event)~1, data=CH)
summary(M1)


###
### PLOTTING SURVIVAL CURVE: 

ggsurvplot(M1, legend="none", 
           xlab="Week (since Feb. 13th)",
           risk.table = T, 
           tables.height = 0.2,
           risk.table.fontsize=3,
           tables.theme = theme_cleantable(),
           tables.y.text=FALSE,
           ggtheme = theme_bw(), 
           palette = c("#E7B800"))

### Testing the models; 
M2_a <- coxph(Surv(entry, exit, event)~1, data=CH)	 # A first model of overall hazard, not accounting for age, sex and differences in proportional hazards between seasons; 									
M2_b <- coxph(Surv(entry, exit, event)~as.factor(Year), data=CH)	 # A first model of overall hazard, not accounting for age, sex and differences in proportional hazards between seasons; 									
M2_c <- coxph(Surv(entry, exit, event)~as.factor(sex), data=CH)	 # A first model of overall hazard, not accounting for age, sex and differences in proportional hazards between seasons; 									
M2_d <- coxph(Surv(entry, exit, event)~as.factor(age), data=CH)	 # A first model of overall hazard, not accounting for age, sex and differences in proportional hazards between seasons; 									

###Assessing fit; 

prop_test1 <- cox.zph(M2_b)
prop_test2 <- cox.zph(M2_c)
prop_test3 <- cox.zph(M2_d)

#### Non-harvest mortality; 

CH2 <- CH %>% mutate(event=if_else(cause_spec==1, 0, event))

M3 <- survfit(Surv(entry, exit, event)~1, data=CH2)
summary(M3)


#########################################################
## Harvest mort in february; 
## Dataset for harvest mort in February; 
Censor_week_h <- week(dmy("01-03-2012"))-Begin_week

CH_h <- CH %>% mutate(exit2=if_else(exit>Censor_week_h, Censor_week_h, exit), event2=if_else(cause_temp2=="dead (harvest)", 1, 0)) %>%
  filter(entry<Censor_week_h)

## S-analysis
M_h <- survfit(Surv(entry, exit2, event2)~1, data=CH_h)
summary(M_h)

########################################################
## survival from Feb 1st - April 30th
Censor_week_e <- week(dmy("30-04-2012"))-Begin_week

CH_e <- CH %>% mutate(exit2=if_else(exit>Censor_week_e, Censor_week_e, exit), event2=if_else(exit>Censor_week_e, 0,event)) 

## S-analysis
M_e <- survfit(Surv(entry, exit2, event2)~1, data=CH_e)
summary(M_e)

########################################################
## survival from April13th - July 13th
Begin_week_l <- week(dmy("01-05-2012"))-Begin_week

CH_l <- CH %>% mutate(entry2=Begin_week_l) %>%
  filter(exit>Begin_week_l)

## S-analysis
M_l <- survfit(Surv(entry2, exit, event)~1, data=CH_l)
summary(M_l)

##################################################################################
##################################################################################
### MOVEMENT AND DISPLACEMENT

library(geosphere)


Cap_sites <- d %>% filter(event=="capture") %>% 
             select(organismName, decimalLatitude, decimalLongitude) %>%
             rename(capLat=decimalLatitude, capLong=decimalLongitude)

d <- full_join(d, Cap_sites)
d <- d %>% mutate(replaceDist=distGeo(matrix(c(d$decimalLongitude, d$decimalLatitude), ncol=2), 
                matrix(c(d$capLong, d$capLat), ncol=2))/1000) 

d_Displace <- d %>% filter(replaceDist<10) %>%
              mutate(Track_month=month(eventDate)) %>% 
              select(organismName, state, event, sex, lifeStage, Track_month, eventDate, replaceDist)              

d_max_each <- d %>% filter(state=="alive") %>%
              mutate(Track_month=month(eventDate)) %>% 
              mutate(New_period=cut(Track_month, breaks=c(0,3.9 , 4.1, 7, 12))) %>%
              group_by(organismName, New_period) %>%
              summarise(max_each=max(replaceDist)) %>%
              group_by(New_period) %>%
              summarise(Mean_rep=mean(max_each), sd_rep=sd(max_each), nn=n())    



temp <- d_Displace %>% group_by(organismName, state) %>%
        filter(Track_month<9) %>%
        summarise(last=max(Track_month)) %>%
        mutate(state=if_else(state=="alive", "alive", "dead")) %>%
        spread(key=state, value=last, fill=9)

########################
##### Plotting; 

temp1 <- d %>% group_by(organismName) %>%
         filter(event=="capture") %>%
         mutate(cap_year=year(eventDate)) %>%
         summarise(cap_year=min(cap_year)) 

d <- full_join(d, temp1)

d2 <- d %>% mutate(Track_week=week(eventDate), Track_year=year(eventDate)) %>%
      mutate(temp_year=Track_year-cap_year) %>%
      filter(temp_year==0 & Track_week<40) %>%
      mutate(state2=if_else(state=="alive", "alive", "dead")) %>%
      mutate(sex=recode(sex, "M"="Male", "F"= "Female"))


      p2 <- ggplot(data=d2, aes(x=Track_week, y=replaceDist, group=organismName)) +
      geom_line(linetype=2, col="#E7B800") + 
      geom_point(colour="dark red", shape=3, data=(filter(d2, state2=="dead")))

p2 + theme (legend.position = "top") + 
facet_wrap(~sex, labeller="label_both") +
xlab("Week of the year") +
ylab("Distance to capture location (km)") 

############################################
### FINDING MAX REPLACEMENT FOR EACH BIRD

Cap_age_sex <- d %>% filter(event=="capture") %>%
            select(organismName, sex, lifeStage)

Max_repl <- d %>% group_by(organismName) %>%
            summarise(max_repl=max(replaceDist)) %>%
            left_join(., Cap_age_sex)


######################################
#### Map of study area; 
library(ggmap)
library(maps)

study_site <- d %>% summarise(long=mean(decimalLongitude), lat=mean(decimalLatitude))

ggplot() + 
  xlim(-20, 50)+
  ylim(40, 80) +
  borders("world", fill="grey60") +
  borders(regions = "Norway(?!:Svalbard)", colour = "gray50", fill = "orange") +
  geom_point(data=study_site, aes(x=long, y=lat), col="dark red") + 

  theme_bw()

#########################################################################
### Adding locations and captures; 

library(ggmap)
library(mapproj)
map <- get_map(location = 'Norway', zoom = 4)
ggmap(map)

Study_area <- c(left = 12.5, bottom = 64, right = 15, top = 64.8)
get_stamenmap(Study_area, zoom = 10, maptype = "terrain") %>% 

ggmap() +
  geom_path(data=temp, aes(x=long, y=lat), size=1.2, type=2, col="dark red")+
  geom_point(data=d, aes(x=decimalLongitude, y=decimalLatitude, colour=state))+
  theme(legend.position = "top", legend.title = element_blank())

  
  borders("world", colour = "dark red", fill = NA) 
 

temp <- map_data("world") %>%
        filter(region=="Norway")
  
  
#####





































