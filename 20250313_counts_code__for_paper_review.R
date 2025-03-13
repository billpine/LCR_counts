############
## GLM Trends - Yearly Analysis of Trend Data
## Jennifer Moore
## for Publication
## January 2024
############

#import libraries 
library(waterData)
library(reshape)
library(ggplot2) #creates ggplots
library(cowplot) #format for ggplots
library(brms) #runs brm models
library(tidybayes)
library(performance)
library(modelr)
library(dplyr)

################
## Import Data and Source Oyster functions code
################

#import data
d <- read.csv("data/transect_data.csv", header = TRUE)
source("script/oyster_functions.R")

###################
## Question 1: Are overall trends in oyster counts over time similar across the
## entire inshore, nearshore, and offshore oyster reefs?
#####################

# Organize the full data set
d2 <- organizeData(d)
d3 <- calculateCountsDensity(d, d2)

# Make sure the counts are integers
d3$count_live <- as.integer(d3$count_live)

# Add harvest and strata columns
d3$fish <- substr(as.character(d3$strata), 1, 1)
d3$rockBin <- ifelse(d3$rocks %in% c("LG", "PI", "SM"), "Y", "N")
d3$strata4 <- paste(d3$harvest, "_", d3$rockBin)

# Subset for Lone Cabbage
dLC <- subset(d3, d3$locality == "LC")

dLC$timerock <- ifelse(dLC$rockBin == "Y" & dLC$period > 17, "afterY", 
                       ifelse(dLC$rockBin == "Y" & dLC$period < 18, "beforeY", dLC$rockBin))

dLC$period_fct <- as.factor(dLC$period)

# Convert period to factor for modeling
dLC$p_2 <- as.factor(dLC$period)
dLC$site_rock <- paste0(dLC$site, "_", dLC$rockBin)

dLC$period_res = paste0(dLC$period, "_", dLC$rockBin)


msite2 <- brm(count_live ~ 0 + site_rock + offset(log(tran_length)) + (1|period_res), 
              data = dLC, family = negbinomial, cores = 4)

summary(msite2)

# Create Figure 2

# Plot pre-period
datPre <- dLC[dLC$period < 10, ] %>% 
  select(period_res, site_rock) %>% 
  distinct() %>% 
  mutate(tran_length = 1)
datPre <- droplevels(datPre)

catPre <- datPre %>% add_epred_draws(msite2) %>%
  group_by(period_res, site_rock, .row) %>% 
  summarise(mu = mean(.epred),
            up = quantile(.epred, .975),
            lw = quantile(.epred, .025))

# Extract the number from period_res
catPre <- catPre %>%
  mutate(period_num = as.numeric(gsub("_.*", "", period_res)))

# Plot with the extracted number on the x-axis
prePlot <- ggplot(catPre, aes(x = period_num, y = mu, color = site_rock)) +
  geom_pointrange(aes(ymin = lw, ymax = up, fill = site_rock), size = 1.5,
                  position = position_dodge(width = 0.25), shape = 22, linewidth = 1)+
  ylim(0, 125) +
  labs(x = "Period", y = "Live oyster counts per m of transect", fill = "Site") +
  theme_bw(base_size = 20) +
  guides(color = "none") +
  scale_x_continuous(n.breaks = 7) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a")) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a"))

# Plot pilot-period
datPilot <- dLC[dLC$period %in% c(10, 11, 16), ] %>% 
  select(period_res, site_rock) %>% 
  distinct() %>% 
  mutate(tran_length = 1)
datPilot <- droplevels(datPilot)

catPilot <- datPilot %>% add_epred_draws(msite2) %>%
  group_by(period_res, site_rock, .row) %>% 
  summarise(mu = mean(.epred),
            up = quantile(.epred, .975),
            lw = quantile(.epred, .025))

# Extract the number from period_res
catPilot <- catPilot %>%
  mutate(period_num = as.numeric(gsub("_.*", "", period_res)))

pilotPlot <- ggplot(catPilot, aes(x = period_num, y = mu, color = site_rock)) +
  geom_pointrange(aes(ymin = lw, ymax = up, fill = site_rock), size = 1.5,
                  position = position_dodge(width = 0.25), shape = 22, linewidth = 1) +
  ylim(0, 125) +
  labs(x = "Period", y = "", fill = "Site") +
  theme_gray(base_size = 20) +
  guides(color = "none") +
  scale_x_continuous(n.breaks = 7) +
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c"))

# Plot post-period
datPost <- dLC[dLC$period > 17 & dLC$period != 19, ] %>% 
  select(period_res, site_rock) %>% 
  distinct() %>% 
  mutate(tran_length = 1)
datPost <- droplevels(datPost)

catPost <- datPost %>% add_epred_draws(msite2) %>%
  group_by(period_res, site_rock, .row) %>% 
  summarise(mu = mean(.epred),
            up = quantile(.epred, .95),
            lw = quantile(.epred, .05))

catPost <- catPost %>%
  mutate(period_num = as.numeric(gsub("_.*", "", period_res)))


postPlot <- ggplot(catPost, aes(x = period_num, y = mu, color = site_rock)) +
  geom_pointrange(aes(ymin = lw, ymax = up, fill = site_rock), size = 1.5,
                  position = position_dodge(width = 0.25), shape = 21, linewidth = 1) +
  ylim(0, 125) +
  labs(x = "Period", y = "", fill = "Site") +
  theme_gray(base_size = 20) +
  guides(color = "none") +
  scale_x_continuous(n.breaks = 9) +
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#33a02c")) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#33a02c"))

# Combine plots
f2<-plot_grid(prePlot, pilotPlot, postPlot, nrow = 1)

##

# Save the plot as a JPEG file
ggsave("D:\\f2.jpeg", plot = f2, width = 12, height = 6, units = "in", dpi = 300)

# Combine dataframes
combined_cat <- bind_rows(catPre, catPilot, catPost)

# Add 'site' and 'rock' columns
# Combine dataframes and retain a single column heading
combined_cat <- bind_rows(catPre, catPilot, catPost)

# Add 'site', 'rock' columns and calculate m2 columns
combined_cat <- combined_cat %>%
  mutate(
    site_rock = as.character(site_rock),
    site = case_when(
      grepl("^I", site_rock) ~ "Inshore",
      grepl("^N", site_rock) ~ "Nearshore",
      grepl("^O", site_rock) ~ "Offshore",
      TRUE ~ NA_character_
    ),
    rock = case_when(
      grepl("N$", site_rock) ~ "No",
      grepl("Y$", site_rock) ~ "Yes",
      TRUE ~ NA_character_
    ),
    # Add m2 columns with proper rounding
    mu_m2 = round(mu / 0.1524, 2),
    lw_m2 = round(lw / 0.1524, 2),
    up_m2 = round(up / 0.1524, 2)
  ) %>%
  # Round all numeric columns to two decimal places
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  # Select only specified columns in the desired order
  select(site, rock, period_num, mu, lw, up, mu_m2, lw_m2, up_m2)

# Write selected columns to a CSV file
write.csv(combined_cat, "D:\\table1.csv", row.names = FALSE)


##################
### Question 2: At the repeated measures sites, how did counts of oysters 
### change over time after restoration?
##################

#start with the original dataset (d) so that only transects 1-3 are kept for the
#4 repeated measures sites

#period 17 is 0 oysters
#data from period 17-26 only 

#subset for 4 stations
d$station <- as.factor(d$station) #make station a factor
t3d <- subset(d, d$station == "LCO9C" | d$station == "LCO10A" | d$station == "LCO11B" | d$station == "LCO12")
t3d$station <- droplevels(t3d$station)
#subset just for transects 1-3
t3d2 <- subset(t3d, t3d$transect < 4)

t3d3 <- organizeData(t3d2)
t3d4 <- calculateCountsDensity(t3d2,t3d3)

t3d4 <- subset(t3d4, t3d4$period > 16)

#check data is correct
table(t3d4$station, t3d4$period)

##      18 20 22 24 26
#LCO10A 1  1  1  1  1
#LCO11B 1  1  1  1  1
#LCO12  1  1  1  1  1
#LCO9C  1  1  1  1  1

#make sure counts are an integer
t3d4$count_live <- as.integer(t3d4$count_live)

# explore period as random effect here
t3d4$period_fct = as.factor(t3d4$period)

mm_fit_1 <- brm(count_live ~ (1 | period_fct) + offset(log(tran_length)), data = t3d4, family = negbinomial(), chains = 4, 
                cores = 4, control = list(adapt_delta = 0.99))

summary(mm_fit_1)

new_dat = select(t3d4, period_fct, treatment) %>% 
  distinct() %>% 
  mutate(tran_length = 1)

cat = new_dat %>% add_epred_draws(mm_fit_1) %>% 
  group_by(period_fct, treatment, .row) %>% 
  summarise(mu = mean(.epred),
            up = quantile(.epred, .975),
            lw = quantile(.epred, .025))  # note: only 90% interval here

#Figure 3
f3<-ggplot(cat, aes(x = as.numeric(as.character(period_fct)), y = mu)) +
  geom_pointrange(aes(ymin = lw, ymax = up), color = "red", size = 1.5, 
                  position = position_dodge(width = .25),
                  fill = "red", shape = 21, linewidth = 1) +
  geom_jitter(data = t3d4, aes(y = count_live / tran_length, size = 1.5), color = "black", 
              width = .1) +
  labs(y = "Live oyster counts per m of transect", x = "Period")+
  theme_bw(base_size = 20) +
  guides(color = "none", size = "none") +
  scale_x_continuous(n.breaks = 9)

# Save the plot as a JPEG file
ggsave("D:\\f3.jpeg", plot = f3, width = 10, height = 6, units = "in", dpi = 300)


##MAKE Q2 TABLE
# Extract coefficients for all models except coef_mc3_13 (base) and mc3_14b2
coef_q2 <- extract_coefficients(mm_fit_1, "mm_fit_1")

# Exponentiate and prepare data for plotting
library(dplyr)

all_q2 <- coef_q2 %>%
  mutate(exp_mean = exp(mean),
         exp_lower = exp(lower_95),
         exp_upper = exp(upper_95),
         mean_m2 = (exp_mean / 0.1524),
         lower_m2 = (exp_lower / 0.1524),
         upper_m2 = (exp_upper / 0.1524)) 
# Round numeric columns to 2 decimal places
q2_rounded <- all_q2 %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save to CSV
write.csv(q2_rounded, "d:\\q2_rounded.csv", row.names = FALSE)


##########################
### Question 3: Across Suwannee Sound (Lone Cabbage reef), are temporal trends in
### live oyster counts similar between restored and unrestored oyster reefs open
### and closed to harvest?
###########################

#organize the full data set
d2 <- organizeData(d)
d3 <- calculateCountsDensity(d,d2)

#make sure the counts are integers
d3$count_live <- as.integer(d3$count_live)

#add in harvest column
#fishing (yes/no) based on strata column
d3$fish <- substr(as.character(d3$strata), 1, 1)
d3$rockBin <- "Y"
d3$rockBin[d3$rocks == "NA"] <- "N"

#create strata column for just 4 strata
d3$strata4 <- paste(d3$harvest, "_", d3$rockBin)

#subset just for Lone Cabbage
dLC <- subset(d3, d3$locality == "LC")

dLC$timerock <- dLC$rockBin
dLC$timerock[dLC$rockBin == "Y" & dLC$period > 17] <- "afterY"
dLC$timerock[dLC$rockBin == "Y" & dLC$period < 18] <- "beforeY"

dLC$period_fct = as.factor(dLC$period)
# # 

library(stringr)
library(tidyverse)
library(ggthemes)

dLC$site_rock <- paste0(dLC$site, "_", dLC$rockBin)
table(dLC$site_rock)

# new RE levels with only those periods & restoration combos observed
dLC$period_res = paste0(dLC$period, "_", dLC$rockBin)
table(dLC$period_res)

m_test <- brm(count_live ~ 0 + site_rock + offset(log(tran_length)) + 
                (1|period_res), data = dLC, family = negbinomial, cores = 4)

summary(m_test)

variables(m_test)  # what are the brms parameter names

tat = gather_draws(m_test, b_site_rockI_N, b_site_rockN_N, b_site_rockO_N, b_site_rockO_Y) %>% 
  summarise(mu = mean(exp(.value)),
            up = quantile(exp(.value), .975),
            lw = quantile(exp(.value), .025)) %>% 
  mutate(label = substr(.variable, 12, 14), 
         site = case_when(label == "I_N" ~ "Inshore",
                          label == "N_N" ~ "Nearshore",
                          label == "O_N" ~ "Offshore",
                          label == "O_Y" ~ "Offshore"),
         restored = case_when(str_detect(label, "N") ~ "No",
                              str_detect(label, "Y") ~ "Yes"))

#Figure 4
p<-ggplot(tat, aes(x = label, y = mu)) + 
  labs(y = "Live oyster counts per m of transect", x = "", color = "Restored", fill = "Restored") +
  geom_jitter(data = dLC, aes(x = site_rock, y = count_live/tran_length, size = 1), color = "black",
              width = .1, alpha = .3) +
  geom_pointrange(aes(ymax = up, ymin = lw, color = restored, fill = restored), size = 1.5, linewidth = 1, shape = 21) +
  scale_x_discrete(labels = tat$site) +
  scale_color_solarized(accent = "blue") +
  scale_fill_solarized(accent = "blue") +
  theme_bw(base_size = 15) +
  guides(size = "none") +
  theme(legend.position = c(.7,.85))

# Save the plot as a JPEG file
ggsave("D:\\oyster_counts.jpeg", plot = p, width = 10, height = 6, units = "in", dpi = 300)

# look at the random effects
bat = gather_draws(m_test, r_period_res[period,]) %>%
  summarise(mu = mean(.value),
            up = quantile(.value, .975),
            lw = quantile(.value, .025)) %>%
  mutate(restored_lab = if_else(str_detect(period, "N"), "No", "Yes"),
         period_lab = parse_number(period))

q <- ggplot(bat, aes(x = period_lab, y = mu, ymax = up, ymin = lw)) +
  geom_pointrange(aes(color = restored_lab), position = position_dodge(width = 0.5)) +
  labs(y = "Random effects estimates", x = "Period", color = "Restored", fill = "Restored") +
  scale_x_continuous(breaks = 1:26) +
  scale_color_solarized(accent = "blue") +
  scale_fill_solarized(accent = "blue") +
  theme_bw(base_size = 20) +
  theme(legend.position = c(0.9, 0.175), legend.text = element_text(size = 12), legend.title = element_text(size = 12))

# Save the plot as a JPEG file
ggsave("D:\\random_effects.jpeg", plot = q, width = 10, height = 6, units = "in", dpi = 300)


###
####
# Extract coefficients for all models except coef_mc3_13 (base) and mc3_14b2
coef_q3 <- extract_coefficients(m_test, "m_test")

# Exponentiate and prepare data for plotting
library(dplyr)

all_q3 <- coef_q3 %>%
  mutate(exp_mean = exp(mean),
         exp_lower = exp(lower_95),
         exp_upper = exp(upper_95),
         mean_m2 = (exp_mean / 0.1524),
         lower_m2 = (exp_lower / 0.1524),
         upper_m2 = (exp_upper / 0.1524)) 
# Round numeric columns to 2 decimal places
q3_rounded <- all_q3 %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save to CSV
write.csv(q3_rounded, "d:\\q3_rounded.csv", row.names = FALSE)


################################
## Question 4: How do covariates of management interest explain variation in 
## oyster counts over time between restored and unrestored reefs?
#################################

########
## organize covariate information

#pull in covariates
#river discharge data - suwannee
#station to analyze
station = '02323500'   
#get site name to use in plot titles and such
stinfo  = siteInfo(station)
#read entire time series
dis = importDVs(staid=station,code='00060',stat='00003', sdate= "1941-01-01") 
dis$year = as.numeric(strftime(dis$dates,format="%Y"))
dis$month = as.numeric(strftime(dis$dates, format="%m"))
#Naming columns, using the Diver sensors, collects date, pressure, temp, conductivity
colnames(dis) <- c("StaID", "Discharge", "Date", "QualCode", "Year", "Month")

#add period to the data frame based on the month/year
dis$Period <- NA
firstyear <- 2008
endyear <- 2022
years <- sort(rep(firstyear:endyear, times = 1, each = 2))
for(i in unique(years)){
  y <- i #year
  p <- which(years == i) #period number - 2010 = 1 and 2, 2011 = 3 and 4, and so forth.
  for(j in 1:nrow(dis)){
    if(dis$Year[j] == y & dis$Month[j] > 3 & dis$Month[j] < 10) dis$Period[j] = p[1] #year i months 4-9
    if(dis$Year[j] == y & dis$Month[j] > 9) dis$Period[j] = p[2] #year i months 10-12
    if(dis$Year[j] == y+1 & dis$Month[j] < 4) dis$Period[j] = p[2] #year i+1 months 1-3
  }
}

#remove data w/o a period (before sampling started)
dis2 <- subset(dis, !is.na(dis$Period))

#organize water quality data
dis_mean_year <- aggregate(Discharge ~ Period, dis, mean)
dis_mean_year$Period <- dis_mean_year$Period - 4
dis_mean_year$YearLag1 <- dis_mean_year$Period + 2
dis_mean_year$YearLag2 <- dis_mean_year$Period + 4
dis_total_year <- aggregate(Discharge ~ Period, dis, sum)
dis_total_year$Period <- dis_total_year$Period - 4
dis_total_year$YearLag1 <- dis_total_year$Period + 2
dis_total_year$YearLag2 <- dis_total_year$Period + 4

#add covariate data to original data set
d4 <- merge(d3, dis_mean_year[,c("Period", "Discharge")], 
            by.x = "period", by.y = "Period", all.x = TRUE)
names(d4)[length(names(d4))] <- "AD"

d4 <- merge(d4, dis_mean_year[,c("YearLag1", "Discharge")], 
            by.x = "period", by.y = "YearLag1", all.x = TRUE)
names(d4)[length(names(d4))] <- "ADLag1"

d4 <- merge(d4, dis_mean_year[,c("YearLag2", "Discharge")], 
            by.x = "period", by.y = "YearLag2", all.x = TRUE)
names(d4)[length(names(d4))] <- "ADLag2"

d4 <- merge(d4, dis_total_year[,c("Period", "Discharge")], 
            by.x = "period", by.y = "Period", all.x = TRUE)
names(d4)[length(names(d4))] <- "TD"

d4 <- merge(d4, dis_total_year[,c("YearLag1", "Discharge")], 
            by.x = "period", by.y = "YearLag1", all.x = TRUE)
names(d4)[length(names(d4))] <- "TDLag1"

d4 <- merge(d4, dis_total_year[,c("YearLag2", "Discharge")], 
            by.x = "period", by.y = "YearLag2", all.x = TRUE)
names(d4)[length(names(d4))] <- "TDLag2"

#
# Calculate the correlation matrix for columns 20 to 25
cor_matrix <- cor(d4[, 20:25], use = "complete.obs")

# Plot the correlation matrix
library(corrplot)
corrplot(cor_matrix, method = "color", addCoef.col = "black", tl.cex = 0.7)


d4$ADsc <- scale(d4$AD)
d4$ADLag1sc <- scale(d4$ADLag1)
d4$ADLag2sc <- scale(d4$ADLag2)
d4$TDsc <- scale(d4$TD)
d4$TDLag1sc <- scale(d4$TDLag1)
d4$TDLag2sc <- scale(d4$TDLag2)

#######
#oyster landing data
ol <- read.csv("data/landings_data_updated2024.csv", header = T)

#add period to the data frame based on the month/year
ol$Period <- NA
firstyear <- 2008
endyear <- 2022
years <- sort(rep(firstyear:endyear, times = 1, each = 2))
for(i in unique(years)){
  y <- i #year
  p <- which(years == i) #period number - 2010 = 1 and 2, 2011 = 3 and 4, and so forth.
  for(j in 1:nrow(ol)){
    if(ol$Year[j] == y & ol$Month[j] > 3 & ol$Month[j] < 10) ol$Period[j] = p[1] #year i months 4-9
    if(ol$Year[j] == y & ol$Month[j] > 9) ol$Period[j] = p[2] #year i months 10-12
    if(ol$Year[j] == y+1 & ol$Month[j] < 4) ol$Period[j] = p[2] #year i+1 months 1-3
  }
}

#remove data w/o a period (before/after sampling)
ol2 <- subset(ol, !is.na(ol$Period))

#organize landings and trips data
landing<- aggregate(as.numeric(Pounds..meat.weight.) ~ Period, ol2, sum)
colnames(landing) <- c("Period", "Value")
landing$Period <- landing$Period - 4 #line the periods back up with period 1 = 2010
landing$YearLag1 <- landing$Period + 2
landing$YearLag2 <- landing$Period + 4
trips <- aggregate(Trips ~ Period, ol2, sum)
colnames(trips) <- c("Period", "Value")
trips$Period <- trips$Period - 4
trips$YearLag1 <- trips$Period + 2
trips$YearLag2 <- trips$Period + 4

#add to dataset
d5 <- merge(d4, landing[,c("Period", "Value")], 
            by.x = "period", by.y = "Period", all.x = TRUE)
names(d5)[length(names(d5))] <- "Landings"
d5 <- merge(d5, landing[,c("YearLag1", "Value")], 
            by.x = "period", by.y = "YearLag1", all.x = TRUE)
names(d5)[length(names(d5))] <- "LandingsLag1"
d5 <- merge(d5, landing[,c("YearLag2", "Value")], 
            by.x = "period", by.y = "YearLag2", all.x = TRUE)
names(d5)[length(names(d5))] <- "LandingsLag2"

d5 <- merge(d5, trips[,c("Period", "Value")], 
            by.x = "period", by.y = "Period", all.x = TRUE)
names(d5)[length(names(d5))] <- "Trips"
d5 <- merge(d5, trips[,c("YearLag1", "Value")], 
            by.x = "period", by.y = "YearLag1", all.x = TRUE)
names(d5)[length(names(d5))] <- "TripsLag1"
d5 <- merge(d5, trips[,c("YearLag2", "Value")], 
            by.x = "period", by.y = "YearLag2", all.x = TRUE)
names(d5)[length(names(d5))] <- "TripsLag2"

d5$Landingssc <- scale(d5$Landings)
d5$LandingsLag1sc <- scale(d5$LandingsLag1)
d5$LandingsLag2sc <- scale(d5$LandingsLag2)
d5$Tripssc <- scale(d5$Trips)
d5$TripsLag1sc <- scale(d5$TripsLag1)
d5$TripsLag2sc <- scale(d5$TripsLag2)


#subset for lone cabbage
dLC <- subset(d5, d5$locality == "LC")

##########################################
#GLM w/ covariates
#look at trend in covariate and then run model

# fit model set with period as a random effect 
dLC$p_2 = as.factor(dLC$period) # convert to factor

dLC$period_res = paste0(dLC$period, "_", dLC$rockBin) #rock bin is only rocks = yes

dLC$site_rock <- paste0(dLC$site, "_", dLC$rockBin) #rock bin is only rocks = yes

#####model set

#river discharge year of sampling (yearly mean)
mc3_1b2 <- brm(count_live ~  0+site_rock * ADsc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
##river discharge 1 year lag (yearly mean)
mc3_2b2 <- brm(count_live ~  0+site_rock * ADLag1sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#river discharge 2 year lag (yearly mean)
mc3_3b2 <- brm(count_live ~  0+site_rock * ADLag2sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#river discharge year of sampling (yearly total)
mc3_4b2 <- brm(count_live ~  0+site_rock * TDsc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
##river discharge 1 year lag (yearly total)
mc3_5b2 <- brm(count_live ~  0+site_rock * TDLag1sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#river discharge 2 year lag (yearly total)
mc3_6b2 <- brm(count_live ~  0+site_rock * TDLag2sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#oyster fishing trips year of sampling
mc3_7b2 <- brm(count_live ~  0+site_rock * Tripssc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#oyster fishing trips 1 year lag
mc3_8b2 <- brm(count_live ~  0+site_rock * TripsLag1sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#oyster fishing trips 2 year lag
mc3_9b2 <- brm(count_live ~  0+site_rock * TripsLag2sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#oyster meat weight (landings) year of sampling
mc3_10b2 <- brm(count_live ~  0+site_rock * Landingssc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
#oyster meat weight 1 year lag
mc3_11b2 <- brm(count_live ~  0+site_rock * LandingsLag1sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
##oyster meat weight 2 year lag
mc3_12b2 <- brm(count_live ~  0+site_rock * LandingsLag2sc + offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 
##base model
mc3_13b2 <- brm(count_live ~  0+site_rock * offset(log(tran_length)) + (1|period_res), data = dLC, family = negbinomial, cores = 4) 


##making sure we are comparing offshore_restored that is
##fished and not to the other sites, which are not restored

library(tidyr)
##open/closed to harvest
# Summarize observations from dLC
table_fish_site_rock <- dLC %>%
  group_by(site_rock, fish) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = fish, values_from = count, values_fill = 0)

#this is inshore not restored, then no_fishing, yes_fishing
#then nearshore not restored, then no_fishing, yes_fishing
#then offshore not restored, then no_fishing, yes_fishing
#then offshore restored, then no_fishing, yes_fishing

#now I want to combine all of the not restored, keeping no_fishing, yes_fishing
#and compare to restored, keeping no_fishing, yes_fishing

filtered_dLC <- dLC %>%
  mutate(
    site_rock = ifelse(site_rock == "O_Y", "O_Y", "all_N")
  )

names(filtered_dLC)

table_fish_site_rock2 <- filtered_dLC %>%
  group_by(site_rock, fish) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = fish, values_from = count, values_fill = 0)

#coding in this way gets around the holes, below is collapsed version
#see how it is calling filtered_dLC

filtered_dLC$rock_fish <- paste0(filtered_dLC$site_rock, "_", filtered_dLC$fish) #rock bin is only rocks = yes

#this is uncollapsed (uncollapsed meaning the site effects are still there)
#see how it is calling dLC

dLC$rock_fish <- paste0(dLC$site_rock, "_", dLC$fish) #rock bin is only rocks = yes


#so the first row is offshore, yes to restoration and the number
#of transects no fishing and yes fishing
#then we have all other locations, not restored, no fishing, yes fishing


# Fit the model with unfiltered  - has sites
#uses rock_fish and not site_rock

#

mc3_14b2.2 <- brm(
  count_live ~ 0 + rock_fish + offset(log(tran_length)) + (1 | period_res),
  data = dLC,
  family = negbinomial,
  cores = 4
)
summary(mc3_14b2.2)

#ok this mode gives you I, N, O and then Y, N restored/not restored, and
#then yes/no fishing.  But it gives them to you as individual
#factors. So the parameter estimate is the effect
#there isn't a slope b/c it isn't continuous

#so you can compare rock_fishI_N_N to rock_fishI_N_Y for 
#example and that is comparing inshore, not restored areas
#that are closed or open to fishing

#write.csv(dLC, "D:/dLC.csv", row.names = FALSE)

summary(mc3_12b2)


# Print model summary
print(summary(mc3_14b2.2))


variables(mc3_14b2.2) 



#####################################
##now the LOO model comparison stuff
#####################################


#loo comparisons below. looking to see if some models
#ELPO <2

loo_1b2 = loo(mc3_1b2)
loo_2b2 = loo(mc3_2b2)
loo_3b2 = loo(mc3_3b2)
loo_4b2 = loo(mc3_4b2)
loo_5b2 = loo(mc3_5b2)
loo_6b2 = loo(mc3_6b2)
loo_7b2 = loo(mc3_7b2)
loo_8b2 = loo(mc3_8b2)
loo_9b2 = loo(mc3_9b2)
loo_10b2 = loo(mc3_10b2)
loo_11b2 = loo(mc3_11b2)
loo_12b2 = loo(mc3_12b2)
loo_13b2 = loo(mc3_13b2)
loo_14b2.2 = loo(mc3_14b2.2)

loo_tblb2 = loo_compare(loo_1b2, loo_2b2, loo_3b2, loo_4b2, loo_5b2, loo_6b2, loo_7b2,
                        loo_8b2, loo_9b2, loo_10b2, loo_11b2, loo_12b2, loo_13b2, loo_14b2.2)

# Convert the loo_compare() result to a data frame
loo_tblb2_df <- as.data.frame(loo_tblb2)
loo_tblb2_df$Model <- rownames(loo_tblb2_df) # Add model names as a column
rownames(loo_tblb2_df) <- NULL               # Remove row names


#need to fix this again so references model names in words

mod_list2 = list("mc3_1" = mc3_1b2,
                 "mc3_2" = mc3_2b2,
                 "mc3_3" = mc3_3b2,
                 "mc3_4" = mc3_4b2,
                 "mc3_5" = mc3_5b2,
                 "mc3_6" = mc3_6b2,
                 "mc3_7" = mc3_7b2,
                 "mc3_8" = mc3_8b2,
                 "mc3_9" = mc3_9b2,
                 "mc3_10" = mc3_10b2, 
                 "mc3_11" = mc3_11b2,
                 "mc3_12" = mc3_12b2,
                 "mc3_13" = mc3_13b2,
                 "mc3_14" = mc3_14b2.2)

modnames = c("AD", "ADLag1", "ADLag2", "TD", "TDLag1", "TDLag2",
             "Trips", "TripsLag1", "TripsLag2", "Landings", "LandingsLag1", 
             "LandingsLag2", "Base", "Fish")

loo_tbl_22 = as.data.frame(loo_tblb2)
loo_tbl_22$rank = 1:14
loo_tbl_22$model = rownames(loo_tbl_22)


elpo_plot<-ggplot(loo_tblb2_df, aes(x = reorder(Model, -elpd_diff), y = elpd_diff)) +
  geom_pointrange(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff), color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Model Comparison (LOO)",
    x = "Model",
    y = "ELPD Difference from Best"
  )


# Save the plot as a JPEG file
ggsave("D:\\elpo.jpeg", plot = elpo_plot, width = 10, height = 6, units = "in", dpi = 300)


##posterior predictive checks

library(bayesplot)

# PPC for top 3 models
ppc_12 <- pp_check(mc3_12b2, ndraws = 100) + ggtitle("PPC: mc3_12b2")
ppc_9 <- pp_check(mc3_9b2, ndraws = 100) + ggtitle("PPC: mc3_9b2")
ppc_11 <- pp_check(mc3_11b2, ndraws = 100) + ggtitle("PPC: mc3_11b2")

# Display all PPC plots together
library(patchwork)
ppc_12 + ppc_9 + ppc_11

#squiggles

###


# Compute multi level R-squared for each model
r2_12 <- bayes_R2(mc3_12b2)
#> r2_12
#Estimate  Est.Error      Q2.5     Q97.5
#R2 0.786005 0.03584226 0.7019954 0.8384891

r2_9  <- bayes_R2(mc3_9b2)
#> r2_9
#Estimate  Est.Error      Q2.5     Q97.5
#R2 0.7868211 0.03648014 0.6991481 0.8401687

r2_11 <- bayes_R2(mc3_11b2)
#> r2_11
#Estimate  Est.Error     Q2.5    Q97.5
#R2 0.7861224 0.03650104 0.698571 0.839404

##multi level R2 suggests models perform equally well

##########


##below start extracting predictors

#top 3 mc3_12b2, mc3_9b2, mc3_11b2
#special mc3_14b2.2

variables(mc3_12b2) #want 1:8
variables(mc3_9b2) #want 1:8
variables(mc3_11b2) #want 1:8
#variables(mc3_14b2.2) #want 1:6 (no interaction)

###
# Collect coefficients individually for all models, then mc3_14b2 separately

extract_coefficients <- function(model, model_name, positions = c(1:8)) {
  as_draws_df(model)[, positions] %>%
    summarise(
      variable_name = get_variables(model)[positions],
      mean = colMeans(., na.rm = TRUE),
      lower_95 = apply(., 2, quantile, 0.025, na.rm = TRUE),
      upper_95 = apply(., 2, quantile, 0.975, na.rm = TRUE),
      model_name = model_name
    )
}

# Extract coefficients for all models except coef_mc3_13 (base) and mc3_14b2
coef_mc3_9 <- extract_coefficients(mc3_9b2, "mc3_9")
coef_mc3_11 <- extract_coefficients(mc3_11b2, "mc3_11")
coef_mc3_12 <- extract_coefficients(mc3_12b2, "mc3_12")

# Combine results, adding mc3_14b2 last and no mc3_13 (base)
all_coef_outputs <- bind_rows(
  coef_mc3_9, coef_mc3_11, coef_mc3_12
)

# Print results
print(all_coef_outputs)

# Round numeric columns to 2 decimal places
all_coef_outputs_rounded <- all_coef_outputs %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save to CSV
write.csv(all_coef_outputs_rounded, "d:\\all_coef_outputs.csv", row.names = FALSE)



# Exponentiate and prepare data for plotting
all_effects_top3 <- all_coef_outputs %>%
  mutate(exp_mean = exp(mean),
         exp_lower = exp(lower_95),
         exp_upper = exp(upper_95),
         mean_m2 = (exp_mean / 0.1524),
         lower_m2 = (exp_lower / 0.1524),
         upper_m2 = (exp_upper / 0.1524))

#####
# Round numeric columns to 2 decimal places
all_effects_rounded_top3 <- all_effects_top3 %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save to CSV
write.csv(all_effects_rounded_top3, "d:\\all_effects_top3.csv", row.names = FALSE)



# Plotting the effect of all covariates
ggplot(all_effects, aes(x = variable_name, y = mean, ymin = lower_95, ymax = upper_95, color = model_name)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Top 3 LOO Models - All Covariates with 95% Credible Intervals",
       x = "Covariate",
       y = "Effecy") +
  facet_wrap(~ model_name, scales = "free_y") +
  coord_flip()

# Filter results where 95% CI does not include zero
significant_results <- all_coef_outputs %>%
  filter(lower_95 > 0 | upper_95 < 0)

# Display the filtered table
print(significant_results)
# 

#######################
##contrasts top LOO models
# Extract correct parameters for O_Y - O_N contrasts
#######################

# mc3_9b2: Trips 2-Year Lag
int_9 <- as_draws_df(mc3_9b2)[, c("b_site_rockO_Y:TripsLag2sc", "b_site_rockO_N:TripsLag2sc")] %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  mutate(model = "mc3_9b2 (TripsLag2sc)")

# mc3_11b2: Landings 1-Year Lag
int_11 <- as_draws_df(mc3_11b2)[, c("b_site_rockO_Y:LandingsLag1sc", "b_site_rockO_N:LandingsLag1sc")] %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  mutate(model = "mc3_11b2 (LandingsLag1sc)")

# mc3_12b2: Landings 2-Year Lag
int_12 <- as_draws_df(mc3_12b2)[, c("b_site_rockO_Y:LandingsLag2sc", "b_site_rockO_N:LandingsLag2sc")] %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  mutate(model = "mc3_12b2 (LandingsLag2sc)")

# Combine all extracted results
int_all <- bind_rows(int_9, int_11, int_12)

# Pivot Wider for Comparison
int_wide <- int_all %>%
  pivot_wider(
    names_from = variable,
    values_from = value,
    values_fn = list
  ) %>%
  unnest(cols = c(`b_site_rockO_Y:TripsLag2sc`, `b_site_rockO_N:TripsLag2sc`,
                  `b_site_rockO_Y:LandingsLag1sc`, `b_site_rockO_N:LandingsLag1sc`,
                  `b_site_rockO_Y:LandingsLag2sc`, `b_site_rockO_N:LandingsLag2sc`))

# Compute Contrasts: O_Y - O_N (Rock Effect)
int_contrasts <- int_wide %>%
  mutate(
    contrast_LandingsLag1 = `b_site_rockO_Y:LandingsLag1sc` - `b_site_rockO_N:LandingsLag1sc`,
    contrast_LandingsLag2 = `b_site_rockO_Y:LandingsLag2sc` - `b_site_rockO_N:LandingsLag2sc`,
    contrast_TripsLag2    = `b_site_rockO_Y:TripsLag2sc` - `b_site_rockO_N:TripsLag2sc`
  )

# Summarize and Rank Predictors

int_contrast_summary <- int_contrasts %>%
  summarise(
    mean_contrast_Lag1 = mean(contrast_LandingsLag1, na.rm = TRUE),
    lower_95_Lag1 = quantile(contrast_LandingsLag1, 0.025, na.rm = TRUE),
    upper_95_Lag1 = quantile(contrast_LandingsLag1, 0.975, na.rm = TRUE),
    
    mean_contrast_Lag2 = mean(contrast_LandingsLag2, na.rm = TRUE),
    lower_95_Lag2 = quantile(contrast_LandingsLag2, 0.025, na.rm = TRUE),
    upper_95_Lag2 = quantile(contrast_LandingsLag2, 0.975, na.rm = TRUE),
    
    mean_contrast_Trips = mean(contrast_TripsLag2, na.rm = TRUE),
    lower_95_Trips = quantile(contrast_TripsLag2, 0.025, na.rm = TRUE),
    upper_95_Trips = quantile(contrast_TripsLag2, 0.975, na.rm = TRUE)
  )

# Print summary of contrasts with 95% credible intervals
print(int_contrast_summary)

# Rank Predictors by Mean Contrast
ranked_contrasts <- int_contrast_summary %>%
  pivot_longer(cols = everything(), names_to = "contrast_type", values_to = "value") %>%
  filter(grepl("mean_contrast", contrast_type)) %>%
  arrange(desc(value)) %>%
  mutate(rank = row_number())

# Display Ranked Predictors
print(ranked_contrasts)


# Visualize O_Y - O_N Contrasts 
#think about how to interpret this

int_contrasts_long <- int_contrasts %>%
  pivot_longer(cols = starts_with("contrast_"), names_to = "predictor", values_to = "contrast_value")

ggplot(int_contrasts_long, aes(x = predictor, y = contrast_value, fill = predictor)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = "Comparison of O_Y - O_N Contrasts Across Predictors",
    subtitle = "Rock Effect by Predictor (Trips, LandingsLag1, LandingsLag2)",
    x = "Predictor",
    y = "O_Y - O_N Contrast (with Posterior Draws)",
    fill = "Predictor"
  ) +
  theme(legend.position = "top")


##now predict the contrasts as number of oysters

##convert contrasts to number of oysters

# Define a Baseline Oyster Count 
baseline_oysters <- 4.4  # from table 1, unrestored per m, 16

# Convert Contrasts to Oyster Counts 

# Compute Oyster Increase from O_Y - O_N Contrast
int_contrasts_oysters <- int_contrasts %>%
  mutate(
    oysters_Lag1 = baseline_oysters * exp(contrast_LandingsLag1),
    oysters_Lag2 = baseline_oysters * exp(contrast_LandingsLag2),
    oysters_Trips = baseline_oysters * exp(contrast_TripsLag2)
  )

# Summarize 
int_contrast_summary_oysters <- int_contrasts_oysters %>%
  summarise(
    mean_oysters_Lag1 = mean(oysters_Lag1, na.rm = TRUE),
    lower_95_Lag1 = quantile(oysters_Lag1, 0.025, na.rm = TRUE),
    upper_95_Lag1 = quantile(oysters_Lag1, 0.975, na.rm = TRUE),
    
    mean_oysters_Lag2 = mean(oysters_Lag2, na.rm = TRUE),
    lower_95_Lag2 = quantile(oysters_Lag2, 0.025, na.rm = TRUE),
    upper_95_Lag2 = quantile(oysters_Lag2, 0.975, na.rm = TRUE),
    
    mean_oysters_Trips = mean(oysters_Trips, na.rm = TRUE),
    lower_95_Trips = quantile(oysters_Trips, 0.025, na.rm = TRUE),
    upper_95_Trips = quantile(oysters_Trips, 0.975, na.rm = TRUE)
  )

# Results 
print(int_contrast_summary_oysters)

write.csv(int_contrast_summary_oysters, "D:/int_contrast_summary_oysters.csv", row.names = FALSE)

#This means that having rock (O_Y) increases oyster counts approximately 3 times (from 100 to 293 oysters) compared to no rock (O_N) 
#under the impact of 1-year lagged harvest.




########
##ROCKFISH, this is the "simplified" fishing model
#below is where i extract and work with the posteriors
########

summary(mc3_14b2.2)
variables(mc3_14b2.2)

table_fish_site_rock3 <- dLC %>%
  group_by(rock_fish) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = rock_fish, values_from = count, values_fill = 0)



# Extract Posterior Samples for Restoration and Fishing Effects
rockfish_posterior <- as_draws_df(mc3_14b2.2)[, c(
  "b_rock_fishI_N_N", "b_rock_fishI_N_Y",
  "b_rock_fishN_N_Y", "b_rock_fishO_N_N",
  "b_rock_fishO_Y_N", "b_rock_fishO_Y_Y"
)]

# Compute Contrasts of rockfish model
rockfish_posterior <- rockfish_posterior %>%
  mutate(
    inshore_fishing = b_rock_fishI_N_Y - b_rock_fishI_N_N,  # Inshore fishing effect, no restoration fishing yes - no restoration fishing no
    offshore_fishing = b_rock_fishO_Y_Y - b_rock_fishO_Y_N,  # Offshore fishing effect, no restoration fishing yes - no restoration fishing no
    offshore_restoration = b_rock_fishO_Y_Y - b_rock_fishO_N_N  # Offshore restoration effect, restoration fishing no - no restoration fishing no
  )

####
extract_coefficients <- function(model, model_name, positions = c(5,6,7,8)) {
  as_draws_df(model)[, positions] %>%
    summarise(
      variable_name = get_variables(model)[positions],
      mean = colMeans(., na.rm = TRUE),
      lower_95 = apply(., 2, quantile, 0.025, na.rm = TRUE),
      upper_95 = apply(., 2, quantile, 0.975, na.rm = TRUE),
      model_name = model_name
    )
}
####



# Extract coefficient
coef_rockfish <- extract_coefficients(mc3_14b2.2, "mc3_14b2.2")

# Exponentiate and prepare data for plotting
library(dplyr)

all_rockfish <- coef_rockfish %>%
  mutate(exp_mean = exp(mean),
         exp_lower = exp(lower_95),
         exp_upper = exp(upper_95),
         mean_m2 = (exp_mean / 0.1524),
         lower_m2 = (exp_lower / 0.1524),
         upper_m2 = (exp_upper / 0.1524)) 
# Round numeric columns to 2 decimal places
all_rockfish_rounded <- all_rockfish %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Save to CSV
write.csv(all_rockfish_rounded, "d:\\all_rockfish.csv", row.names = FALSE)



###


# Summarize Results
summary_contrasts <- rockfish_posterior %>%
  summarise(
    mean_inshore_fishing = mean(inshore_fishing),
    lower_95_inshore_fishing = quantile(inshore_fishing, 0.025),
    upper_95_inshore_fishing = quantile(inshore_fishing, 0.975),
    
    mean_offshore_fishing = mean(offshore_fishing),
    lower_95_offshore_fishing = quantile(offshore_fishing, 0.025),
    upper_95_offshore_fishing = quantile(offshore_fishing, 0.975),
    
    mean_offshore_restoration = mean(offshore_restoration),
    lower_95_offshore_restoration = quantile(offshore_restoration, 0.025),
    upper_95_offshore_restoration = quantile(offshore_restoration, 0.975)
  )

# Print Summary
print(summary_contrasts)

# Store effects
comparison_effects <- rockfish_posterior %>%
  select(inshore_fishing, offshore_fishing, offshore_restoration)

# Reshape for Plotting
comparison_long <- comparison_effects %>%
  pivot_longer(cols = everything(), names_to = "Effect", values_to = "Value")

# Plot Comparison - Density Plot
ggplot(comparison_long, aes(x = Value, fill = Effect)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "Comparison of Restoration vs Fishing Effects",
    x = "Effect Size (Log Scale)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Boxplot to Compare Effect Sizes
ggplot(comparison_long, aes(x = Effect, y = Value, fill = Effect)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "Effect Size Comparison",
    x = "Effect",
    y = "Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

####

##mike's half eye plot for paper 

# theme better suited to pubs
theme_pub <- function(){
  
  theme_minimal() %+replace%
    
    theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"),
          panel.spacing.y = unit(1, "lines"),
          strip.text = element_text(size = 14, hjust = 0.02, vjust = 1),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          title = element_text(size = 14),
          # axis.title = element_text(size = 14),
          axis.title.x = element_text(size = 14, vjust = -1),
          # axis.title.y = element_text(size = 14, angle = 90, vjust = 1.5)
          axis.title.y = element_text(size = 14, angle = 90, vjust = 1.5),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(size = 12, colour = "black"))
}



# key that relates rock_fish to other vars (unpacks rock_fish)
key = select(dLC, rock_fish, site, rockBin, harvest) %>% 
  distinct() %>% 
  mutate(site_lab = case_when(site == "I" ~ "Inshore",  # add some better plot labels. 
                              site == "N" ~ "Nearshore",
                              site == "O" ~ "Offshore"),
         harvest_lab = if_else(harvest == "N", "No", "Yes"),
         res_lab = if_else(rockBin == "N", "No", "Yes"))

variables(mc3_14b2.2)

# show_col(solarized_pal()(2))  # these are the colors I picked, change as you like 
# cols = solarized_pal()(2)

dat = gather_draws(mc3_14b2.2, b_rock_fishI_N_N, b_rock_fishI_N_Y, b_rock_fishN_N_Y,
                   b_rock_fishO_N_N, b_rock_fishO_Y_N, b_rock_fishO_Y_Y) %>% 
  mutate(rock_fish = substr(.variable, 12, 16)) %>% 
  left_join(key, by = "rock_fish")

ggplot(dat, aes(x = harvest_lab, y = exp(.value))) +
  stat_halfeye(aes(color = res_lab), position = position_dodge(width = .5),
               slab_fill = "grey", slab_alpha = .5, shape = 21,
               point_color = "black") + 
  scale_color_solarized(accent = "blue") +
  labs(y = "Oysters / m", x = "Harvest", color = "Restoration") +
  coord_cartesian(ylim = c(0, 60)) +
  geom_segment(x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black") +
  geom_segment(x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black") +
  facet_wrap(~site_lab, scales = "free_x") +
  theme_pub() +  # see below for this function 
  theme(legend.position = c(.925,.875),
        panel.grid.major.y = element_line(color = "grey70", size = 0.5),
        panel.grid.minor.y = element_line(color = "grey85", size = 0.3)) +
  scale_y_continuous(breaks = seq(0, 60, by = 10), minor_breaks = seq(0, 60, by = 5)) +
  guides(color = guide_legend(override.aes = list(size = 8)))


ggsave("d:\\halfeye_plot.jpeg", width = 8, height = 6, dpi = 300)


####
##predictions against baseline
###


# Define Baseline Oyster Count in Unrestored, Fished Areas
baseline_oysters <- 4.4  # per m transect (period 16, offshore unrestored)

# Extract Posterior Samples for Offshore Restoration Effects
post_samples <- as_draws_df(mc3_14b2.2)[, c("b_rock_fishO_N_N", "b_rock_fishO_Y_N",
                                            "b_rock_fishO_Y_Y","b_rock_fishI_N_Y",
                                            "b_rock_fishI_N_N", "b_rock_fishN_N_Y")]

# Compute Predicted Oyster Counts
post_samples <- post_samples %>%
  mutate(
    oysters_unrestored = baseline_oysters * exp(b_rock_fishO_N_N - b_rock_fishO_N_N),  # Baseline (30)
    oysters_restoredunfished = baseline_oysters * exp(b_rock_fishO_Y_N - b_rock_fishO_N_N),  # Effect of Restoration unfished reefs
    oysters_restoredfished = baseline_oysters * exp(b_rock_fishO_Y_Y - b_rock_fishO_N_N),  # Effect of Restoration+fished reefs
    oysters_off_fished = baseline_oysters * exp(b_rock_fishO_Y_Y - b_rock_fishO_Y_N), #effect of fishing offshore restored
    oysters_in_fished = baseline_oysters * exp(b_rock_fishI_N_Y - b_rock_fishI_N_N) # effect of fishing inshore
  )

colnames(post_samples)

# Summarize Results

# Summarize Results
oyster_summary <- post_samples %>%
  summarise(
    mean_unrestored = mean(oysters_unrestored),
    lower_95_unrestored = quantile(oysters_unrestored, 0.025),
    upper_95_unrestored = quantile(oysters_unrestored, 0.975),
    
    mean_restoredunfished = mean(oysters_restoredunfished),
    lower_95_restoredunfished = quantile(oysters_restoredunfished, 0.025),
    upper_95_restoredunfished = quantile(oysters_restoredunfished, 0.975),
    
    
    mean_restoredfished = mean(oysters_restoredfished),
    lower_95_restoredfished = quantile(oysters_restoredfished, 0.025),
    upper_95_restoredfished = quantile(oysters_restoredfished, 0.975),
    
    mean_fishedoff = mean(oysters_off_fished),
    lower_95_fishedoff = quantile(oysters_off_fished, 0.025),
    upper_95_fishedoff = quantile(oysters_off_fished, 0.975),
    
    mean_fishedin = mean(oysters_in_fished),
    lower_95_fishedin = quantile(oysters_in_fished, 0.025),
    upper_95_fishedin = quantile(oysters_in_fished, 0.975)
  )

# Print Summary
print(oyster_summary)

# Convert oyster_summary to long format correctly
oyster_summary_long <- oyster_summary %>%
  pivot_longer(
    cols = everything(),
    names_to = c("stat", "condition"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  pivot_wider(names_from = stat, values_from = value)

# Print to verify structure
print(oyster_summary_long)

# Create the plot
# Filter out 'unrestored' condition before plotting
oyster_summary_long_filtered <- oyster_summary_long %>%
  filter(condition != "unrestored")

# Create a named vector of custom labels
custom_labels <- c(
  "restoredunfished" = "Rock, Offshore, No fishing",
  "restoredfished" = "Rock, Offshore, Fishing",
  "fishedoff" = "No rock, Offshore, Fishing",
  "fishedin" = "No rock, Inshore, Fishing"
)

# Create the plot
oyster_plot <- ggplot(oyster_summary_long_filtered, aes(x = condition, y = mean, fill = condition)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2) +
  geom_hline(yintercept = 4.4, linetype = "dashed", color = "black") +  # Baseline reference line (black)
  theme_minimal() +
  labs(
    title = "Predicted Oyster Counts Across Condition",
    x = "Condition",
    y = "Estimated Oyster Count per m Transect"
  ) +
  scale_x_discrete(labels = custom_labels) +  # Apply custom x-axis labels
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")

# Print the plot
print(oyster_plot)

# Save the specific plot
#
ggsave("D:/oyster_predicted.png", plot = oyster_plot, width = 8, height = 6, units = "in", dpi = 300)

##



##ranked predictors from top 3 models
##this is included in appendix


# Summarize and Rank Predictors 

int_contrast_summary <- int_contrasts %>%
  summarise(
    mean_contrast_Lag1 = mean(contrast_LandingsLag1, na.rm = TRUE),
    lower_95_Lag1 = quantile(contrast_LandingsLag1, 0.025, na.rm = TRUE),
    upper_95_Lag1 = quantile(contrast_LandingsLag1, 0.975, na.rm = TRUE),
    
    mean_contrast_Lag2 = mean(contrast_LandingsLag2, na.rm = TRUE),
    lower_95_Lag2 = quantile(contrast_LandingsLag2, 0.025, na.rm = TRUE),
    upper_95_Lag2 = quantile(contrast_LandingsLag2, 0.975, na.rm = TRUE),
    
    mean_contrast_Trips = mean(contrast_TripsLag2, na.rm = TRUE),
    lower_95_Trips = quantile(contrast_TripsLag2, 0.025, na.rm = TRUE),
    upper_95_Trips = quantile(contrast_TripsLag2, 0.975, na.rm = TRUE)
  )

# Rank Predictors by Mean Contrast
ranked_contrasts <- int_contrast_summary %>%
  pivot_longer(cols = everything(), names_to = "contrast_type", values_to = "value") %>%
  filter(grepl("mean_contrast", contrast_type)) %>%
  arrange(desc(value)) %>%
  mutate(rank = row_number())

# Display Ranked Predictors
print(ranked_contrasts)


##
# Extract Posterior Samples from top 3 models
library(tidyverse)
library(posterior)

# Ensure the variables exist - print variable names if needed
variables(mc3_9b2)
variables(mc3_11b2)
variables(mc3_12b2)



###


# Visualize O_Y - O_N Contrasts 
#think about how to interpret this

int_contrasts_long <- int_contrasts %>%
  pivot_longer(cols = starts_with("contrast_"), names_to = "predictor", values_to = "contrast_value")

ggplot(int_contrasts_long, aes(x = predictor, y = contrast_value, fill = predictor)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = "Comparison of O_Y - O_N Contrasts Across Predictors",
    subtitle = "Rock Effect by Predictor (Trips, LandingsLag1, LandingsLag2)",
    x = "Predictor",
    y = "O_Y - O_N Contrast (with Posterior Draws)",
    fill = "Predictor"
  ) +
  theme(legend.position = "top")


##now predict the contrasts as number of oysters

##convert contrasts to number of oysters

# Define a Baseline Oyster Count 
baseline_oysters <- 4.4  # from table 1, unrestored per m2, last period

# Convert Contrasts to Oyster Counts 

# Compute Oyster Increase from O_Y - O_N Contrast
int_contrasts_oysters <- int_contrasts %>%
  mutate(
    oysters_Lag1 = baseline_oysters * exp(contrast_LandingsLag1),
    oysters_Lag2 = baseline_oysters * exp(contrast_LandingsLag2),
    oysters_Trips = baseline_oysters * exp(contrast_TripsLag2)
  )

# Summarize 
int_contrast_summary_oysters <- int_contrasts_oysters %>%
  summarise(
    mean_oysters_Lag1 = mean(oysters_Lag1, na.rm = TRUE),
    lower_95_Lag1 = quantile(oysters_Lag1, 0.025, na.rm = TRUE),
    upper_95_Lag1 = quantile(oysters_Lag1, 0.975, na.rm = TRUE),
    
    mean_oysters_Lag2 = mean(oysters_Lag2, na.rm = TRUE),
    lower_95_Lag2 = quantile(oysters_Lag2, 0.025, na.rm = TRUE),
    upper_95_Lag2 = quantile(oysters_Lag2, 0.975, na.rm = TRUE),
    
    mean_oysters_Trips = mean(oysters_Trips, na.rm = TRUE),
    lower_95_Trips = quantile(oysters_Trips, 0.025, na.rm = TRUE),
    upper_95_Trips = quantile(oysters_Trips, 0.975, na.rm = TRUE)
  )

# Results 
print(int_contrast_summary_oysters)

write.csv(int_contrast_summary_oysters, "D:/int_contrast_summary_oysters.csv", row.names = FALSE)

#This means that having rock (O_Y) increases oyster counts approximately 3 times (from 100 to 293 oysters) compared to no rock (O_N) 
#under the impact of 1-year lagged harvest.


##the below may not be needed
##because it goes through and extracts
## the slopes from all of the models
##instead of just focusing on the top 3

##it may be
##a good idea to include all of the results
##then you can see in the plot the ones that
##are different compared across models

##then that is also supported by the LOO results

#then you do the contrasts for the top 3


################################################################
# Updated code to extract slopes using provided variable names
# using tidybayes

# Collect coefficients individually for all models, then mc3_14b2 separately

extract_coefficients <- function(model, model_name, positions = c(5,6,7,8)) {
  as_draws_df(model)[, positions] %>%
    summarise(
      variable_name = get_variables(model)[positions],
      mean = colMeans(., na.rm = TRUE),
      lower_95 = apply(., 2, quantile, 0.025, na.rm = TRUE),
      upper_95 = apply(., 2, quantile, 0.975, na.rm = TRUE),
      model_name = model_name
    )
}

# Extract coefficients for all models except coef_mc3_13 (base) and mc3_14b2
coef_mc3_1 <- extract_coefficients(mc3_1b2, "mc3_1")
coef_mc3_2 <- extract_coefficients(mc3_2b2, "mc3_2")
coef_mc3_3 <- extract_coefficients(mc3_3b2, "mc3_3")
coef_mc3_4 <- extract_coefficients(mc3_4b2, "mc3_4")
coef_mc3_5 <- extract_coefficients(mc3_5b2, "mc3_5")
coef_mc3_6 <- extract_coefficients(mc3_6b2, "mc3_6")
coef_mc3_7 <- extract_coefficients(mc3_7b2, "mc3_7")
coef_mc3_8 <- extract_coefficients(mc3_8b2, "mc3_8")
coef_mc3_9 <- extract_coefficients(mc3_9b2, "mc3_9")
coef_mc3_10 <- extract_coefficients(mc3_10b2, "mc3_10")
coef_mc3_11 <- extract_coefficients(mc3_11b2, "mc3_11")
coef_mc3_12 <- extract_coefficients(mc3_12b2, "mc3_12")
#coef_mc3_13 <- extract_coefficients(mc3_13b2, "mc3_13")
#13 is base model, no covariates

# Special extraction for mc3_14b2 (only positions 3 and 4)
coef_mc3_14_special <- as_draws_df(mc3_14b2)[, c(3,4)] %>%
  summarise(
    variable_name = get_variables(mc3_14b2)[c(3,4)],
    mean = colMeans(., na.rm = TRUE),
    lower_95 = apply(., 2, quantile, 0.025, na.rm = TRUE),
    upper_95 = apply(., 2, quantile, 0.975, na.rm = TRUE),
    model_name = "mc3_14_special"
  )

# Combine results, adding mc3_14b2 last and no mc3_13 (base)
all_coef_outputs <- bind_rows(
  coef_mc3_1, coef_mc3_2, coef_mc3_3, coef_mc3_4,
  coef_mc3_5, coef_mc3_6, coef_mc3_7, coef_mc3_8,
  coef_mc3_9, coef_mc3_10, coef_mc3_11, coef_mc3_12,
  coef_mc3_14_special
)

# Print results
print(all_coef_outputs)

# Add new exponentiated columns to the data frame
all_coef_outputs <- all_coef_outputs %>%
  mutate(
    exp_mean = exp(mean),
    exp_lower_95 = exp(lower_95),
    exp_upper_95 = exp(upper_95)
  )

# Display the updated table
print(all_coef_outputs)

# Exponentiate and prepare data for plotting
all_effects <- all_coef_outputs %>%
  mutate(exp_mean = exp(mean),
         exp_lower = exp(lower_95),
         exp_upper = exp(upper_95))

# Plotting the effect of all covariates
ggplot(all_effects, aes(x = variable_name, y = mean, ymin = lower_95, ymax = upper_95, color = model_name)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "All Covariates with 95% Credible Intervals",
       x = "Covariate",
       y = "Effecy") +
  facet_wrap(~ model_name, scales = "free_y") +
  coord_flip()


# Filter results where 95% CI does not include zero
significant_results <- all_coef_outputs %>%
  filter(lower_95 > 0 | upper_95 < 0)

# Display the filtered table
print(significant_results)

# Plot coefficients with 95% credible intervals
ggplot(all_coef_outputs, aes(x = variable_name, y = mean, color = model_name)) +
  stat_halfeye(aes(ymin = lower_95, ymax = upper_95), alpha = 0.6) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.2) +
  theme_minimal() +
  coord_flip() +
  labs(title = "Model Coefficients with 95% Credible Intervals",
       x = "Coefficient",
       y = "Estimate")



####
###Appendix figures

################################
#### heatmap for discharge

#get monthly sum, mean, sd, and var
#discharge
dis.period = subset(dis, dis$Period < 27) #only want data for period 1-26
dis.period.agg  = aggregate(Discharge~Period,data=dis.period,FUN = function(x) c(mean(x,na.rm=T),sd(x,na.rm=T),var(x,na.rm=T),sum(x)))
dis.period.agg  = do.call('data.frame',dis.period.agg)
names(dis.period.agg)[2:5] = c('avg','sd','var','sumflow') 

#calculate how far each month average is from the overall average (average over period of record) - absolute difference
avgPOR <- mean(dis$Discharge, na.rm = T) #9671.727

dis.period.agg$diff <- dis.period.agg$avg - avgPOR
dis.period.agg$diff_percent <- ((dis.period.agg$avg - avgPOR)/(avgPOR))*100

#Figure A2
par(mfrow=c(1,1))
plot(dis.period.agg$Period, dis.period.agg$diff_percent, xlab = "Period/Year", 
     ylab = "Percent difference from 1941-2024 average", pch = 16, 
     ylim = c(-100, 100), cex = 2, xaxt = 'n', yaxt = 'n', cex.lab = 1.2)
axis(side = 1, at = 1:26, labels = 1:26)
axis(side = 1, at = seq(2.5,24.5,2), labels = seq(2011, 2022), line = 1, tick = FALSE)
axis(side = 2, las = 2)
abline(h=0, col = 'red', lwd = 2)



######standard discharge heat map by month and year
#get monthly sum, mean, sd, and var
#discharge
disE.mo  = aggregate(Discharge~Month+Year,data=dis,FUN = function(x) c(mean(x,na.rm=T),sd(x,na.rm=T),var(x,na.rm=T),sum(x)))
disE.mo  = do.call('data.frame',disE.mo)
names(disE.mo)[3:6] = c('avg','sd','var','sumflow') 
disE.mo$yrmo = disE.mo$Year+(disE.mo$Month-0.5)/12       

#value just by month
disE.month  = aggregate(Discharge~Month,data=dis,FUN = function(x) c(mean(x,na.rm=T),sd(x,na.rm=T),var(x,na.rm=T),sum(x)))
disE.month  = do.call('data.frame',disE.month)
names(disE.month)[2:5] = c('avg','sd','var','sumflow')

#calculate how far each month average is from the overall average - absolute difference
for(i in 1:nrow(disE.mo)){
  m = disE.mo$Month[i]
  disE.mo$diff[i] = disE.mo$avg[i] - disE.month$avg[disE.month$Month == m]
  disE.mo$diff_percent[i] = ((disE.mo$avg[i] - disE.month$avg[disE.month$Month == m])/disE.month$avg[disE.month$Month == m])*100
}

library(reshape2)
dis_mat = dcast(disE.mo, Year ~ Month, value.var = "diff")
dis_mat2 = as.matrix(dis_mat[,2:ncol(dis_mat)], dimnames = list(rownames = dis_mat$Year, colnames = colnames(dis_mat)))
rownames(dis_mat2) <- dis_mat$Year
#create matrix for percentage
dis_per = dcast(disE.mo, Year ~ Month, value.var = "diff_percent")
dis_per2 = as.matrix(dis_per[,2:ncol(dis_per)], dimnames = list(rownames = dis_mat$Year, colnames = colnames(dis_mat)))
rownames(dis_per2) <- dis_per$Year
library(RColorBrewer)
col2 <- brewer.pal(9, "Blues") #for the +%s
col3 <- brewer.pal(11, "RdBu") #red is - blue is +
colAll <- c(col3[1], col3[2], col3[3], col3[6], col2[4:7], col3[10],col2[8:9], col3[11], "#042333ff") 
require("plot.matrix")
include_years <- as.character(as.numeric(2010:2023))
dis_per20102023 <- dis_per2[include_years, ]
colnames(dis_per20102023) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
                               "Oct", "Nov", "Dec")
#Figure A3
par(oma=c(0.5,0.5,0.5,3))
plot(dis_per20102023, xlab = "", ylab = "", main = "",
     breaks = c(-100, -50, -25, -10, 10, 25, 50, 100, 150, 200, 250, 300, 350, 400),
     axis.row = list(las = 2),
     col = colAll, na.col = "black", xaxt = 'n')


