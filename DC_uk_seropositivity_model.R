# Replication code: UK Seropositivity model of infections

# Note:
# IFR = Infection Fatality Ratio

# Load packages
library(readr)
library(ggplot2)
library(lubridate)
library(readxl)

# Load data from ons
temp <- tempfile(fileext = ".xlsx")
# Also saved in folder as "publishedweek512020corrected.xlsx" 
download.file("https://www.ons.gov.uk/file?uri=%2fpeoplepopulationandcommunity%2fbirthsdeathsandmarriages%2fdeaths%2fdatasets%2fweeklyprovisionalfiguresondeathsregisteredinenglandandwales%2f2020/publishedweek532020.xlsx", destfile=temp, mode='wb')

# Read in relevant sheet
ons <- data.frame(readxl::read_excel(temp, sheet = 7, skip = 4, col_types = "text"))

# Get date row and deaths by sex and total deaths rows
ons$sex <- NA
ons$sex[which(ons[, 2] == "Persons 4"):nrow(ons)] <- "all"
ons$sex[which(ons[, 2] == "Males 4"):nrow(ons)] <- "m"
ons$sex[which(ons[, 2] == "Females 4"):nrow(ons)] <- "f"
ons$sex[which(ons[, 2] == "Deaths by region of usual residence 5"):nrow(ons)] <- NA

colnames(ons)[2] <- "age"
ons <- ons[, c("sex", "age", paste0("X", 1:53))]
ons <- ons[!is.na(ons$sex), ]
ons <- ons[!ons$age %in% c("Deaths by age group",
                           "Persons 4",
                           "Males 4",
                           "Females 4"), ]

# From wide to long format
library(tidyr)
ons <- unique(ons) %>% gather(week, deaths, -c(sex, age))
ons$deaths <- as.numeric(ons$deaths)
ons$week <- as.numeric(unlist(gsub("X", "", ons$week)))

# Change to weeks since Jan 2020, leveraging that they come in order
ons$week <- ons$week + cumsum(ons$week < c(0, ons$week[-length(ons$week)]))*52

# Getting date
dates <- data.frame(readxl::read_excel(temp, sheet = 7, col_types = "text"))[5, ]
dates <- cbind.data.frame(date = as.Date(as.numeric(unlist(dates[, 3:c(ncol(dates)-1)])),  origin = "1899-12-30"), week = unique(ons$week))
ons <- merge(ons, dates, by = "week", all.x = T)

# Inspect to ensure they roughly match, and match other stats (e.g. gov.uk)
library(ggplot2)
ggplot(ons, aes(x=date, y=deaths, 
                group = interaction(sex, age)))+
  geom_area(data = ons[ons$sex == "all", ], aes(fill = "Total"))+
  geom_area(data = ons[ons$sex != "all", ], aes(fill = "Sum"), 
            alpha = 0.5)+
  theme_minimal()+xlab("")+geom_vline(aes(xintercept = unique(ons$date[ons$week == 52])[1]))

# We also check manually that they match:
sum(ons$deaths[ons$sex == "f"] + ons$deaths[ons$sex == "m"], na.rm = T) == sum(ons$deaths[ons$sex == "all"], na.rm = T)

# We consider men and women separately as they have different IFRs, and exclude the combined category because all deaths are categorized by sex
ons <- ons[ons$sex != "all", ]

# We next load estimated IFRs by gender and age from a study in Nature (link is available in csv):
  ifr <- read.csv("ifr_by_age_nature.csv")
  ifr$age <- ifr[, 1]
  ifr$ifr_raw <- ifr$ifr
  ifr$ifr <- as.numeric(unlist(lapply(ifr$ifr_raw, FUN = function(x) 
    unlist(strsplit(x, " "))[1])))
  ifr$ifr_low <- as.numeric(unlist(lapply(ifr$ifr_raw, FUN = function(x){
    x <- gsub("\\(", "-", x)
    x <- gsub("\\)", "", x)
    unlist(strsplit(x, "-"))[2]})))
  ifr$ifr_high <- as.numeric(unlist(lapply(ifr$ifr_raw, FUN = function(x){
    x <- gsub("\\(", "-", x)
    x <- gsub("\\)", "", x)
    unlist(strsplit(x, "-"))[3]})))
  ifr <- ifr[, c("age", "sex", "ifr_raw", "ifr", "ifr_low", "ifr_high", "source")]
  ifr <- ifr[ifr$age != "", ]
  ifr$age[ifr$age == "80"] <- "80+"


# We then combine deaths in the under 1 category and the 1-4 category because no source provides an IFR for them disaggregated, and same with 80+
ons$age[ons$age %in% c("<1", "1-4")] <- "0-4" 
ons$age[ons$age %in% c("80-84", "85-89", "90+")] <- "80+" 
ons$deaths <- ave(ons$deaths, interaction(ons$age, ons$week, ons$sex), FUN = function(x) sum(x))

ons <- unique(ons) # Avoid double-counting of categories
  
# Make sure categories match
unique(ons$age) == unique(ifr$age) 
unique(ons$sex) == unique(ifr$sex)

# Merge with test that all interactions matched
pre_nrow = nrow(ons)
ons <- merge(ons, ifr[, c("age", "sex", "ifr", 
                          "ifr_low", "ifr_high")], 
             by = c("age", "sex"))
pre_nrow == nrow(ons)

# Make IFR non-zero (using half of lowest observed ifr, relevant for deaths of the very young)
half_of_min <- function(x){min(na.omit(x[x > 0])/2)}
ons$ifr[ons$ifr == 0] <- half_of_min(ons$ifr)
ons$ifr_low[ons$ifr_low == 0] <- half_of_min(ons$ifr_low)
ons$ifr_high[ons$ifr_high == 0] <- half_of_min(ons$ifr_high)

# Construct estimates of IFR-implied infections
ons$sero_estimate <- ons$deaths*(100/ons$ifr)
ons$sero_estimate_low <- ons$deaths*(100/ons$ifr_high) # opposite necessary here
ons$sero_estimate_high <- ons$deaths*(100/ons$ifr_low) # opposite necessary here

# Get official uk cases:
# Load cases and deaths
temp <- tempfile(fileext = ".csv")
download.file("https://coronavirus.data.gov.uk/api/v1/data?filters=areaType=nation&structure=%7B%22areaType%22:%22areaType%22,%22areaName%22:%22areaName%22,%22areaCode%22:%22areaCode%22,%22date%22:%22date%22,%22newCasesByPublishDate%22:%22newCasesByPublishDate%22,%22cumCasesByPublishDate%22:%22cumCasesByPublishDate%22%7D&format=csv", destfile=temp, mode='wb')

# Also saved in folder as "publishedweek512020corrected.xlsx" 

# Read in relevant sheet
dat <- data.frame(read.csv(temp))

# Restrict to England and Wales and drop last 5 days
dat <- dat[dat$areaName %in% c("England", "Wales") & dat$date <= as.Date(Sys.Date()-5), ]
dat$new_cases <- ave(dat$newCasesByPublishDate, dat$date, FUN = function(x) sum(x, na.rm = T))
dat$confirmed <- ave(dat$cumCasesByPublishDate, dat$date, FUN = function(x) sum(x, na.rm = T))
dat$date <- as.Date(dat$date)
dat <- dat[!duplicated(dat$date), ]
for(i in 1:nrow(dat)){
dat$new_cases_7dma[i] <- mean(dat$new_cases[max(c(1,i-6)):i])
} 
dat$pop.2020 <- 59439840
ggplot(dat, aes(x=date, y=new_cases_7dma))+geom_line()

# Get date and week in order
dat$date <- as.Date(dat$date)
library(lubridate)
dat$week <- week(dat$date)
dat$week[year(dat$date) == "2021" & dat$week <= 53] <- dat$week[year(dat$date) == "2021"  & dat$week <= 53] + max(dat$week) # ensure weeks since jan 2020 rather than week of year

# Check to ensure correct merging:
max(table(dat$week)) == 7

# Define functions to get weekly sums of cases
weekly_sum <- function(x){
  x <- na.omit(x)
  if(length(x) == 7){
    return(sum(x))
  } else {
    print(x)
    return(7*sum(x)/length(x))
  }
}
end_of_week <- function(x){
  rev(na.omit(x))[1]*7
}

# We then make and rename a few variables for convenience
dat$new_cases_reported <- ave(dat$new_cases, dat$week, FUN = weekly_sum)
# dat$new_deaths_reported <- ave(dat$new_deaths, dat$week, FUN = weekly_sum)

dat$new_cases_reported_7dma <- ave(dat$new_cases_7dma, dat$week, FUN = end_of_week)
# dat$new_deaths_reported_7dma <- ave(dat$new_deaths_7dma, dat$week, FUN = end_of_week)

# dat$total_deaths <- dat$deaths
dat$total_cases <- dat$confirmed
dat <- rev(dat)

# Add empty rows to the ONS to facilitate interpolation of recent weeks
for(i in (max(ons$week)+1):max(df$week)){
  temp <- ons[ons$week == max(ons$week), ]
  temp[, c("deaths", "sero_estimate", "sero_estimate_low", "sero_estimate_high")] <- NA
  temp$week <- i
  ons <- rbind(ons, temp)
}

# Merging official and estimated cases
df <- merge(ons, dat[!duplicated(dat$week), c( "total_cases", "new_cases_reported",  "new_cases_reported_7dma", "pop.2020", "week")], by = "week", all.y= T)

# Get date for df:
for(i in 52:max(df$week)){
  df$date[df$week == i] <- unique(as.Date(df$date[df$week == i - 1] + 7))
}

# Shift to account for lag between infections and death:
# Time lag from cases to death:
# Source: https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
death_interval <- 14
pred <- ons
pred$date <- pred$date - death_interval
pred$week <- pred$week - round(death_interval/7, 0)
pred$pred_infections <- pred$sero_estimate
pred$pred_infections_low <- pred$sero_estimate_low
pred$pred_infections_high <- pred$sero_estimate_high

df <- merge(df, pred[, c("week", "age", "sex", 
                         "pred_infections", "pred_infections_low", 
                         "pred_infections_high")], all.x=T)

# Plot 1, no correction for missing data for recent weeks:
ggplot(df, aes(x=as.Date(date), y=pred_infections))+
  geom_col(aes(fill = age, col = age, x=as.Date(date)))+
  geom_line(aes(y=new_cases_reported_7dma))+
#  geom_line(aes(y=-100*new_deaths_reported_7dma, x = as.Date(date)))+
  theme_minimal()

# We next augment the model by extending data for recent week based on reported cases
last_week_of_good_data <- 52

for(t in c("pred_infections", "pred_infections_low", "pred_infections_high")){
  df$temp <- df[, t]
  
  df$pred_infections_total_in_week <- NA
  df$pred_infections_total_in_week <- with(df, ave(temp, week, FUN = sum))
 df$pred_infections_total_in_week <- ave(df$pred_infections_total_in_week, df$week, FUN = function(x) mean(x, na.rm = T))
 
 # Three alternatives here, included to inflate confidence intervals appropriately, as the ratio of cases to underlying is uncertain
 
 one <- df[df$week %in% c((last_week_of_good_data-1):last_week_of_good_data-round(death_interval/7, 0)), ]
 two <- df[df$week %in% c((last_week_of_good_data-2):last_week_of_good_data-round(death_interval/7, 0)), ]
 three <- df[df$week %in% c((last_week_of_good_data-3):last_week_of_good_data-round(death_interval/7, 0)), ]
 
 fun_multiplier <- function(x){
   sum(x$pred_infections_total_in_week, na.rm = T)/sum(x$new_cases_reported_7dma, na.rm = T)
 }
 
 temp <- df[df$week %in% c((last_week_of_good_data-3):last_week_of_good_data-round(death_interval/7, 0)), ]
  
  multiplier_with_good_data <- fun_multiplier(temp)
  
  if(t == "pred_infections_low"){
    multiplier_with_good_data <- min(c(fun_multiplier(one), 
                                       fun_multiplier(two),
                                       fun_multiplier(three)))
  }
  if(t == "pred_infections_high"){
    multiplier_with_good_data <- max(c(fun_multiplier(one), 
                                       fun_multiplier(two),
                                       fun_multiplier(three)))
  }
  
  # Last week of good data:
  for(j in (last_week_of_good_data-2):max(df$week)){
    total_pred_cases_target <- mean(multiplier_with_good_data*df$new_cases_reported_7dma[df$week == j])
    
    for(i in interaction(df$age, df$sex)){
      df$temp[as.character(interaction(df$age, df$sex)) == i 
                         & df$week == j] <- 
        mean(df$temp[df$week %in% (last_week_of_good_data-4):last_week_of_good_data &
                                  interaction(df$age, df$sex) == i], na.rm = T)
    }
    
    pred_cases <- sum(df$temp[df$week == j])
    multiplier <- total_pred_cases_target/pred_cases
    
    df$temp[df$week == j] <- df$temp[df$week == j]*multiplier
  }
  
  df[, t] <- df$temp
  }

# Recalculate totals post-extrapolation
df$pred_infections_total_in_week <- NA
df$pred_infections_total_in_week <- with(df, ave(pred_infections, week, FUN = sum))
df$pred_infections_total_in_week <- ave(df$pred_infections_total_in_week, df$week, FUN = function(x) mean(x, na.rm = T))

df$pred_infections_total_in_week_low <- NA
df$pred_infections_total_in_week_low <- with(df, ave(pred_infections_low, week, FUN = sum))
df$pred_infections_total_in_week_low <- ave(df$pred_infections_total_in_week_low, df$week, FUN = function(x) mean(x, na.rm = T))

df$pred_infections_total_in_week_high <- NA
df$pred_infections_total_in_week_high <- with(df, ave(pred_infections_high, week, FUN = sum))
df$pred_infections_total_in_week_high <- ave(df$pred_infections_total_in_week_high, df$week, FUN = function(x) mean(x, na.rm = T))

df <- df[df$week >= 10, ]

# Plot 2:
# Shows infections with age in color
ggplot(df, aes(x=as.Date(date), y=pred_infections))+
  geom_col(aes(fill = age, col = age, x=as.Date(date)))+
  geom_line(aes(y=new_cases_reported_7dma))+
  theme_minimal()+xlab("")+ylab("")#+geom_line(aes(y=-100*new_deaths_reported_7dma))

# Get cumulative infection rates:
df <- df[order(df$date), ]
df$pred_infections_total <- ave(df$pred_infections, interaction(df$age, df$sex), FUN= function(x) cumsum(x))

# Plot 3: (Panel A)
write.csv(df, "uk_covid_model.csv")
df <- read.csv("uk_covid_model.csv")
df$date <- as.Date(df$date)

ggplot(df[!duplicated(df$week), ], aes(x=as.Date(date)))+
  geom_area(aes(y=pred_infections_total_in_week, fill = "estimated"))+
  geom_line(aes(y=pred_infections_total_in_week_low), col = "black")+
  geom_line(aes(y=pred_infections_total_in_week_high), col = "black")+
  geom_area(aes(y=new_cases_reported_7dma, fill = "reported"))+
  theme_minimal()+xlab("")+ylab("")+scale_y_continuous(labels = scales::comma)+theme(legend.position = "bottom", legend.title = element_blank())+ggtitle("Britain's Pandemic\nWeekly totals, England and Wales\n(blue) Reported covid-19 cases, \n(red) Estimated Infections based on The Economist's seropositivity model \n(with 95% confidence interval)")
ggsave("panel_a.png", width = 10, height = 10)


# Plot 4: (Panel B)
ggplot(df, aes(x=date, fill=age, group=age))+
  geom_col(aes(y=100*pred_infections_total/(1000*pop.2020[1])))+ylab("%")+theme_minimal()+ggtitle("Estimated cumulative cases by age group, weekly totals")
ggsave("panel_b.png", width = 10, height = 10)


