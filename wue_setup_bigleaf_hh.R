library(tidyverse)
library(bigleaf)
library(lubridate)
`%notin%` <- Negate(`%in%`)

# Read in half hourly FLUXNET format data
# Filtering:

# Summer data: June, July, August only
# Middle of daytime only, time of day: 9am - 3pm (3:30 if hh?), ie 09:00 to 15:00

# Screen for valid GPP (NEE_VUT_REF_QC) and LE (LE_F_MDS_QC) data, 
# pick high threshold options and see how many samples we get (start at 0 only)
# NEE_VUT_REF_QC: 0 = measured; 1 = good quality gapfill; 2 = medium; 3 = poor

# and screen out for any days where there is rain or two days after
# P_F > 0, get date, then remove that day, that day +1, and that day +2

# Make sure all measurements are real QC?

#####
# testing on Harvard Forest - US-Ha1, and Niwot Ridge 'US-NR1'

#site_id <- 'US-Ha1'
site_id <- 'US-NR1'

# some sites will be hourly (e.g. US-Ha1), others will be half hourly
# note: read.csv() is slow for big files, will need to use faster function
#flx <- read.csv("/Users/davidmiller/dlm_files/keenan_postdoc/fluxnet2015/FLUXNET2015-latest/FLX_US-Ha1_FLUXNET2015_FULLSET_1991-2012_1-4/FLX_US-Ha1_FLUXNET2015_FULLSET_HR_1991-2012_1-4.csv")
flx <- read.csv("/Users/davidmiller/dlm_files/keenan_postdoc/fluxnet2015/FLUXNET2015-latest/FLX_US-NR1_FLUXNET2015_FULLSET_1998-2014_1-4/FLX_US-NR1_FLUXNET2015_FULLSET_HH_1998-2014_1-4.csv")

# TIMESTAMP_START is formatted as YYYYMMDDHHMM
flx$DateTimeStart <- parse_date_time(flx$TIMESTAMP_START, c("ymd HM"))

# get within a certain time of day, let's do 9-3 for "daytime", could change or expand later
# and filter for summertime months
# and filter for measured data only (= 0) from NEE (for GPP) and LE

# these are all ranges to KEEP
hour_range <- seq(9,15)
month_range <- seq(6,8)
nee_qc_range <- c(0)
le_qc_range <- c(0)
time_inds <- which(hour(flx$DateTimeStart) %in% hour_range & 
                     month(flx$DateTimeStart) %in% month_range &
                     flx$NEE_VUT_REF_QC %in% nee_qc_range &
                     flx$LE_F_MDS_QC %in% le_qc_range)

# this is a range to EXCLUDE
p_inds <- which(flx$P_F > 0)
rain_dates <- flx$DateTimeStart[p_inds] %>% date() %>% unique()
rain_dates_plus2 <- c(rain_dates, rain_dates + days(1), rain_dates + days(2)) %>% unique()
clear_good_date_inds <- which(date(flx$DateTimeStart) %notin% rain_dates_plus2) # then invert with %notin% to keep

# Intersect!
keep_inds <- intersect(time_inds, clear_good_date_inds)

# strict range, all measured P_F, T_F, VPD_F
met_range <- c(0)
met_good_inds <- which(flx$P_F_QC %in% met_range & flx$TA_F_MDS_QC %in% met_range & flx$VPD_F_MDS_QC %in% met_range)

keep_inds_strict_met <- intersect(keep_inds, met_good_inds)

flx_sub <- flx[keep_inds,]
flx_sub_strict <- flx[keep_inds_strict_met,] # note: this still screens out days(+2) that P_ERA thought it rained, even if it wasn't measured at the site

ggplot(flx_sub_strict, aes(DateTimeStart, GPP_NT_VUT_REF)) + geom_point()
# using flx_sub_strict since this is the "better, more real" data to use

# use bigleaf to convert le to et
flx_sub_strict$ET <- LE.to.ET(flx_sub_strict$LE_F_MDS, flx_sub_strict$TA_F)
# LE	Latent heat flux (W m-2)
# Tair	Air temperature (deg C)
# ET	Evapotranspiration (kg H2O m-2 s-1)

# GPP in hh data is umolCO2 m-2 s-1, convert to g C
#flx_sub_strict$GPP_NT_VUT_REF_gC <- umolCO2.to.gC(flx_sub_strict$GPP_NT_VUT_REF, constants = bigleaf.constants()) # this seems wrong!
flx_sub_strict$GPP_NT_VUT_REF_gC <- flx_sub_strict$GPP_NT_VUT_REF * 10e-6 * (12/44)

# this does the whole column at once...
# x <- WUE.metrics(
#   data = flx_sub_strict,
#   GPP = "GPP_NT_VUT_REF",
#   NEE = "NEE_VUT_REF",
#   LE = "LE_F_MDS",
#   VPD = "VPD_F_MDS",
#   Tair = "TA_F_MDS",
#   constants = bigleaf.constants()
# )

calc_WUE <- function(GPP, ET){
  # calculates WUE with GPP and ET in the same units
  WUE <- GPP/ET
  return(WUE)
}

calc_uWUE <- function(GPP, ET, VPD){
  # calculates underlying WUE with GPP as g C, ET as kg H2O, and VPD as kPa
  uWUE <- (GPP * sqrt(VPD)) / ET
  return(uWUE)
}

flx_sub_strict$VPD_F_MDS_kPa <- flx_sub_strict$VPD_F_MDS * 0.1

flx_sub_strict$WUE <- calc_WUE(flx_sub_strict$GPP_NT_VUT_REF_gC, flx_sub_strict$ET)
flx_sub_strict$uWUE <- calc_uWUE(flx_sub_strict$GPP_NT_VUT_REF_gC, flx_sub_strict$ET, flx_sub_strict$VPD_F_MDS_kPa)

#ggplot(flx_sub_strict, aes(DateTimeStart, WUE)) + geom_point()
#ggplot(flx_sub_strict, aes(DateTimeStart, uWUE)) + geom_point()

# shows variability in monthly uWUE across years
ggplot(flx_sub_strict) +
  geom_boxplot(aes(x = as.factor(month(DateTimeStart)), y = uWUE, color = as.factor(year(DateTimeStart)))) + 
  labs(x = "Month", color = "Year") +
  theme_classic()

flx_sub_strict$Date_01 <- ymd(paste(year(flx_sub_strict$DateTimeStart), "-", month(flx_sub_strict$DateTimeStart), "-01", sep = ""))

# aggregate for each month
WUE_mean <- aggregate(subset(flx_sub_strict, select = c(WUE, uWUE, GPP_NT_VUT_REF_gC, ET, TA_F_MDS, VPD_F_MDS_kPa)), 
                      by=list(YearMonth = flx_sub_strict$Date_01), FUN = mean)

WUE_sd <- aggregate(subset(flx_sub_strict, select = c(WUE, uWUE, GPP_NT_VUT_REF_gC, ET, TA_F_MDS, VPD_F_MDS_kPa)), 
                      by=list(YearMonth = flx_sub_strict$Date_01), FUN = sd)

WUE_length <- aggregate(subset(flx_sub_strict, select = c(WUE)), 
                    by=list(YearMonth = flx_sub_strict$Date_01), FUN = length)

colnames(WUE_sd)[2:length(colnames(WUE_sd))] <- paste(colnames(WUE_sd)[2:length(colnames(WUE_sd))], "_sd", sep = "")
colnames(WUE_length)[2] <- "n_samples"

WUE_summary <- merge(WUE_mean, WUE_sd, by = "YearMonth")
WUE_summary <- merge(WUE_summary, WUE_length, by = "YearMonth")
WUE_summary <- cbind.data.frame(rep(site_id, nrow(WUE_summary)), WUE_summary)
colnames(WUE_summary)[1] <- "SITE_ID"

ggplot(WUE_summary) +
  geom_point(aes(x = YearMonth, y = uWUE)) + 
  geom_errorbar(aes(x = YearMonth, ymin = uWUE - uWUE_sd, ymax = uWUE + uWUE_sd)) +
  labs(x = "Year + Month") +
  theme_classic()

monthly_list <- WUE_summary$YearMonth # unique list of all the months

# PDSI
site_pdsi <- read.csv("/Users/davidmiller/dlm_files/keenan_postdoc/terraclimate/processed_data/flx_site_combinedsource_terraclimate_PDSI.csv")
pdsi_months <- seq(ymd("19910101"), ymd("20211201"), by = "months") # this is for all columns in the PDSI data, note this is hard coded to match

pdsi_flx_months_cols <- which(pdsi_months %in% WUE_summary$YearMonth) + 7 # pad +7 for number of prep columns in PDSI dataframe

WUE_summary$PDSI <- NA
WUE_summary$PDSI_lag1 <- NA

for(i in 1:nrow(WUE_summary)){
  pdsi_col_index <- which(pdsi_months %in% WUE_summary$YearMonth[i]) + 7 # pad +7 for number of prep columns in PDSI dataframe
  WUE_summary$PDSI[i] <- as.numeric(site_pdsi[which(site_pdsi$SITE_ID == site_id), pdsi_col_index])
  WUE_summary$PDSI_lag1[i] <- as.numeric(site_pdsi[which(site_pdsi$SITE_ID == site_id), pdsi_col_index-1]) # lag by one month in PDSI
}

# OUTPUT THIS FILE
setwd("/Users/davidmiller/dlm_files/keenan_postdoc/wue_drought/outputs")
write.csv(WUE_summary, paste(site_id, "wue_summary.csv", sep = "_"), row.names = F)


# Will need to repeat the above for selection of flux sites
