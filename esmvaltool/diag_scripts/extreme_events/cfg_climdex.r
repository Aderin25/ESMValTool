# #############################################################################
# climdex_df.r
#
# Authors: Christian W. Mohr (CICERO, Norway)
#
#########################################################################################################
# Description
#    This script provides a data frame of all the possible indices (core and user-defined),
#    which can be used with the climdex.pcic.ncdf R-package. The data frame created serves as an index
#    to find which indices have been processed, and what title and labels the plots should have.
#
# Modification history
#    2018 0725-A_cwmohr: Created
#########################################################################################################

idx_df <- data.frame(
  idx_etccdi = c(
    "altcdd", "altcsdi", "altcwd", "altwsdi", "cdd", "csdi", "cwd",
    "dtr", "dtr", "fd", "gsl", "id", "prcptot", "r10mm", "r1mm",
    "r20mm", "r95p", "r99p", "rx1day", "rx1day", "rx5day", "rx5day",
    "sdii", "su", "tn10p", "tn10p", "tn90p", "tn90p", "tnn", "tnn",
    "tnx", "tnx", "tr", "tx10p", "tx10p", "tx90p", "tx90p", "txn",
    "txn", "txx", "txx", "wsdi"
  ),
  time = c(
    "yr", "yr", "yr", "yr", "yr", "yr", "yr", "mon", "yr", "yr",
    "yr", "yr", "yr", "yr", "yr", "yr", "yr", "yr", "mon", "yr",
    "mon", "yr", "yr", "yr", "mon", "yr", "mon", "yr", "mon", "yr",
    "mon", "yr", "yr", "mon", "yr", "mon", "yr", "mon", "yr", "mon",
    "yr", "yr"
  ),
  unit = c(
    "days", "days", "days", "days", "days", "days", "days",
    "deg C", "deg C", "days", "days", "days",
    "mm", "days", "days", "days", "mm", "mm", "mm", "mm", "mm", "mm", "mm/day",
    "days", "Exceedance rate, %", "Exceedance rate, %",
    "Exceedance rate, %", "Exceedance rate, %",
    "deg C", "deg C", "deg C", "deg C", "days",
    "Exceedance rate, %", "Exceedance rate, %",
    "Exceedance rate, %", "Exceedance rate, %",
    "deg C", "deg C", "deg C", "deg C", "days"
  ),
  name = c(
    "Consecutive Dry Days per Year (altCDD)",
    "Cold Spell Duration Index Spanning Years (altCSDI)",
    "Consecutive Wet Days per Year (altCWD)",
    "Warm Spell Duration Index Spanning Years (altWSDI)",
    "Consecutive Dry Days (CDD)", "Cold Spell Duration Index (CSDI)",
    "Consecutive Wet Days (CWD)",
    "Monthly Diurnal Temperature Range (DTR)",
    "Annual Diurnal Temperature Range (DTR)", "Frost Days (FD)",
    "Growing Season Length (GSL)", "Icing Days (ID)",
    "Annual Total Wet-Day Precipitation (PRCPTOT)",
    "Heavy Precipitation Days (R10)", "Precipitation Days (R1)",
    "Very Heavy Precipitation Days (R20)",
    "Very Wet Days (R95p)", "Extremely Wet Days (R99p)",
    "Monthly Max 1-day Precipitation (RX1day)",
    "Annual Max 1-day Precipitation (RX1day)",
    "Monthly Max 5-day Precipitation (RX5day)",
    "Annual Max 5-day Precipitation (RX5day)",
    "Simple Daily Intensity Index (SDII)", "Summer Days (SD)",
    "Monthly Cold Nights (TN10p)", "Annual Cold Nights (TN10p)",
    "Monthly Warm Nights (TN90p)", "Annual Warm Nights (TN90p)",
    "Monthly Minimum Tmin (TNn)", "Annual Minimum Tmin (TNn)",
    "Monthly Maximum Tmin (TNx)", "Annual Maximum Tmin (TNx)",
    "Tropical Nights (TR)", "Monthly Cool Days (TX10p)",
    "Annual Cool Days (TX10p)",
    "Monthly Warm Days (TX90p)", "Annual Warm Days (TX90p)",
    "Monthly Minimum Tmax (TXn)", "Annual Minimum Tmax (TXn)",
    "Monthly Maximum Tmax (TXn)", "Annual Maximum Tmax (TXn)",
    "Warm Spell Duration Index (WSDI)"
  ),
  stringsAsFactors = FALSE
)

idx_df$idx_etccdi_time <- paste(idx_df$idx_etccdi, "ETCCDI_",
                               idx_df$time, sep = "")

# Unfortunatley expressions cannot be added to dataframes.
# These expreesion are required for the timeseries.
idx_ylab <- c(expression(
  "days", "days", "days", "days", "days", "days", "days",
  paste(degree, "C"), paste(degree, "C"), "days", "days", "days",
  "mm", "days", "days", "days", "mm", "mm", "mm",
  "mm", "mm", "mm", "mm day^-1",
  "days", "Exceedance rate, %", "Exceedance rate, %",
  "Exceedance rate, %", "Exceedance rate, %",
  paste(degree, "C"), paste(degree, "C"),
  paste(degree, "C"), paste(degree, "C"),
  "days", "Exceedance rate, %", "Exceedance rate, %",
  "Exceedance rate, %", "Exceedance rate, %",
  paste(degree, "C"), paste(degree, "C"),
  paste(degree, "C"), paste(degree, "C"), "days"
))