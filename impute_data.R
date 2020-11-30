# script to impute missing data
rm(list = ls())
library(ggplot2)
library(mice) 
library(imputeTS)
library(data.table)
library(dplyr)
# last update november 26.. filling data from October 10 and onwards. 
raw_file = fread("deaths_with_missing_data.csv")

# mice dosn't work since it dosn't interpolate 

# states to fill: 
# GA, we have true daily incidence, this needs to be split to confirmed/prob
# ..interpolate deathProbable and deathConfirmed... then divide them to get % 
# ..use that % to split the CUMULATIVE into the number of probable cases 
# ..then reported - probable gives confirmed (and the incidences all add up)
# ..EASIER WAY.. just interpolate the probable.. its the same
# MN ... leave probable cases at 53. Adjust reported death so that reported = confirmed plus prob. 
# .. data gets messed up on october 10
# MN row 40, manual change based on incdeicne of confirmed cases
# NV, since historically they have no probable cases, reported deaths = confirmed deaths
# TX, same as NV.. no prob cases
# VT, same as above, no prob cases -- also easier since report + confirm is given
# WA, same as above

mdat = data.frame(date=raw_file$date, t=seq(304, 1), dd = raw_file$deathConfirmed_MN)
md.pattern(mdat) # from mice, gives us which data is missing
intdat = round(na.interpolation(mdat$dd, option = "linear"))
na_idx = which(is.na(mdat$dd))
cat(intdat)
View(intdat[na_idx])


# use regression? Not a good option
lmHeight = lm(dd~1 + t, data = mdat) #Create the linear regression
summary(lmHeight) #Review the results
predict(lmHeight, mdat)

## to paste data back into excel 
## dosn't work on hpc 
tabout <- function(output){
  print(output)
  capture.output(output, file = "clipboard", append = FALSE, 
                 type = "output", split = FALSE)
  lines <- readClipboard()
  for(i in 1 : 5) {lines <- gsub("  ", " ", lines, fixed=TRUE)}
  lines <- gsub(" ", "\t", lines, fixed=TRUE)
  writeClipboard(lines)
}

# dosn't work on hpc
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

