#2015 categorical data

#####METADATA#####
#loads categorical data about land holding:
#BR = Brecksville Reservation, HI = Hinckley Reservation, MS = Mill Stream Reservation, 
#NC = North Chagrin Reservation, RR = Rocky River Reservation, and WC = West Creek Reservation)

#forest canopy type:
#BM = Beech-Maple
#FP = Floodplain
#M = Mixed
#OAK = Oak

#drainage type (WELL OR POORLY drained) and loam content (1-4) were not included in this analysis

#2015 categorical data.csv

#outputs: a table of land holding and forest type, and a data frame with categorical data for each plot.

#####specify working directory#####
setwd("~/Documents/Research/PCAP/PCAPdata")
#####read in data#####
categorical<-read.csv('2015 categorical data.csv', header = TRUE, sep = ',')
#####view relevant categorical data#####
table(categorical[,3:4])

#############################################################################################
