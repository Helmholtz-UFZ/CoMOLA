##########################################################################################################
#
#                                  Overall Pareto-Frontier Analysis
#
# Description: Identification of the overall Pareto-frontier for multiple optimization runs with CoMOLA,
#              identification of single maxima and compromise solution, plot of results. This is an
#              example for n runs and 4 objectives which can be adapted to n objectives, respectively.
#
# Author: Andrea Kaim & Victor Steffens
# Date: 20-09-2023
#
##########################################################################################################

# Install required packages
# install.packages("dplyr")
# install.packages("plot3D")
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("emoa")
# install.packages("here")
# install.packages("purrr")
# load packages

library(plyr)
library(dplyr)
library(ggplot2)
library(plot3D)
library(viridis)
library(emoa)
library(purrr)
library(here)

# set prerequisites
modelnames <- c("HabStruct", "SAR", "SYM", "WYLD")

#function for extracting fitness-scores for 4 objectives in given dataframe
extract_fitness <- function(df){
  df <- cbind.data.frame(as.numeric(gsub('\\[', '', df[,ncol(df)-3])),
                         as.numeric(df[,ncol(df)-2]),
                         as.numeric(df[,ncol(df)-1]),
                         as.numeric(gsub('\\]', '', df[,ncol(df)])))
}

############################################################
#                 Read best solutions
############################################################

# Adapt paths to your best solutions csv files from each run
# In the best solutions csv file, fitness values for the three objectives are stored in the following order: HabStruct, HI, WQ
# This code creates data frames with fitness values and their corresponding solution IDs for each run

#look for all files in defined folder named 'best_solutions.csv' and parse them as a list of dataframes into files-variable
filenames <- list.files(here('output'), pattern = 'best_solutions.csv')
filenames <- paste(here('output'),filenames, sep ='/')
files <- lapply(filenames, read.csv, h=F, skip=1, as.is=T)

#process list of files
files <- lapply(files, extract_fitness)
files <- lapply(files, setNames, modelnames)
#apply ID numbers
files <- imap(files, ~.x %>% mutate(ID = paste(.y, row_number(), sep = ".")))


#for (i in 1:length(files)){
#  for (j in 1:nrow(files[[i]])){
#    files[[i]]$ID[j] <- paste(i,j, sep='.')
#  }
#}


############################################################
#              Overall Pareto-Frontier
############################################################


# Create list with Pareto-optimal solutions of all optimization runs
S <- subset(do.call(rbind, files), select = -ID)

# transpose S
S_trans <- t(S)

# Identify non-dominated set of solutions (multiplication by -1 because nondominated_points uses minimization as default but here, we maximize)
# Check documentation of emoa package for details
P <- nondominated_points(as.matrix(S_trans * -1)) * -1

# Returns TRUE if a solution is dominated and FALSE if not
# is_dominated(as.matrix(S_trans * -1))

P <- as.data.frame(t(P))

############################################################
#                        Get IDs
############################################################

# Identify IDs from all runs
IDs <- lapply(files, function(df) {
  joined_df <- left_join(P,df)
  return(joined_df)
})

# Join IDs to best solutions (one solution can correspond to multiple IDs if there are solutions with the same fitness values but from different runs - !!! still, maps can look different)
overall_best_solutions <-cbind(IDs[[1]])
for (i in 2:length(IDs)){
  overall_best_solutions <- cbind(overall_best_solutions, IDs[[i]]$ID)
}

#prepare colnames and assign them to overall_best_solutions
modelnumber <- c(1:length(filenames))
IDnames <- paste(rep("ID"),modelnumber, sep='_')
colnames(overall_best_solutions) <- c(modelnames,IDnames)


############################################################
#                   Analyze solutions
############################################################

## Single maxima

# HabStruct
sol_max_HabStruct <- overall_best_solutions[which.max(overall_best_solutions$HabStruct),c(1:4)] # returns fitness values of maximum solution for HabStruct
ID_max_HabStruct  <- overall_best_solutions[which.max(overall_best_solutions$HabStruct),c(5:(5+length(IDs)-1))] # returns map IDs of maximum solution for HabStruct

# SAR
sol_max_SAR <- overall_best_solutions[which.max(overall_best_solutions$SAR),c(1:4)] # returns fitness values of maximum solution for SAR
ID_max_SAR  <- overall_best_solutions[which.max(overall_best_solutions$SAR),c(5:(5+length(IDs)-1))] # returns map IDs of maximum solution for SAR

# SYM
sol_max_SYM <- overall_best_solutions[which.max(overall_best_solutions$SYM),c(1:4)] # returns fitness values of maximum solution for SYM
ID_max_SYM  <- overall_best_solutions[which.max(overall_best_solutions$SYM),c(5:(5+length(IDs)-1))] # returns map IDs of maximum solution for SYM

# WYLD
sol_max_WYLD <- overall_best_solutions[which.max(overall_best_solutions$WYLD),c(1:4)] # returns fitness values of maximum solution for WYLD
ID_max_WYLD  <- overall_best_solutions[which.max(overall_best_solutions$WYLD),c(5:(5+length(IDs)-1))] # returns map IDs of maximum solution for WYLD

## Compromise solution (sum of deviations from single mean values is minimal)
which.mean <- function(w,x,y,z) which.min(abs(w - mean(w))+abs(x - mean(x))+abs(y - mean(y))+abs(z - mean(z)))
sol_compromise <- overall_best_solutions[which.mean(overall_best_solutions$HabStruct,overall_best_solutions$SAR,overall_best_solutions$SYM,overall_best_solutions$WYLD),c(1:4)]
ID_compromise <- overall_best_solutions[which.mean(overall_best_solutions$HabStruct,overall_best_solutions$SAR,overall_best_solutions$SYM,overall_best_solutions$WYLD),c(5:(5+length(IDs)-1))]

## How to find the maps: Check the IDs, if the result looks, e.g. like this:
# ID_1 ID_2 ID_3
# 38 <NA> 2.12  3.5
# there are two maps corresponding to the same fitness values. One from run 2 (solution ID 12) and one from run 3 (solution ID 5). In the output folder of the respective optimization run,
# search for the shape file with the pattern e.g. 15-09-2021_14-09-07_best_ascii_map5.shp (here, output of run 3)


############################################################
#                   3D Plot
############################################################

sol <- overall_best_solutions[,1:4]

# You can change the axes on which objectives should be displed, e.g. by setting x=WQ, y=HabStruct, and fill=HabStruct - also adapt axis labels respectively

ggplot(sol, aes(x=HabStruct, y=SAR)) +
  geom_point(aes(size = WYLD, fill=SYM), shape=21, col="black")+
  
  #higlight the compromised solution
  geom_point(data = sol_compromise, aes(x=HabStruct, y = SAR),col='gold')+
  
  #highlight the maximised solutions
  geom_point(data = sol_max_HabStruct, aes(x=HabStruct, y = SAR), col='white')+
  geom_point(data = sol_max_SAR, aes(x=HabStruct, y=SAR), col='white')+
  geom_point(data = sol_max_SYM, aes(x=HabStruct, y=SAR), col='white')+
  geom_point(data = sol_max_WYLD, aes(x=HabStruct, y=SAR), col ='white')+
  
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(fill = "Agricultural Yield", size = "Water Yield")+
  xlab("Habitat quality based on landscape structure")+
  ylab("Species Area relationship")

