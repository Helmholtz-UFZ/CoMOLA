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
library(raster)

# set prerequisites
modelnames <- c("HabStruct", "SAR", "SYM", "WYLD")
pal <- c('lightgoldenrod','goldenrod1','goldenrod2','goldenrod','goldenrod3','darkgreen','lightgreen','grey')
rasternames <- c('Arable land 1','Arable land 2','Arable land 3','Arable land 4','Arable land 5','Forest','Pasture','Urban Area')

#function for extracting fitness-scores for 4 objectives in given dataframe, can be adjusted to 3 objectives respectively
extract_fitness <- function(df){
  df <- cbind.data.frame(as.numeric(gsub('\\[', '', df[,ncol(df)-3])), #for 3 objectives, change to '-2'
                         as.numeric(df[,ncol(df)-2]), #for 3 objectives, delete line
                         as.numeric(df[,ncol(df)-1]),
                         as.numeric(gsub('\\]', '', df[,ncol(df)])))
}

############################################################
#                 Read best solutions
############################################################

# Adapt paths to your best solutions csv files from each run
# In the best solutions csv file, fitness values for the three objectives are stored in the following order: HabStruct, HI, WQ
# This code creates data frames with fitness values and their corresponding solution IDs for each run

# Parse all best_solution.csv-files from the output-folder into a list called files
filenames <- list.files(here('output'), pattern = 'best_solutions.csv')
filenames <- paste(here('output'),filenames, sep ='/')
files <- lapply(filenames, read.csv, h=F, skip=1, as.is=T)

#process list of files
files <- lapply(files, extract_fitness)
files <- lapply(files, setNames, modelnames)
#apply ID numbers
files <- imap(files, ~.x %>% mutate(ID = paste(.y, row_number(), sep = ".")))


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
#Store the single maxima in a dataframe
max_Solutions <- overall_best_solutions[apply(overall_best_solutions[1:length(modelnames)], MARGIN=2, which.max),]

## Compromise solution (sum of deviations from single mean values is minimal)
which.mean <- function(w,x,y,z) which.min(abs(w - mean(w))+abs(x - mean(x))+abs(y - mean(y))+abs(z - mean(z)))
sol_compromise <- overall_best_solutions[which.mean(overall_best_solutions$HabStruct,overall_best_solutions$SAR,overall_best_solutions$SYM,overall_best_solutions$WYLD),c(1:4)]

#Extract ID's and process the list for easier access to plot, first max_IDs
max_IDs <- max_Solutions[, (length(modelnames)+1):ncol(max_Solutions)]
max_IDs <- lapply(data.frame(t(max_IDs)), function(x) {x[!is.na(x)]})
names(max_IDs) <- paste('Max_', modelnames, sep='')
max_IDs <- lapply(max_IDs, "[", seq(min(lengths(max_IDs))))
#second compromised ID
ID_compromise <- overall_best_solutions[which.mean(overall_best_solutions$HabStruct,overall_best_solutions$SAR,overall_best_solutions$SYM,overall_best_solutions$WYLD),c(5:(5+length(IDs)-1))]
ID_compromise <- ID_compromise[, colSums(is.na(ID_compromise))==0]
names(ID_compromise) <- 'Compromise'
IDs <- append(max_IDs, ID_compromise)
IDs <- lapply(IDs, function(x) {as.integer(unlist(strsplit(x, '[.]')))})

#remove IDs in max_Solutions and rename rownames
max_Solutions <- max_Solutions[,(1:length(modelnames))]
rownames(max_Solutions) <- paste('Max_', modelnames, sep="")

#build the path for max- and compromised solutions and store in list
mapnames <- lapply(IDs, function(x){paste(paste(head(strsplit(tail(strsplit(filenames[x[1]], '/')[[1]],1),'_')[[1]],2),collapse='_'),'_best_ascii_map',x[2],'.asc',sep='')})
mapnames <- paste(here('output'),mapnames, sep='/')
maps <- lapply(mapnames, raster)
names(maps) <- names(IDs)

############################################################
#                  I Plot
############################################################

sol <- overall_best_solutions[,(1:length(modelnames))]

# You can change the axes on which objectives should be displayed - don't forget to adapt axis labels respectively
ggplot(sol, aes(x=HabStruct, y=SAR)) +
  geom_point(aes(size = WYLD, fill=SYM), shape=21, col="black")+
  
  #higlight the compromised solution
  geom_point(data = sol_compromise, aes(x=HabStruct, y = SAR),col='gold')+
  
  #highlight the maximised solutions
  geom_point(data = max_Solutions['Max_HabStruct',], aes(x=HabStruct, y = SAR), col='white')+
  geom_point(data = max_Solutions['Max_SAR',], aes(x=HabStruct, y=SAR), col='white')+
  geom_point(data = max_Solutions['Max_SYM',], aes(x=HabStruct, y=SAR), col='white')+
  geom_point(data = max_Solutions['Max_WYLD',], aes(x=HabStruct, y=SAR), col ='white')+
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(fill = "Agricultural Yield", size = "Water Yield")+
  xlab("Habitat quality based on landscape structure")+
  ylab("Species Area relationship")

#plot every map in maps-list
for (i in 1:length(maps)){
  plot(maps[[i]], 
       col=pal,
       main=names(maps)[i],
       legend = F,
       axes = F,
       box = F)
  legend(x='topright', 
         inset = c(-0.25, 0),
         fill=pal, 
         legend=rasternames, 
         border = F, 
         bty='n',
         xpd=T)
}