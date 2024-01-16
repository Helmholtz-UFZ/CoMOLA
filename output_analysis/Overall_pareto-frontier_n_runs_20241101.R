##########################################################################################################
#
#                                  Overall Pareto-Frontier Analysis
#
# Description: Identification of the overall Pareto-frontier for multiple optimization runs with CoMOLA,
#              identification of single maxima and compromise solution, plot of results. This is an
#              example for n runs and 4 objectives which can be adapted to n objectives, respectively.
#
#              !!! If you run this code with the best_solutions file of only one optimization run, exclude
#                  lines 106-108 and 135
#
# Authors: Andrea Kaim & Victor Steffens
# Date: 10-10-2023
#
##########################################################################################################

# Install required packages
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("emoa")
# install.packages("here")
# install.packages("purrr")
# install.packages("reshape")
# load packages
here::i_am('output_analysis/Overall_pareto-frontier_n_runs_20241101.R')

library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(emoa)
library(purrr)
library(here)
library(raster)
library(plotly)
library(reshape)
library(GGally)

#Set prerequisites
#The order of modelnames corresponds to the order in the config.ini-file from the optimisation-process
modelnames <- c("HabStruct", "SYM", "WYLD", "SAR")
pal <- c('#FFF9C4', '#FFF176', '#FBC02D','goldenrod', 'darkgoldenrod', '#8BC34A', 'darkgreen', 'red' )
rasternames <- c('Cropland 1','Cropland 2','Cropland 3','Cropland 4','Cropland 5','Pasture','Forest','Urban')
location <- here('output') #paste the folder, where you put the results of the optimization process
#location <- "C:/Users/kaim/Nextcloud/Cloud/CoMOLA_Blockseminar/WS2324/Output_virtual_study/output" # use this if here does not work for you, paste directory of your output folder

#old function for extracting fitness values into new dataframe
#function for extracting fitness-scores for 4 objectives in given data frame, can be adjusted to 3 objectives respectively
#extract_fitness <- function(df){
#  df <- cbind.data.frame(as.numeric(gsub('\\[', '', df[,ncol(df)-3])), #for 3 objectives, change to '-2'
#                         as.numeric(df[,ncol(df)-2]), #for 3 objectives, delete line
#                         as.numeric(df[,ncol(df)-1]),
#                         as.numeric(gsub('\\]', '', df[,ncol(df)])))
#}

#new function to extract fitness values into new dataframe -> no need to adjust number of variables anymore
extract_fitness <- function(bs){
  df <- data.frame()
  for (x in 1:nrow(bs)){
    df <- rbind(df, as.numeric(unlist(strsplit(as.vector(regmatches(toString(bs[x,]), gregexpr("(?<=\\[).*?(?=\\])", toString(bs[x,]), perl=T))[[1]][2]),','))))
  }
  bs <- df
}
############################################################
#                 Read best solutions
############################################################

# This code creates data frames with fitness values and their corresponding solution IDs for each run

# Parse all best_solution.csv-files from the output-folder into a list called files
filenames <- list.files(location, pattern = 'best_solutions.csv')
#check if only the best solutions you want to analyze are stored inside filenames, if not, then subset filenames: filenames <- filenames[c(x,x,x)]
print(filenames)
filenames <- paste(location,filenames, sep ='/')
files <- lapply(filenames, read.csv, h=F, skip=1, as.is=T)

#process list of files, the extract_fitness-function could take a while
files <- lapply(files, extract_fitness)
files <- lapply(files, setNames, modelnames)
#apply ID numbers
files <- lapply(1:length(files), \(x) transform(files[[x]], ID = paste0(x, ".", 1:nrow(files[[x]]))))


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
#Duplikate noch rausholen


Pdist <- distinct(P)


############################################################
#                        Get IDs
############################################################

# Identify IDs from all runs
#IDs <- lapply(files, function(df) {
#  joined_df <- left_join(Pdist,df, multiple = "all") # It can happen that within a best_solutions data frame, one solution corresponds to multiple IDs. In this case, select first one (multiple = "all" would be better)
#  return(joined_df)
#})
# Join IDs to best solutions (one solution can correspond to multiple IDs if there are solutions with the same fitness values but from different runs - !!! still, maps can look different)
#overall_best_solutions <- IDs[[1]]
#for (i in 2:length(IDs)){
#  overall_best_solutions <- cbind(overall_best_solutions, IDs[[i]]$ID)
#}
joinlist <- append(list(Pdist), files, 1)
overall_best_solutions <- reduce(joinlist, left_join, by = modelnames)

#prepare colnames and assign them to overall_best_solutions
runnumber <- c(1:length(filenames))
IDnames <- paste(rep("ID"),runnumber, sep='_')
colnames(overall_best_solutions) <- c(modelnames,IDnames)


############################################################
#                   Analyze solutions
############################################################

## Single maxima
#Store the single maxima in a dataframe
max_Solutions <- overall_best_solutions[apply(overall_best_solutions[1:length(modelnames)], MARGIN=2, which.max),]

## Compromise solution (sum of deviations from single mean values is minimal)
which.mean <- function(w,x,y,z) which.min(abs(w - mean(w))+abs(x - mean(x))+abs(y - mean(y))+abs(z - mean(z)))
sol_compromise <- overall_best_solutions[which.mean(overall_best_solutions$HabStruct,overall_best_solutions$SAR,overall_best_solutions$SYM,overall_best_solutions$WYLD),]

#Extract ID's and process the list for easier access to plot, first max_IDs
max_IDs <- max_Solutions[, (length(modelnames)+1):ncol(max_Solutions)]
max_IDs <- lapply(data.frame(t(max_IDs)), function(x) {x[!is.na(x)]})
names(max_IDs) <- paste('Max_', modelnames, sep='')
max_IDs <- lapply(max_IDs, "[", seq(min(lengths(max_IDs))))
#second compromised ID
ID_compromise <- sol_compromise[,(length(modelnames)+1):c(length(modelnames)+length(runnumber))]
ID_compromise <- ID_compromise[, colSums(is.na(ID_compromise))==0]
names(ID_compromise) <- 'Compromise'
IDs <- append(max_IDs, ID_compromise)
IDs <- lapply(IDs, function(x) {as.integer(unlist(strsplit(x, '[.]')))})

#remove IDs in max_Solutions and rename rownames
max_Solutions <- max_Solutions[,(1:length(modelnames))]
rownames(max_Solutions) <- paste('Max_', modelnames, sep="")

#build the path for max- and compromised solutions and store in list
mapnames <- lapply(IDs, function(x){paste(paste(head(strsplit(tail(strsplit(filenames[x[1]], '/')[[1]],1),'_')[[1]],2),collapse='_'),'_best_ascii_map',x[2],'.asc',sep='')})
mapnames <- paste(location,mapnames, sep='/')
maps <- lapply(mapnames, raster)
names(maps) <- names(IDs)

############################################################
#                  I Plot
############################################################

sol <- overall_best_solutions[,(1:length(modelnames))]

#max_Solutions['Current',] <- c(23.86667, 7.924466, 89.67191, 95.39818) # add fitness values of current land-use map

# You can change the axes on which objectives should be displayed - don't forget to adapt axis labels respectively

# Plot with x-axis HabStruct, y-axis SAR, dots WYLD, color SYM
ggplot(sol, aes(x=HabStruct, y=SAR)) +
  geom_point(aes(size = WYLD, fill=SYM), shape=21, col="black")+
  
  #highlight the maximised solutions
  geom_point(data = max_Solutions['Max_HabStruct',], aes(x=HabStruct, y = SAR), col='red', size = 5, shape = 17)+
  geom_point(data = max_Solutions['Max_SAR',], aes(x=HabStruct, y = SAR), col='red', size = 5, shape = 15)+
  geom_point(data = max_Solutions['Max_SYM',], aes(x=HabStruct, y = SAR), col='red', size = 5, shape = 18)+
  geom_point(data = max_Solutions['Max_WYLD',], aes(x=HabStruct, y = SAR), col='red', size = 5, shape = 19)+
  geom_point(data = sol_compromise, aes(x=HabStruct, y = SAR), col='red', size = 7, shape = 8, stroke = 3)+
  #geom_point(data = max_Solutions['Current',], aes(x=HabStruct, y = SAR), col='white', size = 7, shape = 8, stroke = 3)+
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(fill = "Agricultural Yield", size = "Water Yield")+
  xlab("Habitat heterogeneity")+
  ylab("Species richness")


# Plot with x-axis SYM, y-axis WYLD, dots SAR, color HabStruct
ggplot(sol, aes(x=SYM, y=WYLD)) +
  geom_point(aes(size = SAR, fill=HabStruct), shape=21, col="black")+
  
  #highlight the maximised solutions
  geom_point(data = max_Solutions['Max_HabStruct',], aes(x=SYM, y = WYLD), col='red', size = 5, shape = 17)+
  geom_point(data = max_Solutions['Max_SAR',], aes(x=SYM, y = WYLD), col='red', size = 5, shape = 15)+
  geom_point(data = max_Solutions['Max_SYM',], aes(x=SYM, y = WYLD), col='red', size = 5, shape = 18)+
  geom_point(data = max_Solutions['Max_WYLD',], aes(x=SYM, y = WYLD), col='red', size = 5, shape = 19)+
  geom_point(data = sol_compromise, aes(x=SYM, y = WYLD), col='red', size = 7, shape = 8, stroke = 3)+
  #geom_point(data = max_Solutions['Current',], aes(x=SYM, y = WYLD), col='white', size = 7, shape = 8, stroke = 3)+
  
  scale_fill_gradientn(colours = viridis(100))+
  labs(fill = "Habitat heterogeneity", size = "Species richness")+
  xlab("Agricultural Yield")+
  ylab("Water Yield")


# Plot with x-axis SYM, y-axis WYLD, dots SAR, color HabStruct and legend for maximum solutions
ggplot(sol, aes(x = SYM, y = WYLD)) +
  geom_point(aes(size = SAR, fill = HabStruct), shape = 21, col = "black") +
  
  # Highlight the maximized solutions
  geom_point(data = max_Solutions['Max_HabStruct',], aes(x = SYM, y = WYLD, shape = "Max_HabStruct"), col='red', size = 5) +
  geom_point(data = max_Solutions['Max_SAR',], aes(x = SYM, y = WYLD, shape = "Max_SAR"), col='red', size = 5) +
  geom_point(data = max_Solutions['Max_SYM',], aes(x = SYM, y = WYLD, shape = "Max_SYM"), col='red', size = 5) +
  geom_point(data = max_Solutions['Max_WYLD',], aes(x = SYM, y = WYLD, shape = "Max_WYLD"), col='red', size = 5) +
  geom_point(data = sol_compromise, aes(x = SYM, y = WYLD, shape = "Compromise"), col='red', size = 7, stroke = 3) +
  #geom_point(data = max_Solutions['Current',], aes(x=SYM, y = WYLD, shape = "Current"), col='red', size = 7, stroke = 3)+
  
  scale_fill_gradientn(colours = viridis(100)) +
  labs(fill = "Habitat heterogeneity", size = "Species richness") +
  xlab("Agricultural Yield") +
  ylab("Water Yield") +
  
  scale_shape_manual(name = "Shapes",
                     values = c("Compromise" = 8,"Current" = 4, "Max_HabStruct" = 17, "Max_SAR" = 15, "Max_SYM" = 18, "Max_WYLD" = 19),
                     labels = c("Compromise", "Current", "Max_HabStruct", "Max_SAR", "Max_SYM", "Max_WYLD")) +
  
  guides(size = guide_legend(title = "Species Richness"),
         fill = guide_legend(title = "Habitat Heterogeneity"),
         shape = guide_legend(title = "Shapes"))

#do a Pairwise plot
ggpairs(sol)

#plot every map in maps-list
plots <- lapply(names(maps), function(x){
  plot(maps[[x]],
       col= pal,
       main = x,
       legend = F,
       axes = F,
       box = F)
  legend(x ='topright',
         inset = c(-0.2,0.25),
         fill = pal,
         legend = rasternames,
         border =F,
         bty ='n',
         xpd = T)
}
)

# Parallel coordinates plot
# save as web page

parcord <- overall_best_solutions %>%
  plot_ly(type = 'parcoords',
          #line = list(color = 'blue'),
          labelfont=list(size=20), # font size axes
          #width = 1500, height = 800,
          dimensions = list(
            list(range = c(min(overall_best_solutions$HabStruct),max(overall_best_solutions$HabStruct)),
                 label = 'Habitat heterogeneity', values = ~HabStruct),
            list(range = c(min(overall_best_solutions$SAR),max(overall_best_solutions$SAR)),
                 label = 'Forest species richness', values = ~SAR),
            list(range = c(min(overall_best_solutions$SYM),max(overall_best_solutions$SYM)),
                 label = 'Agricultural yield', values = ~SYM),
            list(range = c(min(overall_best_solutions$WYLD),max(overall_best_solutions$WYLD)),
                 label = 'Water yield', values = ~WYLD)
          )
  )
parcord
