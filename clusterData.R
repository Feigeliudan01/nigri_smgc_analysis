#!/usr/bin/env Rscript
library(optparse)

library(igraph)
library(parallel)
library(ggtree)

getwd()
# setwd("/home/seth/aspSM/secMet")

# print("Works on smurf and not smurf_papa!")

log_con <- file("clustering.log")

args<-commandArgs(TRUE)


print("Working on file")
print(args[1])
filename = args[1]
# filename = "sm_data_nigri_test.tsv"
print("Changing directory to")
print(args[2])
setwd(args[2])
# setwd("nigri_test")



# if(!file.exists("clusterBlastAll.csv")){
# option_list = list(
#   make_option(c("-s", "--setName"), action="store", default=NA, type='character',help="Provide set name"),
#   make_option(c("-f", "--fileName"), action="store", default=NA, type='character',help="Provide file name")
#   make_option(c("-t", "--treefile"), action="store", default=NA, type='character',help="Provide file name")
# )
# opt = parse_args(OptionParser(option_list=option_list))


aspDf <- read.table(filename, header = TRUE,  stringsAsFactors = FALSE, sep = '\t', quote = "\n")


clusterBlastAll <- read.csv("clusterBlastAll.csv", header = TRUE, stringsAsFactors = FALSE)
head(clusterBlastAll)

# In older version it was WHERE pident >= 40 on both smurf_bidir_hits
# write.csv( clusterBlastAll, file = "clusterBlastAll.csv", row.names = FALSE )

clusterBlast <- clusterBlastAll
head(clusterBlastAll)
# rm(clusterBlastAll)
str(clusterBlast)
clusterBlast <- clusterBlast[clusterBlast$pident_score != 0, c('q_clust_id', 'h_clust_id', 'pident_score')]

names(clusterBlast)[3] <- 'weight'

# Subsetting clusterBlast to the relevant ones
# Do i have to do this? giving it to vertices might achieve the same


cat("Missing gene clusters in clusterblast query", file = log_con, append = TRUE)
cat(paste(setdiff(clusterBlast$q_clust_id, aspDf$cluster_id), collapse = ",") , file = log_con, append = TRUE)
setdiff(clusterBlast$q_clust_id, aspDf$cluster_id)
cat("Missing gene clusters in downloaded data frame", file = log_con, append = TRUE)
cat(paste(setdiff(aspDf$cluster_id, clusterBlast$q_clust_id), collapse = ",") , file = log_con, append = TRUE)
missingOnes <- setdiff(aspDf$cluster_id, clusterBlast$q_clust_id)


cat("Missing gene clusters in clusterblast query", file = log_con, append = TRUE)
cat(paste(setdiff(clusterBlast$h_clust_id, aspDf$cluster_id), collapse = ","), file = log_con, append = TRUE)

clusterBlast <- clusterBlast[(clusterBlast$q_clust_id %in% aspDf$cluster_id & clusterBlast$h_clust_id %in% aspDf$cluster_id),]
print("Checking clusterBlast")


clusterDups <- clusterBlast[duplicated(clusterBlast[,c('q_clust_id', 'h_clust_id')]),c('q_clust_id', 'h_clust_id')]
cat(paste(clusterDups, collapse = ","), file = log_con, append = TRUE)

sum(clusterBlast$q_clust_id == clusterBlast$h_clust_id)
clusterDups

write.csv( clusterBlast, file = "clusterBlast.csv", row.names = FALSE )
# }else{
#   clusterBlast <- read.csv("clusterBlast.csv")
# }
# setwd(file.path(mainDir, ))



# Creating a whole network, subsetting later with dataframe which is already created


set.seed(123)

cat("Info about clusterBLast", file = log_con, append = TRUE)

# clusterBlastSub <- clusterBlast[clusterBlast$weight >20,]
cat("Minimum clusterBlast weight", file = log_con, append = TRUE)
cat(min(clusterBlast$weight), file = log_con, append = TRUE)

cat("ClusterBlast weight over 100", file = log_con, append = TRUE)
# cat(clusterBlastAll[clusterBlastAll$pident_score > 100,], file = log_con, append = TRUE)

# clusterBlast[clusterBlast$q_clust_id == "29_251500_26" & clusterBlast$h_clust_id == "25_350427_26",]
# clusterBlast[clusterBlast$q_clust_id == "5_494945_7" & clusterBlast$h_clust_id == "38_413327_17",]

print("Creating Network")
# Remove later
present_clusters <- unique(aspDf$cluster_id[aspDf$cluster_id %in% clusterBlast$q_clust_id | aspDf$cluster_id %in% clusterBlast$h_clust_id])
# unique(aspDf$cluster_id)
g.raw <- graph_from_data_frame(d=clusterBlast, vertices= present_clusters, directed=FALSE) # Decided on undirected graph because I it's a bidirectional relationship for all anyway
g <- igraph::simplify(g.raw, edge.attr.comb = list(weight = "mean")) #remove.multiple = TRUE)
walks <- cluster_walktrap(as.undirected(g), steps = 1)
walks
# sizes(walks)
# gw <- induced.subgraph(g, vids = walks[[16]])
# walkMore <- cluster_walktrap(gw)

# fg <- cluster_fast_greedy(g)

if(max(as.numeric(sizes(walks))) > length(unique(aspDf$org_id)) ){
    print("Network families are too large, needs more clustering")
    print("There are too many gene clusters per family, since they are exceeding the number of orgs (we assume that duplications will form another family, hence another round of clustering)")

    # print("Complete this section")
   superNetworks <- as.numeric(which(sizes(walks) > length(unique(aspDf$org_id))))

     firstNetworks <- lapply(names(groups(walks))[!(names(groups(walks)) %in% superNetworks)],function(x){
         sub <- induced.subgraph(g, vids = walks[[x]])
         })

     secondNetwork <- lapply(names(groups(walks))[names(groups(walks)) %in% superNetworks],function(x){
         sub <- induced.subgraph(g, vids = walks[[x]])
         walkMore <- cluster_walktrap(sub, steps = 1)
         lapply(names(groups(walkMore)), function(x){
             induced.subgraph(g, vids = walkMore[[x]])
             } )
         })

     testg <- induced_subgraph(g, vids = walks[[superNetworks[[1]]]])
     # plot(testg)
     # plot(secondNetwork[[1]][[4]])
     # secondNetwork

     secondNetwork <- unlist(secondNetwork, recursive = FALSE)
     cNetwork <- c(firstNetworks, secondNetwork) # complete network
     sizes <- sapply(cNetwork, function(x){length(V(x))})
          counter = 0
 clusterFamDf <- lapply(cNetwork, function(x){
     counter <<- counter + 1
     clusters <- V(x)$name
     data.frame(cluster_id = clusters, clusterFam = rep(counter, length(clusters)),
     stringsAsFactors = FALSE)
     })

 clusterFamDf

 clusterFamDf <- do.call(rbind, clusterFamDf)
 clusterFamDf$cluster_id <- as.character(clusterFamDf$cluster_id)

 # Filling up the ones with no hits
famLimit <-  max(clusterFamDf$clusterFam)
names(clusterFamDf)
missingClusters <- setdiff(unique(aspDf$cluster_id), clusterFamDf$cluster_id)
missingDf <- data.frame(cluster_id = missingClusters, clusterFam = seq(from = famLimit+1, to = famLimit+length(missingClusters), by = 1))

clusterFamDf <- rbind(clusterFamDf, missingDf)

viFams <- merge(aspDf, clusterFamDf, by = c("cluster_id"), all.x = TRUE)

viFams
# write.table(viFams, "smurfOrgIprMibigFams.tsv", row.names = TRUE)
# WRITE IT TO DISK
 # ggplot(as.data.frame(table(sizes))) + geom_bar(aes(sizes, Freq), stat = 'identity') + geom_text(aes(x = sizes, y = Freq, label = Freq) ,position=position_dodge(width=0.9), vjust=-0.25, size =  8 )+ labs( x = "Number of SMGC in family", y = "Family counts" , size = 2) + theme(axis.text = element_text(size = 14) , axis.title=element_text(size=12,face="bold"))

}else{

    print("No further clustering")
    coms <- communities(walks)
    tmpFam <- lapply(names(coms), function(x){
        data.frame(cluster_id = coms[[x]], clusterFam = (x), stringsAsFactors = FALSE)
        })

    famDf <- do.call(rbind, tmpFam)


    viFams <- merge(aspDf, famDf, by = c("cluster_id"), all.x = TRUE)

    # write.table(viFams, "smurfOrgIprMibigFams.tsv", row.names = TRUE)
}

write.table(viFams, paste0(gsub(".tsv","", filename),"_c.tsv"), row.names = FALSE, sep = "\t")

# names(viFams)
# "147_209414_1" %in% missingOnes
# sapply(missingOnes, function(x){
#   x %in% viFams$cluster_id
# })
subset(viFams, cluster_id == "147_209414_1")

