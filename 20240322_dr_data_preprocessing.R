## pre-processing of data for dimensionality reduction

library( EmbedSOM )
library( ConsensusClusterPlus )
library( dplyr )
library( tidyr )
library( data.table )
library( purrr )
library( RColorBrewer )
library( ggridges )
library( readxl )

date.seed <- 20240322

# read csv data

data.dir <- "./scaled_data"

read_plus <- function(flnm) {
  fread(flnm) %>%
    mutate(filename = flnm)
}

do.call_rbind_fread <- function(path, pattern = "\\.csv$") {
  files = list.files(path, pattern, full.names = TRUE)
  do.call(rbind, lapply(files, function(x) read_plus(x)))
}

ln.files <- do.call_rbind_fread( path = data.dir, pattern = "LN")
gut.files <- do.call_rbind_fread( path = data.dir, pattern = "LPL")
flow.data <- rbind(ln.files, gut.files)

# tidy up filenames and set rownames
flow.data <- flow.data %>% separate(col = filename, into = c("X1", "X2", "X3", 
                                                             "tissue", "mouse", "sample", "X4"),
                                    sep = "_" ) %>%
  select( -c("X1", "X2", "X3", "mouse", "X4"))

flow.data$rownames <- paste(flow.data$tissue, flow.data$sample, seq_len(nrow(flow.data)), sep = "_")

rownames(flow.data) <- flow.data$rownames

# sample data, select columns for analysis, keep rownames
dmrd.data <- flow.data %>% group_by( tissue, sample ) %>%
  sample_n(20000)

dmrd.data <- dmrd.data %>% select( -c(rownames, Time, viability, `FJComp-AF-A`, 
                                      `SSC-H`, `SSC-B-H`, `SSC-B-A`, `SSC-A`,
                                      `FSC-H`, `FSC-A`))

# run clustering and assign as a factor in the data
n.clusters <- 50

set.seed(date.seed)

flow.som <- EmbedSOM::SOM(dmrd.data[,1:49], xdim = 24, 
                          ydim = 24, batch = TRUE,
                          parallel = TRUE, threads = 0 )

# get clusters
flow.som.mapping <- flow.som$mapping[ , 1 ]
flow.som.codes <- flow.som$codes

# get clusters from som mapping
consensus.cluster <- ConsensusClusterPlus( t( flow.som.codes ),
                                           maxK = n.clusters, reps = 100, pItem = 0.9, pFeature = 1,
                                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                           distance = "euclidean", 
                                           seed = set.seed(date.seed) )

flow.som.event.cluster <- consensus.cluster[[ n.clusters ]]$
  consensusClass[ flow.som.mapping ]

# reorder clusters from bigger to smaller
flow.som.cluster.rank <- 1 + n.clusters - 
  rank( table( flow.som.event.cluster ), ties.method = "last" )
flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
names( flow.som.event.cluster ) <- NULL

# set clusters as a factor
dmrd.data$Cluster <- factor( flow.som.event.cluster, 
                                  levels = 1 : n.clusters )


color.pool <- c( 
  brewer.pal( 8, "Set1" )[ -6 ], 
  brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
  adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
               red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
  adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
               red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
  adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
               red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
  adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
               red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
  adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
               red.f = 0.3, green.f = 0.3, blue.f = 0.3 ), 
  adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
               red.f = 0.3, green.f = 0.3, blue.f = 0.3 ) )
color.pool.n <- length( color.pool )

cluster.color <- rep( color.pool, 
                          ceiling( n.clusters / color.pool.n ) )[ 1 : n.clusters ]

dmrd.data$Color_clust <- cluster.color[dmrd.data$Cluster]

# add in cluster labels for cell types
dmrd_long <- dmrd.data[,c(1:49,52)] %>%
  pivot_longer( -Cluster, names_to = "Marker", values_to = "Expression")

ggplot(dmrd_long, aes(x = Expression, y = Marker, fill = Expression)) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.001) +
  facet_wrap(~ Cluster) +
  theme_ridges() +
  labs(x = "Expression", y = "Marker")

ggsave("cluster_histograms.jpg", width = 20, height = 50, limitsize = FALSE)

# automated cluster naming

source("prepare_marker_lists.r")
source("flow_cluster_id_score.r")

cell.database <- read_xlsx( "mouse_celltype_database.xlsx" )
cell.database <- dplyr::filter(cell.database, Tissue.restricted == "Immune")

flow.data.cluster.mfi <- dmrd.data[,1:53] %>% 
  ungroup() %>%
  select(-c(tissue, sample, Color_clust)) %>%
  group_by(Cluster) %>%
  summarize_all(mean)

rownames(flow.data.cluster.mfi) <- flow.data.cluster.mfi$Cluster
flow.data.cluster.mfi <- flow.data.cluster.mfi %>% select(-Cluster)

scaled.cluster.mfi <- scale(flow.data.cluster.mfi, scale = FALSE)

colnames(scaled.cluster.mfi)[3] <- "Gr1"
colnames(scaled.cluster.mfi)[5] <- "F480"
colnames(scaled.cluster.mfi)[10] <- "PDCA1"
colnames(scaled.cluster.mfi)[12] <- "Ly6C"
colnames(scaled.cluster.mfi)[16] <- "CTLA4"
colnames(scaled.cluster.mfi)[17] <- "cKit"
colnames(scaled.cluster.mfi)[22] <- "SiglecF"
colnames(scaled.cluster.mfi)[23] <- "TCRbeta"
colnames(scaled.cluster.mfi)[24] <- "PD1"
colnames(scaled.cluster.mfi)[35] <- "GATA3"
colnames(scaled.cluster.mfi)[49] <- "Tbet"

marker.list <- prepare_marker_lists( "mouse_celltype_database.xlsx", 
                                     "Immune", "All" )

cluster.label <- sprintf( "%02d", 1 : n.clusters )


scaled.id.score <- flow_cluster_id_score(t(scaled.cluster.mfi),
                                         marker_pos = marker.list$markers_positive, 
                                         marker_neg = marker.list$markers_negative )

for (cluster in 1:ncol(scaled.id.score)){
  cluster.label[cluster] <- names( which.max(scaled.id.score[,cluster]) )
}

dmrd.data$Cluster_name <- cluster.label[dmrd.data$Cluster]

names(cluster.label) <- ( 1 : n.clusters )

str(dmrd.data)

save.image(file = "preprocessed_data.RData")

