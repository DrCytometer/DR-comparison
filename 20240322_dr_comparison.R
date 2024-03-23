# compare dimensionality reduction tools for visualization of flow cytometry data

# install and start required packages
required.packages <- c("Rtsne", "ggplot2", "EmbedSOM", "scattermore",
                       "uwot", "phateR", "reticulate" )

for (req.package in required.packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    install.packages(req.package, repos='http://cran.us.r-project.org')
  }
}

invisible( lapply( required.packages, library, character.only = TRUE ) )

# set seed, load data
date.seed <- 20240322

load( file = "preprocessed_data.RData")

input.data <- dmrd.data[,1:49]



#### dr methods comparison

# pca----------------------

system.time(
  pca.result <- prcomp( input.data )
)
# user  system elapsed 
# 1.18    0.02    1.36 

dmrd.data <- cbind(dmrd.data, pca.result$x[,1:2])

ggplot(dmrd.data, aes(x = PC1, y = PC2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("pca_plot.jpg", width = 9, height = 6 )



# tSNE (OptSNE-like)----------------------

tsne.learning.rate <- ifelse( nrow(input.data)/4 > 2000, 
                              nrow(input.data)/4, 2000 )
set.seed( date.seed )

system.time(
  tsne.result <- Rtsne( input.data, perplexity = 30, 
                        exaggeration_factor = 4,
                        eta = tsne.learning.rate,
                        max_iter = 750, stop_lying_iter = 75,
                        check_duplicates = FALSE, pca = FALSE, 
                        num_threads = 0 )
)
user  system elapsed 
1446.50    2.36  323.33

tsne.data <- tsne.result$Y

colnames(tsne.data) <- c("tsneX", "tsneY")
dmrd.data <- cbind(dmrd.data, tsne.data)

ggplot(dmrd.data, aes(x = tsneX, y = tsneY, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("tsne_plot.jpg", width = 9, height = 6 )



# tSNE with pca initialization----------------------

set.seed( date.seed )

system.time(
  tsne.pca.init.result <- Rtsne( input.data, perplexity = 30, 
                                 exaggeration_factor = 4,
                                 eta = tsne.learning.rate,
                                 max_iter = 750, stop_lying_iter = 75,
                                 check_duplicates = FALSE, pca = TRUE, 
                                 num_threads = 0 )
)
# user  system elapsed 
# 1370.78    2.10  296.43

tsne.pca.init.data <- tsne.pca.init.result$Y

colnames(tsne.pca.init.data) <- c("tsne_PCA_X", "tsne_PCA_Y")
dmrd.data <- cbind(dmrd.data, tsne.pca.init.data)

ggplot(dmrd.data, aes(x = tsne_PCA_X, y = tsne_PCA_Y, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("tsne_PCA_plot.jpg", width = 9, height = 6 )


# umap via uwot----------------------

set.seed( date.seed )

system.time( 
  umap.result <- uwot::umap( input.data, n_neighbors = 30,
                                        n_epochs = 500, 
                                        n_threads = 0,
                                        n_sgd_threads = 0, 
                                        batch = TRUE, verbose = TRUE )
)
# user  system elapsed 
# 230.62    0.47  309.56 

colnames(umap.result) <- c("UMAP1", "UMAP2")
dmrd.data <- cbind(dmrd.data, umap.result)

ggplot(dmrd.data, aes(x = UMAP1, y = UMAP2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("umap_plot.jpg", width = 9, height = 6 )



# EmbedSOM----------------------

set.seed( date.seed )

system.time(
  flow.som <- EmbedSOM::SOM(input.data, xdim = 24, 
                            ydim = 24, batch = TRUE,
                            parallel = TRUE, threads = 0 )
)
# 4.68s
system.time(
  embed.som <- EmbedSOM::EmbedSOM( data = input.data, map = flow.som, parallel = T )
)
# 1.34s

dmrd.data <- cbind(dmrd.data, embed.som)

ggplot(dmrd.data, aes(x = EmbedSOM1, y = EmbedSOM2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("embedsom_plot.jpg", width = 9, height = 6 )



# embedsom with landmarks from tSNE (OptSNE)----------------------

set.seed( date.seed )

system.time(
  tsne.map <- EmbedSOM::RandomMap( input.data, 2000, 
                                   coords = EmbedSOM::tSNECoords(perplexity = 30, 
                                                                 check_duplicates = FALSE,
                                                                 pca = FALSE, 
                                                                 max_iter = 750, stop_lying_iter = 75,
                                                                 eta = 2000, exaggeration_factor = 4,
                                                                 num_threads = 0 ))
)
#1.68s

system.time(
  embed.tsne <- EmbedSOM( input.data, map = tsne.map,
                          parallel = T, threads = 0 )
)
#4s

colnames(embed.tsne) <- c("EmbedSNE1", "EmbedSNE2")

dmrd.data <- cbind(dmrd.data, embed.tsne)

ggplot(dmrd.data, aes(x = EmbedSNE1, y = EmbedSNE2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("embedsne_plot.jpg", width = 9, height = 6 )



# embedsom with landmarks from umap----------------------

set.seed( date.seed )

system.time(
  umap.map <- EmbedSOM::RandomMap( input.data, 2000, 
                                   coords = EmbedSOM::UMAPCoords() )
)
#6.45s

system.time(
  embed.map <- EmbedSOM( input.data, map = umap.map,
                         parallel = T, threads = 0 )
)
#4.12s

colnames(embed.map) <- c("EmbedMap1", "EmbedMap2")

dmrd.data <- cbind(dmrd.data, embed.map)

ggplot(dmrd.data, aes(x = EmbedMap1, y = EmbedMap2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("embedmap_plot.jpg", width = 9, height = 6 )


# embedsom using GQTSOM----------------------

set.seed( date.seed )

system.time(
  gqt.map <- EmbedSOM::GQTSOM(input.data, target_codes=1000, radius=c(10,.1), rlen=15, parallel=T)
)
#16.58s

system.time(
  embed.gqt <- EmbedSOM( input.data, map = gqt.map,
                         parallel = T, threads = 0 )
)
# 2.14s

colnames(embed.gqt) <- c("EmbedGQT1", "EmbedGQT2")

dmrd.data <- cbind(dmrd.data, embed.gqt)

ggplot(dmrd.data, aes(x = EmbedGQT1, y = EmbedGQT2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()


ggsave("embedgqt_plot.jpg", width = 9, height = 6 )



# phate----------------------

system.time(
  phate.object <- phate(input.data, n.jobs = -1, seed = date.seed)
)

# user  system elapsed 
# 2755.06    7.96  405.25

phate.data <- phate.object$embedding

dmrd.data <- cbind(dmrd.data, phate.data)

ggplot(dmrd.data, aes(x = PHATE1, y = PHATE2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("phate_plot.jpg", width = 9, height = 6 )



# densne----------------------

library(densvis)

set.seed( date.seed )

system.time(
  densne.data <- densne( input.data, perplexity = 30, 
                         exaggeration_factor = 4,
                         eta = tsne.learning.rate,
                         max_iter = 750, stop_lying_iter = 75,
                         check_duplicates = FALSE, pca = FALSE, 
                         dens_frac = 0.5, dens_lambda = 0.5,
                         num_threads = 0 )
)

# user  system elapsed 
# 2010.77    3.06 2193.42 

colnames(densne.data) <- c("denSNEX", "denSNEY")
dmrd.data <- cbind(dmrd.data, densne.data)

ggplot(dmrd.data, aes(x = denSNEX, y = denSNEY, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("densne_plot.jpg", width = 9, height = 6 )



# densmap----------------------

set.seed( date.seed )

system.time(
  densmap.data <- densmap( input.data, 
                           dens_frac = 0.5, dens_lambda = 0.5,
                           n_neighbors = 30L,
                           n_epochs = 500L )
)
# user  system elapsed 
# 0.13    0.00  185.12 

colnames(densmap.data) <- c("densMAP1", "densMAP2")
dmrd.data <- cbind(dmrd.data, densmap.data)

ggplot(dmrd.data, aes(x = densMAP1, y = densMAP2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()

ggsave("densmap_plot.jpg", width = 9, height = 6 )



# run pacmap----------------------

virtualenv_create("dr-comparison")
use_virtualenv("dr-comparison")

py_install("pandas")
#py_install("numpy")
#py_install("annoy")
py_install("pacmap")

python_pandas <- import("pandas")
python_pacmap <- import("pacmap")
python_numpy <- import("numpy")

py_input <- reticulate::r_to_py(input.data)

pacmap.embedding <- python_pacmap$PaCMAP(n_components=2L, n_neighbors=NULL, MN_ratio=0.5, FP_ratio=2.0)

system.time(
  X_transformed <- pacmap.embedding$fit_transform(py_input, init="pca")
)
# user  system elapsed 
# 230.42  129.63  115.49 

pacmap.data <- data.frame(X_transformed)

colnames(pacmap.data) <- c("PaCMAP1", "PaCMAP2")

dmrd.data <- cbind(dmrd.data, pacmap.data)

ggplot(dmrd.data, aes(x = PaCMAP1, y = PaCMAP2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()


ggsave("pacmap_plot.jpg", width = 9, height = 6 )



# run trimap----------------------

py_install("trimap")
python_trimap <- import("trimap")

nparray <- py_input$values

system.time(
  trimap.embedding <- python_trimap$TRIMAP()$fit_transform(nparray)
)
# user  system elapsed 
# 207.06  105.76  180.94

trimap.data <- data.frame(trimap.embedding)

colnames(trimap.data) <- c("TriMap1", "TriMap2")

dmrd.data <- cbind(dmrd.data, trimap.data)

ggplot(dmrd.data, aes(x = TriMap1, y = TriMap2, color = Color_clust)) +
  geom_scattermore() + 
  scale_colour_identity("Cluster", breaks=dmrd.data$Color_clust, 
                        labels=dmrd.data$Cluster_name,
                        guide = "legend") +
  theme_classic()


ggsave("trimap_plot.jpg", width = 9, height = 6 )


# save
save(dmrd.data, file = "dr_comparison_data.RData")


