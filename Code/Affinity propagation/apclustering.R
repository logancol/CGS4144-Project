library(apcluster)

ap_result2_10 <- apcluster(negDistMat(r = 2), top_5000_matrix_numeric[1:10, ], 
                           maxits = 1000, convits = 100, lam = 0.9, q = 0, details = TRUE)
summary(ap_result_10)
png("affinity_propagation_plot_10_genes.png", width = 1200, height = 800)
plot(ap_result_10)
dev.off()

ap_result_100 <- apcluster(negDistMat(r = 2), top_5000_matrix_numeric[1:100, ], 
                           maxits = 1000, convits = 100, lam = 0.9, q = 0, details = TRUE)
summary(ap_result_100)
png("affinity_propagation_plot_100_genes.png", width = 1200, height = 800)
plot(ap_result_100)
dev.off()

ap_result_1000 <- apcluster(negDistMat(r = 2), top_5000_matrix_numeric[1:1000, ], 
                            maxits = 1000, convits = 100, lam = 0.9, q = 0, details = TRUE)
summary(ap_result_1000)
png("affinity_propagation_plot_1000_genes.png", width = 1200, height = 800)
plot(ap_result_1000)
dev.off()

ap_result_5000 <- apcluster(negDistMat(r = 2), top_5000_matrix_numeric[1:5000, ], 
                            maxits = 1000, convits = 100, lam = 0.9, q = 0, details = TRUE)
summary(ap_result_5000)
png("affinity_propagation_plot_5000_genes.png", width = 1200, height = 800)
plot(ap_result_5000)
dev.off()
