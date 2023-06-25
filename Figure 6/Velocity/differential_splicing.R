diff_kinetics <- read.csv(file = "diff_kinetics.csv")

top_genes_diff_splice <- read.csv(file = "top_genes_diff_splice.csv")

GOI <- head(top_genes_diff_splice,10000)

GOI <- right_join(GOI,diff_kinetics, by = "X")

GOI <- head(GOI,10000)

write.csv(GOI, "diff_spliced_genes.csv")

GOI <- GOI[GOI$fit_diff_kinetics %in% c("core","edge","transitory"),]

GOI <- head(GOI,15)

png("diff_spliced_genes.png", units="in", width=6, height=4, res=300, type = "cairo")
ggplot(GOI, aes(x = factor(X, level = c(X)), 
                   y = fit_likelihood, fill = fit_diff_kinetics)) + 
  geom_bar(stat="identity", color = "black", size = 0.25) + theme_minimal() + 
  labs(x = "", y = "Splicing fit likelihood", fill = "Sample") +
  scale_fill_manual( values = c("core" = "#4DBBD5FF", "edge" = "#E64B35FF","transitory" = "#F9E076"))+
  ggpubr::rotate_x_text()+ theme(axis.text = element_text(face="bold", size = 11), legend.text = element_text(face="bold"), 
                                 axis.text.y = element_text(face="bold"), axis.title.y = element_text(face="bold")) + ggpubr::rotate_y_text()
dev.off()
