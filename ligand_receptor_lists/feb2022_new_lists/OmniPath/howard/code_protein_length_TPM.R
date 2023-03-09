
# Downloaded uniprot data from here: Programmatic access - Retrieving entries via queries: https://www.uniprot.org/help/api_queries
uniprot <- read_excel("input/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2022.11.28-17.35.59.13.xlsx")
View(uniprot)
uniprot.size <- uniprot %>%
  dplyr::select(Entry,`Entry Name`,`Protein names`, `Gene Names`, Length, Mass) %>%
  #separate(`Entry Name`, into = c("Gene",NA)) %>%
  separate(`Gene Names`,into = c("Gene1","Gene2","Gene3","Gene4","Gene5",
                                 "Gene6","Gene7","Gene8","Gene9","Gene10",
                                 "Gene11","Gene12","Gene13","Gene14","Gene15",
                                 "Gene16","Gene17","Gene18","Gene19","Gene20")) %>%
  pivot_longer(cols = Gene1:Gene20, names_to = "name", values_to = "Genes") %>%
  filter(!is.na(Genes)) %>%
  dplyr::select(-name)
View(uniprot.size)

write.table(uniprot.size, file = "output/uniprot_size.txt",
            sep = "\t", row.names = F)
uniprot.size <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/uniprot_size.txt",
                            sep = "\t", header = T)
View(uniprot.size)

#(if the protein abundance is in the object "pro", and mean protein abundance is in column "mean"):
pro <- pro %>%
  left_join(uniprot.size,by=c("Gene"="Genes")) %>%
  mutate(mean.signal.length = mean/Length)




## Note: For my correlation purpuse, only proteins that are detected in >=90% samples are kept.

# This is the data with protein abundance normalized to length in column "mean.signal.length"
pro.islet.norm <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/protein_islets_normLength.txt",
                           sep = "\t", header = T)
View(pro.islet.norm)

# This is the same protein data but only contains the 90 common donors in RNAseq
pro.islet.norm.oneBatch <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/protein_islets_normLength_oneBatch.txt",
                             sep = "\t", header = T)
View(pro.islet.norm.oneBatch) 


# This is the RNA TPM data. Mean TPM is in the "mean.TPM" column.
# Note: Only genes that are detected with >=5 raw counts in >=90% samples are kept (and considered to be reliably detected).
rna.islet.tpm.oneBatch <- read.table(file = "ligand_receptor_lists/feb2022_new_lists/OmniPath/howard/RNA.islets.TPM.oneBatch.txt",
                                     sep = "\t", header = T)
View(rna.islet.tpm.oneBatch)


# my plot

## (optional) Some of my plots have SD or CoV as error bars for proteins and RNAs, which could be informative for your purpose as well.

# example code to calculate CoV or SD
df$protein.CoV <- apply(df[,"sample columns here"], 1, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) * 100)
df$protein.SD <- apply(df[,"sample columns here"], 1, function(x) sd(x, na.rm=T))
df <- df %>% 
  mutate(xmin=mean.TPM,
         xmax=mean.TPM+RNA.TPM.SD,
         ymin=mean.signal.length,
         ymax=mean.signal.length+protein.SD)

# Btw, pro.rna.all.gene.oneBatch is the protein and RNA data joined together, not one of the data frame I provided.
ggplot(pro.rna.all.gene.oneBatch,
       aes(
         y=mean.signal.length, 
         x=mean.TPM
         )) +
  
  geom_point(
    size=0.8,
    alpha= 0.4 
  ) +

# For trend line.
  geom_smooth(
    #se=FALSE,
    method=lm
  ) +

# (opitonal) plot error bars
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) + 
  geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +

# (optional) show Pearson R on the plot  
  stat_cor(method = "pearson" #, label.x = 3, label.y = 30
  )+
  
  labs(x= "RNA abundance (mean TPM)", 
       y= "Protein abundance (signal/length)"
  ) +
  
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "bl") +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(colour = "black",size=14),
    axis.title.x = element_text(colour = "black",size=14),
    axis.text.y = element_text(colour = "black",size=12),
    axis.text.x = element_text(colour = "black",size=12),
    #axis.ticks.x = element_blank(),
    #legend.title = element_blank(),
    legend.position="none",
    legend.title=element_text(size=14), 
    legend.text=element_text(size=14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio=1/1) +
  theme(axis.line = element_line()) #+
#  guides(colour = guide_legend(override.aes = list(size=2,alpha=1)))
