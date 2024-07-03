library(gggenomes)

s0 <- read_seqs("xylella_prophages.fasta")
g0 <- read_feats("xylella_prophages.gff") 
l0 <- read_links("xylella_prophages.o6")

colnames(l0)[4] <- "Similarity"

pal = c("lightskyblue1", "lightpink1", "lemonchiffon", "lightsalmon1", "mistyrose2", "thistle", "lightsteelblue3", "indianred1", "bisque3", "grey")

Functional_Modules <- g0$phase
figura <- gggenomes(seqs=s0, links=l0, genes = g0) + geom_seq() + geom_link(aes(alpha=Similarity)) + geom_gene(aes(fill=Functional_Modules)) + 
  geom_bin_label(size=4, nudge_left=0.01, expand_left =0.12 )  + 
  theme(legend.key.size = unit(7, 'mm'), legend.text = element_text(size=10), legend.title = element_text(size=12,face = "bold"), 
        legend.position = "top", legend.spacing.x = unit(0.21,'cm'),legend.box.just = "right",legend.box.margin=margin(-10,0,-20,-5)) + 
  guides(fill = guide_legend(order = 2)) + scale_fill_manual(values=pal, breaks=c('DNA, RNA and nucleotide metabolism', 
                                                                                  'head and packaging',
                                                                                  'integration and excision',
                                                                                  'lysis',
                                                                                  'Moron, auxiliary metabolic gene and host takeover',
                                                                                  'tail',
                                                                                  'rranscription regulation',
                                                                                  'connector',
                                                                                  'other',
                                                                                  'unknown function'))

ggsave("gggenomes.svg", plot =figura, width = 15, height = 5, bg = "transparent")

