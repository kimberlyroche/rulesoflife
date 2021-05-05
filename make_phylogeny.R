library(phangorn)
library(DECIPHER)
library(phyloseq)
library(dada2)

ps <- readRDS(file.path("input", "ps0.rds"))
ftbl <- as(otu_table(ps), "matrix")

seqs <- getSequences(ftbl)
names(seqs) <- seqs # This propagates to the tip labels of the tree

seqs <- seqs[1:100]
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)

phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)

# Two different tree building algorithms—any of these can serve as the starting point for the ML-GTR tree
treeNJ <- NJ(dm) # Note: tip order != sequence order
treeUPGMA <- upgma(dm)

fitNJ = pml(treeNJ, data = phang.align)
fitUPGMA <- pml(treeUPGMA, data = phang.align)

fitNJ_GTR <- update(fitNJ, k=4, inv=0.2)
fitNJ_GTR <- optim.pml(fitNJ_GTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))

is.rooted(fitNJ_GTR$tree) # FALSE
# saveRDS(fitNJ_GTR, "ps_NJtree.RDS")

fitUPGMA_GTR <- update(treeUPGMA, k=4, inv=0.2)
fitUPGMA_GTR <- optim.pml(fitUPGMA_GTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                          rearrangement = "stochastic", control = pml.control(trace = 0))

is.rooted(fitUPGMA_GTR$tree) # FALSE
# saveRDS(fitUPGMA_GTR, "ps_UPGMAtree.RDS")

# ps <- phyloseq(tax_table(ps), sample_data(ps),
#               otu_table(ps), phy_tree(fitGTR$tree))

# is.rooted(phy_tree(ps)) # FALSE

## Root tree (Philr assumes a rooted tree)

# Rooting this way, you select an outgroup taxa (here we choose an Archaea)
phy_tree(ps) <- ape::root(tree_family, taxa_names(ps)[323], resolve.root=TRUE)

# Midpoint rooting in the absence of specific information
phy_tree(ps) <- phangorn::midpoint(phy_tree(ps))
phy_tree(ps) <- ape::multi2di(phy_tree(ps))

is.rooted(phy_tree(ps)) # TRUE
is.binary.tree(phy_tree(ps)) # TRUE

# We name the internal nodes of the tree so they are easier to work with. We prefix the node number with n and thus the root is named n1.
phy_tree(ps) <- ape::makeNodeLabel(phy_tree(ps), method="number", prefix='n')
