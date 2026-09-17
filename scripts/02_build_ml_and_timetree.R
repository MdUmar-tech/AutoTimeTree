#!/usr/bin/env Rscript
# AutoTimeTree - Step 2: Build ML Tree (FastTree), Root with Outgroup, and Compute Time Tree (ape::chronos)

library(ape)

args <- commandArgs(trailingOnly = TRUE)

align_file <- ifelse(length(args) >= 1, args[1], "results/aligned_sequences.fasta")
outgroup_pattern <- ifelse(length(args) >= 2, args[2], "OUTGROUP")
ml_tree_file <- ifelse(length(args) >= 3, args[3], "results/ml_tree.nwk")
rooted_tree_file <- ifelse(length(args) >= 4, args[4], "results/rooted_ml_tree.nwk")
timetree_file <- ifelse(length(args) >= 5, args[5], "results/timetree.nwk")

cat(sprintf("[*] Alignment File: %s\n", align_file))
cat(sprintf("[*] Outgroup Pattern: %s\n", outgroup_pattern))

if (!file.exists(align_file)) {
  stop(sprintf("[Error] Alignment file '%s' does not exist!", align_file))
}

# Ensure directory exists for outputs
out_dir <- dirname(rooted_tree_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 1. Run FastTree (Maximum Likelihood under GTR+CAT model)
cat("[*] Building Maximum Likelihood tree using FastTree (GTR+CAT)...\n")
cmd <- sprintf("fasttree -gtr -nt '%s' > '%s'", align_file, ml_tree_file)
system(cmd)

if (!file.exists(ml_tree_file)) {
  stop("[Error] FastTree failed to generate ML tree file.")
}

# 2. Load ML tree and root with Outgroup
cat("[*] Loading tree and rooting with outgroup...\n")
tree <- read.tree(ml_tree_file)

# Find outgroup tip matching pattern
outgroup_tips <- grep(outgroup_pattern, tree[["tip.label"]], value = TRUE, ignore.case = TRUE)

if (length(outgroup_tips) > 0) {
  cat(sprintf("[+] Found outgroup tip: '%s'. Rooting tree...\n", outgroup_tips[1]))
  tree <- root(tree, outgroup = outgroup_tips[1], resolve.root = TRUE)
} else {
  cat("[!] Warning: Specified outgroup pattern not found in tip labels! Proceeding unrooted or mid-point rooted.\n")
}

write.tree(tree, rooted_tree_file)
cat(sprintf("[+] Saved rooted ML tree to '%s'.\n", rooted_tree_file))

# 3. Compute Time Tree (Chronogram) using ape::chronos
cat("[*] Computing Time Tree (chronogram) using ape::chronos (relaxed molecular clock)...\n")
tree[["node.label"]] <- NULL
# Avoid zero branch lengths
tree[["edge.length"]][tree[["edge.length"]] < 1e-6] <- 1e-6

time_tree <- chronos(tree, lambda = 1, model = "relaxed")
write.tree(time_tree, timetree_file)
cat(sprintf("[+] Saved Time Tree (chronogram) to '%s'.\n", timetree_file))
