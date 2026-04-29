#' @title gennex: from genetic data to nexus phylogeny
#' @description This function pulls mitochondrial DNA from Gen Bank and creates a phylogenetic map
#' @param data a list of accession numbers for mitogenomes in GenBank (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)
#' @param bs how many bootstrap samples do you want it to run
#' @param out one of the species names in the "data" list to be used as the outgroup
#' @param
#' @param
#' @param
#' @keywords phylogeny, mitodna, genbank
#' @export
#' @examples gennex(data = c(Canis_lupus_familiaris = "NC_002008", Lemur_catta = "NC_004025",Galago_moholi = "KC757396", Saimiri_sciureus = "FJ785425", Callithrix_jacchus = "NC_025586"), bs = 100, out = "Canis_lupus_familiaris")

gennex <- function(data = data , bs = bs, out = out) {
  #pulling DNA sequence data from GenBank
  genbank <- read.GenBank(data, species.names = FALSE, as.character = TRUE)

  #create a single DNA sequence for each species
  genbank <- sapply(genbank, paste0, collapse = "")

  #
  dna <- DNAStringSet(genbank)
  names(dna) <- names(data)
  #dna <- call outcome

  #performing a multisequence alignment using {DECIPHER}
  aligned_dna <- AlignSeqs(dna, processors = 1, verbose = FALSE)
  #aligned_dna
  #length(aligned_dna) # 20 -- it's good
  #width(aligned_dna)

  #convert alignment to characters
  aligned_dna <- as.character(aligned_dna)
  #match names
  names(aligned_dna) <- names(data)

  #wrangle alignment into a phylogenetic format
  dnabin <- as.DNAbin(strsplit(aligned_dna, split = ""))

  #use phanghorn::phyDat to create a maximum likelihood phylogeny
  phydat <- phyDat(dnabin, type = "DNA")

  #calculate genetic distance metric
  gen_dist <- dist.ml(phydat)

  #create a starting tree to build the model off of, use Neighbor-Joining
  start_tree <- NJ(gen_dist)

  #combine tree data with aligned DNA using {phanogram}
  pml <- pml(start_tree, data = phydat, k = 4) # divides to 4 rates, default in phylogenetics
  #optimize the tree along maximum likelihood
  optim_pml <- optim.pml(
    pml,
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = pml.control(trace = 0)
  )

  #optim_pml$tree

  tree_rooted <- root(
    optim_pml$tree,
    outgroup = out,
    resolve.root = TRUE
  )

  boot <- bootstrap.pml(
    optim_pml,
    bs = bs,
    optNni = TRUE,
    multicore = FALSE,
    control = pml.control(trace = 0)
  )

  plotBS(tree_rooted, boot, p = 50, main = "Mitogenome ML phylogeny", cex = 0.7, font = 3, label.offset = 0.01)

  nex <- write.nexus(tree_rooted, file = "my_gennex.nex")
  nex <- read.nexus("my_gennex.nex")
  summary(nex)
}
