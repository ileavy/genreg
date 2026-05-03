#' @title nexsimple: a specific workflow for getting a nexus phylogeny to look how i want it
#' @description this function starts with a nexus file and a data frame of variables organized by species, and prunes both to match for easier phylogenetic analysis. It also returns some nice visuals.
#' @param nex a nexus data set
#' @param df a data frame with a column that contains "Genus_species" for multiple species
#' @param var the variable of interest, usually the predictor, a column within df in quotes
#' @param
#' @param
#' @param
#' @keywords phylogeny, nexus, visual
#' @export
#' @examples
#' \dontrun{
#' nex <- read.nexus("gennex_output.nex")
#' df  <- read.csv("VarBySpecies.csv")
#' nexsimple(nex = nex, df = df, var = "CPP")
#' }

nexsimple <- function(nex, df, var) {
  tree <- nex
  df_clean <- df[!is.na(df[[var]]), ]
  tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% df_clean$species])
  df_final <- df_clean[df_clean$species %in% tree_pruned$tip.label, ]
  return(list(tree = tree_pruned, df = df_final))
}



