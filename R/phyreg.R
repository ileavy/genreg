#' @title phyreg: run a single regression on a nexus form
#' @description
#' @param tree a nexus data set
#' @param df dataset with column 1
#' @param pvar predictor variable(s), in the standard format for linear models if multiple
#' @param ovar outcome variable(s)
#' @param
#' @param
#' @keywords
#' @export
#' @examples
#' \dontrun{
#' tree <- read.nexus("gennex_output.nex") #needs to be in class phylo
#' df <- read.csv("VarBySpecies.csv")
#' phyreg(tree = tree, df = df, pvar = "CPP", ovar = "death.prop")
#' }


### NOTE: requires caper, nlme, gt,  modelsummary packages
phyreg <- function(tree = tree, df = df, pvar = pvar, ovar = ovar) {

  #cleaning the data frame
  df <- as.data.frame(df)
  rownames(df) <- NULL

  #create comparative data set
  comp_data <- comparative.data(phy = tree, data = df, names.col = "species", vcv = TRUE)

  #formulas so variables input as strings can be used
  lm_formula   <- as.formula(paste(ovar, "~", pvar))
  pgls_formula <- as.formula(paste(ovar, "~", pvar))

  #models
  lm_model   <- lm(lm_formula, df)
  pgls_model <- pgls(pgls_formula, comp_data, lambda = "ML")

  #pagel's lambda calculation/pull
  pgls_sum   <- summary(pgls_model)
  lambda_val <- pgls_model$param["lambda"]
  lambda_ci  <- pgls_model$param.CI$lambda
  lambda_p   <- pgls_sum$param.CI$lambda$opt

  # Tidy the models into data frames we can bind for comparison
  lm.df <- tidy(lm_model, conf.int = TRUE) |>
    mutate(model = "Linear")

  pgls_coef <- as.data.frame(pgls_sum$coefficients)
  pgls.df <- data.frame(
    term      = rownames(pgls_coef),
    estimate  = pgls_coef[, "Estimate"],
    std.error = pgls_coef[, "Std. Error"],
    statistic = pgls_coef[, "t value"],
    p.value   = pgls_coef[, "Pr(>|t|)"],
    model     = "Phylogenetic",
    conf.low  = pgls_coef[, "Estimate"] - 1.96 * pgls_coef[, "Std. Error"],
    conf.high = pgls_coef[, "Estimate"] + 1.96 * pgls_coef[, "Std. Error"]
  )

  lambda.df <- data.frame(
    model     = "Pagel's Lambda",
    term      = "lambda",
    estimate  = round(lambda_val, 3),
    std.error = NA,
    statistic = NA,
    p.value   = round(lambda_p, 3)
  )

#bind the models for a return
 compare_chart <- bind_rows(lm.df, pgls.df, lambda.df) |>
    dplyr::select(model, term, estimate, std.error, statistic, p.value) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    gt(groupname_col = "model") |>
    cols_label(
      term      = "Term",
      estimate  = "Estimate",
      std.error = "SE",
      statistic = "t",
      p.value   = "p"
    ) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) |>
   tab_spanner(label = "Model comparison", columns = everything()) |>
   fmt_number(decimals = 3) |>
   fmt_missing(columns = everything(), missing_text = "—") |>  # display NA as dash
   tab_options(table.width = pct(100))

 return(compare_chart)

}
