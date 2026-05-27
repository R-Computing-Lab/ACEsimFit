#' Harvest results from MCMCglmm analyses
#' @description A function to extract the A, C, E estimates from the model
#' fitting results of the MCMCglmm analyses. This function is designed to work with the output of \code{Sim_Fit2}, which includes both interval and ordinal model fitting results.
#' @param results_fit A list of model fitting results generated from \code{Sim_Fit2}. This list should have two sub-lists: one for "Interval" model results and one for "Ordinal" model results. Each sub-list should contain the fitting results for each iteration (e.g., "Iteration1", "Iteration2", etc.).
#' @param type A character string specifying the type of model results to extract. Should be
#' either "Interval" for the interval model results or "Ordinal" for the ordinal model results.
#' @return Returns a \code{data.frame} with the following columns:
#' \item{Iteration}{The iteration name (e.g., "Iteration1", "
#' Iteration2", etc.)}
#' \item{Analysis}{The type of analysis ("Interval" or "Ordinal")}
#' \item{VA11}{The estimated additive genetic variance component (A) for the
#' first group of kin pairs}
#' \item{VC11}{The estimated common environmental variance component (C) for the
#' first group of kin pairs}
#' \item{VE11}{The estimated unique environmental variance component (E) for the
#' first group of kin pairs}
#' @export
#' @importFrom dplyr filter select relocate
#' @importFrom tidyr pivot_wider
#' @importFrom purrr imap_dfr

harvest <- function(results_fit, type) {

     # 1. Grab the correct sub-list based on type ("Interval" or "Ordinal")
     target_list <- results_fit[[type]]

     # 2. Loop through every iteration using imap (which tracks the iteration name)
     master_df <- purr::imap_dfr(target_list, function(iteration_data, iter_name) {

          # pull the parameters
          param_table <- as.data.frame(iteration_data[["Results"]][["summary"]][["parameters"]])

          # reshape the data
          param_table |>
                  dplyr::filter(name %in% c("VA11", "VC11", "VE11")) |>
               dplyr::select(name, Estimate) |>
               tidyr::pivot_wider(
                    names_from = name,
                    values_from = Estimate
               ) |>
                 dplyr::mutate(
                    Analysis = type,
                    Iteration = iter_name # Captures "Iteration1", "Iteration2", etc.
               )
     })

     master_df <- master_df |>
             dplyr::relocate(Iteration, Analysis)

     return(master_df)
}
