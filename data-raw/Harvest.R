#takes the results list and type of analysis and makes a df with the iteration number, analysis name, and A,C,E

harvest <- function(results_fit, type) {


     # 1. Grab the correct sub-list based on type ("Interval" or "Ordinal")
     target_list <- results_fit[[type]]

     # 2. Loop through every iteration using imap (which tracks the iteration name)
     master_df <- imap_dfr(target_list, function(iteration_data, iter_name) {

          if (length(iteration_data) == 1 && is.na(iteration_data)) { #skips the iterations that fail
               return(NULL) #might need to be NA - nope NULL gives me what I want
               #also I learned the && trick - that means if the length isn't 1, it doesn't even check the second part of the statement, it just moves along
          }

          # pull the parameters
          param_table <- as.data.frame(iteration_data[["Results"]][["summary"]][["parameters"]])

          #send it to the env so we can see what we are working with
          assign("testtable", param_table, envir = .GlobalEnv)

          # reshape the data
          param_table <- param_table %>%
               dplyr::filter(name %in% c("VA11", "VC11", "VE11")) %>%
               select(name, Estimate) %>%
               pivot_wider(
                    names_from = name,
                    values_from = Estimate
               ) %>%
               mutate(
                    Analysis = type,
                    Iteration = iter_name # Captures "Iteration1", "Iteration2", etc.
               )
     })


     master_df <- master_df %>%
          relocate(Iteration, Analysis)

     return(master_df)
}


