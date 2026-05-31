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

          # reshape the data
          param_table %>%
               filter(name %in% c("VA11", "VC11", "VE11")) %>%
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

#-------I can't remmeber if this stuff at the bottom was a relic of testing or if it is relevant somehow ---------------

table_ord <- harvest(results_fit, "Ordinal")
table_int <- harvest(results_fit, "Interval")

final_results <- rbind(table_ord, table_int)
final_results <- final_results %>% mutate(TV = VA11 + VC11 + VE11,
                                          A = VA11 / TV,
                                          C = VC11 / TV,
                                          E = VE11 / TV) %>%
     select(-c("VA11","VC11","VE11","TV"))
final_results

