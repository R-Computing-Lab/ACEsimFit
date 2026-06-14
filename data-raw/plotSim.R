#making a function for plotting the chart in the way I want it specifically for this simulation
#calling the function plotSim


#its going to take results_fit list from SimFit2
make_parameter_df <- function(results_fit, prop_missing, GroupRel, true_A, true_C, true_E) {

     table_ord <- data.frame()
     table_int <- data.frame()
     final_results <- data.frame()

#use harvest functions to put the parameters into a data frame - repeat twice for ord and int results, and then rbind it

table_ord <- harvest(results_fit, "Ordinal") %>% mutate(Analysis = "Ordinal")
table_int <- harvest(results_fit, "Interval") %>% mutate(Analysis = "Interval")
final_results <- rbind(table_ord, table_int)

#next add in the mutation to get total variance and standardization, censorship columns, and a column with the relatedness values

final_results <- final_results %>% mutate(TV = VA11 + VC11 + VE11,
                                            A = VA11 / TV,
                                            C = VC11 / TV,
                                            E = VE11 / TV,

                                            censored1 = prop_missing[1], #prop missing for the first kin type, where prop_missing is a vector of two values
                                          censored2 = prop_missing[2],

                                          R1 = GroupRel[1], #tells us the relatedness for pair one, where group rel is a vector with two values
                                          R2 = GroupRel[2], #tells us relatedness of pair two

                                          Bias_A = true_A - A, #tells us how far off A was from the true A
                                          Bias_C = true_C - C, #tells us how far off from the true C
                                          Bias_E = true_E - E #tells us how far off from true E
                                        ) %>%

     #this is the part AI thinks I should add on to make my function better - will see if thats actually the way to go or not
     #this will give me a different plot option - I can graph the bias as a function of missingnes (for example) for the two analysis types
     #AI was correct - this addition makes it much better and now gives me more info I need for more interesting figures

     # 3. COLLAPSE THE 10,000 ITERATIONS INTO MEAN + SE
     # Group by your X-axis variable and your line color variable

      group_by(Analysis, censored1, censored2, R1, R2) %>%
     mutate(
          # Mean bias (this forms the continuous line)
          mean_bias_A = mean(Bias_A, na.rm = TRUE),
          mean_bias_C = mean(Bias_C, na.rm = TRUE),
          mean_bias_E = mean(Bias_E, na.rm = TRUE),

          # Standard Error (this forms the shaded confidence band)
          se_bias_A = sd(Bias_A, na.rm = TRUE) / sqrt(n()),
          se_bias_C = sd(Bias_C, na.rm = TRUE) / sqrt(n()),
          se_bias_E = sd(Bias_E, na.rm = TRUE) / sqrt(n())) %>%

ungroup()

return(final_results) #return this data frame and then use it to plot

}


#I am not sure if I will do anything with this - I don't think this particular graph is all that informative for the poster given it can only show one condition at a time

#change this so true_A, true_C, etc are flexible
# true_A <- 0.50
# true_C <- 0.20
# true_E <- 0.30
#
# #plot <-
# ggplot(final_results2) +
#      geom_density(aes(x = A, fill = "Additive Genetic (A)"), alpha = 0.4, color = "#9E7E38") +
#      geom_density(aes(x = C, fill = "Shared Environment (C)"), alpha = 0.4, color = "#D4AF37") +
#      geom_density(aes(x = E, fill = "Error (E)"), alpha = 0.4, color = "#4A4E41") +
#
#      geom_vline(xintercept = true_A, color = "#9E7E38", linetype = "dashed", size = 1) +
#      geom_vline(xintercept = true_C, color = "#D4AF37", linetype = "dashed", size = 1) +
#      geom_vline(xintercept = true_E, color = "#4A4E41", linetype = "dashed", size = 1) +
#      scale_fill_manual(
#           values = c(
#                "Additive Genetic (A)" = "#9E7E38",
#                "Shared Environment (C)" = "#D4AF37",
#                "Error (E)" = "#4A4E41")) +
#
#      facet_wrap(~Analysis) +
#
#      theme_minimal() +
#
#      theme(axis.title.x = element_blank(),
#            text = element_text(family = "Baskerville", color = "black"),#not to self turn the text to white and then its ready for poster (already has transparent background)
#            panel.background = element_rect(fill = "transparent", color = NA),
#            plot.background = element_rect(fill = "transparent", color = NA),
#            panel.grid = element_blank(),
#            axis.ticks.y = element_blank()) +
#
#      labs(title = "Distribution of Simulation Results",
#           subtitle = "MZs and DZs",
#           y = "Density",
#           fill = "Parameter Estimate")
