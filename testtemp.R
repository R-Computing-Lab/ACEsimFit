library(SimFit)
test <- kinsim_double(
     GroupNames = c("tiger","lion"),
     GroupSizes = c(100,300),
     GroupRel = c(.8,.6),
     GroupR_c = c(.99,.95),
     ifComb = TRUE)
tiger <- test[which(test$GroupName=="tiger"), c("y1","y2")]
lion <- test[which(test$GroupName=="lion"), c("y1","y2")]

testResult <- fit_uniACE(tiger, lion, GroupRel = c(.8,.6),GroupR_c = c(.99,.95))
