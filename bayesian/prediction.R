library(spatstat)

woak <- lansing[lansing$marks == "whiteoak"]
x.lower <- 0.1
x.upper <- 0.25
y.lower <- 0.5
y.upper <- 0.7
true.count <- sum(woak$x > x.lower &
                    woak$x < x.upper &
                    woak$y > y.lower &
                    woak$y < y.upper)
vol <- (x.upper - x.lower) * 10 * (y.upper - y.lower) * 10
index <- (int_df$x > x.lower * 10 & 
            int_df$x < x.upper * 10 & 
            int_df$y > y.lower * 10 & 
            int_df$y < y.upper * 10)

pred.count <- mean(int_df$intesity_function[index]) * vol

cat("True count is ", true.count, " and predicted count is ", round(pred.count, 2), ".\n")
