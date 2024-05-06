library(plotrix)

set.seed(2440)
x <- runif(1000, 0, 10)
y <- runif(1000, 0, 10)
grp = sample(size=1000, c("red", "blue", "green", "orange"), replace=T)
grp[3] <- "green"
plot(x, y, col=grp)
draw.circle(3, 5, radius=1)
mtext("Circle of radius 1 centered at (3,5) \n 9 points of each color located inside circle")

in.circle <- (x-3)^2 + (y-5)^2 <= 1
sum(in.circle)
sum(in.circle[grp=="red"])
sum(in.circle[grp=="blue"])
sum(in.circle[grp=="green"])
sum(in.circle[grp=="orange"])

