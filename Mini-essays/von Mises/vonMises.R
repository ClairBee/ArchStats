
setwd("~/Desktop/Dissertation/von Mises")

library(circular)

#===============================================================================
# calculate proportion of vM distribution within certain distance

par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), mgp=c(2.2,0.45,0))
plot(circular(x, rotation = "clock"), dvonmises(x, mu = 0, kappa = 0), lty = l,
     shrink = 2, zero = pi/2, ylim = c(-1,2))

# k = 4: 99% of dist lies within pi/2 radius
points(qvonmises(0.005, mu = circular(0), kappa = 4), zero = pi/2, col = "darkred", pch = 20)
points(qvonmises(0.995, mu = circular(0), kappa = 4), zero = pi/2, col = "darkred", pch = 20)

round(suppressWarnings(pvonmises(pi/2, mu = 0, kappa = 4) - pvonmises(-pi/2, mu = 0, kappa = 4))*100,2)
round(suppressWarnings(pvonmises(pi/2, mu = 0, kappa = 1) - pvonmises(-pi/2, mu = 0, kappa = 1))*100,2)
round(suppressWarnings(pvonmises(pi/2, mu = 0, kappa = 0) - pvonmises(-pi/2, mu = 0, kappa = 0))*100,2)

# k = 1: 78% of dist lies within pi/2 radius
points(qvonmises(0.1, mu = circular(0), kappa = 1), zero = pi/2, col = "darkblue", pch = 20)
points(qvonmises(0.9, mu = circular(0), kappa = 1), zero = pi/2, col = "darkblue", pch = 20)


#===============================================================================
# plot vM distributions with various values of kappa
x <- c(-(101*pi):(101 * pi))/100
breaks_pi <- pretty(range(x/pi), n = 3)

pdf(file = "vM-linear-plot.pdf", height = 5)
# linear plot
par(oma = c(.5, 0, 0, 0), mar = c(1, 0, 0, 0), mgp=c(2.2,0.45,0))
plot(x, dvonmises(x, mu = 0, kappa = 0), type = "l", lty = 3, ylim = c(0,1.3),
     yaxt = "none", xaxt = "none", ylab = "", frame.plot = F, lwd = 2)
axis(1, at = c(-1,-0.5,0,0.5,1) * pi,
     labels = c(expression(paste("-", pi)), expression(paste("-", pi, "/2")),
                "0", expression(paste(pi, "/2")), expression(paste(pi))))
lines(x, dvonmises(x, mu = 0, kappa = 1), lty = 2, col = "darkblue", lwd = 2)
lines(x, dvonmises(x, mu = 0, kappa = 4), lty = 5, col = "darkred", lwd = 2)
lines(x, dvonmises(x, mu = 0, kappa = 10), lty = 4, col = "darkgreen", lwd = 2)
legend(1.5, 1.3, title = expression(kappa), lwd = 2,
       legend = c(0, 1, 4, 10),
       lty = c(3, 2, 5, 4), 
       col = c("black", "darkblue", "darkred", "darkgreen"), 
       bty = "n")
dev.off()

# CREATE A TEMPLATE FOR CIRCULAR PLOTS TO ENSURE CONSISTENCY
pdf(file = "vM-circular-plot.pdf")
# circular plot
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), mgp=c(2.2,0.45,0))
plot(circular(x, rotation = "clock"), dvonmises(x, mu = 0, kappa = 0), lty = l,
     shrink = 2, zero = pi/2, ylim = c(-1,2))
lines(circular(x), dvonmises(x, mu = 0, kappa = 0), lty = 3, lwd = 2, zero = pi/2)
lines(circular(x), dvonmises(x, mu = 0, kappa = 1), lty = 2, lwd = 2, col = "darkblue", zero = pi/2)
lines(circular(x), dvonmises(x, mu = 0, kappa = 4), lty = 5, lwd = 2, col = "darkred", zero = pi/2)
lines(circular(x), dvonmises(x, mu = 0, kappa = 10), lty = 4, lwd = 2, col = "darkgreen", zero = pi/2)
legend(2, 3, title = expression(kappa), lwd = 2,
       legend = c(0, 1, 4, 10),
       lty = c(3, 2, 5, 4), 
       col = c("black", "darkblue", "darkred", "darkgreen"), 
       bty = "n")
dev.off()
