
#######################################################
# ATTEMPT AT DOING A SPECIES ACCUMULATION CURVE
#######################################################

library(rcompanion)
library(dplyr)

data <- read.csv("estimate_mlg/mlg_accumulation_long.csv")
#head(data)

sample <- data$Sample_size
mlg    <- data$MLGs

# plot(sample, mlg)

# linear models
model0 <- lm(mlg~I(sample^1.00)+sample, data=data)
model1 <- lm(mlg~I(sample^1.15)+sample, data=data)
model2 <- lm(mlg~I(sample^1.20)+sample, data=data)
model3 <- lm(mlg~I(sample^1.30)+sample, data=data)
model4 <- lm(mlg~I(sample^1.40)+sample, data=data)
model5 <- lm(mlg~I(sample^1.50)+sample, data=data)
model6 <- lm(mlg~log(sample), data=data)

# compare models
compare_lm <- compareLM(model0,model1,model2,model3,model4,model5,model6)
compare_lm

# model with the lowest AIC?
# mlg ~ I(sample^1.15) + sample
# order(compare_lm$Fit.criteria$AIC)
print("Best model")
compare_lm$Models[order(compare_lm$Fit.criteria$AIC)[1]]

# Build the model
model <- lm(mlg ~ I(sample^1.15) + sample, data=data)
new.speeds <- data.frame(sample = c(seq(0,780,1)))
new <- predict(model, newdata = new.speeds, se.fit=F)
CI <- predict(model, newdata = new.speeds, se.fit=T, interval="prediction")
plot(mlg ~ sample, data = data, xlim=c(0,780), ylim=c(0,200))
lines(new ~ new.speeds$sample, lwd=2)
lines(CI$fit[,3] ~ new.speeds$sample, lty=2, lwd=1.5)
lines(CI$fit[,2] ~ new.speeds$sample, lty=2, lwd=1.5)

# max(new)

max <- max(new)
index <- which(new==max)
upper <- CI$fit[index,3]
lower <- CI$fit[index,2]

print("Model predictions")
cbind(max, lower,upper)

# dev.off()
# pdf("mlg_accumulation_curve.pdf", height = 6, width = 8)
# plot(mlg ~ sample, data = data, xlim=c(0,800), ylim=c(0,220), xlab="Samples (N)", ylab="MLGs (N)", cex=0.5, cex.lab=1.2, cex.axis=1.2)
# lines(new ~ new.speeds$sample, lwd=2, col="navy")
# lines(CI$fit[,3] ~ new.speeds$sample, lty=2, lwd=1.5, col="navy")
# lines(CI$fit[,2] ~ new.speeds$sample, lty=2, lwd=1.5, col="navy")
# abline(h=max, col="red", lty=3, lwd=2)
# legend("right", legend = c("Extrapolation",
#                            "Extrapolation upper/lower",
#                            "Max. MLGs"),
#        col=c("navy", "navy", "red"), lty = c(1,2,3), lwd = 1.5, inset = 0.05)
# dev.off()


dev.off()
png("estimate_mlg/mlg_accumulation_curve.png", height = 8, width = 8, res=600, units = "in")
plot(mlg ~ sample, data = data, xlim=c(0,800), ylim=c(0,220), xlab="Samples (N)", ylab="MLGs (N)", cex=0.5, cex.lab=1.2, cex.axis=1.2)
lines(new ~ new.speeds$sample, lwd=2, col="navy")
lines(CI$fit[,3] ~ new.speeds$sample, lty=2, lwd=1.5, col="navy")
lines(CI$fit[,2] ~ new.speeds$sample, lty=2, lwd=1.5, col="navy")
abline(h=max, col="red", lty=3, lwd=2)
legend("right", legend = c("Extrapolation",
                           "Extrapolation upper/lower",
                           "Max. MLGs"),
       col=c("navy", "navy", "red"), lty = c(1,2,3), lwd = 1.5, inset = 0.05)
dev.off()
