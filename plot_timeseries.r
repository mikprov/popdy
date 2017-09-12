# Plot timeseries
# by: Mikaela Provost
# on: 2017-8-16

# ---
# goals:
# ---
# 1. make a three panel plot with these graphs:
# 2. plot egg production for multiple harvest levels
# 3. plot catch 

#par.old = par()
par(mfcol = c(2,2))

# legend code
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", 
       legend = c("F = 0.00", "F = 0.4", "F = 0.8", "F = 1.2", "F = 1.6", "F = 2","F = 2.4"), 
       xpd = TRUE,
       lwd = 2, col = c(colorRampPalette(c("blue", "red"))(length(output))), title = "Exploitation Rate F")

# 1a) plot egg production for multiple harvest levels
# with multiple harvest levels
for(i in 1:length(output)) output[[i]]$eggs = log(output[[i]]$eggs) # take log of eggs
# plot egg production sensitivity for all fishing levels for a short window at end of simulation
plot(output[[length(output)]]$eggs[(time-100):time], type = "n", 
     main = "FORCING RECRUITMENT",
     #ylim = c(min(output[[length(output)]]$eggs[(time-100):time]), max(output[[1]]$eggs[(time-100):time])), 
     ylim = c(13,30),
     ylab = "Egg production",
     xlab = "Time (yr)")
mtext(letters[1], side = 3, line = -1, adj = 0.9, cex = 0.8)
for(i in 1:(length(output))) {
  lines(output[[i]]$eggs[(time-100):time], 
        col = colorRampPalette(c("blue", "red"))(length(output))[i] )
}


# 1b) plot catch at all fishing levels
for(i in 1:length(output)) output[[i]]$catch = log(output[[i]]$catch)
plot(output[[length(output)]]$catch[(time-100):time], type = "n", 
     #ylim = c(min(output[[length(output)]]$catch[(time-100):time]), max(output[[2]]$catch[(time-100):time])),
     main = "Forcing Recruitment",
     ylim = c(21.25, 28.5),
     ylab = "Catch (kg)",
     xlab = "Time (yr)")
mtext(letters[2], side = 3, line = -1, adj = 0.9, cex = 0.8)
for(i in 2:(length(output))) {
  lines(output[[i]]$catch[(time-100):time], col = colorRampPalette(c("blue", "red"))(length(output))[i])
}

# 1c) plot all fishing levels
for(i in 1:length(output)) output[[i]]$n = log(output[[i]]$n)
plot(output[[length(output)]]$n[(time-100):time], type = "n", 
     #ylim = c(min(output[[length(output)]]$catch[(time-100):time]), max(output[[2]]$catch[(time-100):time])),
     main = "Forcing Recruitment",
     ylim = c(10,22),
     ylab = "Recruitment",
     xlab = "Time (yr)")
mtext(letters[2], side = 3, line = -1, adj = 0.9, cex = 0.8)
for(i in 2:(length(output))) {
  lines(output[[i]]$n[(time-100):time], col = colorRampPalette(c("blue", "red"))(length(output))[i])
}

par(mfcol = c(1,1))




# --- 
# plot egg production sensitivity for all fishing levels
# ---
# for some reason, the periodogram for log(eggs) looks really funky, not sure why yet
# here i take the inverse of the log
for(i in 1:length(output)) output[[i]]$eggs = exp(output[[i]]$eggs) # take inverse of log
# set up spans, the vector of integers giving the widths of smoothers
tmp <- ceiling(sqrt(length(1:(time-100))))
if (tmp %% 2 == 0) {m <- tmp+1} else {m <- tmp}
# Looks better if smoother, set m higher
# m = (2 * m) + 1
m = 3 * m
m = 1 * m
low_f = spec.pgram(output[[1]]$eggs[2:time], spans = c(m,m), plot = FALSE)
high_f = spec.pgram(output[[length(output)]]$eggs[2:time], spans = c(m,m), plot = FALSE)
spec.pgram(output[[length(output)]]$eggs[2:time], spans = c(m,m), type = "n", 
           main = paste("", "\ "),
           ylim = c(min(high_f$spec), max(low_f$spec)), sub = "",
           ylab = expression("Power spectral density "~(eggs^2)),
           #ylab = "",
           xlab = "Frequency (1/yr)")
#mtext(letters[1], side = 3, line = -1.1, adj = 0.98, cex = 0.8)
for(i in 1:(length(output))) {
  spec.pgram(output[[i]]$eggs[2:time], 
             col = colorRampPalette(c("blue", "red"))(length(output))[i], 
             spans = c(m,m), add = TRUE)
}
