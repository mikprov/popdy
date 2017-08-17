

par(mfrow = c(2,1))#, mar = c(5, 4, 4, 2) + 0.3)

# with multiple harvest levels
for(i in 1:length(output)) output[[i]]$eggs = log(output[[i]]$eggs)
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


# plot catch sensitivity for all fishing levels
for(i in 1:length(output)) output[[i]]$catch = log(output[[i]]$catch)
plot(output[[length(output)]]$catch[(time-100):time], type = "n", 
     #ylim = c(min(output[[length(output)]]$catch[(time-100):time]), max(output[[2]]$catch[(time-100):time])),
     ylim = c(21.25, 28.5),
     ylab = "Catch (kg)",
     xlab = "Time (yr)")
mtext(letters[2], side = 3, line = -1, adj = 0.9, cex = 0.8)
for(i in 2:(length(output))) {
  lines(output[[i]]$catch[(time-100):time], col = colorRampPalette(c("blue", "red"))(length(output))[i])
}
par(mfrow = c(1,1))



# legend code
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", 
       legend = c("F = 0.00", "F = 0.4", "F = 0.8", "F = 1.2", "F = 1.6", "F = 2","F = 2.4"), 
       xpd = TRUE,
       lwd = 2, col = c(colorRampPalette(c("blue", "red"))(length(output))), title = "Exploitation Rate F")

# ------------ code below needs to be fixed ---------------- #
# plot egg production sensitivity for all fishing levels
low_f = spec.pgram(output[[1]][[1]]$eggs, spans = c(m,m), plot = FALSE)
high_f = spec.pgram(output[[length(output)]][[1]]$eggs, spans = c(m,m), plot = FALSE)
spec.pgram(output[[length(output)]][[1]]$eggs, spans = c(m,m), type = "n", 
           main = paste("TOP-DOWN FORCING (M)", "\nEgg production"),
           ylim = c(min(high_f$spec), max(low_f$spec)), sub = "",
           ylab = expression("Power spectral density "~(eggs^2)),
           xlab = "Frequency (1/yr)")
mtext(letters[1], side = 3, line = -1.1, adj = 0.98, cex = 0.8)
for(i in 1:(length(output))) {
  spec.pgram(output[[i]][[1]]$eggs, col = colorRampPalette(c("blue", "red"))(length(output))[i], spans = c(m,m), add = TRUE)
}
