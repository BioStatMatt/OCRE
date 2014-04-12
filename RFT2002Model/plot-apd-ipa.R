dat_ipa <- read.csv("ipa.csv")
dat_apd <- read.csv("apd.csv")
layout(matrix(1:6,3,2,byrow=TRUE))
par(mar = c(4.1,4.1,1.1,2.1))
plot(dat_ipa$time_ms/1000/60,
     log10(dat_ipa$taudiff), 
     xlab="Time (min)",
     ylab=expression(log[10](taudiff)),
     type='s')
plot(dat_ipa$time_ms/1000/60,
     dat_ipa$fatpfactor, 
     xlab="Time (min)",
     ylab="fatpfactor",
     type='s')
plot(dat_ipa$time_ms/1000/60,
     dat_ipa$finhib, 
     xlab="Time (min)",
     ylab="finhib",
     type='s')
plot(dat_ipa$time_ms/1000/60,
     dat_ipa$vcleft, 
     xlab="Time (min)",
     ylab="vcleft",
     type='s')
plot(dat_ipa$time_ms/1000/60,
     dat_ipa$inasfinal, 
     xlab="Time (min)",
     ylab="inasfinal",
     type='s')
plot(dat_apd$time_ms/1000/60,
     dat_apd$APD_ms,
     xlab="Time (min)",
     ylab=expression(APD[90]),
     type="l")

