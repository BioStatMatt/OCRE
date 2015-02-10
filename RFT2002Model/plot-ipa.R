dat <- read.csv("ipa.csv")
layout(matrix(1:6,3,2,byrow=TRUE))
plot(dat$time_ms/1000/60,
     log10(dat$taudiff), 
     xlab="Time (min)",
     ylab=expression(log[10](taudiff)),
     type='s')
plot(dat$time_ms/1000/60,
     dat$fatpfactor, 
     xlab="Time (min)",
     ylab="fatpfactor",
     type='s')
plot(dat$time_ms/1000/60,
     dat$finhib, 
     xlab="Time (min)",
     ylab="finhib",
     type='s')
plot(dat$time_ms/1000/60,
     dat$vcleft, 
     xlab="Time (min)",
     ylab="vcleft",
     type='s')
plot(dat$time_ms/1000/60,
     dat$inasfinal, 
     xlab="Time (min)",
     ylab="inasfinal",
     type='s')
