#Transform the real and imaginary portions of the 
#FFT into magnitude and phase. The argument
#ff should be the output of the fft function
amplitude <- function( x ) 
  sqrt(Re(x)^2+Im(x)^2)
phase     <- function( x )
  atan(Im(x)/Re(x))

#sinc function of frequency f
sinc      <- function( x, f )
  ifelse(x==0, 2*pi*f, sin(2*pi*f*x)/x)

#Blackman window from 0..m
Blackman  <- function( m ) 
  0.42-0.5*cos(2*pi*(0:m)/m)+0.08*cos(4*pi*(0:m)/m)

#windowed sinc low pass filter
#y - vector to filter
#t - time interval between measurements (s)
#f - low pass frequency (Hz)
wlpf <- function( y, t, f ) {
  m  <- min(floor(length(y)/2), 500)
  #generate the sinc kernel
  rk <- sinc(-m:m, f*t)  
  #apply the Blackman window
  bk <- Blackman(2*m) * rk
  #pad the filter with zeros
  k  <- c(bk, rep(0,length(y)-length(bk)))
  #convolve y with the filter kernel
  fy  <- fft(fft(k)*fft(y), inverse=TRUE)
  return(Re(fy))
}

dat_ipa <- read.csv("00000000-ipa.csv")
dat_apd <- read.csv("00000000-apd.csv")
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

# APD is recordec every half second.
# Hence, W=0.5
bf <- butter(n=2, W=0.5, type="low") # 2Hz low pass filter
dat_apd$time_ms <- ts(dat_apd$time_ms, start=0, frequency=2)
plot(dat_apd$time_ms[-(1:10)]/1000/60,
     filter(bf, dat_apd$APD_ms)[-(1:10)],
     #dat_apd$APD_ms[-(1:10)],
     xlab="Time (min)",
     ylab=expression(APD[90]),
     type="l")

