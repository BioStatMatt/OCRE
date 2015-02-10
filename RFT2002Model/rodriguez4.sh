date > standard.out
./rodriguez4 <<EOF
&control
   simname="standard",
   simlen=4320000,
   isch1s=1200000,
   isch1e=1560000,
   isch2s=2760000,
   isch2e=3120000,
   pintvl=500,
   rectype="fullexpo"
/
EOF
date >> standard.out
R CMD BATCH plot-apd-ipa.R
cat plot-apd-ipa.Rout >> standard.out
mutt -s "Job Complete" -a standard.pdf -- matt.shotwell@vanderbilt.edu < standard.out
