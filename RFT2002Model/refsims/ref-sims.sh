date > standard.out
../rodriguez4 < ref-flex.cfg
../rodriguez4 < ref-mirr.cfg
../rodriguez4 < ref-expo.cfg
../rodriguez4 < ref-inst.cfg
date >> standard.out
mutt -s "Job Complete" matt.shotwell@vanderbilt.edu < standard.out

