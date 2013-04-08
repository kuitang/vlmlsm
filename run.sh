
#for N in 2 4 8 16 32; do
# Must run in serial to conserve memory!
for N in 24 32; do
  LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /opt/MATLAB/R2011b/bin/matlab -nodisplay -r "nNodes = $N ; basic_bbp_test" > run_$N.txt
  stty echo
done
  
