library(tidyverse)
library(SparseTransitions)

## Modify below call to use your own Julia installation
#jl_pkg_setup('/Applications/Julia-1.2.app/Contents/Resources/julia/bin')

nf=8
nams=c("E1",	"E2",	"E3",	"pEMT1",	"pEMT2",	"pEMT3",	"M","pMET")

#  Data Set 1, with rowsum normalized to 10,0000
r0 = c(4417, 2240, 2080, 645, 137, 333, 105, 43)
r1 = c(158, 1730, 1145, 2932, 939, 2268, 670, 157)
r2 = c(8, 162, 128, 842, 1835, 2092, 4653, 280)
r3 = c(2, 32, 12, 283, 909, 1147, 7285, 330)
r4 = c(27, 111, 9, 838, 1532, 165, 3826, 3492)
r5 = c(216, 1107, 28, 1365, 1043, 301, 4156, 1785)
r6 = c(1085, 2470, 53, 1027, 783, 198, 2240, 2144)


rcount1 = t(rbind(r0,r1,r2,r3,r4,r5,r6))

dimnames(rcount1)=list(nams,c("0","2","6","10","w2","w6","w10"))

rprop1=scale(rcount1,F,colSums(rcount1))  #make colsums eq to 1






for (i in 1:80){
  print(paste0("boot:",i))
  boot_res_tmp <- boot_run(rcount1, i)
  saveRDS(boot_res_tmp, file=paste0("boot_samples/boot_res_tmp_",i,".Rds"))
}
