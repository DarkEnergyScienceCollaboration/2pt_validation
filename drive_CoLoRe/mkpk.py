import pyccl as ccl
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv)!=9 :
    print "Usage : mkpk.py Om Ob h s8 ns pk_type tf_type fname_out"
    exit(1)
om=float(sys.argv[1])
ob=float(sys.argv[2])
hh=float(sys.argv[3])
s8=float(sys.argv[4])
ns=float(sys.argv[5])
pkt=sys.argv[6]
tft=sys.argv[7]
fno=sys.argv[8]
lk0=-4.
lkf=2.
nk=512
karr=10.**np.linspace(lk0,lkf,nk)

print "Generating power spectrum for params: "
print " Om=",om
print " Ob=",ob
print " hh=",hh
print " s8=",s8
print " ns=",ns
print " pk_type:",pkt
print " tf_type:",tft
print " "
cosmo=ccl.Cosmology(ccl.Parameters(Omega_c=om-ob,Omega_b=ob,h=hh,sigma8=s8,n_s=ns,),
                    matter_power_spectrum=pkt,transfer_function=tft)
pkarr=ccl.linear_matter_power(cosmo,karr*hh,1.)*hh**3

np.savetxt(fno,np.transpose([karr,pkarr]))
