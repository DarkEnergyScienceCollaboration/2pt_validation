#!/usr/bin/env python
import sys
sys.path+=['.','..']
import lsswf

colores=[lsswf.colore([('global/seed',seed)]) for seed in range(10)]
f=[c.run() for c in colores]
print(f)

