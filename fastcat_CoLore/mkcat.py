#!/usr/bin/env python
import sys
from subs import *



def main():
    o,args = setupOptions()
    initMPI(o)
    process(o)
    


if __name__=="__main__":
    main()
    
