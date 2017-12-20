import os

rootDir=os.path.join(os.environ['SCRATCH'],'LSSWF')

def maybeMakeDir(name):
    dname=os.path.join(rootDir,name)
    if not os.path.exists(dname):
        os.makedirs(dname)
    return dname
