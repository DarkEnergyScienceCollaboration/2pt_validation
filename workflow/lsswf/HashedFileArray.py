#
# Abtraction of numpy file arrays ordered by hash
#
from io import BytesIO 
import numpy as np
import hashlib
import os

class HashedFileArray:
    def __init__(self,arr,directory):
        self.directory=directory
        _file=BytesIO()
        np.savetxt(_file,arr)
        self.fcontents=_file.getvalue()
        _file.close()
        m = hashlib.md5()
        m.update(self.fcontents)
        self.name=m.hexdigest()

    def full_filename(self):
        return os.path.join(self.directory,self.name)
        
    def filename(self):
        return self.name

    def genFile(self):
        if not os.path.exists(self.full_filename()):
            open(self.full_filename(),'wb').write(self.fcontents)

