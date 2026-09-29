import os, glob
for p in glob.glob('*'):
 print(p, os.path.getsize(p))