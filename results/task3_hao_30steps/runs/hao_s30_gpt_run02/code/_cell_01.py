import os, glob
for f in glob.glob('*'):
 print(f, os.path.getsize(f))