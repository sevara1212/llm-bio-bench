import os, glob
print(os.getcwd())
for f in glob.glob('*'):
 print(f, os.path.getsize(f))