import os, glob
print(os.getcwd())
for x in glob.glob('*'):
 print(x, os.path.getsize(x))