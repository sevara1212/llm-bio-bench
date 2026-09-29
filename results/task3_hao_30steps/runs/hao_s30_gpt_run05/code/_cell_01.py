import os, glob
print(os.getcwd())
for f in glob.glob('*')+glob.glob('**/*',recursive=True):
 if os.path.isfile(f): print(f,os.path.getsize(f))