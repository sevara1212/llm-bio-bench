import os,glob
for f in glob.glob('*'):print(f,os.path.getsize(f))
try: os.remove('processed.h5ad')
except:pass
print('removed',os.statvfs('.').f_bavail*os.statvfs('.').f_frsize/1e9)