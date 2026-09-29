import os,shutil
for f in ['qc_processed.h5ad']:
 if os.path.exists(f): os.remove(f)
print(shutil.disk_usage('.'))