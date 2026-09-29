import os
for f in os.listdir('.'):
 print(f,os.path.getsize(f)//1024//1024,'MB')
if os.path.exists('processed_markers.h5ad'): os.remove('processed_markers.h5ad')