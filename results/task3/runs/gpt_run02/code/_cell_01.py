import os, glob
for root,dirs,files in os.walk('.'):
    level=root.count(os.sep)
    if level<3:
        print(root, files[:20])