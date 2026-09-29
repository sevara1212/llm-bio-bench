import os
import glob

print("Files in working directory:")
for f in os.listdir('.'):
    print(f, os.path.getsize(f) if os.path.isfile(f) else "dir")