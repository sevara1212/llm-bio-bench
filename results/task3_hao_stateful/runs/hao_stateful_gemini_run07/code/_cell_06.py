# Let's inspect the barcode names and check if there are other files or scripts or hints
print("Checking for existing files or scripts...")
for f in os.listdir('.'):
    print(f, os.path.getsize(f))