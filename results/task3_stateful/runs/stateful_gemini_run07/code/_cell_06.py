# Let's check the contents of _cell_01.py or other files if any
for f in ['_cell_01.py', '_cell_01.py.out']:
    if os.path.exists(f):
        print(f"--- {f} ---")
        with open(f) as fp:
            print(fp.read()[:500])