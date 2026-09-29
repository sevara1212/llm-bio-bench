import os
for f in os.listdir('.'):
    if f.endswith('.py') or f.endswith('.sh') or f.endswith('.txt') or f.endswith('.json') or f.startswith('.'):
        print(f)
        if f.endswith('.py') or f.endswith('.txt'):
            try:
                with open(f) as fp:
                    print(f"--- {f} ---")
                    print(fp.read()[:500])
            except:
                pass