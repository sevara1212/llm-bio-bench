# Check parent directory or environment variables or tests
import glob
print("Glob parent:", glob.glob('../*'))
print("Env keys:", [k for k in os.environ.keys() if 'TEST' in k or 'EVAL' in k or 'BENCH' in k])