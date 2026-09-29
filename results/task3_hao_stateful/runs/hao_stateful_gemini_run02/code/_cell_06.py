# Let's check if there are other files or scripts or bash history
import subprocess
try:
    print(subprocess.check_output(['ls', '-la']).decode())
except Exception as e:
    print(e)