print(open('_cell_01.py').read() if os.path.exists('_cell_01.py') else 'no _cell_01.py')
import glob
print(glob.glob('.tmp/*'))