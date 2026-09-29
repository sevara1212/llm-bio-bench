import importlib.util
for x in ['scrublet','bbknn','harmonypy'] : print(x,importlib.util.find_spec(x))