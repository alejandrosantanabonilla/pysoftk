import numpy as np
import os

cwd=os.path.join(os.getcwd(), 'X_embedded.npy')
cyclic_y_pred=np.load(cwd, allow_pickle=True)
print (cyclic_y_pred)
