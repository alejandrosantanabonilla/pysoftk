import numpy as np
import os
import cloudpickle

def load_array(file_path):
    """Loads a NumPy array potentially containing Numba-compiled code using cloudpickle."""
    with open(file_path, 'rb') as f:
        try:
            # Try loading with NumPy first
            return np.load(f, allow_pickle=True)
        except (ValueError, TypeError):
            # If NumPy fails, try cloudpickle
            f.seek(0)  # Reset file pointer
            return cloudpickle.load(f)

#cwd = os.path.join(os.getcwd(), 'X_embedded.npy')
#x_emb = load_array(cwd)
#print(x_emb)

#cyclic_y_pred = load_array("cyclic_y_pred.npy")
#print(cyclic_y_pred)

#cyclic_x_emb = load_array("cyclic_x_embedding.npy")
#print(cyclic_x_emb)

cyclic_av_rep = load_array("cyclic_av_rep.npy")
print(cyclic_av_rep)

#all_pos = load_array("all_pos.npy")
#print(all_pos)


