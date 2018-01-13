import numpy as np
from sklearn.utils import check_random_state


# create noisy helix
amplitude_a = 3
scale_a = 4
amplitude_b = 2
scale_b = 4
num_points = 1000
t_time = np.linspace(0, np.pi, num_points)
x_data = amplitude_a * np.cos(scale_a * t_time)
y_data = amplitude_b * np.sin(scale_b * t_time)
z_data = amplitude_b * t_time

helix_data = np.vstack((x_data, y_data, z_data))

# add noise
