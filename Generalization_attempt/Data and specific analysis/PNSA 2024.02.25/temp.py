import matplotlib.pyplot as plt
import numpy as np

x = np.load('Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/input/sequences/Pulses_location_in_seq.npz', allow_pickle=True)['arr_0']

print(x)
# plt.plot(x)
# plt.show()
