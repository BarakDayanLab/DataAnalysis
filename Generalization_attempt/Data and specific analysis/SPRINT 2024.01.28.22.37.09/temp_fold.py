import matplotlib.pyplot as plt
import numpy as np
from Generalization_attempt.util.Helpers import merge

### Experiment raw data
experiment_data_file_paths = [
    'Data/222217_Photon_TimeTags/output/Bright(1,2)/Bright_timetags.npz',
    'Data/222217_Photon_TimeTags/output/Dark(3,4)/Dark_timetags.npz',
    'Data/222217_Photon_TimeTags/output/North(8)/North_timetags.npz',
    'Data/222217_Photon_TimeTags/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/222217_Photon_TimeTags/output/South(5)/South_timetags.npz',
]
NORTH_INDICES = [0, 1, 2]
SOUTH_INDICES = [3, 4]

north_sequence_data_file_path = 'Data/222217_Photon_TimeTags/input/sequences/North_sequence_vector.npz'
north_sequence_datum_name = 'arr_0'
SEQUENCE_LENGTH = len(np.load(north_sequence_data_file_path)[north_sequence_datum_name])

experiment_datum_names = ['arr_0'] * len(experiment_data_file_paths)
experiment_data = [
    np.load(experiment_data_file_path, allow_pickle=True)[experiment_datum_name]
    for experiment_data_file_path, experiment_datum_name in zip(experiment_data_file_paths, experiment_datum_names)
]

experiment_data = [
    [
        merge(*(experiment_data[idx][i] for idx in NORTH_INDICES)) for i in range(len(experiment_data[0]))  # North
    ],
    [
        merge(*(experiment_data[idx][i] for idx in SOUTH_INDICES)) for i in range(len(experiment_data[0]))  # South
    ],
]

folded_data = [[val % SEQUENCE_LENGTH for arr in experiment_data[i] for val in arr] for i in
               range(len(experiment_data))]
i0 = 480
i1 = SEQUENCE_LENGTH
h, _ = np.histogram(folded_data[1], bins=np.arange(i0, i1 + 1))
first_half = np.zeros(h.shape, dtype=bool)
first_half[:len(h) // 2] = 1

first_mid = np.median(np.where((h[:-1] <= 0.5 * np.max(h)) & (h[1:] >= 0.5 * np.max(h)) & (first_half[:-1])))
second_mid = np.median(np.where((h[1:] <= 0.5 * np.max(h)) & (h[:-1] >= 0.5 * np.max(h)) & (~first_half[:-1])))
print(i0 + first_mid, i0 + second_mid)

plt.hist(folded_data[0], bins=np.arange(SEQUENCE_LENGTH + 1), alpha=0.5, label='North')
plt.hist(folded_data[1], bins=np.arange(SEQUENCE_LENGTH + 1), alpha=0.5, label='South')
plt.legend(loc='upper right')
plt.show()
