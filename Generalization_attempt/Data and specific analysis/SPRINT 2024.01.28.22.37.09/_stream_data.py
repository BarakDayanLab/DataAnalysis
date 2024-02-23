from time import sleep

import pandas as pd
import numpy as np
from Generalization_attempt.util.Helpers import merge

### Experiment raw data
experiment_data_file_paths = [
    'Data/223709_Photon_TimeTags/output/Bright(1,2)/Bright_timetags.npz',
    'Data/223709_Photon_TimeTags/output/Dark(3,4)/Dark_timetags.npz',
    'Data/223709_Photon_TimeTags/output/North(8)/North_timetags.npz',
    'Data/223709_Photon_TimeTags/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/223709_Photon_TimeTags/output/South(5)/South_timetags.npz',
]
NORTH_INDICES = [0, 1, 2]
SOUTH_INDICES = [3, 4]

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

north_sequence_data_file_path = 'Data/223709_Photon_TimeTags/input/sequences/North_sequence_vector.npz'
north_sequence_datum_name = 'arr_0'
SEQUENCE_LENGTH = len(np.load(north_sequence_data_file_path)[north_sequence_datum_name])

pulses_location_data_file_path = 'Data/223709_Photon_TimeTags/input/sequences/Pulses_location_in_seq.npz'
pulses_location_datum_name = 'arr_0'
SEQUENCE_RULE = np.load(pulses_location_data_file_path)[pulses_location_datum_name]

### Experiment analyzed data
df = pd.read_csv('output_summary.csv')

### manually analyze data
print(SEQUENCE_RULE)
for i in range(len(df)):
    cycle_index = df['cycle_index'].iloc[i]
    sequence_index = df['sequence_index'].iloc[i]
    north_pulses = experiment_data[0][cycle_index]
    north_pulses = north_pulses[
        (north_pulses < (sequence_index + 1) * SEQUENCE_LENGTH) & (north_pulses >= sequence_index * SEQUENCE_LENGTH)]
    south_pulses = experiment_data[1][cycle_index]
    south_pulses = south_pulses[
        (south_pulses < (sequence_index + 1) * SEQUENCE_LENGTH) & (south_pulses >= sequence_index * SEQUENCE_LENGTH)]

    print(cycle_index, sequence_index, north_pulses%SEQUENCE_LENGTH, south_pulses%SEQUENCE_LENGTH)
    sleep(1)
