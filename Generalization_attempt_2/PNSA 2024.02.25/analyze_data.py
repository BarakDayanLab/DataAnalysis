from collections import OrderedDict
from functools import reduce

import numpy as np
import pandas as pd

from Generalization_attempt_2.Analyze import bin_and_reshape

# DEFINE CONSTANTS
transit_conditions = [[1, 2], [2, 1], [2, 0, 2], [1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1]]
transit_conditions = [np.array(x) for x in transit_conditions]
PULSE_TYPES = {
    'N': 0,
    'S': 1,
    'n': 2,
    's': 3,
}
experiment_data_file_paths = [
    'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/output/Bright(1,2)/Bright_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/output/Dark(3,4)/Dark_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/output/North(8)/North_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/output/South(5)/South_timetags.npz',
]
normalization_data_file_paths = [
    'Data/173049_Photon_TimeTags/Iter_1_Seq_1__Without Atoms/output/Bright(1,2)/Bright_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_1__Without Atoms/output/Dark(3,4)/Dark_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_1__Without Atoms/output/North(8)/North_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_1__Without Atoms/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/173049_Photon_TimeTags/Iter_1_Seq_1__Without Atoms/output/South(5)/South_timetags.npz',
]
BRIGHT_INDEX = 0
DARK_INDEX = 1
NORTH_INDICES = [0, 1, 2]
SOUTH_INDICES = [3, 4]

experiment_datum_names = ['arr_0'] * len(experiment_data_file_paths)
normalization_datum_names = ['arr_0'] * len(normalization_data_file_paths)
pulses_location_data_file_path = 'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/input/sequences/Pulses_location_in_seq.npz'
pulses_location_datum_name = 'arr_0'
north_sequence_data_file_path = 'Data/173049_Photon_TimeTags/Iter_1_Seq_2__With Atoms/input/sequences/North_sequence_vector.npz'
north_sequence_datum_name = 'arr_0'

SEQUENCE_LENGTH = len(np.load(north_sequence_data_file_path)[north_sequence_datum_name])
CYCLE_LENGTH = int(8e6)  # [ns], = 8ms
CYCLE_LENGTH = int((CYCLE_LENGTH // SEQUENCE_LENGTH) * SEQUENCE_LENGTH)
seq_rule_temp = np.load(pulses_location_data_file_path)[pulses_location_datum_name]
SEQUENCE_RULE = np.concatenate(
    (seq_rule_temp[:, :2].astype(int), np.array([PULSE_TYPES[x] for x in seq_rule_temp[:, 2]]).reshape((-1, 1))),
    axis=1)

experiment_data = [
    np.load(experiment_data_file_path, allow_pickle=True)[experiment_datum_name]
    for experiment_data_file_path, experiment_datum_name in zip(experiment_data_file_paths, experiment_datum_names)
]

# ANALYZE DATA
bins = SEQUENCE_RULE[:, :2].flatten()
binned_experiment_data = bin_and_reshape(experiment_data, sequence_length=SEQUENCE_LENGTH, cycle_length=CYCLE_LENGTH,
                                         bins=bins)

# binned_experiment_data has shape [number of streams, # of cycles, # of seq per cycle, # bins per seq]

# general binning
N2S_by_seq = np.sum(binned_experiment_data[NORTH_INDICES][:, :, :, SEQUENCE_RULE[:, 2] == PULSE_TYPES['S']], axis=0)
S2N_by_seq = np.sum(binned_experiment_data[SOUTH_INDICES][:, :, :, SEQUENCE_RULE[:, 2] == PULSE_TYPES['N']], axis=0)
N2S = np.sum(N2S_by_seq, axis=2)
S2N = np.sum(S2N_by_seq, axis=2)
reflections = N2S + S2N

# search for transits
transit_condition_met = np.zeros((len(transit_conditions), *reflections.shape), dtype=bool)
relevant_for_condition = np.ones((len(transit_conditions), *reflections.shape), dtype=bool)
for i, transit_condition in enumerate(transit_conditions):
    window_size = len(transit_condition)
    # relevant_for_condition[i] = np.ones(reflections.shape, dtype=bool)
    relevant_for_condition[i, :, :window_size - 1] = 0
    relevant_for_condition[i, :, -window_size + 1:] = 0

    reshaped = np.lib.stride_tricks.sliding_window_view(reflections, (1, window_size))
    temp = np.all(reshaped >= transit_condition, axis=(2, 3))
    temp = np.concatenate((
        np.zeros((temp.shape[0], window_size - 1), dtype=bool),
        temp,
        np.zeros((temp.shape[0], window_size - 1), dtype=bool),
    ), axis=1)
    reshaped = np.lib.stride_tricks.sliding_window_view(temp, (1, window_size))
    transit_condition_met[i] = np.any(reshaped, axis=(2, 3))

LAST_DET_PULSE_INDEX = 5  # TODO: this with code
last_det_pulse_is_north = SEQUENCE_RULE[LAST_DET_PULSE_INDEX, 2] == PULSE_TYPES['N']
if last_det_pulse_is_north:
    last_detection_pulse_reflection = (N2S_by_seq[:, :, -1] > 0)
else:
    last_detection_pulse_reflection = (S2N_by_seq[:, :, -1] > 0)

SECOND_SPRINT_PULSE_INDEX = -1
FIRST_SPRINT_PULSE_INDEX = -2
total_sprint_single_photon = (
        np.sum(binned_experiment_data[:, :, :, (FIRST_SPRINT_PULSE_INDEX, SECOND_SPRINT_PULSE_INDEX)],
               axis=(0, 3)) == 1)

second_sprint_pulse_is_north = SEQUENCE_RULE[SECOND_SPRINT_PULSE_INDEX, 2] == PULSE_TYPES['n']
second_sprint_north = np.sum(binned_experiment_data[NORTH_INDICES, :, :, SECOND_SPRINT_PULSE_INDEX], axis=0)
second_sprint_south = np.sum(binned_experiment_data[SOUTH_INDICES, :, :, SECOND_SPRINT_PULSE_INDEX], axis=0)
second_sprint_bright = binned_experiment_data[BRIGHT_INDEX, :, :, SECOND_SPRINT_PULSE_INDEX]
second_sprint_dark = binned_experiment_data[DARK_INDEX, :, :, SECOND_SPRINT_PULSE_INDEX]
first_sprint_north = np.sum(binned_experiment_data[NORTH_INDICES, :, :, FIRST_SPRINT_PULSE_INDEX], axis=0)
first_sprint_south = np.sum(binned_experiment_data[SOUTH_INDICES, :, :, FIRST_SPRINT_PULSE_INDEX], axis=0)

sprint_reflections, sprint_transmissions = first_sprint_south + second_sprint_south, first_sprint_north + second_sprint_north
if not second_sprint_pulse_is_north:
    # assume first and second sprint pulses are the same
    sprint_transmissions, sprint_reflections = sprint_reflections, sprint_transmissions

# save the data
condition_summary = []
subtitles = [
    'transits',
    'detection_reflection',
    'single_photon_sprint',
    'north_single_photon',
    'south_single_photon',
    'bright_photons',
    'dark_photons',
    'reflected_sprint_photons',
]
titles = []

for i, transit_condition in enumerate(transit_conditions):
    titles.append(f'condition#{i}')
    cond = relevant_for_condition[i] & transit_condition_met[i]
    condition_summary.append([])
    condition_summary[-1].extend(
        [
            np.sum(cond),
            np.sum(cond & last_detection_pulse_reflection),
            np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon),
            np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon & (
                    first_sprint_north | second_sprint_north)),
            np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon & (
                    first_sprint_south | second_sprint_south)),
            np.sum(cond & last_detection_pulse_reflection & second_sprint_bright),
            np.sum(cond & last_detection_pulse_reflection & second_sprint_dark),
            np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon & sprint_reflections),
        ]
    )
cond = reduce(np.logical_or,
              (relevant_for_condition[i] & transit_condition_met[i] for i in range(len(transit_conditions))))
titles.append('condition#any')
condition_summary.append([])
condition_summary[-1].extend([
    np.sum(cond),
    np.sum(cond & last_detection_pulse_reflection),
    np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon),
    np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon & (
            first_sprint_north | second_sprint_north)),
    np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon & (
            first_sprint_south | second_sprint_south)),
    np.sum(cond & last_detection_pulse_reflection & second_sprint_bright),
    np.sum(cond & last_detection_pulse_reflection & second_sprint_dark),
    np.sum(cond & last_detection_pulse_reflection & total_sprint_single_photon & sprint_reflections),
])
# print(summary)

df = pd.DataFrame(data=np.array(condition_summary).T, index=subtitles, columns=titles)
df.to_csv('output_summary.csv')
