import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Generalization_attempt.util.Experiment import PlannedSequence, Experiment, Cycle
from Generalization_attempt.util.Helpers import merge

PULSE_TYPES = {
    'N': 0,
    'S': 1,
    'n': 2,
    's': 3,
}
experiment_data_file_paths = [
    'Data/223709_Photon_TimeTags/output/Bright(1,2)/Bright_timetags.npz',
    'Data/223709_Photon_TimeTags/output/Dark(3,4)/Dark_timetags.npz',
    'Data/223709_Photon_TimeTags/output/North(8)/North_timetags.npz',
    'Data/223709_Photon_TimeTags/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/223709_Photon_TimeTags/output/South(5)/South_timetags.npz',
]
normalization_data_file_paths = [
    'Data/222217_Photon_TimeTags/output/Bright(1,2)/Bright_timetags.npz',
    'Data/222217_Photon_TimeTags/output/Dark(3,4)/Dark_timetags.npz',
    'Data/222217_Photon_TimeTags/output/North(8)/North_timetags.npz',
    'Data/222217_Photon_TimeTags/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/222217_Photon_TimeTags/output/South(5)/South_timetags.npz',
]
NORTH_INDICES = [0, 1, 2]
SOUTH_INDICES = [3, 4]

experiment_datum_names = ['arr_0'] * len(experiment_data_file_paths)
normalization_datum_names = ['arr_0'] * len(normalization_data_file_paths)
pulses_location_data_file_path = 'Data/223709_Photon_TimeTags/input/sequences/Pulses_location_in_seq.npz'
pulses_location_datum_name = 'arr_0'
north_sequence_data_file_path = 'Data/223709_Photon_TimeTags/input/sequences/North_sequence_vector.npz'
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

experiment_data = [
    [
        merge(*(experiment_data[idx][i] for idx in NORTH_INDICES)) for i in range(len(experiment_data[0]))  # North
    ],
    [
        merge(*(experiment_data[idx][i] for idx in SOUTH_INDICES)) for i in range(len(experiment_data[0]))  # South
    ],
]

normalization_data = [
    np.load(normalization_data_file_path, allow_pickle=True)[normalization_datum_name]
    for normalization_data_file_path, normalization_datum_name in
    zip(normalization_data_file_paths, normalization_datum_names)
]

normalization_data = [
    [
        merge(*(normalization_data[idx][i] for idx in NORTH_INDICES)) for i in range(len(normalization_data[0]))
        # North
    ],
    [
        merge(*(normalization_data[idx][i] for idx in SOUTH_INDICES)) for i in range(len(normalization_data[0]))
        # South
    ],
]


class SprintCycle(Cycle):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, length=CYCLE_LENGTH)


class SprintPlannedSequence(PlannedSequence):
    length = SEQUENCE_LENGTH

    def bin(self):
        self.bins = np.zeros((len(SEQUENCE_RULE), 2), dtype=int)
        # self.types = [s[2] for s in SEQUENCE_RULE]

        i0_N = np.searchsorted(SEQUENCE_RULE[:, 0], self.timestamps[0] % self.length)
        i1_N = np.searchsorted(SEQUENCE_RULE[:, 1], self.timestamps[0] % self.length)
        i0_S = np.searchsorted(SEQUENCE_RULE[:, 0], self.timestamps[1] % self.length)
        i1_S = np.searchsorted(SEQUENCE_RULE[:, 1], self.timestamps[1] % self.length)
        relevant_N = (i0_N == i1_N + 1)
        relevant_S = (i0_S == i1_S + 1)
        # NORTH
        for i in i1_N[relevant_N]:
            self.bins[i, 0] += 1
        # SOUTH
        for i in i1_S[relevant_S]:
            self.bins[i, 1] += 1

        # a photon measured in the north, induced by a transmission
        # self.north_transmission = self.bins[SEQUENCE_RULE[:, 2] == 'N', 0]
        # a photon measured in the north, induced by a reflection
        self.north_reflections = self.bins[SEQUENCE_RULE[:, 2] == PULSE_TYPES['S'], 0]
        # a photon measured in the south, induced by a transmission
        # self.south_transmission = self.bins[SEQUENCE_RULE[:, 2] == 'S', 1]
        # a photon measured in the south, induced by a reflection
        self.south_reflections = self.bins[SEQUENCE_RULE[:, 2] == PULSE_TYPES['N'], 1]

        self.sprint_transmissions = self.bins[-1, 0] if SEQUENCE_RULE[-1, 2] == PULSE_TYPES['n'] else self.bins[-1, 1]
        self.sprint_reflections = self.bins[-1, 0] if SEQUENCE_RULE[-1, 2] == PULSE_TYPES['s'] else self.bins[-1, 1]

    def calculate_result_vector(self, others=None, index=None):
        self.x_vals = np.array(
            [
                [f'condition#{i}transit', f'condition#{i}transmission', f'condition#{i}reflection']
                for i in range(2)
            ]
        )
        self.result_vector = np.zeros(self.x_vals.shape)
        self.num_data_points = np.zeros(self.x_vals.shape)

        if index == len(others) - 1:
            return
        next = others[index + 1]

        reflections = np.sum([self.north_reflections, self.south_reflections])
        next_reflections = np.sum([next.north_reflections, next.south_reflections])
        last_detection_pulse_reflected_photons = self.bins[-2, 0] if SEQUENCE_RULE[-2, 2] == PULSE_TYPES['S'] else \
            self.bins[-2, 1]
        sprint_photons = np.sum([self.sprint_transmissions, self.sprint_reflections])

        # print(f'index: {index}, reflections: {reflections}, sprint_photons: {sprint_photons}')

        condition = 0
        if (last_detection_pulse_reflected_photons >= 1) and (reflections >= 1) and (next_reflections >= 1):
            self.num_data_points[condition][0] += 1
            self.result_vector[condition][0] += 1
            if sprint_photons == 1:
                self.result_vector[condition][1] += self.sprint_transmissions
                self.num_data_points[condition][1] += self.sprint_transmissions
                self.result_vector[condition][2] += self.sprint_reflections
                self.num_data_points[condition][2] += self.sprint_reflections
        condition = 1
        if (last_detection_pulse_reflected_photons >= 1) and \
                (((reflections >= 2) and (next_reflections >= 1)) or ((reflections >= 1) and (next_reflections >= 2))):
            self.num_data_points[condition][0] += 1
            self.result_vector[condition][0] += 1
            if sprint_photons == 1:
                self.result_vector[condition][1] += self.sprint_transmissions
                self.num_data_points[condition][1] += self.sprint_transmissions
                self.result_vector[condition][2] += self.sprint_reflections
                self.num_data_points[condition][2] += self.sprint_reflections


experiment = Experiment(experiment_data=experiment_data, cycle_class=SprintCycle,
                        planned_sequence_class=SprintPlannedSequence)
experiment.analyze()
normalization = Experiment(experiment_data=normalization_data, cycle_class=SprintCycle,
                           planned_sequence_class=SprintPlannedSequence)
normalization.analyze()

df_experiment_summary = pd.DataFrame({'x_vals': experiment.x_vals.flatten(),
                                      'experiment_counts': experiment.num_data_points.flatten(),
                                      'experiment_sums': experiment.sum_results.flatten(),
                                      })
df_normalization_summary = pd.DataFrame({'x_vals': normalization.x_vals.flatten(),
                                         'experiment_counts': normalization.num_data_points.flatten(),
                                         'experiment_sums': normalization.sum_results.flatten(),
                                         })
# df_total = pd.DataFrame({
#     'cycle_index': [experiment.all_points[i][0] for i in range(len(experiment.all_points))],
#     'sequence_index': [experiment.all_points[i][1] for i in range(len(experiment.all_points))],
#     'num_data_points': [experiment.all_points[i][2][0] for i in range(len(experiment.all_points))],
#     'transmissions': [experiment.all_points[i][3][0] for i in range(len(experiment.all_points))],
#     'reflections': [experiment.all_points[i][3][1] for i in range(len(experiment.all_points))]
# })
# df_total.to_csv('output_total.csv', index=False)

df_experiment_summary.to_csv('output_experiment_summary.csv', index=False)
df_normalization_summary.to_csv('output_normalization_summary.csv', index=False)

print(df_experiment_summary)
print(df_normalization_summary)
