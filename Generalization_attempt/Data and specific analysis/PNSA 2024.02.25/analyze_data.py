import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Generalization_attempt.util.Experiment import PlannedSequence, Experiment, Cycle
from Generalization_attempt.util.Helpers import merge

transit_conditions = [[1, 2], [2, 1], [2, 0, 2], [1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1]]
transit_conditions = [np.array(x) for x in transit_conditions]
max_length_condition = max([len(x) for x in transit_conditions])
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

experiment_data = [
    [
        merge(*(experiment_data[idx][i] for idx in NORTH_INDICES)) for i in range(len(experiment_data[0]))  # North
    ],
    [
        merge(*(experiment_data[idx][i] for idx in SOUTH_INDICES)) for i in range(len(experiment_data[0]))  # South
    ],
    experiment_data[0],  # Bright
    experiment_data[1],  # Dark
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
    normalization_data[0],  # Bright
    normalization_data[1],  # Dark
]


class SprintCycle(Cycle):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, length=CYCLE_LENGTH)


class SprintPlannedSequence(PlannedSequence):
    length = SEQUENCE_LENGTH

    def bin(self):
        self.bins = np.zeros((len(SEQUENCE_RULE), 4), dtype=int)
        # self.types = [s[2] for s in SEQUENCE_RULE]

        i0_N = np.searchsorted(SEQUENCE_RULE[:, 0], self.timestamps[0] % self.length)
        i1_N = np.searchsorted(SEQUENCE_RULE[:, 1], self.timestamps[0] % self.length)
        i0_S = np.searchsorted(SEQUENCE_RULE[:, 0], self.timestamps[1] % self.length)
        i1_S = np.searchsorted(SEQUENCE_RULE[:, 1], self.timestamps[1] % self.length)
        i0_B = np.searchsorted(SEQUENCE_RULE[:, 0], self.timestamps[2] % self.length)
        i1_B = np.searchsorted(SEQUENCE_RULE[:, 1], self.timestamps[2] % self.length)
        i0_D = np.searchsorted(SEQUENCE_RULE[:, 0], self.timestamps[3] % self.length)
        i1_D = np.searchsorted(SEQUENCE_RULE[:, 1], self.timestamps[3] % self.length)

        relevant_N = (i0_N == i1_N + 1)
        relevant_S = (i0_S == i1_S + 1)
        relevant_B = (i0_B == i1_B + 1)
        relevant_D = (i0_D == i1_D + 1)
        # NORTH
        for i in i1_N[relevant_N]:
            self.bins[i, 0] += 1
        # SOUTH
        for i in i1_S[relevant_S]:
            self.bins[i, 1] += 1
        # BRIGHT
        for i in i1_S[relevant_B]:
            self.bins[i, 2] += 1
        # DARK
        for i in i1_S[relevant_D]:
            self.bins[i, 3] += 1

        # a photon measured in the north, induced by a transmission
        # self.north_transmission = self.bins[SEQUENCE_RULE[:, 2] == 'N', 0]
        # a photon measured in the north, induced by a reflection
        self.north_reflections = self.bins[SEQUENCE_RULE[:, 2] == PULSE_TYPES['S'], 0]
        # a photon measured in the south, induced by a transmission
        # self.south_transmission = self.bins[SEQUENCE_RULE[:, 2] == 'S', 1]
        # a photon measured in the south, induced by a reflection
        self.south_reflections = self.bins[SEQUENCE_RULE[:, 2] == PULSE_TYPES['N'], 1]

        self.detection_reflections = np.sum([self.north_reflections, self.south_reflections])

        self.sprint_transmissions = self.bins[-1, 0] if SEQUENCE_RULE[-1, 2] == PULSE_TYPES['n'] else self.bins[-1, 1]
        self.sprint_reflections = self.bins[-1, 0] if SEQUENCE_RULE[-1, 2] == PULSE_TYPES['s'] else self.bins[-1, 1]
        self.bright_photons = self.bins[-1, 2]
        self.dark_photons = self.bins[-1, 3]

    def calculate_result_vector(self, others=None, index=None):
        self.x_vals = np.array(
            [
                [
                    f'condition#{i}transit',
                    f'condition#{i}datapoint',
                    f'condition#{i}transmission',
                    f'condition#{i}reflection',
                    f'condition#{i}Bright',
                    f'condition#{i}Dark',
                ]
                for i in range(len(transit_conditions))
            ]
        )
        self.result_vector = np.zeros(self.x_vals.shape)
        self.num_data_points = np.zeros(self.x_vals.shape)

        relevant_index = (index + 1 >= max_length_condition) & (index + max_length_condition <= len(others))
        if not relevant_index:
            return

        last_detection_pulse_reflected_photons = self.bins[-2, 0] if SEQUENCE_RULE[-2, 2] == PULSE_TYPES['S'] else \
            self.bins[-2, 1]
        sprint_photons = np.sum([self.sprint_transmissions, self.sprint_reflections])

        # print(f'index: {index}, reflections: {reflections}, sprint_photons: {sprint_photons}')

        for condition_index in range(len(transit_conditions)):
            condition = transit_conditions[condition_index]
            reflections = np.array(
                [others[i].detection_reflections for i in range(index - len(condition) + 1, index + len(condition))])
            reflections = np.lib.stride_tricks.sliding_window_view(reflections, len(condition))
            transit_condition_met = np.any([np.all(reflections[i] >= condition) for i in range(len(reflections))])
            # transit_condition_met = None # there was a transit based on this condition

            if not transit_condition_met:
                return
            self.num_data_points[condition_index][0] += 1
            self.result_vector[condition_index][0] += 1
            if not (last_detection_pulse_reflected_photons >= 1):
                return
            self.num_data_points[condition_index][1] += 1
            self.result_vector[condition_index][1] += 1
            if not (sprint_photons >= 1):
                return
            self.num_data_points[condition_index][4] += ((self.bright_photons + self.dark_photons) > 0)
            self.num_data_points[condition_index][5] += ((self.bright_photons + self.dark_photons) > 0)
            self.result_vector[condition_index][4] += self.bright_photons
            self.result_vector[condition_index][5] += self.dark_photons
            if not (sprint_photons == 1):
                self.result_vector[condition_index][2] += self.sprint_transmissions
                self.num_data_points[condition_index][2] += self.sprint_transmissions
                self.result_vector[condition_index][3] += self.sprint_reflections
                self.num_data_points[condition_index][3] += self.sprint_reflections


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
