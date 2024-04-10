import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Generalization_attempt.util.Experiment import PlannedSequence, Experiment, Cycle

experiment_data_file_path = 'Data/Experiment/South_timetags.npz'
experiment_datum_name = 'arr_0'
normalization_data_file_path = 'Data/Normalization/South_timetags.npz'
normalization_datum_name = 'arr_0'
PULSE_TIME = 500  # [ns]
NUM_REPEATS_PER_DETUNING = 74
DETUNING_FREQUENCIES = np.arange(-12, 67, 3)  # [MHz]
NUM_DETUNINGS = len(DETUNING_FREQUENCIES)
NUM_CYCLES = 5
SPECTRUM_CYCLE_LENGTH = PULSE_TIME * 2 * NUM_REPEATS_PER_DETUNING * NUM_DETUNINGS * NUM_CYCLES


class SpectrumPlannedSequence(PlannedSequence):
    length = PULSE_TIME * 2 * NUM_REPEATS_PER_DETUNING * NUM_DETUNINGS  # [ns]

    def bin(self):
        self.frequencies_indeces = -np.ones(2 * NUM_REPEATS_PER_DETUNING * NUM_DETUNINGS + 1, dtype=int)
        self.frequencies_indeces[1::2] = np.repeat(range(len(DETUNING_FREQUENCIES)), NUM_REPEATS_PER_DETUNING)
        self.bins = np.zeros(self.frequencies_indeces.shape)
        for i in range(len(self.timestamps[0])):
            bin_index = (self.timestamps[0][i] % self.length) // PULSE_TIME
            self.bins[bin_index] += 1

    def calculate_result_vector(self, others=None, index=None):
        # 3 conditions
        self.x_vals = np.array([DETUNING_FREQUENCIES for _ in range(4)])
        self.num_data_points = np.zeros(self.x_vals.shape)
        self.result_vector = np.zeros(self.x_vals.shape)

        for i, freq_index in enumerate(self.frequencies_indeces):
            if freq_index == -1:
                continue
            # freq = self.frequencies[i]
            # freq_index = ((i - 1) // 2) // NUM_REPEATS_PER_DETUNING
            y_prev = self.bins[i - 1]
            y_detuned = self.bins[i]
            y_next = self.bins[i + 1]

            # condition 0
            cond = 0
            if y_prev >= 0 and y_next >= 0:
                self.num_data_points[cond, freq_index] += 1
                self.result_vector[cond, freq_index] += y_detuned

            # condition 1
            cond = 1
            if y_prev >= 1 and y_next >= 1:
                self.num_data_points[cond, freq_index] += 1
                self.result_vector[cond, freq_index] += y_detuned

            # condition 2
            cond = 2
            if (y_prev >= 2 and y_next >= 1) or (y_prev >= 1 and y_next >= 2):
                self.num_data_points[cond, freq_index] += 1
                self.result_vector[cond, freq_index] += y_detuned

            # condition 3
            cond = 3
            if y_prev >= 2 and y_next >= 2:
                self.num_data_points[cond, freq_index] += 1
                self.result_vector[cond, freq_index] += y_detuned


# class SpectrumNormalizationPlannedSequence(SpectrumPlannedSequence):
#     def calculate_result_vector(self, others=None, index=None):
#         self.x_vals = DETUNING_FREQUENCIES
#         self.num_data_points = np.zeros(self.x_vals.shape)
#         self.result_vector = np.zeros(self.x_vals.shape)
#
#         for i, freq_index in enumerate(self.frequencies_indeces):
#             if freq_index == -1:
#                 continue
#             # freq = self.frequencies[i]
#             # freq_index = ((i - 1) // 2) // NUM_REPEATS_PER_DETUNING
#             y_prev = self.bins[i - 1]
#             y_detuned = self.bins[i]
#             y_next = self.bins[i + 1]
#
#             self.num_data_points[freq_index] += 1
#             self.result_vector[freq_index] += y_detuned


class SpectrumCycle(Cycle):
    def __init__(self, *args, **kwargs):
        return super().__init__(*args, **kwargs, length=SPECTRUM_CYCLE_LENGTH)


experiment_data = np.load(experiment_data_file_path, allow_pickle=True)[experiment_datum_name]
experiment = Experiment(experiment_data=[experiment_data], cycle_class=SpectrumCycle,
                        planned_sequence_class=SpectrumPlannedSequence)
experiment.analyze()
normalization_data = np.load(normalization_data_file_path, allow_pickle=True)[normalization_datum_name]
normalization_experiment = Experiment(experiment_data=[normalization_data], cycle_class=SpectrumCycle,
                                      planned_sequence_class=SpectrumPlannedSequence)
normalization_experiment.analyze()

save_data = {'detunings': experiment.x_vals[0]}
for i in range(experiment.x_vals.shape[0]):
    save_data.update({
        f'normalization_counts_{i}': normalization_experiment.num_data_points[i],
        f'normalization_sums_{i}': normalization_experiment.sum_results[i],
        f'experiment_counts_{i}': experiment.num_data_points[i],
        f'experiment_sums_{i}': experiment.sum_results[i],
    })
pd.DataFrame(save_data).to_csv('Data/output_multiple.csv')

fig, axs = plt.subplots(2, 3, squeeze=False)
axs = axs.flatten()
axs[0].plot(experiment.x_vals.T, experiment.sum_results.T / experiment.num_data_points.T)
axs[0].set_xlabel('Detuning frequency [Mhz]')
axs[0].set_ylabel('Meaned counts after test')
axs[0].set_ylim(0, axs[0].get_ylim()[1] * 1.1)
axs[0].legend([f'condition {i}' for i in range(experiment.x_vals.shape[0])])
axs[1].plot(normalization_experiment.x_vals[0].T,
            normalization_experiment.sum_results[0].T / normalization_experiment.num_data_points[0].T)
axs[1].set_xlabel('Detuning frequency [Mhz]')
axs[1].set_ylabel('Total counts no atoms (normalization)')
axs[1].set_ylim(0, axs[1].get_ylim()[1] * 1.1)
axs[2].plot(experiment.x_vals.T, experiment.num_data_points.T)
axs[2].set_xlabel('Detuning frequency [Mhz]')
axs[2].set_ylabel('Number of data points')
axs[2].legend([f'condition {i}' for i in range(experiment.x_vals.shape[0])])
axs[2].set_ylim(0, axs[2].get_ylim()[1] * 1.1)
axs[3].plot(normalization_experiment.x_vals[0].T, (experiment.sum_results.T / experiment.num_data_points.T) /
            (normalization_experiment.sum_results[0].T / normalization_experiment.num_data_points[0].T).reshape((-1, 1)))
axs[3].set_xlabel('Detuning frequency [Mhz]')
axs[3].set_ylabel('Normalized meaned counts after test')
axs[3].legend([f'condition {i}' for i in range(experiment.x_vals.shape[0])])
axs[3].set_ylim(0, axs[3].get_ylim()[1] * 1.1)

axs[4].plot(experiment.x_vals.T, experiment.num_data_points.T / normalization_experiment.num_data_points.T)
axs[4].set_xlabel('Detuning frequency [Mhz]')
axs[4].set_ylabel('SNR')
axs[4].legend([f'condition {i}' for i in range(experiment.x_vals.shape[0])])
axs[4].set_ylim(0, axs[4].get_ylim()[1] * 1.1)

plt.show()
