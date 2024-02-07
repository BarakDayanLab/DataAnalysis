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
        self.x_vals = DETUNING_FREQUENCIES
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

            if y_prev > 0 and y_next > 0:
                self.num_data_points[freq_index] += 1
                self.result_vector[freq_index] += y_detuned


class SpectrumNormalizationPlannedSequence(SpectrumPlannedSequence):
    def calculate_result_vector(self, others=None, index=None):
        self.x_vals = DETUNING_FREQUENCIES
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

            self.num_data_points[freq_index] += 1
            self.result_vector[freq_index] += y_detuned


class SpectrumCycle(Cycle):
    def __init__(self, *args, **kwargs):
        return super().__init__(*args, **kwargs, length=SPECTRUM_CYCLE_LENGTH)


experiment_data = np.load(experiment_data_file_path, allow_pickle=True)[experiment_datum_name]
experiment = Experiment(experiment_data=[experiment_data],
                        planned_sequence_class=SpectrumPlannedSequence)
experiment.analyze()
normalization_data = np.load(normalization_data_file_path, allow_pickle=True)[normalization_datum_name]
normalization_experiment = Experiment(experiment_data=[normalization_data],
                                      planned_sequence_class=SpectrumNormalizationPlannedSequence)
normalization_experiment.analyze()

pd.DataFrame({'detunings': experiment.x_vals,
              'experiment_counts': experiment.num_data_points,
              'experiment_sums': experiment.sum_results,
              'normalization_counts': normalization_experiment.num_data_points,
              'normalization_sums': normalization_experiment.sum_results}).to_csv(
    'Data/output.csv')

fig, axs = plt.subplots(1, 2, squeeze=False)
axs = axs.flatten()
axs[0].plot(experiment.x_vals, experiment.sum_results / experiment.num_data_points)
axs[0].set_xlabel('Detuning frequency [Mhz]')
axs[0].set_ylabel('Meaned counts after test')
axs[0].set_ylim(0, 1.1 * np.max(experiment.sum_results / experiment.num_data_points))
axs[1].plot(normalization_experiment.x_vals,
            normalization_experiment.sum_results / normalization_experiment.num_data_points)
axs[1].set_xlabel('Detuning frequency [Mhz]')
axs[1].set_ylabel('Total counts no atoms (normalization)')
axs[1].set_ylim(0, 1.1 * np.max(normalization_experiment.sum_results / normalization_experiment.num_data_points))

plt.show()
