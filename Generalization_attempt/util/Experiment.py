import numpy as np
from tqdm import tqdm


class Cycle:
    def __init__(self, timestamps, planned_sequence_class, length=10e6):
        self.timestamps = timestamps
        self.planned_sequence_class = planned_sequence_class
        self.length = length
        self.initialize_sequences()

    def initialize_sequences(self):
        n = int(np.floor(self.length / self.planned_sequence_class.length))
        self.sequences = []
        i_0 = 0
        i_1 = 0
        for i in range(n):
            sub_timestamps = []
            for timestamp in self.timestamps:
                i_0 = np.searchsorted(timestamp, i * self.planned_sequence_class.length)
                i_1 = np.searchsorted(timestamp, (i + 1) * self.planned_sequence_class.length)
                sub_timestamps.append(timestamp[i_0:i_1])
            self.sequences.append(self.planned_sequence_class(sub_timestamps))

    def analyze(self):
        for sequence in self.sequences:
            sequence.bin()
        for i, sequence in enumerate(self.sequences):
            sequence.calculate_result_vector(others=self.sequences, index=i)
        self.x_vals = self.sequences[0].x_vals
        self.num_data_points = np.sum([sequence.num_data_points for sequence in self.sequences], axis=0)
        self.sum_results = np.sum([sequence.result_vector for sequence in self.sequences], axis=0)


class PlannedSequence:
    length = int(2e6)  # [ns]

    def __init__(self, timestamps):
        self.timestamps = timestamps

    def bin(self):
        return

    def calculate_result_vector(self, others=None, index=None):
        # should also calculate x_vals
        return


class Experiment:
    def __init__(self, experiment_data, cycle_class=Cycle, planned_sequence_class=PlannedSequence):
        self.experiment_data = experiment_data  # an array of size [N_SPDs, N_Cycles] of timetag arrays
        self.planned_sequence_class = planned_sequence_class
        self.cycle_class = cycle_class
        self.convert_to_cycle_objects()

    def convert_to_cycle_objects(self):
        self.cycles = [
            self.cycle_class(
                timestamps=timestamps,
                planned_sequence_class=self.planned_sequence_class
            )
            for timestamps in tqdm(zip(*self.experiment_data), desc='Converting to cycle objects', total=len(self.experiment_data[0]))
        ]

    def analyze(self):
        for cycle in tqdm(self.cycles, desc='Analyzing experiments'):
            cycle.analyze()
        self.x_vals = self.cycles[0].x_vals
        self.num_data_points = np.sum([experiment.num_data_points for experiment in self.cycles], axis=0)
        self.sum_results = np.sum([experiment.sum_results for experiment in self.cycles], axis=0)
        # self.meaned_results = self.sum_results / self.num_data_points
