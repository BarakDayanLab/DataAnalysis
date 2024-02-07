import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Generalization_attempt.util.Experiment import PlannedSequence, Experiment, Cycle

# from sortednp import merge

experiment_data_file_paths = [
    'Data/223709_Photon_TimeTags/output/Bright(1,2)/Bright_timetags.npz',
    'Data/223709_Photon_TimeTags/output/Dark(3,4)/Dark_timetags.npz',
    'Data/223709_Photon_TimeTags/output/North(8)/North_timetags.npz',
    'Data/223709_Photon_TimeTags/output/FastSwitch(6,7)/FS_timetags.npz',
    'Data/223709_Photon_TimeTags/output/South(5)/South_timetags.npz',
]
experiment_datum_names = ['arr_0'] * len(experiment_data_file_paths)
pulses_location_data_file_path = 'Data/223709_Photon_TimeTags/input/sequences/Pulses_location_in_seq.npz'
pulses_location_datum_name = 'arr_0'
north_sequence_data_file_path = 'Data/223709_Photon_TimeTags/input/sequences/North_sequence_vector.npz'
north_sequence_datum_name = 'arr_0'

SEQUENCE_LENGTH = len(np.load(north_sequence_data_file_path)[north_sequence_datum_name])
CYCLE_LENGTH = int((8e6 // SEQUENCE_LENGTH) * SEQUENCE_LENGTH)  # 8 [ms]
SEQUENCE_RULE = np.load(pulses_location_data_file_path)[pulses_location_datum_name]

class SprintCycle(Cycle):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs, length=CYCLE_LENGTH)

class SprintPlannedSequence(PlannedSequence):
    def bin(self):
        pass

# experiment_data = [
#     np.load(experiment_data_file_path, allow_pickle=True)[experiment_datum_name]
#     for experiment_data_file_path, experiment_datum_name in zip(experiment_data_file_paths, experiment_datum_names)
# ]
# experiment_data = [
#     merge(experiment_data[:3]),  # North
#     merge(experiment_data[3:]),  # South
# ]
