import Experiment_data_load
import numpy as np
import os
import ntpath
import sys
import json
import math
import pymsgbox
import pathlib
from tkinter import *
from tkinter import messagebox
from tkinter.filedialog import askdirectory
from Utilities.Utils import Utils
from Utilities.BDLogger import BDLogger
from Utilities.BDBatch import BDBatch
from Utilities.BDResults import BDResults
from tqdm import tqdm
import matplotlib
from matplotlib.container import ErrorbarContainer
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.widgets import RangeSlider
from matplotlib.widgets import TextBox
from matplotlib.widgets import Button
import numpy as np
import scipy.stats as stat
from scipy.optimize import fsolve

class experiment_data_analysis:
    def popupbox_inquiry(self, message=''):

        if len(message) > 0:
            pymsgbox.alert(text=message, title='Error!')

        experiment = pymsgbox.prompt(text='Please enter the type of experiment for which you want to load the '
                                          'data\nCapital letters!',
                                     title='Experiment data loading', default='QRAM')
        while True:
            if experiment is None:
                break
            elif not experiment.isupper():
                experiment = pymsgbox.prompt(text='Please enter the type of experiment for which you want to load the '
                                                  'data\nall characters must be Capital letters!',
                                             title='Experiment data loading', default='QRAM')
            else:
                break

        if experiment is None:
            return None, None, None

        date = pymsgbox.prompt(text='Please enter the date of the experiment', title='Experiment data loading',
                               default='YYYYMMDD')
        while True:
            date = date.strip()
            if date is None:
                break
            elif not date.isnumeric():
                date = pymsgbox.prompt(
                    text='Please enter the date of the experiment \nall characters must be integers!',
                    title='Experiment data loading', default='For example: 20230719')
            elif len(date) != 8:
                date = pymsgbox.prompt(text='Please enter the date of the experiment \nmust be 8 characters!',
                                       title='Experiment data loading', default='For example: 20230719')
            else:
                break

        if date is None:
            return None, None, None

        time = pymsgbox.prompt(text='Please enter the time of the experiment', title='Experiment data loading',
                               default='HHMMSS')
        while True:
            time = time.strip()
            if time is None:
                break
            elif not time.isnumeric():
                time = pymsgbox.prompt(
                    text='Please enter the time of the experiment \nall characters must be integers!',
                    title='Experiment data loading', default='For example: 115944')
            elif len(time) != 6:
                time = pymsgbox.prompt(text='Please enter the time of the experiment \nmust be 6 characters!',
                                       title='Experiment data loading', default='For example: 115944')
            else:
                break

        if time is None:
            return None, None, None

        return experiment, date, time

    def open_folder_to_dictionary(self):

        # Open folder and load to dictionary
        root = Tk()
        self.exp_data_path = '{}'.format(askdirectory(title='Experiment folder', initialdir=r'U:\Lab_2023\Experiment_results'))
        # self.exp_data_path.replace('/', os.sep)
        data = Experiment_data_load.DictionaryBuilder()
        dictionary = data.load_files_to_dict(self.exp_data_path)
        messagebox.showinfo(title='Success!', message='Experiment data is ready to use.')
        root.destroy()

        return dictionary

    def number_of_pulses_per_seq(self):
        '''

        :return:
        '''
        detection_pulses_per_seq = 0
        SPRINT_pulses_per_seq = 0
        for tup in self.Pulses_location_in_seq:
            if (tup[2] == 'N') or (tup[2] == 'S'):
                detection_pulses_per_seq += 1
            elif (tup[2] == 'n') or (tup[2] == 's'):
                SPRINT_pulses_per_seq += 1
        return detection_pulses_per_seq, SPRINT_pulses_per_seq

    def get_pulses_location_in_seq(self, delay, seq, smearing):
        '''
        :description
        A function that takes a pulse sequence, (a) applies a delay to it, (b) detects the pulses inside and
        (c) widens them pulses by the "smearing" factor, to build a filter (done to match the filter to the performance
        of the physical system, e.g. AOMs).

        :param delay: Between the actual sequence pulse location to the location of the folded data
        :param seq: The sequence of pulses from which the filter is generated.
        :param smearing: The value that is added to the filter before and after each pulse in the sequence.

        :return: the pulses location in the sequence (indices tuples) and the filter signal itself
        '''

        # Create an array consisting of zeros/ones that acts as filter, based on sequence values.
        seq_filter = (np.array(seq) > 0).astype(int)
        # Shift all values ("roll") to the right, based on the delay parameter
        seq_filter = np.roll(seq_filter, delay)
        # Find signal boundries indices - comparing each element with its neighbor
        indices = np.where(seq_filter[:-1] != seq_filter[1:])[0] + 1
        # Create the pulses boundry tuples - with the smearing widening
        pulses_loc = [(indices[i] - smearing, indices[i + 1] + smearing) for i in range(0, len(indices), 2)]
        # Recreate the signal by filling it with 1's in the relevant places
        seq_filter_with_smearing = np.zeros(seq_filter.shape[0])
        for (start, end) in pulses_loc: np.put_along_axis(seq_filter_with_smearing, np.arange(start, end), 1, axis=0)
        return pulses_loc, seq_filter_with_smearing

    def init_params_for_experiment(self, dict):

        self.sequence_len = len(dict['input']['sequences']['South_sequence_vector'])
        config_values_key = list(key for key in dict['meta_data'].keys() if 'config' in key)[0]
        self.M_window = dict['meta_data'][config_values_key]['M_window']
        self.M_time = self.M_window
        self.Pulses_location_in_seq = dict['input']['sequences']['Pulses_location_in_seq'].tolist()
        for index, pulse in enumerate(self.Pulses_location_in_seq):
            self.Pulses_location_in_seq[index][0] = int(pulse[0])
            self.Pulses_location_in_seq[index][1] = int(pulse[1])
        self.MZ_delay = 250  # [ns]

        self.number_of_sequences = math.ceil(self.M_window / self.sequence_len)

        self.number_of_detection_pulses_per_seq, self.number_of_SPRINT_pulses_per_seq = self.number_of_pulses_per_seq()

        self.end_of_det_pulse_in_seq = int(
            self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq - 1][1]) + 6

        # ----------------------------------------------------------
        # Prepare pulses location of South, North and Ancilla
        # ----------------------------------------------------------

        self.filter_delay = [0, 0]
        self.pulses_location_in_seq_S, self.filter_S = self.get_pulses_location_in_seq(self.filter_delay[0],
                                                                                       dict['input']['sequences']['South_sequence_vector'],
                                                                                       smearing=0)  # smearing=int(Config.num_between_zeros/2))
        self.pulses_location_in_seq_N, self.filter_N = self.get_pulses_location_in_seq(self.filter_delay[1],
                                                                                       dict['input']['sequences']['North_sequence_vector'],
                                                                                       smearing=0)  # smearing=int(Config.num_between_zeros/2))

        # define empty variables

        self.tt_measure = []
        self.tt_S_measure = []

        self.folded_transmission = np.zeros(self.sequence_len)
        self.folded_reflection = np.zeros(self.sequence_len)

        # Initiate folded cumulative tt "N", "S", BP, DP, FS
        self.folded_tt_S_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_N_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_BP_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_DP_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_FS_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_S_directional_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_N_directional_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_BP_timebins_cumulative_avg = np.zeros(self.sequence_len, dtype=int)
        self.folded_tt_DP_timebins_cumulative_avg = np.zeros(self.sequence_len, dtype=int)

        self.tt_S_binning = np.zeros(self.number_of_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sequence_len)
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_sequences)
        self.num_of_det_transmissions_per_seq_accumulated = np.zeros(self.number_of_sequences)

        num_of_seq_per_count = 50
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_sequences // num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_sequences // num_of_seq_per_count)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_sequences // num_of_seq_per_count)

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_sequences)

        self.all_transits_seq_indx_per_cond, \
        self.all_transits_length_per_cond = [
            [
                [] for _ in range(len(self.transit_conditions))
            ] for _ in range(2)
        ]

        self.seq_with_data_points, \
        self.reflection_SPRINT_data, \
        self.transmission_SPRINT_data, \
        self.BP_counts_SPRINT_data, \
        self.DP_counts_SPRINT_data = [
            [
                [] for _ in range(len(self.transit_conditions))
            ] for _ in range(5)
        ]

        self.reflection_SPRINT_data_without_transits, \
        self.transmission_SPRINT_data_without_transits, \
        self.BP_counts_SPRINT_data_without_transits, \
        self.DP_counts_SPRINT_data_without_transits = [
            [
                [] for _ in range(len(self.transit_conditions))
            ] for _ in range(4)
        ]

    def ingest_time_tags(self, dict, cycle, delay):
        """
        Takes all the raw results we got from the streams and does some processing on them - preparing "measures"
        """

        # ------------------------------------------
        # Unify detectors and windows within detectors and create vector of tt's for each direction (bright port, dark port, north, south and from FS) and sort them
        # ------------------------------------------
        self.tt_BP_measure = np.array(dict['output']['Bright(1,2)']['Bright_timetags'][cycle]) + delay
        self.tt_DP_measure = np.array(dict['output']['Dark(3,4)']['Dark_timetags'][cycle]) + delay
        self.tt_N_measure = np.array(dict['output']['North(8)']['North_timetags'][cycle]) + delay
        self.tt_S_measure = np.array(dict['output']['South(5)']['South_timetags'][cycle]) + delay
        self.tt_FS_measure = np.array(dict['output']['FastSwitch(6,7)']['FS_timetags'][cycle]) + delay

    def fold_tt_histogram(self, exp_sequence_len):

        self.folded_tt_S = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_FS = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_S_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP_timebins = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP_timebins = np.zeros(exp_sequence_len, dtype=int)

        # for x in self.Exp_dict['output']['South(5)']['South_timetags'][cycle]:
        for x in self.tt_S_measure:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_S[x % exp_sequence_len] += 1
        # for x in self.Exp_dict['output']['North(8)']['North_timetags'][cycle]:
        for x in self.tt_N_measure:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_N[x % exp_sequence_len] += 1
        # for x in self.Exp_dict['output']['Bright(1,2)']['Bright_timetags'][cycle]:
        for x in self.tt_BP_measure:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_BP[x % exp_sequence_len] += 1
        # for x in self.Exp_dict['output']['Dark(3,4)']['Dark_timetags'][cycle]:
        for x in self.tt_DP_measure:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_DP[x % exp_sequence_len] += 1
        # for x in self.Exp_dict['output']['FastSwitch(6,7)']['FS_timetags'][cycle]:
        for x in self.tt_FS_measure:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_FS[x % exp_sequence_len] += 1

        self.folded_tt_S_directional = (np.array(self.folded_tt_S) + np.array(self.folded_tt_FS))
        # self.folded_tt_N_directional = self.folded_tt_N
        # TODO: ask dor -  switched folded_tt_N with folded_tt_N_directional
        self.folded_tt_N_directional[:self.end_of_det_pulse_in_seq] = \
            np.array(self.folded_tt_N[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_BP[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_DP[:self.end_of_det_pulse_in_seq])

        self.folded_tt_BP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_BP[self.end_of_det_pulse_in_seq:]
        self.folded_tt_DP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_DP[self.end_of_det_pulse_in_seq:]
        # if self.pulses_location_in_seq_A or ((Config.sprint_pulse_amp_N[0] > 0) & (len(Config.sprint_pulse_amp_N) > 1)):
        #     self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #         (np.array(
        #             self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
        #          + (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_BP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          + (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_DP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          )
        #     self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #         (np.array(
        #             self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
        #          - (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_BP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          )
        #     self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #         (np.array(
        #             self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
        #          - (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_DP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          )

    def fold_all_cycle_tt_histogram(self, exp_sequence_len):

        self.folded_tt_S = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_FS = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_S_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_N_directional = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_BP_timebins = np.zeros(exp_sequence_len, dtype=int)
        self.folded_tt_DP_timebins = np.zeros(exp_sequence_len, dtype=int)

        for x in [elem for lst in self.Exp_dict['output']['South(5)']['South_timetags'] for elem in lst]:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_S[x % exp_sequence_len] += 1
            # self.folded_tt_S[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['North(8)']['North_timetags'] for elem in lst]:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_N[x % exp_sequence_len] += 1
            # self.folded_tt_N[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['Bright(1,2)']['Bright_timetags'] for elem in lst]:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_BP[x % exp_sequence_len] += 1
            # self.folded_tt_BP[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['Dark(3,4)']['Dark_timetags'] for elem in lst]:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_DP[x % exp_sequence_len] += 1
            # self.folded_tt_DP[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['FastSwitch(6,7)']['FS_timetags'] for elem in lst]:
            if (x > int(0.6e6)) and (x < int(self.M_time - 0.4e6)):
                self.folded_tt_FS[x % exp_sequence_len] += 1
            # self.folded_tt_FS[x % exp_sequence_len] += 1

        self.folded_tt_S_directional = (np.array(self.folded_tt_S) + np.array(self.folded_tt_FS))
        # self.folded_tt_N_directional = self.folded_tt_N
        # TODO: ask dor -  switched folded_tt_N with folded_tt_N_directional
        self.folded_tt_N_directional[:self.end_of_det_pulse_in_seq] = \
            np.array(self.folded_tt_N[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_BP[:self.end_of_det_pulse_in_seq]) \
            + np.array(self.folded_tt_DP[:self.end_of_det_pulse_in_seq])

        self.folded_tt_BP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_BP[self.end_of_det_pulse_in_seq:]
        self.folded_tt_DP_timebins[self.end_of_det_pulse_in_seq:] = self.folded_tt_DP[self.end_of_det_pulse_in_seq:]
        # if self.pulses_location_in_seq_A or ((Config.sprint_pulse_amp_N[0] > 0) & (len(Config.sprint_pulse_amp_N) > 1)):
        #     self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #         (np.array(
        #             self.folded_tt_N_directional[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
        #          + (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_BP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          + (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_DP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          )
        #     self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #         (np.array(
        #             self.folded_tt_BP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
        #          - (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_BP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          )
        #     self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]] = \
        #         (np.array(
        #             self.folded_tt_DP_timebins[self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]])
        #          - (Config.sprint_pulse_amp_Early[1]
        #             * np.array(self.folded_tt_DP[
        #                        self.pulses_location_in_seq[-2][0]:self.pulses_location_in_seq[-2][1]]))
        #          )

    def divide_tt_to_reflection_trans_extended(self):
        '''
        A function designed to count the number of photons reflected or transmitted for each sequence, such that,
        for the detection pulses the number of photons will be accumulated for each sequence and for the SPRINT
        pulses there will be number of reflected or transmitted photons for each SPRINT pulse.
        :param sprint_pulse_len: the length in [ns] of the SPRINT pulses in the sequence.
        :param num_of_det_pulses: the number of detection pulses in the sequence.
        :return:
        '''

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_sequences)

        self.num_of_det_reflections_per_seq_S_, \
        self.num_of_det_reflections_per_seq_N_, \
        self.num_of_det_reflections_per_seq_, \
        self.num_of_det_transmissions_per_seq_S_, \
        self.num_of_det_transmissions_per_seq_N_, \
        self.num_of_det_transmissions_per_seq_ = [
            [
                [
                    [] for _ in range(self.number_of_detection_pulses_per_seq)
                ]
                for _ in range(self.number_of_sequences)
            ] for _ in range(6)
        ]

        self.num_of_SPRINT_reflections_per_seq_S_, \
        self.num_of_SPRINT_reflections_per_seq_N_, \
        self.num_of_SPRINT_reflections_per_seq_, \
        self.num_of_SPRINT_transmissions_per_seq_S_, \
        self.num_of_SPRINT_transmissions_per_seq_N_, \
        self.num_of_SPRINT_transmissions_per_seq_, \
        self.num_of_BP_counts_per_seq_in_SPRINT_pulse, \
        self.num_of_DP_counts_per_seq_in_SPRINT_pulse = [
            [
                [
                    [] for _ in range(self.number_of_SPRINT_pulses_per_seq)
                ]
                for _ in range(self.number_of_sequences)
            ] for _ in range(8)
        ]

        # N direction:
        # self.N_tt = np.array(self.tt_N_measure)
        self.N_tt = self.tt_N_measure
        tt_inseq_ = self.N_tt % self.sequence_len
        seq_num_ = self.N_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.N_tt, tt_inseq_, seq_num_):
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?
            #  TODO: Q: what is the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_window - 0.4e6)):
                for indx, tup in enumerate(self.Pulses_location_in_seq):
                    # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
                    # if (Config.PNSA_Exp_Square_samples_Early[tup[0]] == 0) or \
                    #         (self.number_of_SPRINT_pulses_per_seq > 1):
                    start_pulse_time_in_seq = tup[0]
                    end_pulse_time_in_seq = tup[1]
                    # else:
                    #     start_pulse_time_in_seq = tup[0] + self.MZ_delay
                    #     end_pulse_time_in_seq = tup[1] + self.MZ_delay
                    if (tt_inseq >= start_pulse_time_in_seq) and (tt_inseq <= end_pulse_time_in_seq):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(int(element))
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(int(element))
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(int(element))
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(int(element))
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(int(element))
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(int(element))

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # Bright port direction:
        # self.BP_tt = np.array(self.tt_BP_measure)
        self.BP_tt = self.tt_BP_measure
        tt_inseq_ = self.BP_tt % self.sequence_len
        seq_num_ = self.BP_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.BP_tt, tt_inseq_, seq_num_):
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?
            #  TODO: Q: what is the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_window - 0.4e6)):
                for indx, tup in enumerate(self.Pulses_location_in_seq):
                    # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
                    if (self.Exp_dict['input']['sequences']['Early_sequence_vector'][tup[0]] == 0) or \
                            (self.number_of_SPRINT_pulses_per_seq > 1):
                        start_pulse_time_in_seq = tup[0]
                        end_pulse_time_in_seq = tup[1]
                    else:
                        start_pulse_time_in_seq = tup[0] + self.MZ_delay
                        end_pulse_time_in_seq = tup[1] + self.MZ_delay
                    if (tt_inseq >= start_pulse_time_in_seq) and (tt_inseq <= end_pulse_time_in_seq):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(int(element))
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(int(element))
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(int(element))
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(int(element))
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(int(element))
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(int(element))
                            self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_num][ind].append(int(element))

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # Dark port direction:
        # self.DP_tt = np.array(self.tt_DP_measure)
        self.DP_tt = self.tt_DP_measure
        tt_inseq_ = self.DP_tt % self.sequence_len
        seq_num_ = self.DP_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.DP_tt, tt_inseq_, seq_num_):
            # TODO: Q: I assume this is a time-window to ignore some time on start/end, how did we decide on this time?
            # TODO: Q: what is the meaning of this specific time-window?
            if (element > int(0.6e6)) and (element < int(self.M_window - 0.4e6)):
                for indx, tup in enumerate(self.Pulses_location_in_seq):
                    # if all(x == 0 for x in Config.PNSA_Exp_Square_samples_Early[tup[0]:tup[1]]):
                    if (self.Exp_dict['input']['sequences']['Early_sequence_vector'][tup[0]] == 0) or \
                            (self.number_of_SPRINT_pulses_per_seq > 1):
                        start_pulse_time_in_seq = tup[0]
                        end_pulse_time_in_seq = tup[1]
                    else:
                        start_pulse_time_in_seq = tup[0] + self.MZ_delay
                        end_pulse_time_in_seq = tup[1] + self.MZ_delay
                    if (tt_inseq >= start_pulse_time_in_seq) and (tt_inseq <= end_pulse_time_in_seq):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(int(element))
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(int(element))
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(int(element))
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(int(element))
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(int(element))
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(int(element))
                            self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_num][ind].append(int(element))

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # self.S_tt = np.array(self.tt_S_measure + self.tt_FS_measure)
        self.S_tt = np.concatenate((self.tt_S_measure, self.tt_FS_measure))
        tt_inseq_ = self.S_tt % self.sequence_len
        seq_num_ = self.S_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.S_tt, tt_inseq_, seq_num_):
            if (element > int(0.6e6)) and (element < int(self.M_window - 0.4e6)):
                for indx, tup in enumerate(self.Pulses_location_in_seq):
                    if (tt_inseq >= tup[0]) and (tt_inseq <= tup[1]):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_reflections_per_seq_N_[seq_num][indx].append(int(element))
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(int(element))
                            if tup[2] == 'S':
                                self.num_of_det_transmissions_per_seq_S_[seq_num][indx].append(int(element))
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(int(element))
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_reflections_per_seq_N_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(int(element))
                            if tup[2] == 's':
                                self.num_of_SPRINT_transmissions_per_seq_S_[seq_num][ind].append(int(element))
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(int(element))

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.sequence_len
                    self.num_of_det_reflections_per_seq_N[seq_num] += self.filter_N[tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[seq_num] += self.filter_S[tt_inseq]

    def experiment_calculations(self):

        self.divide_tt_to_reflection_trans_extended()

        self.num_of_det_reflections_per_seq = self.num_of_det_reflections_per_seq_S \
                                              + self.num_of_det_reflections_per_seq_N
        self.num_of_det_transmissions_per_seq = self.num_of_det_transmissions_per_seq_S \
                                                + self.num_of_det_transmissions_per_seq_N

        pass

    def calculate_running_averages(self, cycle_number):
        self.folded_tt_S_cumulative_avg = Utils.running_average(self.folded_tt_S_cumulative_avg, self.folded_tt_S, cycle_number)
        self.folded_tt_N_cumulative_avg = Utils.running_average(self.folded_tt_N_cumulative_avg, self.folded_tt_N, cycle_number)
        self.folded_tt_BP_cumulative_avg = Utils.running_average(self.folded_tt_BP_cumulative_avg, self.folded_tt_BP, cycle_number)
        self.folded_tt_DP_cumulative_avg = Utils.running_average(self.folded_tt_DP_cumulative_avg, self.folded_tt_DP, cycle_number)
        self.folded_tt_FS_cumulative_avg = Utils.running_average(self.folded_tt_FS_cumulative_avg, self.folded_tt_FS, cycle_number)
        self.folded_tt_S_directional_cumulative_avg = Utils.running_average(self.folded_tt_S_directional_cumulative_avg, self.folded_tt_S_directional, cycle_number)
        self.folded_tt_N_directional_cumulative_avg = Utils.running_average(self.folded_tt_N_directional_cumulative_avg, self.folded_tt_N_directional, cycle_number)
        self.folded_tt_BP_timebins_cumulative_avg = Utils.running_average(self.folded_tt_BP_timebins_cumulative_avg, self.folded_tt_BP_timebins, cycle_number)
        self.folded_tt_DP_timebins_cumulative_avg = Utils.running_average(self.folded_tt_DP_timebins_cumulative_avg, self.folded_tt_DP_timebins, cycle_number)
        pass

    def find_transit_events(self, cond_list=None, minimum_number_of_seq_detected=2):
        '''
        Find transits of atoms by searching events that satisfy the number of reflected photons per sequence required at
        each cond place with minimum number of conditions needed to be satisfied, defined by the minimum_number_of_seq_detected.
        For example:
        given cond=[2,1,2] and minimum_number_of_seq_detected=2, if in 3 consecutive sequences we get either 2-1-0 or
        0-2-1 or 2-0-2, the condition is satisfied and it is defined as a transit.
        :param cond: The condition that need to be met for number of reflections per detection pulses in sequence.
        :param minimum_number_of_seq_detected: The number of terms needed to be satisfied per cond vector.
        :return:
        '''

        if cond_list is None:
            cond_list = [[2, 2]]

        all_transits_seq_indx = []

        # Find all tansits for each of the conditions given in the list:
        for cond in cond_list:
            current_transit = []
            for i in range(self.number_of_sequences - len(cond) + 1):
                cond_check = (self.num_of_det_reflections_per_seq[i:(i + len(cond))] >= cond).astype(int)
                if sum(cond_check) >= minimum_number_of_seq_detected:
                    current_transit = np.unique(
                        current_transit + [*range(i + np.where(cond_check != 0)[0][0], (i + len(cond)))]).tolist()
                elif len(current_transit) > 1:
                    current_transit = current_transit[
                                      :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][
                                           -1] + 1]
                    if len(all_transits_seq_indx) > 0:
                        if bool(set(current_transit) & set(all_transits_seq_indx[-1])):
                            current_transit = all_transits_seq_indx[-1] + current_transit[1:]
                            all_transits_seq_indx = all_transits_seq_indx[:-1]
                    all_transits_seq_indx.append(current_transit)
                    current_transit = []
            if len(current_transit) > 1:
                current_transit = current_transit[
                                  :np.where(self.num_of_det_reflections_per_seq[current_transit] >= min(cond))[0][-1] + 1]
                if all_transits_seq_indx:
                    if bool(set(current_transit) & set(all_transits_seq_indx[-1])):
                        current_transit = all_transits_seq_indx[-1] + current_transit[1:]
                        all_transits_seq_indx = all_transits_seq_indx[:-1]
                all_transits_seq_indx.append(current_transit)

        # Remove duplicate values from the list of transits and merge consecutive transits:
        if len(all_transits_seq_indx) > 0:
            # remove duplicate lists:
            all_transits_seq_indx = sorted(list(set(tuple(transit) for transit in all_transits_seq_indx)), key=lambda x: x[0])
            all_transits_seq_indx = [sorted(list(x_)) for x_ in all_transits_seq_indx]
            # merge all lists that are contained or are consecutive to another list:
            for i in range(len(all_transits_seq_indx) - 1, 0, -1):
                if all_transits_seq_indx[i][0] <= all_transits_seq_indx[i - 1][-1] + 1:
                    all_transits_seq_indx[i - 1] = sorted(list(set(all_transits_seq_indx[i - 1] + all_transits_seq_indx[i])))
                    all_transits_seq_indx.pop(i)
        return all_transits_seq_indx

    def analyze_SPRINT_data_points(self, all_transits_seq_indx, SPRINT_pulse_number=[1],
                                   number_of_reflection_in_preparation_pulse=1, background=False):
        '''
        For a given vector of sequence indexes, find the relevant data points for SPRINT pulse and analyze results.
        :param all_transits_seq_indx: The vector of indexes of sequences that need to be analyzed.
        :param SPRINT_pulse_number: The SPRINT pulse number for which we want to check the results.
        :param background: If the check is for background data and not for actual transits, the "potential data" check,
                           for which we condition the relevance of the data if we get at least 1 reflection from the
                           last detection pulse, is irrelevant.
        '''
        reflection_SPRINT_data = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        transmission_SPRINT_data = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.
        seq_with_data_points = []
        BP_counts_SPRINT_data = []
        DP_counts_SPRINT_data = []

        if len(self.Pulses_location_in_seq) > 0:
            for transit in all_transits_seq_indx:
                for seq_indx in transit[:-1]:
                    # Checking potential for data point by looking for a single photon at the reflection of the last
                    # detection pulse:
                    potential_data = False if not background else True
                    if self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq-1][2] == 'N' and \
                       len(self.num_of_det_reflections_per_seq_N_[seq_indx][-1]) >= number_of_reflection_in_preparation_pulse:
                        potential_data = True
                        # self.potential_data += 1
                    elif self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq-1][2] == 'S' and \
                            len(self.num_of_det_reflections_per_seq_S_[seq_indx][-1]) >= number_of_reflection_in_preparation_pulse:
                        potential_data = True
                        # self.potential_data += 1
                    # Getting SPRINT data if the SPRINT pulse has one photon in the reflection or transmission
                    if potential_data and len(self.Pulses_location_in_seq) > self.number_of_detection_pulses_per_seq:
                        if self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq-1+SPRINT_pulse_number[0]][2] == 'n':
                            transmissions = 0
                            reflections = 0
                            for sprint_pulse in SPRINT_pulse_number:
                                transmissions += len(self.num_of_SPRINT_transmissions_per_seq_N_[seq_indx][sprint_pulse-1])
                                reflections += len(self.num_of_SPRINT_reflections_per_seq_N_[seq_indx][sprint_pulse-1])
                            if (transmissions + reflections) == 1:
                                seq_with_data_points.append(seq_indx)
                                reflection_SPRINT_data.append(reflections)
                                transmission_SPRINT_data.append(transmissions)
                            if (transmissions + reflections) > 0:
                                BP_counts_SPRINT_data.append(len(self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
                                DP_counts_SPRINT_data.append(len(self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
                        elif self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq-1+SPRINT_pulse_number[0]][2] == 's':
                            transmissions = 0
                            reflections = 0
                            for sprint_pulse in SPRINT_pulse_number:
                                transmissions += len(self.num_of_SPRINT_transmissions_per_seq_S_[seq_indx][sprint_pulse-1])
                                reflections += len(self.num_of_SPRINT_reflections_per_seq_S_[seq_indx][sprint_pulse-1])
                            if (transmissions + reflections) == 1:
                                seq_with_data_points.append(seq_indx)
                                reflection_SPRINT_data.append(reflections)
                                transmission_SPRINT_data.append(transmissions)
                            if (transmissions + reflections) > 0:
                                BP_counts_SPRINT_data.append(len(self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
                                DP_counts_SPRINT_data.append(len(self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_indx][SPRINT_pulse_number[-1]-1]))
        return seq_with_data_points, reflection_SPRINT_data, transmission_SPRINT_data, BP_counts_SPRINT_data, DP_counts_SPRINT_data

    def get_transit_data(self, transit_condition, cond_num, reflections_in_preparation_cond):
        '''
        Use the data loaded from the experiment to find atomic transits under different conditions and some data to
        compare the conditions, such as, number of transits, transit length distribution and more...
        :param transit_condition: List of lists, so that each condition is of the structure [[1,2],[2,1],....]
               and so on...
        :return:
        '''

        self.all_transits_seq_indx_per_cond[cond_num] = self.find_transit_events(cond_list=transit_condition, minimum_number_of_seq_detected=2)
        self.all_transits_length_per_cond[cond_num] = self.get_transits_length(self.all_transits_seq_indx_per_cond[cond_num])

        # Analyze SPRINT data during transits:
        (self.seq_with_data_points[cond_num], self.reflection_SPRINT_data[cond_num], self.transmission_SPRINT_data[cond_num],
         self.BP_counts_SPRINT_data[cond_num], self.DP_counts_SPRINT_data[cond_num]) = \
            self.analyze_SPRINT_data_points(self.all_transits_seq_indx_per_cond[cond_num], SPRINT_pulse_number=[1, 2],
                                            number_of_reflection_in_preparation_pulse=reflections_in_preparation_cond,
                                            background=False)  # Enter the index of the SPRINT pulse for which the data should be analyzed
        # print(self.potential_data)
        # Analyze SPRINT data when no transit occur:
        self.all_seq_without_transits = [
            np.delete(np.arange(0, self.number_of_sequences, 1, dtype='int'),
                      sum(self.all_transits_seq_indx_per_cond[cond_num], [])).tolist()
        ]
        (_, self.reflection_SPRINT_data_without_transits[cond_num], self.transmission_SPRINT_data_without_transits[cond_num],
         self.BP_counts_SPRINT_data_without_transits[cond_num], self.DP_counts_SPRINT_data_without_transits[cond_num]) = \
            self.analyze_SPRINT_data_points(self.all_seq_without_transits, SPRINT_pulse_number=[1, 2],
                                            number_of_reflection_in_preparation_pulse=reflections_in_preparation_cond,
                                            background=True)  # Enter the index of the SPRINT pulse for which the data should be analyzed

        # self.number_of_transits_live = len(self.all_transits_seq_indx_per_cond)
        # self.number_of_transits_total = len(
        #     [vec for lst in self.batcher['all_transits_seq_indx_batch'] for vec in lst])
        #
        # self.num_of_total_SPRINT_reflections = sum(self.reflection_SPRINT_data_without_transits)
        # self.num_of_total_SPRINT_transmissions = sum(self.transmission_SPRINT_data_without_transits)
        # self.num_of_total_SPRINT_BP_counts = sum(self.BP_counts_SPRINT_data_without_transits)
        # self.num_of_total_SPRINT_DP_counts = sum(self.DP_counts_SPRINT_data_without_transits)

    def results_to_dictionary(self, is_exp=False):
        '''
        Take all the results from the batcher and summarize them into a dictionary.
        :return:
        '''

        dictionary = dict()

        ## General data for background experiment ##
        # folded data:
        dictionary['folded_tt_N'] = self.folded_tt_N
        dictionary['folded_tt_BP'] = self.folded_tt_BP
        dictionary['folded_tt_DP'] = self.folded_tt_DP
        dictionary['folded_tt_S'] = self.folded_tt_S
        dictionary['folded_tt_FS'] = self.folded_tt_FS

        dictionary['filter_N'] = self.filter_N
        dictionary['filter_S'] = self.filter_S
        dictionary['Pulses_location_in_seq'] = self.Pulses_location_in_seq
        # if is_exp:
        #     config_values_key = list(key for key in self.Exp_dict['meta_data'].keys() if 'config' in key)[0]
        #     dictionary['exp_config_values'] = self.Exp_dict['meta_data'][config_values_key]
        # else:
        #     config_values_key = list(key for key in self.Background_dict['meta_data'].keys() if 'config' in key)[0]
        #     dictionary['exp_config_values'] = self.Background_dict['meta_data'][config_values_key]
        # #
        # # transmission data in detection pulses in sequences:
        # dictionary['transmission_data_in_detection_pulses_per_seq_per_cycle'] = self.batcher[
        #      'num_of_det_transmissions_per_seq']
        # dictionary['transmission_data_in_detection_pulses_per_seq_per_cycle_N'] = self.batcher[
        #      'num_of_det_transmissions_per_seq_N']
        # dictionary['transmission_data_in_detection_pulses_per_seq_per_cycle_S'] = self.batcher[
        #      'num_of_det_transmissions_per_seq_S']
        #
        # reflection data in detection pulses in sequences:
        # dictionary['reflection_data_in_detection_pulses_per_seq_per_cycle'] = self.batcher[
        #      'num_of_det_reflections_per_seq']
        # dictionary['reflection_data_in_detection_pulses_per_seq_per_cycle_N'] = self.batcher[
        #      'num_of_det_reflections_per_seq_N']
        # dictionary['reflection_data_in_detection_pulses_per_seq_per_cycle_S'] = self.batcher[
        #      'num_of_det_reflections_per_seq_S']
        #
        # # transmission data in SPRINT pulses in sequences:
        # dictionary['transmission_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
        #     'num_of_SPRINT_transmissions_per_seq']
        # dictionary['transmission_data_in_SPRINT_pulses_per_seq_per_cycle_N'] = self.batcher[
        #     'num_of_SPRINT_transmissions_per_seq_N']
        # dictionary['transmission_data_in_SPRINT_pulses_per_seq_per_cycle_S'] = self.batcher[
        #     'num_of_SPRINT_transmissions_per_seq_S']
        #
        # # reflection data in SPRINT pulses in sequences:
        # dictionary['reflection_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
        #     'num_of_SPRINT_reflections_per_seq']
        # dictionary['reflection_data_in_SPRINT_pulses_per_seq_per_cycle_N'] = self.batcher[
        #     'num_of_SPRINT_reflections_per_seq_N']
        # dictionary['reflection_data_in_SPRINT_pulses_per_seq_per_cycle_S'] = self.batcher[
        #     'num_of_SPRINT_reflections_per_seq_S']
        #
        # # Bright port data in SPRINT pulses in sequence:
        # dictionary['bright_port_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
        #     'num_of_BP_counts_per_seq_in_SPRINT_pulse']
        #
        # # Dark port data in SPRINT pulses in sequence:
        # dictionary['dark_port_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
        #     'num_of_DP_counts_per_seq_in_SPRINT_pulse']

        ## Sorting and summarizing the data analysis for each condition ##
        dictionary['SNR'] = []
        for indx, condition in enumerate(self.transit_conditions):
            dictionary[str(condition)] = dict()
            # Handle transits related data:
            dictionary[str(condition)]['number_of_transits'] = 0
            dictionary[str(condition)]['number_of_sequences_with_transits'] = 0
            dictionary[str(condition)]['indices_of_all_sequences_with_transits_per_cycle'] = []
            dictionary[str(condition)]['all_transits_length_per_cond'] = []
            for i, lst in enumerate(self.batcher['all_transits_seq_indx_per_cond']):
                for vec in lst[indx]:
                    # Get number of transits and number of sequences with transits
                    dictionary[str(condition)]['number_of_transits'] += 1
                    dictionary[str(condition)]['number_of_sequences_with_transits'] += len(vec)
                # Save index of sequences with transits per cycle:
                dictionary[str(condition)]['indices_of_all_sequences_with_transits_per_cycle'].append(lst[indx])
                dictionary[str(condition)]['all_transits_length_per_cond'].append(self.batcher['all_transits_length_per_cond'][i][indx])
            dictionary[str(condition)]['number_of_sequences_with_transits'] -= dictionary[str(condition)]['number_of_transits']
            if is_exp:
                if self.Background_dict['Analysis_results'][str(condition)]['number_of_transits'] == 0:
                    dictionary[str(condition)]['SNR'] = np.nan
                else:
                    dictionary[str(condition)]['SNR'] = ((dictionary[str(condition)]['number_of_transits'] / self.number_of_cycles) /
                                                         (self.Background_dict['Analysis_results'][str(condition)]['number_of_transits'] /
                                                          len(list(self.Background_dict['output'][list(self.Background_dict['output'].keys())[0]].values())[0])))
                dictionary['SNR'].append([str(condition), dictionary[str(condition)]['SNR']])

            # Handle SPRINT data results for transit events:

            # For each cycle, get a list of all sequence indices with data points, number of photon transmitted or
            # reflected per data point:
            dictionary[str(condition)]['sequence_indices_with_data_points_per_cycle'] = \
                [lst[indx] for lst in self.batcher['seq_with_data_points']]
            dictionary[str(condition)]['transmission_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['transmission_SPRINT_data']]
            dictionary[str(condition)]['reflection_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['reflection_SPRINT_data']]
            # Total number of Transmissions/reflections counts for the entire experiment length (All cycles combined):
            dictionary[str(condition)]['transmission'] = (
                sum(sum(dictionary[str(condition)]['transmission_per_data_point_per_cycle'], [])))
            dictionary[str(condition)]['reflection'] = (
                sum(sum(dictionary[str(condition)]['reflection_per_data_point_per_cycle'], [])))

            # For each cycle, get the number of photon in the Bright Port(BP) or Dark Port(DP) channel of the
            # Mach-Zehnder(MZ) per data point:
            dictionary[str(condition)]['BP_counts_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['BP_counts_SPRINT_data']]
            dictionary[str(condition)]['DP_counts_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['DP_counts_SPRINT_data']]
            # Total number of Transmissions/reflections counts for the entire experiment length (All cycles combined):
            dictionary[str(condition)]['Bright'] = (
                sum(sum(dictionary[str(condition)]['BP_counts_per_data_point_per_cycle'], [])))
            dictionary[str(condition)]['Dark'] = (
                sum(sum(dictionary[str(condition)]['DP_counts_per_data_point_per_cycle'], [])))

            # Handle SPRINT data results for non transit events:

            # For each cycle, get a list of the number of photons transmitted or reflected per data point:
            dictionary[str(condition)]['transmission_without_transits_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['transmission_SPRINT_data_without_transits']]
            dictionary[str(condition)]['reflection_without_transits_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['reflection_SPRINT_data_without_transits']]
            # Total number of Transmissions/reflections counts for the entire experiment length (All cycles combined):
            dictionary[str(condition)]['transmission_without_transits'] = (
                sum(sum(dictionary[str(condition)]['transmission_without_transits_per_data_point_per_cycle'], [])))
            dictionary[str(condition)]['reflection_without_transits'] = (
                sum(sum(dictionary[str(condition)]['reflection_without_transits_per_data_point_per_cycle'], [])))

            # For each cycle, get a list of the number of photons in the Bright Port(BP) or Dark Port(DP) channel of the
            # Mach-Zehnder(MZ) per data point:
            dictionary[str(condition)]['BP_counts_without_transits_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['BP_counts_SPRINT_data_without_transits']]
            dictionary[str(condition)]['DP_counts_without_transits_per_data_point_per_cycle'] = \
                [lst[indx] for lst in self.batcher['DP_counts_SPRINT_data_without_transits']]
            # Total number of Transmissions/reflections counts for the entire experiment length (All cycles combined):
            dictionary[str(condition)]['Bright_without_transits'] = (
                sum(sum(dictionary[str(condition)]['BP_counts_without_transits_per_data_point_per_cycle'], [])))
            dictionary[str(condition)]['Dark_without_transits'] = (
                sum(sum(dictionary[str(condition)]['DP_counts_without_transits_per_data_point_per_cycle'], [])))

        return dictionary

    def get_transits_length(self, all_transits_seq_indx):
        '''
        Build an array with the length of all the transits.
        :return:
        '''
        all_transits_length = []
        for transit in all_transits_seq_indx:
            transit_length = (max(sum(self.num_of_det_reflections_per_seq_[transit[-1]], [])) -
                              min(sum(self.num_of_det_reflections_per_seq_[transit[0]], [])))
            all_transits_length.append(transit_length)

        return all_transits_length

    def plot_results(self):
        '''

        :return:
        '''

        self.f1 = plt.figure('North, Bright and Dark ports with filter lines')
        plt.plot(self.folded_tt_N_directional_cumulative_avg, label='"N" detectors')
        plt.plot(self.folded_tt_BP_timebins_cumulative_avg, label='"BP" detectors')
        plt.plot(self.folded_tt_DP_timebins_cumulative_avg, label='"DP" detectors')
        plt.plot(self.filter_N * max(self.folded_tt_N_directional_cumulative_avg +
                                     self.folded_tt_S_directional_cumulative_avg), '--k', label='Filter "N"')
        plt.legend(loc='upper right')

        self.f2 = plt.figure('South and Fast Switch (FS) ports with filter lines')
        plt.plot(self.folded_tt_S_directional_cumulative_avg, label='"S" detectors')
        plt.plot(self.filter_S * max(self.folded_tt_N_directional_cumulative_avg +
                                     self.folded_tt_S_directional_cumulative_avg), '--k', label='Filter "S"')
        plt.legend(loc='upper right')

        self.transits_fig = []
        self.SPRINT_figure = []
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        for indx, cond in enumerate(self.transit_conditions):
            self.transit_subplots = []
            # create figure for transits information:
            self.transits_fig.append(plt.figure(str(cond) + ' - Transit information'))
            self.transit_subplots.append(plt.subplot2grid((2, 2), (0, 0), colspan=2, rowspan=1))
            self.transit_subplots.append(plt.subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1))

            # plot hist of transit durations:
            all_transit_durations = [elem for vec in self.batcher['all_transits_length_per_cond'] for elem in vec[indx]]
            self.transit_subplots[0].hist(all_transit_durations, range(1, max(all_transit_durations), 50))
            self.transit_subplots[0].set(xlabel='Duration [ns]', ylabel='Counts [# of transits]')

            # plot all sequences with transit events
            seq_with_transit_events = np.zeros(self.number_of_sequences)
            for cycle in range(self.number_of_cycles):
                seq_with_transit_events[[vec for elem in
                                         self.Exp_dict['Analysis_results'][str(cond)]['indices_of_all_sequences_with_transits_per_cycle'][cycle]
                                         for vec in elem]] += 1
            textstr_transit_event_counter = (r'$N_{Transits Total} = %s $' %
                                             (self.Exp_dict['Analysis_results'][str(cond)]['number_of_transits'],) +
                                             '[Counts]')
            self.transit_subplots[1].plot(range(self.number_of_sequences), seq_with_transit_events,
                                          label='Transit events accumulated')
            self.transit_subplots[1].set(xlabel='Sequence [#]', ylabel='Counts [Photons]')
            self.transit_subplots[1].text(0.02, 0.92, textstr_transit_event_counter,
                                          transform=self.transit_subplots[1].transAxes,
                                          fontsize=14, verticalalignment='top', bbox=props)


            # SPRINT results box
            reflections = self.Exp_dict['Analysis_results'][str(cond)]['reflection']
            transmissions = self.Exp_dict['Analysis_results'][str(cond)]['transmission']
            reflections_without_transits = self.Exp_dict['Analysis_results'][str(cond)]['reflection_without_transits']
            transmissions_without_transits = self.Exp_dict['Analysis_results'][str(cond)]['transmission_without_transits']
            if (reflections_without_transits + transmissions_without_transits) > 0:
                SPRINT_reflections_percentage_without_transits = (
                        '%.1f' % ((reflections_without_transits * 100) /
                                  (reflections_without_transits + transmissions_without_transits)))
                SPRINT_transmissions_percentage_without_transits = (
                        '%.1f' % ((transmissions_without_transits * 100) /
                                  (reflections_without_transits + transmissions_without_transits)))
            else:
                SPRINT_reflections_percentage_without_transits = '%.1f' % 0
                SPRINT_transmissions_percentage_without_transits = '%.1f' % 0

            SPRINT_reflections_with_transits = '%d' % reflections
            SPRINT_reflections = f'${SPRINT_reflections_with_transits}_{{({SPRINT_reflections_percentage_without_transits}\%)}}$'
            SPRINT_reflections_text = '$R_{SPRINT}$'
            SPRINT_transmissions_with_transits = '%d' % transmissions
            SPRINT_transmissions = f'${SPRINT_transmissions_with_transits}_{{({SPRINT_transmissions_percentage_without_transits}\%)}}$'
            SPRINT_transmissions_text = '$T_{SPRINT}$'
            SPRINT_Score = f'{SPRINT_reflections} - {SPRINT_transmissions}'
            SPRINT_text = f'{SPRINT_reflections_text} {SPRINT_reflections} - {SPRINT_transmissions} {SPRINT_transmissions_text}'
            props_SPRINT = dict(boxstyle='round', edgecolor='gray', linewidth=2, facecolor='gray', alpha=0.5)

            # Coherence results box
            bright = self.Exp_dict['Analysis_results'][str(cond)]['Bright']
            dark = self.Exp_dict['Analysis_results'][str(cond)]['Dark']
            bright_without_transits = self.Exp_dict['Analysis_results'][str(cond)]['Bright_without_transits']
            dark_without_transits = self.Exp_dict['Analysis_results'][str(cond)]['Dark_without_transits']
            if (bright_without_transits + dark_without_transits) > 0:
                SPRINT_BP_percentage_without_transits = (
                        '%.1f' % ((bright_without_transits * 100) / (bright_without_transits + dark_without_transits)))
                SPRINT_DP_percentage_without_transits = (
                        '%.1f' % ((dark_without_transits * 100) / (bright_without_transits + dark_without_transits)))
            else:
                SPRINT_BP_percentage_without_transits = '%.1f' % 0
                SPRINT_DP_percentage_without_transits = '%.1f' % 0

            SPRINT_BP_counts_with_transits = '%d' % bright
            SPRINT_BP_counts = f'${SPRINT_BP_counts_with_transits}_{{({SPRINT_BP_percentage_without_transits}\%)}}$'
            SPRINT_BP_counts_text = '$BP_{SPRINT}$'
            SPRINT_DP_counts_with_transits = '%d' % dark
            SPRINT_DP_counts = f'${SPRINT_DP_counts_with_transits}_{{({SPRINT_DP_percentage_without_transits}\%)}}$'
            SPRINT_DP_counts_text = '$DP_{SPRINT}$'
            SPRINT_Coherence_Score = f'{SPRINT_BP_counts} - {SPRINT_DP_counts}'
            SPRINT_Coherence_text = f'{SPRINT_BP_counts_text} {SPRINT_BP_counts} - {SPRINT_DP_counts} {SPRINT_DP_counts_text}'
            props_SPRINT_coherence = dict(boxstyle='round', edgecolor='gray', linewidth=2, facecolor='gray', alpha=0.5)

            # create figure for SPRINT information:
            self.SPRINT_figure = plt.figure(str(cond) + ' - SPRINT information')
            ax = plt.gca()
            seq_with_SPRINT_data = np.zeros(self.number_of_sequences)
            seq_with_SPRINT_data[[vec for elem in
                                  self.Exp_dict['Analysis_results'][str(cond)]['sequence_indices_with_data_points_per_cycle']
                                  for vec in elem]] += 1
            plt.plot(range(self.number_of_sequences), seq_with_SPRINT_data,
                       label='All transit events with data')
            plt.xlabel('Sequence [#]')
            plt.ylabel('Counts [Photons]')
            # plt.legend(loc='upper right')
            plt.text(0.5, 1.2, SPRINT_text, transform=ax.transAxes, fontsize=24, verticalalignment='top',
                     horizontalalignment='center', bbox=props_SPRINT)
            plt.text(0.5, 0.9, SPRINT_Coherence_text, transform=ax.transAxes, fontsize=24, verticalalignment='top',
                     horizontalalignment='center', bbox=props_SPRINT_coherence)

        plt.show(block=True)

    def save_dict_as_folders_and_variables(self, data_dict, base_path):
        """
        Save a dictionary to folders with .npz files.

        Parameters:
            data_dict (dict): Dictionary to be saved.
            base_path (str): Base path where folders and .npz files will be created.
        """
        # Ensure the base path exists
        if not os.path.exists(base_path):
            os.makedirs(base_path)
        self.save_to_folder(base_path, data_dict)

    def save_to_folder(self, folder_path, data):
        """
        Recursively save dictionary data into folders and .npz files.

        Parameters:
            folder_path (str): Path of the current folder.
            data (dict): Dictionary data to be saved.
        """
        for key, value in data.items():
            current_path = os.path.join(folder_path, key)

            if isinstance(value, dict):
                # If value is a dictionary, create a new folder
                if not os.path.exists(current_path):
                    os.makedirs(current_path)
                # Recursively process this dictionary
                self.save_to_folder(current_path, value)
            else:
                # If value is not a dictionary, save it as a .npz file
                if isinstance(value, np.ndarray):
                    np.savez(os.path.join(current_path + '.npz'), data=value)
                # elif isinstance(value, list):
                else:
                    with open(os.path.join(current_path + '.json'), 'w') as file:
                        try:
                            json.dump(value, file)
                        except Exception as error:
                            print("An error occurred:", error)

    def load_all_analysis_data(self):
        '''

        :return:
        '''
        # Open folder and load to dictionary
        root = Tk()
        self.analysis_data_path = '{}'.format(askdirectory(title='Experiment folder', initialdir=r'U:\Lab_2023\Experiment_results'))
        # Create a list of paths for each sub-experiment analysis results
        list_of_analysis_data_dirs = []
        list_of_exp_config_values_path = []
        for path, dirs, files in os.walk(self.analysis_data_path):
            if 'Analysis_results' in path.split("\\")[-1]:
                list_of_analysis_data_dirs.append(path)
            if 'meta_data' in path.split("\\")[-1]:
                list_of_exp_config_values_path.append(path)
        data = Experiment_data_load.DictionaryBuilder()
        # Create a list of dictionaries for each sub-experiment analysis results
        dictionary = []
        for indx, path in enumerate(list_of_analysis_data_dirs):
            dictionary.append(data.load_files_to_dict(path))
        for indx, path in enumerate(list_of_exp_config_values_path):
            if 'without' in path.lower():
                exp_type = 'background'
            elif 'with' in path.lower():
                exp_type = 'experiment'
            temp_dict = data.load_files_to_dict(path)
            config_values_key = list(key for key in temp_dict.keys() if 'config' in key)[0]
            dictionary[indx//2][exp_type]['Exp_config_values'] = data.load_files_to_dict(path)[config_values_key]
        messagebox.showinfo(title='Success!', message='Analysis results are ready!')
        root.destroy()

        return dictionary

    def analyze_all_results(self, north_efficiency=0.357, south_efficiency=0.465375, min_transit_len=1000,
                            max_transit_len=10000):
        '''

        :return:
        '''
        self.number_of_photons_per_pulse = []
        self.fidelity, self.qber, self.qber_err, self.SNR, self.qber_bg, self.qber_err_bg = [
            [
                [] for _ in range(len(self.list_of_all_analysis_dictionaries))
            ] for _ in range(6)
        ]
        self.position_of_data = [
            [
                [] for _ in range(len(self.list_of_all_analysis_dictionaries))
            ] for _ in range(len(self.conditions))
        ]
        for indx, dicionary in enumerate(self.list_of_all_analysis_dictionaries):
            (self.number_of_photons_per_pulse.
             append(sum(self.calc_average_photons_per_SPRINT_pulse(dicionary, north_efficiency=0.357,
                                                                   south_efficiency=0.465375))))
            for cond_index, cond in enumerate(self.conditions):
                # Calculate fidelity for each condition
                if self.experiment_type == 'T':
                    self.fidelity[indx].append(dicionary['experiment'][cond]['transmission'] /
                                               (dicionary['experiment'][cond]['transmission'] +
                                                dicionary['experiment'][cond]['reflection']))
                elif self.experiment_type == 'R':
                    self.fidelity[indx].append(dicionary['experiment'][cond]['reflection'] /
                                               (dicionary['experiment'][cond]['transmission'] +
                                                dicionary['experiment'][cond]['reflection']))

                # QBER for each condition:

                ## get all indices where there are no more then 1 count in "bright" + "dark" port
                b = []
                d = []
                for cyc_index, lst_of_transits in enumerate(dicionary['experiment'][cond]['all_transits_length_per_cond']):
                    for transit_indx, transit_len in enumerate(lst_of_transits):
                        if (transit_len <= min_transit_len) or (transit_len >= max_transit_len):
                            continue
                        for i, seq_indx in enumerate(dicionary['experiment'][cond][
                                    'sequence_indices_with_data_points_per_cycle'][cyc_index]):
                            if seq_indx in dicionary['experiment'][cond][
                                'indices_of_all_sequences_with_transits_per_cycle'][cyc_index][transit_indx]:
                                b.append(dicionary['experiment'][cond][
                                             'BP_counts_per_data_point_per_cycle'][cyc_index][i])
                                d.append(dicionary['experiment'][cond][
                                             'DP_counts_per_data_point_per_cycle'][cyc_index][i])
                                self.position_of_data[cond_index][indx].append([cyc_index, seq_indx])

                # bright_counts_per_data_point = np.array(sum(dicionary['experiment'][cond]
                #                                             ['BP_counts_per_data_point_per_cycle'], []))
                # dark_counts_per_data_point = np.array(sum(dicionary['experiment'][cond]
                #                                           ['DP_counts_per_data_point_per_cycle'], []))
                bright_counts_per_data_point = np.array(b)
                dark_counts_per_data_point = np.array(d)

                # rel_indices = np.where((bright_counts_per_data_point + dark_counts_per_data_point) < 2)
                rel_indices = np.where((bright_counts_per_data_point + dark_counts_per_data_point) > 0)
                bright_counts_per_data_point[rel_indices] = (bright_counts_per_data_point[rel_indices] /
                                                             (bright_counts_per_data_point[rel_indices] +
                                                              dark_counts_per_data_point[rel_indices]))
                dark_counts_per_data_point[rel_indices] = (dark_counts_per_data_point[rel_indices] /
                                                           (bright_counts_per_data_point[rel_indices] +
                                                            dark_counts_per_data_point[rel_indices]))

                ## calculate qber
                Bright_counts = sum(bright_counts_per_data_point[rel_indices])
                Dark_counts = sum(dark_counts_per_data_point[rel_indices])
                if (Bright_counts + Dark_counts) == 0:
                    self.qber[indx].append(0)
                    self.qber_err[indx].append(0)
                else:
                    self.qber[indx].append(Dark_counts / (Dark_counts + Bright_counts))
                    qber_std = np.sqrt((Dark_counts / (Dark_counts + Bright_counts) ** 2 *
                                        np.sqrt(Bright_counts)) ** 2 + (Bright_counts / (Dark_counts + Bright_counts) ** 2 *
                                                                        np.sqrt(Dark_counts)) ** 2)
                    self.qber_err[indx].append(qber_std)
                self.SNR[indx].append(dicionary['experiment'][cond]['SNR'])

                # QBER without atoms:

                ## get all indices where there are no more then 1 count in "bright" + "dark" port
                bright_counts_per_data_point_bg = np.array(sum(dicionary['background'][cond]
                                                               ['BP_counts_without_transits_per_data_point_per_cycle'],
                                                               []))
                dark_counts_per_data_point_bg = np.array(sum(dicionary['background'][cond]
                                                             ['DP_counts_without_transits_per_data_point_per_cycle'],
                                                             []))

                # rel_indices = np.where((bright_counts_per_data_point + dark_counts_per_data_point) < 2)
                rel_indices_bg = np.where((bright_counts_per_data_point_bg + dark_counts_per_data_point_bg) > 0)
                bright_counts_per_data_point_bg[rel_indices_bg] = (bright_counts_per_data_point_bg[rel_indices_bg] /
                                                                   (bright_counts_per_data_point_bg[rel_indices_bg] +
                                                                    dark_counts_per_data_point_bg[rel_indices_bg]))
                dark_counts_per_data_point_bg[rel_indices_bg] = (dark_counts_per_data_point_bg[rel_indices_bg] /
                                                                 (bright_counts_per_data_point_bg[rel_indices_bg] +
                                                                  dark_counts_per_data_point_bg[rel_indices_bg]))

                ## calculate qber
                Bright_counts_bg = sum(bright_counts_per_data_point_bg[rel_indices_bg])
                Dark_counts_bg = sum(dark_counts_per_data_point_bg[rel_indices_bg])
                self.qber_bg[indx].append(Dark_counts_bg / (Dark_counts_bg + Bright_counts_bg))
                qber_std = np.sqrt((Dark_counts_bg / (Dark_counts_bg + Bright_counts_bg) ** 2 *
                                    np.sqrt(Bright_counts_bg)) ** 2 +
                                   (Bright_counts_bg / (Dark_counts_bg + Bright_counts_bg) ** 2 *
                                    np.sqrt(Dark_counts_bg)) ** 2)
                self.qber_err_bg[indx].append(qber_std)

                # ## Calculate qber
                # self.qber[indx].append(dicionary['experiment'][cond]['Dark'] /
                #                        (dicionary['experiment'][cond]['Dark'] +
                #                         dicionary['experiment'][cond]['Bright']))
                # qber_std = np.sqrt((dicionary['experiment'][cond]['Dark'] / (dicionary['experiment'][cond]['Dark'] +
                #                                                              dicionary['experiment'][cond][
                #                                                                  'Bright']) ** 2 *
                #                     np.sqrt(dicionary['experiment'][cond]['Bright']))**2 +
                #                    (dicionary['experiment'][cond]['Bright'] / (dicionary['experiment'][cond]['Dark'] +
                #                                                                dicionary['experiment'][cond][
                #                                                                  'Bright']) ** 2 *
                #                     np.sqrt(dicionary['experiment'][cond]['Dark'])) ** 2)
                # self.qber_err[indx].append(qber_std)
                # self.SNR[indx].append(dicionary['experiment'][cond]['SNR'])

    def calc_average_photons_per_SPRINT_pulse(self, dictionary, coupling_tranmission=0.49, north_efficiency=0.357,
                                              south_efficiency=0.465375):
        '''

        :param dictionary:
        :param north_efficiency:
        :param south_efficiency:
        :param coupling_tranmission:
        :return:
        '''
        clicks_in_north_per_seq = (dictionary['background']['folded_tt_N'] + dictionary['background']['folded_tt_BP']
                                   + dictionary['background']['folded_tt_DP'])
        clicks_in_south_per_seq = dictionary['background']['folded_tt_S'] + dictionary['background']['folded_tt_FS']
        # plt.plot(clicks_in_north_per_seq)
        # plt.plot(clicks_in_south_per_seq)
        total_clicks_in_SPRINT_pulse_north = []
        total_clicks_in_SPRINT_pulse_south = []
        total_clicks_in_SPRINT_pulse_transmission_considering_efficiency = []
        total_clicks_in_SPRINT_pulse_reflection_considering_efficiency = []
        pulse_duration = []
        for lst in dictionary['background']['Pulses_location_in_seq']:
            if lst[-1] == 'n':
                pulse_duration.append(lst[1]-lst[0])
                total_clicks_in_SPRINT_pulse_north.append(sum(clicks_in_north_per_seq[lst[0]:lst[1]]))
                (total_clicks_in_SPRINT_pulse_transmission_considering_efficiency.
                 append(sum(clicks_in_north_per_seq[lst[0]:lst[1]])/(coupling_tranmission*north_efficiency)))
                total_clicks_in_SPRINT_pulse_south.append(sum(clicks_in_south_per_seq[lst[0]:lst[1]]))
                (total_clicks_in_SPRINT_pulse_reflection_considering_efficiency.
                 append(sum(clicks_in_south_per_seq[lst[0]:lst[1]])/south_efficiency))
            elif lst[-1] == 's':
                pulse_duration.append(lst[1]-lst[0])
                total_clicks_in_SPRINT_pulse_north.append(sum(clicks_in_north_per_seq[lst[0]:lst[1]]))
                (total_clicks_in_SPRINT_pulse_reflection_considering_efficiency.
                 append(sum(clicks_in_north_per_seq[lst[0]:lst[1]])/north_efficiency))
                total_clicks_in_SPRINT_pulse_south.append(sum(clicks_in_south_per_seq[lst[0]:lst[1]]))
                (total_clicks_in_SPRINT_pulse_transmission_considering_efficiency.
                 append(sum(clicks_in_south_per_seq[lst[0]:lst[1]])/(coupling_tranmission*south_efficiency)))

        number_of_seq_in_exp = (math.ceil((dictionary['background']['Exp_config_values']['M_window']-1e6) /
                                          (len(dictionary['background']['filter_N']))))
        avg_photons_per_pulse = (np.array(total_clicks_in_SPRINT_pulse_transmission_considering_efficiency) /
                                 number_of_seq_in_exp)

        return avg_photons_per_pulse

    def reflection_or_transmission_exp(self, dictionary):
        '''
        Go through the pulse sequence and check if the direction of the last detection pulse (the preparation pulse) is
        the same as the direction of the SPRINT pulse
        :param dictionary:
        :return:
        '''

        self.list_of_pulses = [s[-1] for s in dictionary['background']['Pulses_location_in_seq']]
        for i in range(1, len(self.list_of_pulses)):
            if (self.list_of_pulses[i-1] == 'N' and self.list_of_pulses[i] == 'n') or (self.list_of_pulses[i-1] == 'S' and self.list_of_pulses[i] == 's'):
                return 'T'
            if (self.list_of_pulses[i-1] == 'S' and self.list_of_pulses[i] == 'n') or (self.list_of_pulses[i-1] == 'N' and self.list_of_pulses[i] == 's'):
                return 'R'

    def plot_all_QBER_results_bg(self):
        '''

        :return:
        '''
        self.plot_lines_bg = []
        self.f_bg, self.ax_bg = plt.subplots(1, 1, num='QBER without atoms')
        my_handler_map = {ErrorbarContainer: CustomErrorbarHandler(numpoints=1)}
        for indx, cond in enumerate(self.conditions):
            qber_per_cond = [qber_bg[indx] for qber_bg in self.qber_bg]
            qber_err_per_cond = [qber_err_bg[indx] for qber_err_bg in self.qber_err_bg]
            SNR_text = '%.1f' % np.mean([snr[indx] for snr in self.SNR])
            # plot_line, = plt.plot(self.number_of_photons_per_pulse, qber_per_cond, label=('condition: ' + cond))
            plot_line = self.ax_bg.errorbar(self.number_of_photons_per_pulse, qber_per_cond, yerr=qber_err_per_cond, fmt='-o',
                                      capsize=3, label=('condition: ' + cond + ', SNR = ' + SNR_text))
            self.plot_lines_bg.append(plot_line)
            # plt.text(np.array(self.number_of_photons_per_pulse)*1.03, np.array(qber_per_cond)*1.03, SNR_text, fontsize=10)
        self.lgnd = self.ax_bg.legend(loc='upper right', handler_map=my_handler_map)
        self.ax_bg.set_ylabel('QBER', fontsize=24)
        self.ax_bg.set_xlabel('$\mu $[#photons from Alice]', fontsize=24)
        self.ax_bg.grid(visible=True, which='both', axis='both')
        # self.ax.grid(visible=True, which='major', axis='both')
        # self.ax.minorticks_on()
        self.ax_bg.xaxis.set_minor_locator(AutoMinorLocator(2))

        # self.f.canvas.mpl_connect('pick_event', self.onpick)
        # plt.show()

    def plot_all_QBER_results(self):
        '''

        :return:
        '''

        # General plot parameters
        mpl.rcParams['font.family'] = 'Cambria'
        mpl.rcParams['font.size'] = 18
        mpl.rcParams['axes.linewidth'] = 2
        # mpl.rcParams['axes.spines.top'] = False
        # mpl.rcParams['axes.spines.right'] = False
        mpl.rcParams['xtick.major.size'] = 10
        mpl.rcParams['xtick.major.width'] = 2
        mpl.rcParams['ytick.major.size'] = 10
        mpl.rcParams['ytick.major.width'] = 2

        self.plot_lines = []
        # self.f, self.ax = plt.subplots(1, 1, num='QBER')
        self.f = plt.figure(num='QBER', figsize=(10, 8.5))
        self.ax = self.f.add_subplot(111)
        self.f.subplots_adjust(bottom=0.10, top=0.90)
        my_handler_map = {ErrorbarContainer: CustomErrorbarHandler(numpoints=1)}

        # Create axes for sliders
        ax_transit_length = self.f.add_axes([0.3, 0.95, 0.4, 0.05])
        ax_transit_length.spines['top'].set_visible(True)
        ax_transit_length.spines['right'].set_visible(True)

        # Create slider
        self.s_transit_length = RangeSlider(ax=ax_transit_length, label='Transit length', valmin=500, valmax=10000,
                                            valinit=(1000, 10000), valfmt=' %d [ns]', facecolor='#cc7000')

        for indx, cond in enumerate(self.conditions):
            qber_per_cond = [qber[indx] for qber in self.qber]
            qber_err_per_cond = [qber_err[indx] for qber_err in self.qber_err]
            SNR_text = '%.1f' % np.mean([snr[indx] for snr in self.SNR])
            # plot_line, = plt.plot(self.number_of_photons_per_pulse, qber_per_cond, label=('condition: ' + cond))
            plot_line = self.ax.errorbar(self.number_of_photons_per_pulse, qber_per_cond, yerr=qber_err_per_cond, fmt='-o',
                                         capsize=3, label=('condition: ' + cond + ', SNR = ' + SNR_text))
            self.plot_lines.append(plot_line)
            # plt.text(np.array(self.number_of_photons_per_pulse)*1.03, np.array(qber_per_cond)*1.03, SNR_text, fontsize=10)
        self.lgnd = self.ax.legend(loc='upper right', handler_map=my_handler_map)
        self.ax.set_ylabel('QBER', fontsize=24)
        self.ax.set_xlabel('$\mu $[#photons from Alice]', fontsize=24)
        self.ax.grid(visible=True, which='both', axis='both')
        # self.ax.grid(visible=True, which='major', axis='both')
        # self.ax.minorticks_on()
        self.ax.xaxis.set_minor_locator(AutoMinorLocator(2))


        self.lined = dict()
        for legend_line, plot_line in zip(self.lgnd.legend_handles, self.plot_lines):
            legend_line.set_picker(5)
            self.lined[legend_line] = plot_line

        # self.f.canvas.mpl_connect('pick_event', self.onpick)
        # plt.show()

    def onpick(self, event):
        legend_line = event.artist
        plot = self.lined[legend_line]
        for p in plot.get_children():
            vis = not p.get_visible()
            p.set_visible(vis)
        legend_line.set_alpha(1.0 if vis else 0.2)
        self.f.canvas.draw()

    def update(self, val):
        (min_length, max_length) = self.s_transit_length.val
        self.analyze_all_results(north_efficiency=0.357, south_efficiency=0.465375, min_transit_len=min_length,
                                 max_transit_len=max_length)

        for indx, cond in enumerate(self.conditions):
            qber_per_cond = np.array([qber[indx] for qber in self.qber])
            qber_err_per_cond = np.array([qber_err[indx] for qber_err in self.qber_err])
            # self.plot_lines[indx].set_data([self.number_of_photons_per_pulse], [qber_per_cond])
            ln, (err_top, err_bot), (bars, ) = self.plot_lines[indx]

            ln.set_data([self.number_of_photons_per_pulse], [qber_per_cond])
            x_base = self.number_of_photons_per_pulse
            y_base = qber_per_cond

            yerr_top = y_base + qber_err_per_cond
            yerr_bot = y_base - qber_err_per_cond

            err_top.set_ydata(yerr_top)
            err_bot.set_ydata(yerr_bot)

            new_segments = [np.array([[x, yt], [x, yb]]) for
                            x, yt, yb in zip(x_base, yerr_top, yerr_bot)]

            bars.set_segments(new_segments)

        self.f.canvas.draw_idle()

        # Class's constructor
    def __init__(self, analyze_results=False, transit_conditions=None):

        # self.exp_type, self.exp_date, self.exp_time = self.popupbox_inquiry()
        # while True:
        #     if self.exp_type is None or self.exp_date is None or self.exp_time is None:
        #         option = pymsgbox.confirm('Missing an input, do you wish to retry?', 'Experiment data loading',
        #                                   ['Yes please', 'No thanks'])
        #         if option == 'Yes please':
        #             self.exp_type, self.exp_date, self.exp_time = self.popupbox_inquiry()
        #         else:
        #             break
        #     else:
        #         # print(exp_type, exp_date, exp_time)
        #         self.data = Experiment_data_load.DictionaryBuilder(self.exp_type, self.exp_date, self.exp_time)
        #         if self.data is None:
        #             self.exp_type, self.exp_date, self.exp_time = self.popupbox_inquiry()
        #         else:
        #             self.Exp_dict = self.data.load_files_to_dict(self.data.exp_path)
        #             pymsgbox.alert(text='Experiment data is ready to use.', title='Success!')
        #             break
        #
        # self.init_params_for_experiment()
        # self.fold_tt_histogram(self.sequence_len)
        # self.background_noise_N = ((sum(self.folded_tt_N[int(self.Pulses_location_in_seq[-1][1]):300]) +
        #                             sum(self.folded_tt_BP[int(self.Pulses_location_in_seq[-1][1]):300]) +
        #                             sum(self.folded_tt_DP[int(self.Pulses_location_in_seq[-1][1]):300])) /
        #                            (300 - int(self.Pulses_location_in_seq[-1][
        #                                           1])))  # Per [ns] * number_of_sequences_per_cycle * number_of_cycles
        # # ((self.sequence_len - int(self.Pulses_location_in_seq[-1][1])) * len(self.Exp_dict['output']['South(5)']['South_timetags'])))  # Per [ns] * number_of_sequences_per_cycle
        # self.background_noise_S = ((sum(self.folded_tt_S[int(self.Pulses_location_in_seq[-1][1]):300]) +
        #                             sum(self.folded_tt_FS[int(self.Pulses_location_in_seq[-1][1]):300])) /
        #                            (300 - int(self.Pulses_location_in_seq[-1][
        #                                           1])))  # Per [ns] * number_of_sequences_per_cycle * number_of_cycles
        # if self.Pulses_location_in_seq[0][2] == 'N':
        #     print((sum(self.folded_tt_S_directional[
        #                int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])]) -
        #            self.background_noise_S * (
        #                        int(self.Pulses_location_in_seq[0][1]) - int(self.Pulses_location_in_seq[0][0]))) /
        #           sum(self.folded_tt_N_directional[
        #               int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])]))
        # else:
        #     print((sum(self.folded_tt_N_directional[
        #                int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])]) -
        #            self.background_noise_N * (
        #                        int(self.Pulses_location_in_seq[0][1]) - int(self.Pulses_location_in_seq[0][0]))) /
        #           sum(self.folded_tt_S_directional[
        #               int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])]))

        if transit_conditions is None:
            transit_conditions = [[[2, 1, 2]]]
        self.paths_map = {
            "name": __name__,
            "file": __file__,
            "cwd": os.getcwd(),
            "root": str(pathlib.Path(__file__).parent.resolve())
        }

        # Setup console logger. We do this first, so rest of code can use logging functions.
        self.logger = BDLogger()

        # Initialize the BDBatch helper - to serve us when batching experiment samples
        self.batcher = BDBatch(json_map_path=self.paths_map['cwd'])

        # Initialize the BDResults helper - for saving experiment results
        self.bd_results = BDResults(json_map_path=self.paths_map['cwd'], version="0.1", logger=self.logger)
        self.bd_results.create_experiment_run_folder()

        # Check network driver availability
        network_drive_letter = 'U'
        network_drive_available = self.bd_results.is_network_drive_available(f'{network_drive_letter}:\\Lab_2023')
        if not network_drive_available:
            self.logger.error(f'Network drive {network_drive_letter} is not available/connected. PLEASE FIX.')
            sys.exit(1)

        self.transit_conditions = transit_conditions

        if analyze_results:
            self.list_of_all_analysis_dictionaries = self.load_all_analysis_data()
            self.conditions = [key for key in self.list_of_all_analysis_dictionaries[0]['background'].keys() if
                               '[[' in key]
            self.experiment_type = self.reflection_or_transmission_exp(
                self.list_of_all_analysis_dictionaries[0])  # returns 'T' for transmission and 'R' for reflection
            self.analyze_all_results(north_efficiency=0.357, south_efficiency=0.465375, min_transit_len=1000,
                                     max_transit_len=6000)
            self.plot_all_QBER_results_bg()
            self.plot_all_QBER_results()
        else:
            # Open folder and load to dictionary
            self.Exp_dict = self.open_folder_to_dictionary()

            # Divide data to background and experiment data:
            exp_key = list(self.Exp_dict.keys())[0]
            for key in list(self.Exp_dict.keys()):
                if 'without' in key.lower():
                    self.Background_dict = self.Exp_dict[key]
                elif 'with' in key.lower():
                    exp_key = key
            self.Exp_dict = self.Exp_dict[exp_key]

            # background analysis:
            # check number of cycles in experiment:
            self.number_of_cycles = len(list(self.Background_dict['output'][list(self.Background_dict['output'].keys())[0]].values())[0])
            self.init_params_for_experiment(self.Background_dict)
            # Initialize the batcher
            self.batcher.set_batch_size(self.number_of_cycles)
            self.batcher.empty_all()
            for cycle in tqdm(range(self.number_of_cycles)):
                self.ingest_time_tags(self.Background_dict, cycle, delay=0)
                self.experiment_calculations()

                # fold time-tags and accumulate to a running average:
                self.fold_tt_histogram(self.sequence_len)
                self.calculate_running_averages(cycle + 1)

                for condition_number, transit_condition in enumerate(self.transit_conditions):
                    ### Find transits and extract SPRINT data:  ###
                    self.get_transit_data(transit_condition[0], condition_number, transit_condition[1])
                self.batcher.batch_all(self)
            self.Background_dict['Analysis_results'] = self.results_to_dictionary()
            self.save_dict_as_folders_and_variables({'Analysis_results': {
                                                              'background': self.Background_dict['Analysis_results']
                                                              }}, self.exp_data_path)

            # "real" experiment analysis:
            # check number of cycles in experiment:
            self.number_of_cycles = len(list(self.Exp_dict['output'][list(self.Exp_dict['output'].keys())[0]].values())[0])
            self.init_params_for_experiment(self.Exp_dict)

            # Initialize the batcher
            self.batcher.set_batch_size(self.number_of_cycles)
            self.batcher.empty_all()

            for cycle in tqdm(range(self.number_of_cycles)):
                self.ingest_time_tags(self.Exp_dict, cycle, delay=0)
                self.experiment_calculations()

                # fold time-tags and accumulate to a running average:
                self.fold_tt_histogram(self.sequence_len)
                self.calculate_running_averages(cycle + 1)

                for condition_number, transit_condition in enumerate(self.transit_conditions):
                    ### Find transits and extract SPRINT data:  ###
                    self.get_transit_data(transit_condition[0], condition_number, transit_condition[1])
                self.batcher.batch_all(self)

            self.Exp_dict['Analysis_results'] = self.results_to_dictionary(True)
            self.save_dict_as_folders_and_variables({'Analysis_results': {
                                                                    'experiment': self.Exp_dict['Analysis_results']
                                                                }}, self.exp_data_path)
            self.plot_results()
            self.batcher.empty_all()


class CustomErrorbarHandler(matplotlib.legend_handler.HandlerErrorbar):
    """
    Sub-class the standard error-bar handler
    """

    def create_artists(self, *args, **kwargs):
        #  call the parent class function
        a_list = matplotlib.legend_handler.HandlerErrorbar.create_artists(self, *args, **kwargs)
        # re-order the artist list, only the first artist is added to the
        # legend artist list, this is the one that corresponds to the markers
        a_list = a_list[-1:] + a_list[:-1]
        return a_list

if __name__ == '__main__':
    # tansit_conditions is a list of lists. Each list comprises a list of conditions for transits and the number
    # of photons reflected in the "preparation" pulse to take for SPRINT data analysis.
    # self = experiment_data_analysis(transit_conditions=[[[[1, 1, 2]], 0], [[[1, 1, 2]], 1], [[[1, 2, 1]], 0],
    #                                                     [[[1, 2, 1]], 1], [[[1, 2]], 0], [[[1, 2]], 1],
    #                                                     [[[2, 1], [1, 2]], 0], [[[2, 1], [1, 2]], 1],
    #                                                     [[[2, 1, 2]], 0], [[[2, 1, 2]], 1]])
    # self = experiment_data_analysis(transit_conditions=[[[[2, 1, 2]], 1]])
    # self = experiment_data_analysis(transit_conditions=[[[[2, 1]], 0],  [[[2, 1]], 1]])

    self = experiment_data_analysis(analyze_results=True)
    self.s_transit_length.on_changed(self.update)
    self.f.canvas.mpl_connect('pick_event', self.onpick)
    plt.show()