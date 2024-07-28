import Experiment_data_load
import numpy as np
import os
import sys
import json
import math
import pymsgbox
import pathlib
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import messagebox
from tkinter.filedialog import askdirectory
from Utilities.Utils import Utils
from Utilities.BDLogger import BDLogger
from Utilities.BDBatch import BDBatch
from Utilities.BDResults import BDResults
from tqdm import tqdm


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

        self.number_of_PNSA_sequences = math.ceil(self.M_window / self.sequence_len)

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

        self.tt_S_binning = np.zeros(self.number_of_PNSA_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_PNSA_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_PNSA_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sequence_len)
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_accumulated = np.zeros(self.number_of_PNSA_sequences)

        num_of_seq_per_count = 50
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences // num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences // num_of_seq_per_count)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences // num_of_seq_per_count)

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_PNSA_sequences)

        self.all_transits_seq_indx_per_cond = [
            [] for _ in range(len(self.transit_conditions))
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

    def ingest_time_tags(self, dict, cycle):
        """
        Takes all the raw results we got from the streams and does some processing on them - preparing "measures"
        """

        # ------------------------------------------
        # Unify detectors and windows within detectors and create vector of tt's for each direction (bright port, dark port, north, south and from FS) and sort them
        # ------------------------------------------
        self.tt_BP_measure = dict['output']['Bright(1,2)']['Bright_timetags'][cycle]
        self.tt_DP_measure = dict['output']['Dark(3,4)']['Dark_timetags'][cycle]
        self.tt_N_measure = dict['output']['North(8)']['North_timetags'][cycle]
        self.tt_S_measure = dict['output']['South(5)']['South_timetags'][cycle]
        self.tt_FS_measure = dict['output']['FastSwitch(6,7)']['FS_timetags'][cycle]

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

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_PNSA_sequences)

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
                for _ in range(self.number_of_PNSA_sequences)
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
                for _ in range(self.number_of_PNSA_sequences)
            ] for _ in range(8)
        ]

        # N direction:
        self.N_tt = np.array(self.tt_N_measure)
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
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # Bright port direction:
        self.BP_tt = np.array(self.tt_BP_measure)
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
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(element)
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(element)
                            self.num_of_BP_counts_per_seq_in_SPRINT_pulse[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        # Dark port direction:
        self.DP_tt = np.array(self.tt_DP_measure)
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
                                self.num_of_det_transmissions_per_seq_N_[seq_num][indx].append(element)
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(element)
                            self.num_of_DP_counts_per_seq_in_SPRINT_pulse[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    self.num_of_det_reflections_per_seq_S[seq_num] += self.filter_S[tt_inseq]
                    self.num_of_det_transmissions_per_seq_N[seq_num] += self.filter_N[tt_inseq]

        self.S_tt = np.array(self.tt_S_measure + self.tt_FS_measure)
        tt_inseq_ = self.S_tt % self.sequence_len
        seq_num_ = self.S_tt // self.sequence_len
        for (element, tt_inseq, seq_num) in zip(self.S_tt, tt_inseq_, seq_num_):
            if (element > int(0.6e6)) and (element < int(self.M_window - 0.4e6)):
                for indx, tup in enumerate(self.Pulses_location_in_seq):
                    if (tt_inseq >= tup[0]) and (tt_inseq <= tup[1]):
                        if indx < self.number_of_detection_pulses_per_seq:
                            if tup[2] == 'N':
                                self.num_of_det_reflections_per_seq_N_[seq_num][indx].append(element)
                                self.num_of_det_reflections_per_seq_[seq_num][indx].append(element)
                            if tup[2] == 'S':
                                self.num_of_det_transmissions_per_seq_S_[seq_num][indx].append(element)
                                self.num_of_det_transmissions_per_seq_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_reflections_per_seq_N_[seq_num][ind].append(element)
                                self.num_of_SPRINT_reflections_per_seq_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_transmissions_per_seq_S_[seq_num][ind].append(element)
                                self.num_of_SPRINT_transmissions_per_seq_[seq_num][ind].append(element)

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

        ## fold reflections and transmission
        # self.fold_tt_histogram(exp_sequence_len=self.sequence_len)

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
            for i in range(self.number_of_PNSA_sequences - len(cond) + 1):
                cond_check = (self.num_of_det_reflections_per_seq[i:(i + len(cond))] >= cond).astype(int)
                if sum(cond_check) >= minimum_number_of_seq_detected:
                    # TODO: ask dor (08.01.24) - what happens at [0,4,0]? and why including the middle at [2,0,2]?
                    # adding to current transit the indices from first element satisfing the condition to the last element checked.
                    # for example:
                    # if the condition is [1,1,1] and for i=7000 the reflections were [0(i=7000),1 (i=7001),1 (i=7002)]
                    # than current transit would add [7001,7002]
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

    def analyze_SPRINT_data_points(self, all_transits_seq_indx, SPRINT_pulse_number=[1], background=False):
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
                       len(self.num_of_det_reflections_per_seq_N_[seq_indx][-1]) >= 1:
                        potential_data = True
                        # self.potential_data += 1
                    elif self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq-1][2] == 'S' and \
                            len(self.num_of_det_reflections_per_seq_S_[seq_indx][-1]) >= 1:
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

    def get_transit_data(self, transit_condition, cond_num):
        '''
        Use the data loaded from the experiment to find atomic transits under different conditions and some data to
        compare the conditions, such as, number of transits, transit length distribution and more...
        :param transit_condition: List of lists, so that each condition is of the structure [[1,2],[2,1],....]
               and so on...
        :return:
        '''

        self.all_transits_seq_indx_per_cond[cond_num] = self.find_transit_events(cond_list=transit_condition, minimum_number_of_seq_detected=2)

        # Analyze SPRINT data during transits:
        (self.seq_with_data_points[cond_num], self.reflection_SPRINT_data[cond_num], self.transmission_SPRINT_data[cond_num],
         self.BP_counts_SPRINT_data[cond_num], self.DP_counts_SPRINT_data[cond_num]) = \
            self.analyze_SPRINT_data_points(self.all_transits_seq_indx_per_cond[cond_num], SPRINT_pulse_number=[1, 2],
                                            background=False)  # Enter the index of the SPRINT pulse for which the data should be analyzed
        # print(self.potential_data)
        # Analyze SPRINT data when no transit occur:
        self.all_seq_without_transits = [
            np.delete(np.arange(0, self.number_of_PNSA_sequences, 1, dtype='int'),
                      sum(self.all_transits_seq_indx_per_cond[cond_num], [])).tolist()
        ]
        (_, self.reflection_SPRINT_data_without_transits[cond_num], self.transmission_SPRINT_data_without_transits[cond_num],
         self.BP_counts_SPRINT_data_without_transits[cond_num], self.DP_counts_SPRINT_data_without_transits[cond_num]) = \
            self.analyze_SPRINT_data_points(self.all_seq_without_transits, SPRINT_pulse_number=[1, 2],
                                            background=True)  # Enter the index of the SPRINT pulse for which the data should be analyzed

        # self.number_of_transits_live = len(self.all_transits_seq_indx_per_cond)
        # self.number_of_transits_total = len(
        #     [vec for lst in self.batcher['all_transits_seq_indx_batch'] for vec in lst])
        #
        # self.num_of_total_SPRINT_reflections = sum(self.reflection_SPRINT_data_without_transits)
        # self.num_of_total_SPRINT_transmissions = sum(self.transmission_SPRINT_data_without_transits)
        # self.num_of_total_SPRINT_BP_counts = sum(self.BP_counts_SPRINT_data_without_transits)
        # self.num_of_total_SPRINT_DP_counts = sum(self.DP_counts_SPRINT_data_without_transits)

    def save_experiment_results(self, experiment_comment, daily_experiment_comments):
        """
        Responsible for saving all the results gathered in the experiment and required by analysis
        Return a flag telling if save should occur at all.
        TODO: add 'save_experiment_results' to BaseExperiment, Make 'run' method invoke it at the end
        """

        # Save all other files
        results = {
            "folded_tt_S_cumulative_avg": self.folded_tt_S_cumulative_avg,
            "folded_tt_N_cumulative_avg": self.folded_tt_N_cumulative_avg,
            "folded_tt_BP_cumulative_avg": self.folded_tt_BP_cumulative_avg,
            "folded_tt_DP_cumulative_avg": self.folded_tt_DP_cumulative_avg,
            "folded_tt_FS_cumulative_avg": self.folded_tt_FS_cumulative_avg,
            "folded_tt_S_directional_cumulative_avg": self.folded_tt_S_directional_cumulative_avg,
            "folded_tt_N_directional_cumulative_avg": self.folded_tt_N_directional_cumulative_avg,
            "folded_tt_BP_timebins_cumulative_avg": self.folded_tt_BP_timebins_cumulative_avg,
            "folded_tt_DP_timebins_cumulative_avg": self.folded_tt_DP_timebins_cumulative_avg,

            "MZ_BP_counts_balancing_batch": self.batcher['MZ_BP_counts_balancing_batch'],
            "MZ_BP_counts_balancing_check_batch": self.batcher['MZ_BP_counts_balancing_check_batch'],
            "MZ_DP_counts_balancing_batch": self.batcher['MZ_DP_counts_balancing_batch'],
            "MZ_DP_counts_balancing_check_batch": self.batcher['MZ_DP_counts_balancing_check_batch'],
            "Phase_Correction_vec_batch": self.batcher['Phase_Correction_vec_batch'],
            "Phase_Correction_min_vec_batch": self.batcher['Phase_Correction_min_vec_batch'],
            "Phase_Correction_value": self.batcher['Phase_Correction_value'],
            "MZ_S_tot_counts": self.batcher['MZ_S_tot_counts'],

            "Index_of_Sequences_with_data_points": self.batcher['seq_with_data_points_batch'],
            "Reflections_per_data_point": self.batcher['reflection_SPRINT_data_batch'],
            "Transmissions_per_data_point": self.batcher['transmission_SPRINT_data_batch'],
            "Bright_port_counts_per_data_point": self.batcher['BP_counts_SPRINT_data_batch'],
            "Dark_port_counts_per_data_point": self.batcher['DP_counts_SPRINT_data_batch'],

            "FLR_measurement": self.batcher['flr_batch'],
            "lock_error": self.batcher['lock_err_batch'],
            "k_ex": self.batcher['k_ex_batch'],
            "interference_error": self.batcher['interference_error_batch'],
            "exp_timestr": experiment_comment,

            "exp_comment": f'transit condition: {self.transit_condition}; reflection threshold: {self.reflection_threshold} @ {int(self.reflection_threshold_time / 1e6)} ms',
            "daily_experiment_comments": daily_experiment_comments,

            "experiment_config_values": self.Exp_Values,

            "run_parameters": self.run_parameters
        }

        # Save the results
        if self.playback['active']:
            resolved_path = self.playback['save_results_path']
        else:
            resolved_path = self.bd_results.get_sequence_folder(sequence_definitions)
        self.bd_results.save_results(results, resolved_path)

    def results_to_dictionary(self):
        '''
        Take all the results from the batcher and summarize them into a dictionary.
        :return:
        '''

        dict = {}
        ## General data for background experiment ##
        # transmission data in detection pulses in sequences:
        dict['transmission_data_in_detection_pulses_per_seq_per_cycle'] = self.batcher[
             'num_of_det_transmissions_per_seq']
        dict['transmission_data_in_detection_pulses_per_seq_per_cycle_N'] = self.batcher[
             'num_of_det_transmissions_per_seq_N']
        dict['transmission_data_in_detection_pulses_per_seq_per_cycle_S'] = self.batcher[
             'num_of_det_transmissions_per_seq_S']
        # reflection data in detection pulses in sequences:
        dict['reflection_data_in_detection_pulses_per_seq_per_cycle'] = self.batcher[
             'num_of_det_reflections_per_seq']
        dict['reflection_data_in_detection_pulses_per_seq_per_cycle_N'] = self.batcher[
             'num_of_det_reflections_per_seq_N']
        dict['reflection_data_in_detection_pulses_per_seq_per_cycle_S'] = self.batcher[
             'num_of_det_reflections_per_seq_S']
        # transmission data in SPRINT pulses in sequences:
        dict['transmission_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
            'num_of_SPRINT_transmissions_per_seq']
        dict['transmission_data_in_SPRINT_pulses_per_seq_per_cycle_N'] = self.batcher[
            'num_of_SPRINT_transmissions_per_seq_N']
        dict['transmission_data_in_SPRINT_pulses_per_seq_per_cycle_S'] = self.batcher[
            'num_of_SPRINT_transmissions_per_seq_S']
        # reflection data in SPRINT pulses in sequences:
        dict['reflection_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
            'num_of_SPRINT_reflections_per_seq']
        dict['reflection_data_in_SPRINT_pulses_per_seq_per_cycle_N'] = self.batcher[
            'num_of_SPRINT_reflections_per_seq_N']
        dict['reflection_data_in_SPRINT_pulses_per_seq_per_cycle_S'] = self.batcher[
            'num_of_SPRINT_reflections_per_seq_S']
        # Bright port data in SPRINT pulses in sequence:
        dict['bright_port_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
            'num_of_BP_counts_per_seq_in_SPRINT_pulse']
        # Dark port data in SPRINT pulses in sequence:
        dict['dark_port_data_in_SPRINT_pulses_per_seq_per_cycle'] = self.batcher[
            'num_of_DP_counts_per_seq_in_SPRINT_pulse']

        ## Sorting and summarizing the data analysis for each condition ##
        for indx, condition in enumerate(self.transit_conditions):
            # Handle transits related data:
            dict[str(condition)]['number_of_transits'] = 0
            dict[str(condition)]['number_of_sequences_with_transits'] = 0
            dict[str(condition)]['indices_of_all_sequences_with_transits_per_cycle'] = []
            for lst in self.batcher['all_transits_seq_indx_per_cond']:
                for vec in lst[indx]:
                    # Get number of transits and number of sequences with transits
                    dict[str(condition)]['number_of_transits'] += 1
                    dict[str(condition)]['number_of_sequences_with_transits'] += len(vec)
                # Save index of sequences with transits per cycle:
                dict[str(condition)]['indices_of_all_sequences_with_transits_per_cycle'].append(lst[indx])
            dict[str(condition)]['number_of_sequences_with_transits'] -= dict[str(condition)]['number_of_transits']



        return dict


    # Class's constructor
    def __init__(self, exp_type='QRAM', exp_date='20230719', exp_time=None, transit_conditions=[[[2, 1, 2]]]):

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

        # Open folder and load to dictionary
        root = Tk()
        self.exp_data_path = '{}'.format(askdirectory(title='Experiment folder', initialdir=r'U:\Lab_2023\Experiment_results'))
        self.data = Experiment_data_load.DictionaryBuilder()
        self.Exp_dict = self.data.load_files_to_dict(self.exp_data_path)
        messagebox.showinfo(title='Success!', message='Experiment data is ready to use.')
        root.destroy()

        # Divide data to background and experiment data:
        exp_key = list(self.Exp_dict.keys())[0]
        for key in list(self.Exp_dict.keys()):
            if 'without' in key.lower():
                self.Background_dict = self.Exp_dict[key]
            else:
                exp_key = key
        self.Exp_dict = self.Exp_dict[exp_key]

        # background analysis:
        # check number of cycles in experiment:
        number_of_cycles = len(list(self.Background_dict['output'][list(self.Background_dict['output'].keys())[0]].values())[0])
        self.init_params_for_experiment(self.Background_dict)
        # Initialize the batcher
        self.batcher.set_batch_size(number_of_cycles)
        self.batcher.empty_all()
        # self.Background_dict['Analysis_results'] = {}
        for cycle in tqdm(range(number_of_cycles)):
            self.ingest_time_tags(self.Background_dict, cycle)
            self.experiment_calculations()
            # self.calculate_running_averages(cycle+1)
            for condition_number, transit_condition in enumerate(self.transit_conditions):
                ### Find transits and extract SPRINT data:  ###
                self.get_transit_data(transit_condition, condition_number)
            self.batcher.batch_all(self)
            self.Background_dict['Analysis_results'] = self.results_to_dictionary()

if __name__ == '__main__':
    self = experiment_data_analysis(transit_conditions=[[[2, 1], [1, 2]], [[2, 1, 2]]])