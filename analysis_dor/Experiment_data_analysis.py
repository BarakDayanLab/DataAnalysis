import Experiment_data_load
import numpy as np
import os
import math
import pymsgbox
import matplotlib.pyplot as plt

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
                date = pymsgbox.prompt(text='Please enter the date of the experiment \nall characters must be integers!',
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
                time = pymsgbox.prompt(text='Please enter the time of the experiment \nall characters must be integers!',
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

    def init_params_for_experiment(self):

        self.sequence_len = len(self.Exp_dict['input']['sequences']['South_sequence_vector'])
        self.M_window = self.Exp_dict['meta_data']['Exp_config_values']['M_window']
        self.Pulses_location_in_seq = self.Exp_dict['input']['sequences']['Pulses_location_in_seq']
        self.MZ_delay = 250  # [ns]

        self.number_of_PNSA_sequences = math.ceil(self.M_window / self.sequence_len)

        self.number_of_detection_pulses_per_seq, self.number_of_SPRINT_pulses_per_seq = self.number_of_pulses_per_seq()

        self.end_of_det_pulse_in_seq = int(self.Pulses_location_in_seq[self.number_of_detection_pulses_per_seq-1][1]) + 6
        # define empty variables

        self.tt_measure = []
        self.tt_S_measure = []

        self.folded_transmission = np.zeros(self.sequence_len)
        self.folded_reflection = np.zeros(self.sequence_len)

        self.tt_S_binning = np.zeros(self.number_of_PNSA_sequences + 1)
        self.seq_transit_events_live = np.zeros(self.number_of_PNSA_sequences)
        self.seq_transit_events_batched = np.zeros(self.number_of_PNSA_sequences)
        self.tt_S_SPRINT_events = np.zeros(self.sequence_len)
        self.tt_S_SPRINT_events_batch = np.zeros(self.sequence_len)
        self.num_of_det_reflections_per_seq_accumulated = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_accumulated = np.zeros(self.number_of_PNSA_sequences)

        num_of_seq_per_count = 50
        self.num_of_BP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences//num_of_seq_per_count)
        self.num_of_DP_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences//num_of_seq_per_count)
        self.num_of_S_counts_per_n_sequences = np.zeros(self.number_of_PNSA_sequences//num_of_seq_per_count)

        self.num_of_det_reflections_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_reflections_per_seq_N = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_S = np.zeros(self.number_of_PNSA_sequences)
        self.num_of_det_transmissions_per_seq_N = np.zeros(self.number_of_PNSA_sequences)

        self.seq_with_data_points = []
        self.reflection_SPRINT_data = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        self.transmission_SPRINT_data = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.
        self.BP_counts_SPRINT_data = []
        self.DP_counts_SPRINT_data = []

        self.reflection_SPRINT_data_without_transits = []  # Array of vectors with data on the number of reflections per SPRINT pulse in sequence.
        self.transmission_SPRINT_data_without_transits = []  # Array of vectors with data on the number of transmissions per SPRINT pulse in sequence.
        self.BP_counts_SPRINT_data_without_transits = []
        self.DP_counts_SPRINT_data_without_transits = []

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
            self.folded_tt_S[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['North(8)']['North_timetags'] for elem in lst]:
            self.folded_tt_N[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['Bright(1,2)']['Bright_timetags'] for elem in lst]:
            self.folded_tt_BP[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['Dark(3,4)']['Dark_timetags'] for elem in lst]:
            self.folded_tt_DP[x % exp_sequence_len] += 1
        for x in [elem for lst in self.Exp_dict['output']['FastSwitch(6,7)']['FS_timetags'] for elem in lst]:
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

        self.num_of_det_reflections_per_seq_S_,\
        self.num_of_det_reflections_per_seq_N_, \
        self.num_of_det_transmissions_per_seq_S_, \
        self.num_of_det_transmissions_per_seq_N_ = [
            [
                [
                    [] for _ in range(self.number_of_detection_pulses_per_seq)
                ]
                for _ in range(self.number_of_PNSA_sequences)
            ] for _ in range(4)
        ]

        self.num_of_SPRINT_reflections_per_seq_S_, \
        self.num_of_SPRINT_reflections_per_seq_N_, \
        self.num_of_SPRINT_transmissions_per_seq_S_, \
        self.num_of_SPRINT_transmissions_per_seq_N_, \
        self.num_of_BP_counts_per_seq_in_SPRINT_pulse, \
        self.num_of_DP_counts_per_seq_in_SPRINT_pulse = [
            [
                [
                    [] for _ in range(self.number_of_SPRINT_pulses_per_seq)
                ]
                for _ in range(self.number_of_PNSA_sequences)
            ] for _ in range(6)
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
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)

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
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
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
                            if tup[2] == 'S':
                                self.num_of_det_reflections_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_transmissions_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_reflections_per_seq_S_[seq_num][ind].append(element)
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
                            if tup[2] == 'S':
                                self.num_of_det_transmissions_per_seq_S_[seq_num][indx].append(element)
                        else:
                            ind = indx - self.number_of_detection_pulses_per_seq
                            if tup[2] == 'n':
                                self.num_of_SPRINT_reflections_per_seq_N_[seq_num][ind].append(element)
                            if tup[2] == 's':
                                self.num_of_SPRINT_transmissions_per_seq_S_[seq_num][ind].append(element)

                if tt_inseq <= self.end_of_det_pulse_in_seq:  # The part of the detection pulses in the sequence
                    seq_num = (element - 1) // self.sequence_len
                    self.num_of_det_reflections_per_seq_N[seq_num] += self.filter_N[tt_inseq]
                    self.num_of_det_transmissions_per_seq_S[seq_num] += self.filter_S[tt_inseq]


    # Class's constructor
    def __init__(self, exp_type='QRAM', exp_date='20230719', exp_time=None):

        self.exp_type, self.exp_date, self.exp_time = self.popupbox_inquiry()
        while True:
            if self.exp_type is None or self.exp_date is None or self.exp_time is None:
                option = pymsgbox.confirm('Missing an input, do you wish to retry?', 'Experiment data loading',
                                          ['Yes please', 'No thanks'])
                if option == 'Yes please':
                    self.exp_type, self.exp_date, self.exp_time = self.popupbox_inquiry()
                else:
                    break
            else:
                # print(exp_type, exp_date, exp_time)
                self.data = Experiment_data_load.DictionaryBuilder(self.exp_type, self.exp_date, self.exp_time)
                if self.data is None:
                    self.exp_type, self.exp_date, self.exp_time = self.popupbox_inquiry()
                else:
                    self.Exp_dict = self.data.load_files_to_dict(self.data.exp_path)
                    pymsgbox.alert(text='Experiment data is ready to use.', title='Success!')
                    break

        self.init_params_for_experiment()
        self.fold_tt_histogram(self.sequence_len)
        if self.Pulses_location_in_seq[0][2] == 'N':
            print(sum(self.folded_tt_S_directional[int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])])/
                  sum(self.folded_tt_N_directional[int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])]))
        else:
            print(sum(self.folded_tt_N_directional[int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])])/
                  sum(self.folded_tt_S_directional[int(self.Pulses_location_in_seq[0][0]):int(self.Pulses_location_in_seq[0][1])]))


if __name__ == '__main__':

    analize = experiment_data_analysis()