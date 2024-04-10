import Experiment_data_load
import numpy as np
import os
import pymsgbox

def popupbox_inquiry(message=''):

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


if __name__ == '__main__':

    exp_type, exp_date, exp_time = popupbox_inquiry()
    while True:
        if exp_type is None or exp_date is None or exp_time is None:
            option = pymsgbox.confirm('Missing an input, do you wish to retry?', 'Experiment data loading',
                                      ['Yes please', 'No thanks'])
            if option == 'Yes please':
                exp_type, exp_date, exp_time = popupbox_inquiry()
            else:
                break
        else:
            # print(exp_type, exp_date, exp_time)
            data = Experiment_data_load.DictionaryBuilder(exp_type, exp_date, exp_time)
            if data is None:
                exp_type, exp_date, exp_time = popupbox_inquiry()
            else:
                Exp_dict = data.load_files_to_dict(data.exp_path)
                pymsgbox.alert(text='Experiment data is ready to use.', title='Success!')
                break
    #
    # exp_type = 'QRAM'
    # exp_date = '20230719'
    # exp_time = '011453'
    # data = Experiment_data_load.DictionaryBuilder(exp_type, exp_date, exp_time)
    # pass
    exp_sequence_len = len(Exp_dict['input']['sequences']['South_sequence_vector'])

    folded_tt_S = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_N = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_BP = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_DP = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_FS = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_S_directional = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_N_directional = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_BP_timebins = np.zeros(exp_sequence_len, dtype=int)
    folded_tt_DP_timebins = np.zeros(exp_sequence_len, dtype=int)

    for x in [elem for lst in Exp_dict['output']['South(5)']['South_timetags'] for elem in lst]:
        folded_tt_S[x % exp_sequence_len] += 1
    for x in [elem for lst in Exp_dict['output']['North(8)']['North_timetags'] for elem in lst]:
        folded_tt_N[x % exp_sequence_len] += 1
    for x in [elem for lst in Exp_dict['output']['Bright(1,2)']['Bright_timetags'] for elem in lst]:
        folded_tt_BP[x % exp_sequence_len] += 1
    for x in [elem for lst in Exp_dict['output']['Dark(3,4)']['Dark_timetags'] for elem in lst]:
        folded_tt_DP[x % exp_sequence_len] += 1
    for x in [elem for lst in Exp_dict['output']['FastSwitch(6,7)']['FS_timetags'] for elem in lst]:
        folded_tt_FS[x % exp_sequence_len] += 1

    folded_tt_S_directional = (np.array(folded_tt_S) + np.array(folded_tt_FS))
    folded_tt_N_directional = (np.array(folded_tt_N) + np.array(folded_tt_BP) + np.array(folded_tt_DP))
