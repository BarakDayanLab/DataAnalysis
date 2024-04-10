import numpy as np
from tqdm import tqdm


def bin_and_reshape(timetag_streams, sequence_length, cycle_length, bins):
    N_cycles = cycle_length // sequence_length
    bins_ = (np.arange(N_cycles).reshape((-1, 1)) * sequence_length + bins.reshape((1, -1))).reshape(-1)
    bins_ = np.append(bins_, cycle_length)
    reshaped = np.array(
        [
            [
                np.histogram(timetag_array, bins=bins_)[0].reshape((N_cycles, -1)) for timetag_array in
                tqdm(timetags, desc='binning stream', leave=False)
            ]
            for timetags in tqdm(timetag_streams, desc='binning streams')
        ]
    )
    reshaped = reshaped[:, :, :, ::2]
    return reshaped


if __name__ == '__main__':
    x = np.arange(16).reshape((4, 4))
    print(np.lib.stride_tricks.sliding_window_view(x, (1, 2)))

    # print(np.logical_or(np.array(False), np.array(False), np.array(False), np.array(False)))

