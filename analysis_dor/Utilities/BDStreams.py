import os
import re
import time
import json
import base64
import struct
import numpy as np
from io import BytesIO
from Utilities.Utils import Utils


class NumpyEncoder(json.JSONEncoder):
    """
    The class is used upon json.dump and json.dumps functions to encode an object member.
    If it is a numpy object we need to encode, we will do it using numpy 'save', otherwise we return the object to let
    the standard decoder open it

    Encoding: (1) Numpy save is used to create bytes (2) Decode bytes to utf-8 (3) Encode to Base64
    """

    def default(self, obj):
        # If input object is a ndarray it will be converted into a dict holding the nparray bytes
        # TODO: do we want to handle individually np.ndarray/np.void/np.int64/np.float64 ? via isinstance(obj, nd,void)
        if hasattr(obj, 'dtype'):
            np_bytes = BytesIO()
            np.save(np_bytes, obj.data, allow_pickle=True)
            base64_bytes = base64.b64encode(np_bytes.getvalue())
            base64_string = base64_bytes.decode('utf-8')
            return dict(__nparray__=True, __base64__=base64_string)

        # Let the base class default method work on the obj
        return super().default(obj)
        #return json.JSONEncoder.default(self, obj)

def json_numpy_obj_hook(dct):
    """
    Function is used upon json.load and json.loads functions to decode an object member.
    If it is a numpy object we decoded, we will open it using numpy 'load', otherwise we return the object to let
    the standard decoder open it

    Decoding: (1) Encode bytes from utf-8 (2) Decode from Base64 (3) Numpy load will reconstruct it
    """
    # Is this is a numpy object we encoded (e.g. has __nparray__ attribute), we need to open the bytes
    if isinstance(dct, dict) and '__nparray__' in dct:
        base64_bytes = dct['__base64__'].encode('utf-8')
        data = base64.b64decode(base64_bytes)
        loaded_np = np.load(BytesIO(data), allow_pickle=True)
        return loaded_np

    return dct


class BDStreams:

    def __init__(self, streams=None, save_raw_data=False, logger=None):

        # Flag - whether we should save raw data or not
        self.save_raw_data = save_raw_data

        if logger is not None:
            self.logger = logger

        self.streams = streams
        self.number_of_rows_saved = 0

        pass

    def __del__(self):
        pass

    def set_streams_definitions(self, streams):
        self.streams = streams
        pass

    def get_results_from_comm_channels(self, comm_messages):

        for stream in self.streams.values():
            if 'source' in stream and stream['source'] == 'comm':
                comm_name = stream['name'].lower()
                channel = stream['channel']
                value = stream['default']
                if comm_messages != {} and comm_name in comm_messages[channel] and comm_messages[channel][comm_name] is not None:
                    value = comm_messages[channel][comm_name]
                stream['results'] = value
        pass

    def get_results_from_opx_streams(self):
        """
        Read all data from OPX streams
        """
        if self.streams is None:
            return

        for stream in self.streams.values():
            if 'handler' in stream and stream['handler'] is not None:
                stream['handler'].wait_for_values(1)

        for stream in self.streams.values():
            if 'handler' in stream and stream['handler'] is not None:
                stream['results'] = stream['handler'].fetch_all()
            else:
                stream['results'] = None

        pass

    def set_opx_handlers(self, opx_job, playback=False):

        # Get a Handler for all the opx-related streams
        for key, value in self.streams.items():
            if value['source'] == "opx":
                if playback:
                    value['handler'] = "playback: data not loaded yet..."
                else:
                    value['handler'] = opx_job.result_handles.get(key)

        return self.streams

    """
        There are two vectors we get from the OPX/OPD: Counts and Time-Tags
        For each "cycle", we get two outputs:

        +---------------+------------------+
        +  Counts       |  Time-Tags Array |
        +---------------+------------------+

        Counts:
        +------+------+------+------+------+------+
        +  12  |  9   |   1  |  23  |  54  | ...  +
        +------+------+------+------+------+------+

        Timetags:
        +------+------+------+------+------+------+
        +  23  |  79  |  91  | 1023 | 2354 | ...  +
        +------+------+------+------+------+------+

        For example: the photon at position 0 was received at time 23ns on a time-window of 10e7 nano seconds
        The size of this array is dynamic - it holds values as the number of counts received

        NOTE: we have today an issue that the stream may still hold values from the previous detections        
    """

    def normalize_stream(self, stream_data, detector_index, detector_delays, M_window):
        """
        Given a stream data, normalize it:
        - Remove the first element - the count
        - Remove time-tags who are exactly on the M_window
        - Remove junk values coming from the OPX (9999984)
        - Remove time-tags that after adding the delay will be outside the window
        """
        stream_len = stream_data[0]
        delay = detector_delays[detector_index]

        # Iterate starting from the first position after the count and build the normalized stream
        normalized = []
        for tt_index in range(1, stream_len+1):
            data = stream_data[tt_index]
            if data % M_window == 0 or data == 9999984:
                continue

            delayed_time_tag = data + delay
            if delayed_time_tag < M_window:
                normalized.append(delayed_time_tag)

        normalized.sort()
        return normalized

    def clean_streams(self):
        for stream in self.streams.values():
            stream['all_rows'] = []
            stream['results'] = []

    def load_entire_folder(self, playback_files_path):
        """
        Iterates over all playback files in folder and loads them into memory (streams object)
        """

        # Clean all streams
        self.clean_streams()

        # Get all files in playback folder (ignore Python source files and directories)
        playback_files = Utils.get_files_in_path(path=playback_files_path, opt_out_filter='.py')

        # Take time
        load_time_start = time.time()

        # Iterate over all files in folder
        i = 1
        for playback_file in playback_files:
            if i % 100 == 0:
                self.logger.info(f'Loading playback files ({i}/{len(playback_files)})')
            i += 1
            self.load_streams(os.path.join(playback_files_path, playback_file))

        self.logger.info(f'Finished loading {len(playback_files)} playback files. Took {time.time() - load_time_start} secs')
        pass

    def load_streams(self, data_file):

        try:
            # Load the json from file
            with open(data_file, 'r') as file:
                loaded_json = json.load(file, object_hook=json_numpy_obj_hook)

            # Populate streams with the data loaded
            for stream_data in loaded_json['streams_data'].values():
                name = stream_data['name']
                results = stream_data['data']
                stream = self.streams[name]
                if 'all_rows' not in stream:  # If it's the first results we're adding, create a new array
                    stream['all_rows'] = [results]
                else:
                    stream['all_rows'].append(results)

                # TODO: we should have an individual timestamp per stream
                stream['timestamp'] = loaded_json['timestamp']

        except Exception as err:
            self.logger.warn(f'Failed to load playback data file "{data_file}". {err}')

        pass

    def load_streams_DEP(self, data_file):

        # Open the binary file
        with open(data_file, 'rb') as file:

            try:
                bytes = file.read(4)
                bytes_unpacked = struct.unpack('>i', bytes)[0]
                current_time = int(bytes_unpacked)

                number_of_streams = int(struct.unpack('>b', file.read(1))[0])
                for i in range(0, number_of_streams):

                    # Get stream name len and then name
                    name_len = struct.unpack(f'>I', file.read(4))[0]
                    name = struct.unpack(f'>{name_len}s', file.read(name_len))[0]
                    name = name.decode()

                    # Get data len and then data
                    data_len = struct.unpack('>I', file.read(4))[0]
                    data_bytes = file.read(data_len)

                    # Reconstruct the numpy array from bytes
                    load_bytes = BytesIO(data_bytes)
                    loaded_np = np.load(load_bytes, allow_pickle=True)

                    # Append the data we got to all rows
                    stream = self.streams[name]
                    if 'all_rows' not in stream:  # If it's the first results we're adding, create a new array
                        stream['all_rows'] = [loaded_np]
                    else:
                        stream['all_rows'].append(loaded_np)

                    # TODO: we should append these - to create the sequence of time stamps recorded. Now we override.
                    stream['timestamp'] = current_time

            except Exception as err:
                self.logger.warn(f'Failed to load playback data file "{data_file}". {err}')
        pass

    def construct_streams_data_as_json(self, time_formatted):

        # Set the global data for all streams
        data_to_save = {
            "timestamp": time_formatted,
            "streams_data": {}
        }

        # Iterate over all streams and add their name and data to the json
        for name, stream in self.streams.items():
            if 'save_raw' in stream and stream['save_raw']:
                data_to_save['streams_data'][name] = {
                    "name": name,
                    "data": stream['results']
                }

        return data_to_save

    def save_streams(self, playback_files_path):

        # Do we need to save raw data at all?
        if not self.save_raw_data:
            return

        # If there are no streams defined, no playback data to save, ignore.
        if self.streams is None:
            return

        # Format the file name for the playback
        time_formatted = time.strftime("%Y%m%d_%H%M%S")
        save_name = os.path.join(playback_files_path, f'{time_formatted}_streams.json')

        # Construct the data to save
        data_to_save = self.construct_streams_data_as_json(time_formatted)

        try:
            current_time = time.time()

            # Save the file
            with open(save_name, "w") as file:
                json.dump(data_to_save, file, indent=4, cls=NumpyEncoder)

            self.number_of_rows_saved += 1

            total_prep_time = time.time() - current_time
        except Exception as err:
            self.logger.warn(f'Failed to save raw data [{save_name}]: {err}. Skipping.')

    def save_streams_DEP(self, playback_files_path):
        """

                        +-----------------+-------------+
        Stream Name ==> | Stream Name Len | Stream Name |
                        +-----------------+-------------+

                        +-----------------+-------------+
        Stream Data ==> | Stream Data Len | Stream Data +
                        +-----------------+-------------+

        +-----------+-------------+---------------+---------------+-------+---------------+---------------+
        + Timestamp | Num Streams | Stream-1 Name | Stream-1 Data | ..... | Stream-N Name | Stream-N Data |
        +-----------+-------------+---------------+---------------+-------+---------------+---------------+

        """
        if not self.save_raw_data:
            return

        # If there are no streams defined, no playback data to save, ignore.
        if self.streams is None:
            return

        try:
            current_time = time.time()
            bytes_array = bytearray()

            # Get only those streams that their configuration indicates they need to be saved
            streams_to_save = [s for s in self.streams.values() if 'save_raw' in s and s['save_raw']]

            # Start packing the data with 'b' (1-byte=unsigned-char) for number of streams
            bytes_array += struct.pack('>ib', int(current_time), len(streams_to_save))

            # Iterate over all streams and pack their name and data
            for stream in streams_to_save:

                # Pack the name (I=unsigned int for len and then name)
                name = stream['name']
                bytes_array += struct.pack(f'>I{len(name)}s', len(name), bytes(name, 'utf-8'))

                # Pack the data
                data = [] if stream['handler'] is None else stream['results']

                # Save in to BytesIo buffer
                np_bytes = BytesIO()
                np.save(np_bytes, data, allow_pickle=True)

                # get bytes value
                np_bytes = np_bytes.getvalue()
                bytes_array += (struct.pack(f'>I', len(np_bytes)) + np_bytes)

            # Format the file name for the playback
            time_formatted = time.strftime("%Y%m%d_%H%M%S")
            save_name = os.path.join(playback_files_path, f'{time_formatted}_streams.dat')

            # Save the file
            with open(save_name, "wb") as file:
                file.write(bytes_array)

            self.number_of_rows_saved += 1

            total_prep_time = time.time() - current_time

        except Exception as err:
            self.logger.warn(f'Failed to save raw data [{save_name}]: {err}. Skipping.')

        pass

