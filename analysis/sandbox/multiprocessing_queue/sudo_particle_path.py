# module
from numpy import arange
from h5py import File

array_size = 25

def create_pp(print_keys = False): 

    '''
    Create particle path sudo data 
    - set print to True to view keys
    '''

    in_file = 'particle_path_data.hdf5'

    tof_data = []
    for i in range(12):
        arr = arange(array_size)*0 + i + 1 
        tof_data.append(arr)

    names = [
            'T1_prop_dist', 
            'T1_endpoint_dist',
            'T1_prop_times',
            'T1_interactions',
            'T1_part_ids',
            'T1_input_times',
            'T4_prop_dist',
            'T4_endpoint_dist',
            'T4_prop_times',
            'T4_interactions',
            'T4_part_ids',
            'T4_input_times',
            ]

    # print(len(names), len(tof_data))
    with File(in_file, 'w', track_order=True) as f: 
        # loop through the names and data lists 
        # to create a sudo file
        for i in range(12): 
            f.create_dataset(names[i], data = tof_data[i])

        # verify the data was stored 
        if print_keys == True: 
            for key in f.keys(): 
                print(key, f[key][0:5])
    print('order of keys: input sort matches stored sort')

    
    return in_file

def check_contents(in_file): 

    ''' Check contents of hdf5 files '''
    with h5py.File(in_file, 'r') as f: 
        for key in f.keys(): 
            # Print the name and the 1st 5 elements
            print(key, f[key][0:5])