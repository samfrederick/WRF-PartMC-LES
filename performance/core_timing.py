import sys, os
import pandas as pd

if __name__ == '__main__':
    path = sys.argv[1]
    elapsed = pd.read_csv(os.path.join(path, 'elapsed.txt'), header=None)
    elapsed = elapsed.replace(' elapsed seconds' , '', regex=True)
    elapsed['Duration(s)']= elapsed[0].astype('float')
    elapsed = elapsed.drop(columns=0)

    clock = pd.read_csv(os.path.join(path, 'clock.txt'), header=None)
    clock = clock.replace(' on domain', '', regex=True)
    clock['Domain_Time'] = pd.to_datetime(clock[0], format='%Y-%m-%d_%H:%M:%S')
    clock = clock.drop(columns=0)

    data = clock.join(elapsed)
    data.to_csv('core_timing.csv', index=False)