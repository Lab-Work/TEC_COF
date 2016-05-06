__author__ = 'Yanning Li'

# This function simulates the virtual MATLAB to test AIMSUN signals

from os.path import exists
import time
import os
import sys
import numpy as np



files = {}

# two hours simulation
all_signal = np.array([50, 0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550]*10)


# MAIN
def main(argv):

    com_file = 'E:\\AIMSUN_MATLAB_COM\\COM_CONFIG.txt'

    readComConfigFile(com_file)

    readNet()

    initMatlab()

    while not simulationCompleted():

        cur_time = readData()

        end_time = 7200

        new_signal = generateNewSignal(cur_time, end_time)

        writeSignal(new_signal, cur_time)

        time.sleep(0.1)





# read the COMCONFIG file
# Here are the utility functions we defined for properly read control signals.
def readComConfigFile(com_file):

    global files

    f = open(com_file, 'r')

    for line in f:
        if line[0] == '%':
            # comment line, skip
            continue
        else:
            line = line.strip()
            items = line.split(',')

            # There are four keys:
            # data; data_write_done; signal; signal_write_done; matlab_init_done
            files[items[0]] = items[1]

    f.close()

    print 'MATLAB: Finished reading communication configuration file.'




# read the network
def readNet():

    while not exists(files['net_write_done']):

        print 'MATLAB: Waiting for network file...'

        time.sleep(1)

    print 'MATLAB: Finished reading network file.'



# initialize MALTAB
def initMatlab():

    f = open(files['matlab_init_done'],'w')

    f.write('matlab init done')

    f.close()

    print 'Finished initialize MATLAB'



# read data; find out the current time in the simulator
def readData():

    # wait for new data
    while not exists(files['simulation_completed']):

        print 'MALTAB: Waiting for data file...'
        time.sleep(1)

        if exists(files['data_write_done']):

            f = open(files['data'],'r')

            for line in f:

                line = line.strip()

                if line[0] == '%':
                    continue

                items = line.split(',')
                cur_time = float(items[2])
                break


            f.close()

            os.remove(files['data'])
            os.remove(files['data_write_done'])

            return cur_time


# simulation finished file flag
def simulationCompleted():

    if exists(files['simulation_completed']):
        return True


# Generate signal.
# write 10 minutes preset signal to a formatted file
def generateNewSignal(cur_time, end_time):

    # write 10 mins (10 instances every time)
    start_index = int(cur_time)/60 - 1
    end_index = min(cur_time+10, int(end_time)/60)

    signal_to_write = all_signal[start_index:end_index]

    signal_file = []

    for i in range(0,len(signal_to_write)):

        line = '{0},60,{1}\n'.format(cur_time+i*60, signal_to_write[i])

        signal_file.append(line)

    return signal_file



# write signal into files
def writeSignal(new_signal, cur_time):

    if len(new_signal) != 0:

        f = open(files['signal'],'wb')

        for line in new_signal:
            f.write(line)

        f.close()

        time.sleep(0.5)

        # write flag
        f = open(files['signal_write_done'], 'w')
        f.write('signal write done')
        f.close()

        print 'MATLAB: Wrote new signal at {0}'.format(cur_time)





if __name__ == "__main__":
    sys.exit(main(sys.argv))