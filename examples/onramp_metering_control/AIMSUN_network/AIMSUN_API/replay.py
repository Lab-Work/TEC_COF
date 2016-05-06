__author__ = 'Yanning Li'

# This script reads the signal and replay the simulation


from os.path import exists
import time as time_module
import os
import sys
from AAPI import *

SITEPACKAGES = "C:\\Python27\\Lib\\site-packages"
sys.path.append(SITEPACKAGES)
import numpy as np


_debug = False

# import the shared file in this folder
com_file = 'E:\\AIMSUN_MATLAB_COM\\COM_CONFIG.txt'

# The dictionary that saves all the paths for the share communication files
files = {}

# traffic signal array.
signal_array = None

# the detection interval
# This is the hardware configuration of the interval of the detectors
det_interval = AKIDetGetIntervalDetection()

# This is the interval used to write the data
# e.g, detectors collect data at 1 min interval, and output 5 rows of data every 5 min
data_aggregate_interval = None
# the data to be write every data_aggregate_interval; a list of lines
data_to_write = []

# The unit used in this simulation
# BUG: It seems to be a bug. According to the manual, 1-metric; 0-english; However, no mater which unit was used, the value
# is always 0.
unit_AIMSUN = AKIInfNetGetUnits()
MILE2METER = 1609.34

# Save the ids of objects in the network
meter_id = 0
det_dict = {}
link_id = {}

meter_qmax = 1500.0/3600.0

# Manual configuration
# Since we are working on a very simple network example. It is not necessary to make this example code working for any
# general networks. In this case, just manually set some information here if it is not easily done automatically
# 1. the detector ids
# 2. the up/dn stream links and the onramp id
def manualConfig():

    # The following id must confirmed before simulation
    global data_aggregate_interval

    det_dict['us_in'] = 394
    det_dict['us_out'] = 395
    det_dict['ds_in'] = 400
    det_dict['ds_out'] = 403
    det_dict['on_in'] = 397
    det_dict['on_out'] = 396

    link_id['us'] = 329
    link_id['ds'] = 330
    link_id['on'] = 390

    data_aggregate_interval = 60


# This function is called when this module is correctly loaded
def AAPILoad():
    AKIPrintString( "AAPILoad" )
    return 0

# This function is called to perform any initialization
# 1. It sets up the communication between MATLAB
# 2. It exports the network topology to MATLAB.
# 3. It waits for MATLAB to initialize
def AAPIInit():
    AKIPrintString( "AAPIInit" )

    # manually set the detectors
    manualConfig()

    readComConfigFile()

    AKIPrintString( "Completed reading com config file." )

    initializeMeter()

    AKIPrintString( "Completed initializing the meter." )

    # if unit_AIMSUN == 1:
    #     AKIPrintString( "The system unit is metric" )
    # elif unit_AIMSUN == 0:
    #     AKIPrintString( "The system unit is english" )
    #

    # export the network to MATLAB
    exportNetwork()

    AKIPrintString( "Completed exporting network." )

    # read the signal file
    print files['all_signal']
    if not exists(files['all_signal']):

        AKIPrintString( "Error: all_signal.txt not exists..." )

    else:
        readAllSignalFile(0)
        # now start simulation
        AKIPrintString( "Have read all signals" )

    return 0


# This function is called before every time step
# 1. It should first check if a new signal file is written by MATLAB
# 2. If written, update is signal array.
# 3. Set the meter according to the signal array.
def AAPIManage(time, timeSta, timeTrans, acycle):

    if ECIGetTypeMeteringById(meter_id) != 1:
        AKIPrintString( "Warning: the meter type is not Green." )
    else:

        # have signals to apply
        if signal_array is not None:

            # need to update signal
            apply_signal_index = np.where(signal_array[:,0]==time)[0]
            if len(apply_signal_index) != 0:
                # the signal updating time point
                cycle_dur = signal_array[apply_signal_index[0],1]
                new_flow = signal_array[apply_signal_index[0],2]
                green_time = cycle_dur*new_flow/(meter_qmax*3600)

                # we approximate the flow by green time over signal duration cycle.
                ECIChangeParametersGreenMeteringById( meter_id, timeSta, cycle_dur, green_time, 0, cycle_dur,
                                                      0, 0)
                # ECIChangeParametersFlowMeteringById(meter_id, timeSta, 900, new_flow, 0)

                AKIPrintString( "---- Time {0}: applied signal flow {1} veh/hr: {2}/{3}.".format(time, new_flow, green_time, cycle_dur) )

            # printSignalParameters(timeSta)

    return 0


# This function should read the detector data
# 1. Check if the detection interval is reached
# 2. If yes, then read the detector count data.
# 3. Write the data into files.
# The file format [report_time, flow (veh/hr), speed (kph)]
def AAPIPostManage(time, timeSta, timeTrans, acycle):
    # AKIPrintString( "AAPIPostManage" )

    return 0


def AAPIFinish():
    # AKIPrintString( "AAPIFinish" )

    f = open(files['simulation_completed'],'wb')
    f.write('simulation completed')
    f.close()

    return 0

def AAPIUnLoad():
    # AKIPrintString( "AAPIUnLoad" )
    return 0




# This function read the network configuration information and write to a file for MATLAB
# 1. It read the network topology, including the links and juncs
# 2. All read information will be exported to a txt file, which can be read by MATLAB.
# 3. See the net file for the format
def exportNetwork():

    # net is a list of strings
    net = []

    # read link section information
    num_section = AKIInfNetNbSectionsANG()
    if _debug:
        AKIPrintString( '-- The number of sections:{0}\n'.format(num_section))

    for i in range(0, num_section):
        sec_id = AKIInfNetGetSectionANGId(i)

        # read link information
        link_info = AKIInfNetGetSectionANGInf(sec_id)

        if _debug:
            AKIPrintString( '--- Link {0}: speedLimit, {1} kph; capacity {2} veh/h; length {3} meters; num_lanes {4}'.format(
                link_info.id, link_info.speedLimit, link_info.capacity, link_info.length, link_info.nbCentralLanes
            ) )

            for j in range(0, link_info.nbTurnings ):

                AKIPrintString('------ {0} turns to {1}'.format(link_info.nbTurnings,
                                                                AKIInfNetGetIdSectionANGDestinationofTurning(sec_id, j)))

        # We know the link id from the AIMSUN GUI
        if sec_id == link_id['us'] or sec_id == link_id['ds']:

            link_type = 'freeway'

            # all unit in meters and seconds
            vf= 24.6713
            w= -3.5
            kc_pl= 0.0274
            km_pl= 0.2203
            qmax_pl= 0.6767
            v_min= 0
            v_max = 33.5208

        elif sec_id == link_id['on']:

            link_type = 'onramp'

            # all unit in meters and seconds
            # TODO: to be calibrated
            vf= 16.3582
            w= -2.7979
            kc_pl= 0.0322
            km_pl= 0.2204
            qmax_pl= 0.5266
            v_min= 0
            v_max = 33.5208

        else:

            AKIPrintString('\nERROR: check the link ids\n')

        # Make sure AIMSUN is using metric km
        lengthKM = link_info.length/1000.0

        line = 'link;link_id,{0};num_lanes,{1};lengthKM,{2};link_type,{3};'.format(
            sec_id, link_info.nbCentralLanes, lengthKM, link_type) + \
                'vf,{0};w,{1};kc_pl,{2};km_pl,{3};qmax_pl,{4};v_min,{5};v_max,{6}\n'.format(
                     vf, w, kc_pl, km_pl, qmax_pl, v_min, v_max )

        net.append(line)

    # collect the junction info (only one junction)
    junc_id = AKIInfNetGetJunctionId(0)
    ratio = [1,1]
    line = 'junc;junc_id,{0};inlink,{1},{2};outlink,{3};assignment_ratio,{4},{5};junc_type,onrampjunc'.format(
        junc_id, link_id['us'], link_id['on'], link_id['ds'], ratio[0], ratio[1]
    )

    net.append(line)

    AKIPrintString('-- Finished extracting network')

    # export the net to txt file
    f = open(files['net'],'wb')
    for line in net:
        f.write(line)

    f.close()

    AKIPrintString('-- Finished writing network')

    # flag write done
    f = open(files['net_write_done'],'wb')
    f.write('net write done\n')
    f.close()

    AKIPrintString('-- Finished writing net_write_done')


# This function initialize the Meter
# 1. It reads the Meter id and set the global variable
# 2. It verifies the configuration of the meter
def initializeMeter():

    global meter_id

    # read meter information
    num_meter = ECIGetNumberMeterings()
    if _debug:
        AKIPrintString( '-- Number of meters: {0}'.format(num_meter))

    # for each meter. We only have one meter
    for i in range(0, num_meter):
        meter_id = ECIGetMeteringIdByPosition(i)

        # 1-green; 2-flow; 3-delay; 4-ALINEA; 5-greenByLane
        meter_type = ECIGetTypeMeteringById(meter_id)

        if meter_type != 1:
            AKIPrintString('WARNING: Meter {0} is NOT green type'.format(meter_id))

        AKIPrintString('-- Meter {0} found with type {1}\n'.format(meter_id, meter_type))




# Here are the utility functions we defined for properly read control signals.
def readComConfigFile():

    global files

    f = open(com_file, 'r')

    for line in f:
        if line[0] == '%' or line[0] == '#':
            # comment line, skip
            continue
        else:
            line = line.strip()
            items = line.split(',')

            # There are four keys:
            # data; data_write_done; signal; signal_write_done; matlab_init_done
            files[items[0]] = items[1]

    f.close()




# This function
# 1. check if signal file exist
# 2. If no, do nothing.
# 3. If yes, read the signal (time, duration, flow (veh/hr)), and update signal array
# input:
#       signal file: each row, [signal_set_time, duration, flow (veh/hr)]
# output:
#       updated signal_array: flow in veh/hr
def readAllSignalFile(timeSta):

    global signal_array

    if exists(files['all_signal']):
        # read signal
        f_signal = open(files['all_signal'])

        # new signal
        new_signal = []

        for line in f_signal:
            if line[0] == '%' or line[0] == '#':
                continue
            else:
                line = line.strip()
                items = line.split(',')

                new_signal.append([float(items[0]), float(items[1]), float(items[2])])

        f_signal.close()

        new_signal = np.array(new_signal)

        # update signal array
        if signal_array is None:
            signal_array = new_signal
        else:
            # fine the first value that is greater than the first new set time
            first_index = np.argmax(signal_array[:,0] >= new_signal[0,0])

            # update global signal array
            signal_array = np.concatenate( (signal_array[0:first_index-1,:], new_signal))


        AKIPrintString( "-- Time {0} Read new signal\n".format(timeSta) )

    return





# The following functions are nto used
def AAPIPreRouteChoiceCalculation(time, timeSta):
    return 0

def AAPIEnterVehicle(idveh, idsection):
    return 0

def AAPIExitVehicle(idveh, idsection):
    return 0

def AAPIEnterPedestrian(idPedestrian, originCentroid):
    return 0

def AAPIExitPedestrian(idPedestrian, destinationCentroid):
    return 0

def AAPIEnterVehicleSection(idveh, idsection, atime):
    return 0

def AAPIExitVehicleSection(idveh, idsection, atime):
    return 0
