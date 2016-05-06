__author__ = 'Yanning Li'


# This scripts assumes the follows:
# - The network contains only one merge: us freeway and onramp merges to the ds freeway
# - The link id are known from AIMSUN GUI and set here in the manual config
# - Each detector across all lanes at its location
# - Only one meter exist which is at the exit of the onramp
# - The meter is preconfigured as a flow meter, we only need to update is flow as needed.
# - We will use the upstream most sensors on the upstream freeway and the onramp, the downstream most sensor on the
#   downstream freeway, and the sensors at the junction to run the MPC.
# - We fit a triangular fundamental diagram to the data collected in the AIMSUN simulation.
# - MATLAB and AIMSUN communicates in km/hr, veh/km, s; MATLAB internally use meters and seconds; AIMSUN
#   internally use metric or english.

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

    data_aggregate_interval = 5*60


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


    # wait for matlab to initialize
    while not exists(files['matlab_init_done']) and not exists(files['stop_control']):
        # sleep 0.1 second
        AKIPrintString( "Waiting for MATLAB to initialize..." )
        time_module.sleep(1)

    # matlab init done, then remove the flag file
    # os.remove(files['matlab_init_done'])

    # now start simulation
    AKIPrintString( "AAPIInit" )

    return 0


# This function is called before every time step
# 1. It should first check if a new signal file is written by MATLAB
# 2. If written, update is signal array.
# 3. Set the meter according to the signal array.
def AAPIManage(time, timeSta, timeTrans, acycle):

    if ECIGetTypeMeteringById(meter_id) != 2:
        AKIPrintString( "Warning: the meter type is not Flow." )
    else:

        # check the signal file; if exists,
        # the global signal array will be updated.
        readSignalFile(timeSta)

        # have signals to apply
        if signal_array is not None:

            # need to update signal
            apply_signal_index = np.where(signal_array[:,0]==time)[0]
            if len(apply_signal_index) != 0:
                # the signal updating time point
                new_flow = signal_array[apply_signal_index[0],2]

                # by default, we set the f_min = 0, f_max = 900
                ECIChangeParametersFlowMeteringById(meter_id, timeSta, 900, new_flow, 0)

                AKIPrintString( "---- Time {0}: applied signal flow {1} veh/hr.".format(time, new_flow) )

            # printSignalParameters(timeSta)

    return 0


# This function should read the detector data
# 1. Check if the detection interval is reached
# 2. If yes, then read the detector count data.
# 3. Write the data into files.
# The file format [report_time, flow (veh/hr), speed (kph)]
def AAPIPostManage(time, timeSta, timeTrans, acycle):
    # AKIPrintString( "AAPIPostManage" )

    if time%det_interval == 0:
        # interval reached
        extractDetectionData(time)

    if time%data_aggregate_interval == 0:

        if not exists(files['stop_control']):

            # write data interval reached
            writeData(timeSta)

            # if in warm up period, wait for MATLAB to read the data
            s_counter = 0
            while exists(files['data']) and not exists(files['stop_control'])\
                    and exists(files['matlab_init_done']) and not exists(files['start_control']) :
                time_module.sleep(0.1)
                s_counter += 1
                if s_counter == 20:
                    AKIPrintString('------ Warming up: wait for MATLAB to read previous data...')
                    s_counter = 0

            # if in control period, wait for MALTAB signal
            s_counter = 0
            while not exists(files['signal_write_done']) and not exists(files['stop_control'])\
                  and exists(files['start_control']) :
                time_module.sleep(0.1)
                s_counter += 1
                if s_counter == 20:
                    AKIPrintString('------ Waiting for MATLAB signal...')
                    s_counter = 0



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

        if meter_type != 2:
            AKIPrintString('WARNING: Meter {0} is NOT Flow type'.format(meter_id))

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
def readSignalFile(timeSta):

    global signal_array

    if exists(files['signal_write_done']):
        # read signal
        f_signal = open(files['signal'])

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


        # remove signal and flag file
        os.remove(files['signal'])
        os.remove(files['signal_write_done'])

        AKIPrintString( "-- Time {0} Read new signal\n".format(timeSta) )

    return



# print signal control parameters at time
def printSignalParameters(timeSta):

    # read the current control parameters of the meter
        fmax = doublep()
        flow = doublep()
        fmin = doublep()

        status = ECIGetParametersFlowMeteringById(meter_id, timeSta, fmax, flow, fmin)

        # make sure it is correctly red
        if status != 0:
            AKIPrintString( "Error: could not read the control parameters of Meter {0}.".format(meter_id) )
            return 0

        AKIPrintString('-- TIME: {0}'.format(timeSta))
        AKIPrintString('---- Control parameter of Meter {0}: Flow, {1} veh/hr in [{2}, {3}]'.format(
            meter_id, flow.value(), fmin.value(), fmax.value()
        ))
        del fmin
        del fmax
        del flow



# This function extract the detection data, aggregated in certain time interval
# input:
#       now_time: the current time in the simulator
# output:
#       data_to_write: the global list; each row is a string
#       each row: link_id, up/downstream, time, flow (veh/hr), speed (km/hr)
def extractDetectionData(now_time):

    num_det = AKIDetGetNumberDetectors()

    # for each detector
    for i in range(0, num_det):

        det_id = AKIDetGetIdDetector(i)

        # identify this detector
        det_info = AKIDetGetPropertiesDetectorById(det_id)
        sec_id = det_info.IdSection

        # avg flow in the last detection interval, veh/hr
        det_flow = AKIDetGetCounterAggregatedbyId(det_id,0)*3600/det_interval

        # REMARK: Make sure AIMSUN is using metric
        det_speed = AKIDetGetSpeedAggregatedbyId(det_id,0)

        if det_id == det_dict['us_in'] or det_id == det_dict['ds_in'] or det_id == det_dict['on_in']:
            bound = 'upstream'
            line = '{0},{1},{2},{3},{4}\n'.format(sec_id,bound,now_time,det_flow,det_speed)
            data_to_write.append(line)

        elif det_id == det_dict['us_out'] or det_id == det_dict['ds_out'] or det_id == det_dict['on_out']:
            bound = 'downstream'
            line = '{0},{1},{2},{3},{4}\n'.format(sec_id,bound,now_time,det_flow,det_speed)
            data_to_write.append(line)

    AKIPrintString("---- Time {0}: extracted new data.".format(now_time) )



# This function writes the aggregated data every data_aggregate_interval
# input:
#       timeSta: the current time in the stationary period in simulator
# output:
#       data_to_write: the global list; each row is a string
#       each row: link_id, up/downstream, time, flow (veh/hr), speed (km/hr)
def writeData(timeSta):

    global data_to_write

    if len(data_to_write) != 0:

        # # Before writing data, check if the data has been read by MATLAB. (If read, MATLAB would delete it),
        # # so we do not overwrite the data
        # p_counter = 0
        # while exists(files['data']) and not exists(files['stop_control']):
        #
        #     if p_counter == 10:
        #         AKIPrintString('------ Waiting for MATLAB to read previous data...')
        #         p_counter = 0
        #
        #     time_module.sleep(0.1)
        #     p_counter += 1


        # write the data into file
        f = open(files['data'],'wb')
        for line in data_to_write:
            f.write(line)
        f.close()

        # write the data flag file
        f = open(files['data_write_done'],'wb')
        f.write('data write done\n')
        f.close()

        # empty the global variable data_to_write
        del data_to_write[:]

        AKIPrintString("-- Time {0}: Wrote new data.".format(timeSta) )
















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
