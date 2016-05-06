__author__ = 'Carlos and Yanning'

import sys
import os
from os.path import exists
import time
from datetime import datetime
import csv
import random

SITEPACKAGES = "C:\\Python27\\Lib\\site-packages"
sys.path.append(SITEPACKAGES)
import numpy as np
import matplotlib.pyplot as plt

from PyANGBasic import *
from PyANGKernel import *
from PyANGConsole import *
from PyANGAimsun import *

# Files used in this file
# parsFilePath = 'c:/tmp/pars.txt'
# flowsFilePath = 'c:/tmp/flows.txt'
# turnsFilePath = 'c:/tmp/turns.txt'

#=====================================================
# The following are just to remove the error messages
'''
def QString():
    pass
def QVariant():
    pass
def GKScheduleDemandItem():
    pass
GKTimeDuration = None
GKVehicle = None
GKSystem = None
GKSection = None
GKExperiment = None
GKVehicleReactionTimes = None
def GAimsunSimulator():
    pass
GKSimulationTask = None
GKReplication = None
GGetramModule = None
GKTimeSerieIndex = None
GK = None
GKColumnIds = None
Qt = None
QTime = None
'''
#=====================================================



# The following variables need to be shared by functions (init, and update)
cmd_logger = None
opt_logger = None
opt_logger_file = None
# TODO: Windows is not properly configured yet. Hence could not plot figure.
# Figure handle is to plot the optimization steps (fig, ax)
fig_handle = None

# Define the columns for the files to be read
# pars.txt
INDEX_PARS_STATE = 0
INDEX_PARS_VEH = 1
INDEX_PARS_FROMTIME = 2
INDEX_PARS_DUR = 3
# flows.txt
INDEX_FLOWS_STATE = 0
INDEX_FLOWS_SECID = 1
INDEX_FLOWS_INFLOW = 2
# turns.txt
INDEX_TURNS_STATE = 0
INDEX_TURNS_FROM = 1
INDEX_TURNS_TO = 2
INDEX_TURNS_PERC = 3

KMH2MPH = 0.62137

_debug = False
_show_detector_data = False


#=================================================================================================
# The following functions are used to read the flows.txt, turns.txt, and paras.txt files and create states

# Read the demand parameters from file
# input:
#       parsFilePath: the file name with full path. 
#           Each row of file: state name, car type and name, state start time, state duration
#           e.g: cong_state1_car,53 Car,15:30:00,00:05:00
# output:
#       pars_dict: [state name] = [vehicle type, state start time, duration], all strings
def readParsFile(parsFilePath):
    pars_dict = {}
    pars_file = open(parsFilePath, 'r')
    while True:
        line = pars_file.readline()
        if not (bool(line)):
            break
        line = line.strip()
        items = line.split(',')
        tmp_key = items[INDEX_PARS_STATE].strip()
        if tmp_key not in pars_dict.keys():
            pars_dict[tmp_key] = list()
        pars_dict[tmp_key].append(
            (items[INDEX_PARS_VEH].strip(), items[INDEX_PARS_FROMTIME].strip(), items[INDEX_PARS_DUR].strip()))
    pars_file.close()
    return pars_dict

# Read the demand flows from file
# input:
#       flowsFilePath: the file name with full path. 
#           Each row of file: state name, section id, inflow (veh/hr)
#           e.g: cong_state3_car,30109,50.0
# output:
#       flows_dict: [state name] = [section_id, inflow], all strings
def readFlowsFile(flowsFilePath):
    flows_dict = {}
    flows_file = open(flowsFilePath, 'r')
    while True:
        line = flows_file.readline()
        if not (bool(line)):
            break
        line = line.strip()
        items = line.split(',')
        tmp_key = items[INDEX_FLOWS_STATE].strip()
        if tmp_key not in flows_dict.keys():
            flows_dict[tmp_key] = list()
        flows_dict[tmp_key].append((items[INDEX_FLOWS_SECID].strip(), items[INDEX_FLOWS_INFLOW].strip()))
    flows_file.close()
    return flows_dict


# Read the demand turns from file
# input:
#       turnsFilePath: the file name with full path. 
#           Each row of file: state name, from section_id, to section_id, percent
#           e.g: cong_state1_car,21217,23587,2.0
# output:
#       turns_dict: [state name] = [from_sec_id, to_sec_id, percent]
def readTurnsFile(turnsFilePath):
    turns_dict = {}
    turns_file = open(turnsFilePath, 'r')
    while True:
        line = turns_file.readline()
        if not (bool(line)):
            break
        line = line.strip()
        items = line.split(',')
        tmp_key = items[INDEX_TURNS_STATE].strip()
        if tmp_key not in turns_dict.keys():
            turns_dict[tmp_key] = list()
        turns_dict[tmp_key].append(
            (items[INDEX_TURNS_FROM].strip(), items[INDEX_TURNS_TO].strip(), items[INDEX_TURNS_PERC].strip()))
    turns_file.close()
    # print '\n\nturns_dict:{0}\n\n'.format(turns_dict)
    return turns_dict


# This function creates state from read dicts
# It also has a main_entrance_id, which can be discounted.
# input:
#       model: the AIMSUN model object
#       state_name: the statename to be created
#       _dicts: the dicts created by files
#       main_entrance_id: the inflow to be discounted. Set as None if do not want any discount
#       main_entrance_discount_ratio: set as 1 if not discount.
def createState(model, state_name, pars_dict, flows_dict, turns_dict, main_entrance_id, main_entrance_discount_ratio):

    # create new state and set parameters
    state = GKSystem.getSystem().newObject("GKTrafficState", model)
    state.setName(state_name)
    it = QTime.fromString((pars_dict[state_name])[0][1], Qt.ISODate)
    duration = (pars_dict[state_name])[0][2]
    state.setInterval(it, GKTimeDuration.fromString(duration))

    if _debug is True:
        pass
        # print_cmd('AIMSUNFUNCTIONS: state from {0} duration {1}'.format(it.toString(), duration))

    # set the vehicle for this state
    vehicleString = str((pars_dict[state_name])[0][0]).split(' ')

    # a quick fix for the FD calibration
    # vehicleString = str((pars_dict[state_name])[0]).split(' ')

    vehId = int(vehicleString[0])   # make sure this is correct
    vehName = vehicleString[1]
    vehicle = state.getModel().getCatalog().find(vehId)
    if vehicle is None:
        # is wont work since the name is not unique, UserClass has object called car
        vehicle = state.getModel().getCatalog().findByName(vehName)
    state.setVehicle(vehicle)

    # set the inflow of the state
    for entrance in range(0, len(flows_dict[state_name])):
        # print (flows_dict[state_name])[entrance]


        fromSection = findSection(model, (flows_dict[state_name])[entrance][0])

        # discount the main entrance flow
        if fromSection.getId() == main_entrance_id:
            state.setEntranceFlow(fromSection, None,
                                  float((flows_dict[state_name])[entrance][1])*main_entrance_discount_ratio)
        else:
            state.setEntranceFlow(fromSection, None,
                                  float((flows_dict[state_name])[entrance][1]))

    # set the turn percentage of the state
    if state_name in turns_dict.keys():
        for turn in range(0, len(turns_dict[state_name])):
            # print_cmd('For state {0}, has turns: {1}'.format(state_name, turns_dict[state_name]))
            fromSection = findSection(model, (turns_dict[state_name])[turn][0])
            toSection = findSection(model, (turns_dict[state_name])[turn][1])
            state.setTurningPercentage(fromSection, toSection, None, float((turns_dict[state_name])[turn][2]))
    else:
        print_cmd('No Turn information for state {0}'.format(state_name))

    # for testing the aimsun automation
    if _debug:
        # print_cmd('AIMSUNFUNCTION: state turn 330->(340, 341): ({0},{1})'.format(
        #     state.getTurningPercentage(model.getCatalog().find(int(330)),
        #                                model.getCatalog().find(int(340)), None),
        #     state.getTurningPercentage(model.getCatalog().find(int(330)),
        #                                model.getCatalog().find(int(341)), None)))
        pass

    # append the state to the state folder
    folder = getStateFolder(model)
    folder.append(state)

    return state


# return the GKSection object of the entry id.
# input: entry a string '330' which gives the id
# ouput: GKSection
def findSection(model, entry):
    section = model.getCatalog().find(int(entry))
    if section.isA(QString("GKSection")) is False:
        section = None
    return section


# Returns (and creates if needed) the folder for the traffic state
def getStateFolder(model):
    folderName = "GKModel::trafficStates"
    folder = model.getCreateRootFolder().findFolder(folderName)
    if folder is None:
        folder = GKSystem.getSystem().createFolder(
            model.getCreateRootFolder(), folderName)

    # print_cmd('getStateFolder: type: {0}'.format(type(folder))
    # print_cmd('getStateFolder: name: {0}'.format(folder.getName())
    return folder


# Returns (and creates if needed) the folder for the traffic demand
def getDemandFolder(model):
    folderName = "GKModel::trafficDemand"
    folder = model.getCreateRootFolder().findFolder(folderName)
    if folder is None:
        folder = GKSystem.getSystem().createFolder(
            model.getCreateRootFolder(), folderName)
    return folder


# Returns (and creates if needed) the folder for the traffic demand
def getScenariosFolder(model):
    folderName = "GKModel::top::scenarios"
    folder = model.getCreateRootFolder().findFolder(folderName)
    if folder is None:
        folder = GKSystem.getSystem().createFolder(
            model.getCreateRootFolder(), folderName)
    return folder


# Before the states can be added to the demand, each state must be added to a schedule item before being added to demand.
# NOTE: states contains the flow, turn percentage, and vehicle type. However, to add the state to the demand, we need to
# assign the state to a scheduleDemandItem. The simulator uses the fromTime and duration of the scheduleDemandItem, NOT
# the fromTime and duration set in the state!
# input: 
#       state: the state object
# output: 
#       GKScheduleDemandItem object: associated with the state
def createScheduleItem(state):
    schedule = GKScheduleDemandItem()
    schedule.setTrafficDemandItem(state)

    if _debug:
        # print_cmd('schedule.state.duration {0}'.format(schedule.getTrafficDemandItem().getDuration().toString()))
        pass

    hr = state.getFrom().hour()
    minute = state.getFrom().minute()
    sec = state.getFrom().second()
    schedule.setFrom(3600 * hr + 60 * minute + sec)

    if _debug:
        # print_cmd('state.getDuration().toSeconds(): {0}'.format(state.getDuration().toSeconds()[0]))
        # print_cmd(type(state.getDuration().toSeconds()[0]))
        pass

    schedule.setDuration(state.getDuration().toSeconds()[0])

    return schedule


# reads the demand files:
# paras.txt; flow.txt; turn.txt
# input: model, model class of AIMSUN as global variables; and file paths
# output: return GKDemand class object; this demand has been added to the model;
#         the return is not necessary, only for debugging
def read_demand_from_file(model, demand_name, pars_file, flows_file, turns_file, main_entrance_id, main_entrance_discount_ratio):

    # Set up a Traffic Demand and add it to the model
    if demand_name is not None:
        demand = model.getCatalog().findByName(QString(demand_name))
    else:
        demand = None

    if demand is None or demand.isA(QString("GKTrafficDemand")) is False:
        # create new demand
        demand = GKSystem.getSystem().newObject("GKTrafficDemand", model, -1, True)
        demand.setName(QString(demand_name))
        print_cmd('Demand {0} not found. Creating new one...'.format(demand_name))

        # create under demand folder
        folder = getDemandFolder(model)
        folder.append(demand)

    # clear schedule in demand to prevent overlapping
    demand.removeSchedule()

    # Read pars.txt into a dictionary
    pars_dict = readParsFile(pars_file)
    flows_dict = readFlowsFile(flows_file)
    turns_dict = readTurnsFile(turns_file)

    for state_name in pars_dict.keys():
        state = createState(model, state_name, pars_dict, flows_dict, turns_dict, main_entrance_id, main_entrance_discount_ratio)

        if _debug:
            # print_cmd('state from {0}, duration:{1}\n'.format(state.getFrom().toString(), state.getDuration().toString()))
            #
            # print_cmd('state entrance flow {0}'.format(state.getEntranceFlow(model.getCatalog().find(int(330)), None)))
            #
            # print_cmd('state turn 330->(340, 341): ({0},{1})'.format(state.getTurningPercentage(model.getCatalog().find(int(330)),
            #                                                                             model.getCatalog().find(int(340)), None),
            #                                                  state.getTurningPercentage(model.getCatalog().find(int(330)),
            #                                                                             model.getCatalog().find(int(341)), None)))
            pass

        schedule = createScheduleItem(state)

        demand.addToSchedule(schedule)

    # append the demand to the demands folder
    # folder = getDemandFolder(model)
    # folder.append(demand)

    return demand


# load demand and state from .ang
# Only for debugging purpose
def load_demand_from_ang(model, demand_name):
    # find the demand from the model
    demand = model.getCatalog().findByName( QString(demand_name) )
    if demand is None or not demand.isA(QString("GKTrafficDemand")):
        print_cmd('Error: no traffic demand named {0}\n'.format(demand_name))
        return None

    demand.removeSchedule()

    # find state
    stateType = model.getType("GKTrafficState")
    for state in model.getCatalog().getObjectsByType(stateType).itervalues():
        print_cmd('loading state {0}'.format(state.getName()))

        schedule = createScheduleItem(state)
        demand.addToSchedule(schedule)

    return demand



#=================================================================================================
# This code is to read validation data from csv file
# data format should be as follows: Make SURE there is Header
# Date/Time (UTC-06:00) Central Time (US & Canada),RadarSpeed (MPH),RadarVehiclesCount
# 4/26/2015 15:55,,
# 4/26/2015 16:00,72.23863137,51
# 4/26/2015 16:05,71.23257753,91
# input: read data between [start_time, end_time]; if None, then read all rows
#        name_strs is a list of detector names
#        valid_file_path: if '', then current folder same as the script
# output: a dict: validation_data[detector_name] = np.array([[speed],[count]])
# Note: to read data for a simulation from 03:00 to 04:00 (AIMSUN aggragate data at 03:05 + 00:05...)
# need to read 12 lines with start_time and end_time 03:00~to 04:00 (Carlos shifted the time 5 min earlier)
def read_validation_data(start_time_str, end_time_str, name_strs, valid_file_path):

    dict_valid = {}

    # if *_time_str is None, then assume the entire file is the validation data
    if start_time_str is None or end_time_str is None:

        for name in name_strs:

            # first list speed, second list count
            dict_valid[name] = [[], []]

            file_name = valid_file_path + name + '.csv'
            f_handle = open(file_name, 'r')

            data_set = csv.reader(f_handle, delimiter=',')
            # skip header
            next(data_set, None)

            for row in data_set:
                dict_valid[name][0].append(float(row[1]))
                dict_valid[name][1].append(float(row[2]))

            f_handle.close()

    # Otherwise on read data in a specific time period
    else:
        dt_format = '%m/%d/%Y %H:%M'
        t_start = datetime.strptime(start_time_str, dt_format)
        t_end = datetime.strptime(end_time_str, dt_format)

        for name in name_strs:
            # first list speed, second list count
            dict_valid[name] = [[], []]

            file_name = valid_file_path + name + '.csv'
            f_handle = open(file_name, 'r')

            data_set = csv.reader(f_handle, delimiter=',')
            # skip header
            next(data_set, None)

            for row in data_set:

                cur_dt = datetime.strptime(row[0], dt_format)

                if t_start <= cur_dt < t_end:
                    # no need to save the time
                    dict_valid[name][0].append(float(row[1]))
                    dict_valid[name][1].append(float(row[2]))

            f_handle.close()

        # print 'Loaded validation data from time {0} to time {1}, [[speed/mph],[count/5min]]: {2}'.format(start_time_str, end_time_str, dict_valid[name])

    return dict_valid
    # check if correctly read
    # for key in dict_valid:
        # print_cmd('key {0}: {1}'.format(key, dict_valid[key])




#=================================================================================================
# The following functions are creating scenario, experiment, model, replication, simulate and get data...

# set up the scenario
# the scenario could be free flow or congested using different demand
# there are also some parameters we would like to set (but not used later as the decision variable)
def setup_scenario(model, scenario_name, demand):

    print_cmd('\nSetting up scenario...')

    scenario = model.getCatalog().findByName( QString(scenario_name) )
    if scenario is None or not scenario.isA(QString("GKScenario")):
        scenario = GKSystem.getSystem().newObject("GKScenario", model, -1, True)
        scenario.setName(QString(scenario_name))
        print_cmd('Error: no traffic scenario named {0}. Creating new one...\n'.format(scenario_name))

    scenario.setDemand(demand)

    # append the state to the state folder
    # folder = getScenariosFolder(model)
    # folder.append(scenario)

    # set parameters here
    # parameters are set in the ScenarioInput data class
    paras = scenario.getInputData()
    # set the detection and statistical intervals as 5 min
    paras.setDetectionInterval(GKTimeDuration.fromString(QString("00:01:00")))
    paras.setStatisticalInterval(GKTimeDuration.fromString(QString("00:01:00")))

    # There are also five environment parameters:
    # day of week; season; weather; event; methodology
    # TODO: the following parameters should have been set in GUI
    # scenario.setDataValueByID( GKGenericScenario.weekdayAtt, QVariant('Monday') )
    # scenario.setDataValueByID( GKGenericScenario.seasonAtt, QVariant( 'Summer' ) )
    # scenario.setDataValueByID( GKGenericScenario.weatherAtt, QVariant( 'Sunny' ) )
    # scenario.setDataValueByID( GKGenericScenario.eventAtt, QVariant( 'Fair' ) )
    # print_cmd('senario: {0}'.format(scenario.getDataValueByID( GKGenericScenario.seasonAtt )[0].toString() )
    # scenario.setValueForVariable(QString('seasonAtt'), QString('Spring'))
    # print_cmd('setup_scenario: '.format(scenario.getValueForVariable( QString('seasonAtt') ) )

    # var_dict = scenario.getVariables()
    # for key in var_dict:
    #     print_cmd('key {0}: {1}'.format(key, var_dict[key])
    # print_cmd('setup_scenario: get_variables: {0}'.format()

    return scenario


# set up the experiment for I80_EB
# This function reads and setup the preset parameters from preset_paras.txt
def setup_experiment_I80_EB(model, experiment_name, scenario, preset_paras_file):
    print_cmd('\nSetting up experiment...\n')

    experiment = model.getCatalog().findByName( QString(experiment_name) )
    if experiment is None or not experiment.isA(QString("GKExperiment")):
        experiment = GKSystem.getSystem().newObject("GKExperiment", model, -1, True)
        print_cmd('ERROR: No traffic experiment named {0}. Creating a new one...\n'.format(experiment_name))

        # attach the new experiment to folder
        folder = getScenariosFolder(model)
        folder.append(experiment)

    experiment.setScenario(scenario)


    #==================================================================================
    # TODO: some parameters still needs to be set in the .ang file which is easier, including:
    # TODO: the speed limit, simulation step, strategy plan...

    # read preset paras, which contains the initial values
    # Here the preset parameters is limited to a few parameters
    preset_paras = read_preset_paras(preset_paras_file)

    #==================================================================================
    # vehicle class parameters
    car_type = model.getCatalog().find( 53 )
    if not car_type.isA(QString("GKVehicle")):
        print_cmd('Error: Car type is not correct: demand may not correctly loaded\n')
        return None

    # unique truck id is 56
    truck_type = model.getCatalog().find(56)
    if not truck_type.isA(QString("GKVehicle")):
        print_cmd('Error: Truck type is not correct: demand may not correctly loaded\n')
        return None

    print_cmd('Setting preset parameters:')

    for key in preset_paras.keys():

        # key format: car_maxSpeed
        para_name = key.split('_')
        para_value = preset_paras[key]
        if para_name[0] == 'car':
            # car paras
            if para_name[1] == 'maxSpeed':
                car_type.setDataValueByID( GKVehicle.maxSpeedMean, QVariant(para_value[0]))
                car_type.setDataValueByID( GKVehicle.maxSpeedDev, QVariant(para_value[1]))
                car_type.setDataValueByID( GKVehicle.maxSpeedMin, QVariant(para_value[2]))
                car_type.setDataValueByID( GKVehicle.maxSpeedMax, QVariant(para_value[3]))
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'speedAcceptance':
                car_type.setDataValueByID( GKVehicle.speedAcceptanceMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.speedAcceptanceDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.speedAcceptanceMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.speedAcceptanceMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'maxAccel':
                car_type.setDataValueByID( GKVehicle.maxAccelMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.maxAccelDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.maxAccelMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.maxAccelMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'reactionTime':
                # [reaction_time, reaction_stop, reaction_light, reaction_prob]
                car_react = GKVehicleReactionTimes(para_value[0], para_value[1],
                                                   para_value[2], para_value[3])

                car_type.setVariableReactionTimes( [car_react] )
                experiment.setVariableReactionTimesMicro(car_type, [car_react])
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'minDist':
                car_type.setDataValueByID( GKVehicle.minDistMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.minDistDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.minDistMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.minDistMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'sensitivityFactor':
                car_type.setDataValueByID( GKVehicle.sensitivityFactorMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.sensitivityFactorDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.sensitivityFactorMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.sensitivityFactorMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            else:
                print_cmd('\n---- ERROR: Could not recognize preset parameter entry {0}: {1}\n'.format(key, para_value))


        elif key.split('_')[0] == 'truck':
            # truck paras
            if para_name[1] == 'maxSpeed':
                truck_type.setDataValueByID( GKVehicle.maxSpeedMean, QVariant(para_value[0]))
                truck_type.setDataValueByID( GKVehicle.maxSpeedDev, QVariant(para_value[1]))
                truck_type.setDataValueByID( GKVehicle.maxSpeedMin, QVariant(para_value[2]))
                truck_type.setDataValueByID( GKVehicle.maxSpeedMax, QVariant(para_value[3]))
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'speedAcceptance':
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'maxAccel':
                truck_type.setDataValueByID( GKVehicle.maxAccelMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.maxAccelDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.maxAccelMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.maxAccelMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'reactionTime':
                # [reaction_time, reaction_stop, reaction_light, reaction_prob]
                truck_react = GKVehicleReactionTimes(para_value[0], para_value[1],
                                                   para_value[2], para_value[3])

                truck_type.setVariableReactionTimes( [truck_react] )
                experiment.setVariableReactionTimesMicro(truck_type, [truck_react])
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'minDist':
                truck_type.setDataValueByID( GKVehicle.minDistMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.minDistDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.minDistMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.minDistMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'sensitivityFactor':
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            else:
                print_cmd('\n---- ERROR: Could not recognize preset parameter entry {0}: {1}\n'.format(key, para_value))

        else:
            print_cmd('\nERROR: Could not recognize preset parameter entry {0}\n '.format(key))

    return experiment


# setup num_rep replications as well as the average_result
# create num_rep replications
# Single Run in  main
# input: seed_list is the seed of replications. We may want seeds to be static
def setup_replication(model, experiment, num_rep, seed_list):

    print_cmd('\nSetting up replications...')

    if experiment != None and experiment.isA( "GKExperiment" ) \
        and experiment.getSimulatorEngine() == GKExperiment.eMicro:

        # add 10 replications here
        replication_list = experiment.getReplications()
        if len(replication_list) == 0:
            # create 10 replications
            for i in range(0, num_rep):
                replication = GKSystem.getSystem().newObject("GKReplication", model, -1, True)
                replication.setExperiment(experiment)
                replication_list.append(replication)

                if seed_list is not None:
                    replication.setRandomSeed(seed_list[i])
                print_cmd('---- Created replication {0} with seed {1}'.format(replication.getId(),
                                                                 replication.getRandomSeed()))
        else:
            # show replcations:
            print_cmd('---- Reloading {0} replications: {1} \n'.format( len(replication_list),
                                                           [ replication.getId() for replication in replication_list  ]))

        # create the average experiment result
        avg_result = GKSystem.getSystem().newObject( "GKExperimentResult", model )
        avg_result.setName('average_result')
        print_cmd('Created new average replication: {0}'.format(avg_result.getName()))
        # print_cmd('Total number of replications is: {0}',format(len(experiment.getReplications()))

        # set the experiment of this result object
        avg_result.setExperiment( experiment )
        # add replcations to the average
        for replication in replication_list:
            avg_result.addReplication( replication )
            print_cmd('---- Added replication {0} to {1}'.format(replication.getId(), avg_result.getName()))

        # compute the average; add to the experiment.
        experiment.addReplication( avg_result )

        return avg_result


# set random seed
def set_random_seed(avg_result):

    replication_list = avg_result.getReplications()

    print_cmd('\nResetting seeds for replications:')

    i = 0
    for replication in replication_list:
        replication.setRandomSeed(random.randint(0, 10000))
        print_cmd('----Reset replication {0} with seed {1}'.format(replication.getId(),
                                                                 replication.getRandomSeed()))
        # replication.setRandomSeed(int(seed_list[i]))
        i += 1


# Create a simulator which could be used for all simulations
# Single Run in Main
def create_simulator(model):
    simulator = GAimsunSimulator()
    simulator.setModel(model)

    return simulator


# simulate the experiment after the experiment parameters are set
# input: simulator
#      : avg_result returned from setup_replication
def simulate_experiment(simulator, avg_result):

    print_cmd('\nReset replications...')

    # first reset replications
    avg_result.resetReplications()

    # add replications to simulator
    replication_list = avg_result.getReplications()

    for replication in replication_list:
        simulation_task = GKSimulationTask(replication, GKReplication.eBatch, "", "", True)  # Other approach
        simulator.addSimulationTask(simulation_task)
        # print_cmd('Added replication {0} to simulator with status {1}. '.format(replication.getId(),
        #                                                                        replication.getSimulationStatus())
        # print_cmd('pending {0}; done {1}; discarded {2}; loaded {3}'.format(GKGenericExperiment.ePending,
        #                                                                    GKGenericExperiment.eDone,
        #                                                                    GKGenericExperiment.eDiscarded,
        #                                                                    GKGenericExperiment.eLoaded)


    # simulate model
    if not simulator.isBusy():
        print_cmd('Simulating...\n')
        sim_status = simulator.simulate()
    else:
        print_cmd('Simulator is busy\n')

    # make sure correctly simulated
    if sim_status is True:
        print_cmd('Simulation finished\n')
    else:
        print_cmd('ERROR: Simulation failed\n')

    # simulator.postSimulate()


# set new parameters to the models
def set_new_paras(model, experiment, paras):

    # Unique car id is 53
    car_type = model.getCatalog().find(53)
    if not car_type.isA(QString("GKVehicle")):
        print_cmd('Error: Car type is not correct: demand may not correctly loaded\n')
        print_cmd(type(car_type))
        # debug_catalog(model)
        return None

    # unique truck id is 56
    truck_type = model.getCatalog().find(56)
    if not truck_type.isA(QString("GKVehicle")):
        print_cmd('Error: Truck type is not correct: demand may not correctly loaded\n')
        return None

    # note the parameter should be in the same format as the preset parameters
    # even if we are only calibrating the mean
    for key in paras.keys():

        para_name = key.split('_')
        para_value = paras[key]

        if para_name[0] == 'car':
            # car paras
            # calibrate speedAcceptanceMean and speedAcceptanceDev for freeflow
            if para_name[1] == 'speedAcceptance':
                car_type.setDataValueByID( GKVehicle.speedAcceptanceMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.speedAcceptanceDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.speedAcceptanceMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.speedAcceptanceMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            # calibrate maxAccelMean for congflow
            elif para_name[1] == 'maxAccel':
                car_type.setDataValueByID( GKVehicle.maxAccelMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.maxAccelDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.maxAccelMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.maxAccelMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            # calibrate reaction_time for congflow
            elif para_name[1] == 'reactionTime':
                # [reaction_time, reaction_stop, reaction_light, reaction_prob]
                car_react = GKVehicleReactionTimes(para_value[0], para_value[1],
                                                   para_value[2], para_value[3])

                car_type.setVariableReactionTimes( [car_react] )
                experiment.setVariableReactionTimesMicro(car_type, [car_react])
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            # calibrate minDistMean for congflow
            elif para_name[1] == 'minDist':
                car_type.setDataValueByID( GKVehicle.minDistMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.minDistDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.minDistMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.minDistMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            # calibrate sensitivityFactorMean = min = max for congflow
            elif para_name[1] == 'sensitivityFactor':
                car_type.setDataValueByID( GKVehicle.sensitivityFactorMean, QVariant(para_value[0]) )
                car_type.setDataValueByID( GKVehicle.sensitivityFactorDev, QVariant(para_value[1]) )
                car_type.setDataValueByID( GKVehicle.sensitivityFactorMin, QVariant(para_value[2]) )
                car_type.setDataValueByID( GKVehicle.sensitivityFactorMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            # This should already been preset as unrealistically high such that speed can be fully
            # controlled by the speed acceptance
            elif para_name[1] == 'maxSpeed':
                car_type.setDataValueByID( GKVehicle.maxSpeedMean, QVariant(para_value[0]))
                car_type.setDataValueByID( GKVehicle.maxSpeedDev, QVariant(para_value[1]))
                car_type.setDataValueByID( GKVehicle.maxSpeedMin, QVariant(para_value[2]))
                car_type.setDataValueByID( GKVehicle.maxSpeedMax, QVariant(para_value[3]))
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            else:
                print_cmd('\n---- ERROR: Could not recognize preset parameter entry {0}: {1}\n'.format(key, para_value))


        elif para_name[0] == 'truck':
            # truck paras
            if para_name[1] == 'speedAcceptance':
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.speedAcceptanceMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'maxAccel':
                truck_type.setDataValueByID( GKVehicle.maxAccelMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.maxAccelDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.maxAccelMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.maxAccelMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'reactionTime':
                # [reaction_time, reaction_stop, reaction_light, reaction_prob]
                truck_react = GKVehicleReactionTimes(para_value[0], para_value[1],
                                                   para_value[2], para_value[3])

                truck_type.setVariableReactionTimes( [truck_react] )
                experiment.setVariableReactionTimesMicro(truck_type, [truck_react])
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'minDist':
                truck_type.setDataValueByID( GKVehicle.minDistMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.minDistDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.minDistMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.minDistMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'sensitivityFactor':
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorMean, QVariant(para_value[0]) )
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorDev, QVariant(para_value[1]) )
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorMin, QVariant(para_value[2]) )
                truck_type.setDataValueByID( GKVehicle.sensitivityFactorMax, QVariant(para_value[3]) )
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            elif para_name[1] == 'maxSpeed':
                truck_type.setDataValueByID( GKVehicle.maxSpeedMean, QVariant(para_value[0]))
                truck_type.setDataValueByID( GKVehicle.maxSpeedDev, QVariant(para_value[1]))
                truck_type.setDataValueByID( GKVehicle.maxSpeedMin, QVariant(para_value[2]))
                truck_type.setDataValueByID( GKVehicle.maxSpeedMax, QVariant(para_value[3]))
                print_cmd('---- set {0}: {1}'.format(key, para_value))

            else:
                print_cmd('\n---- ERROR: Could not recognize preset parameter entry {0}: {1}\n'.format(key, para_value))

        elif para_name[0] == 'main':
                # TODO: if we want to add ramp flows as decision variables,
                # TODO: we need to write the corresponding ramp adjust functions
                adjust_I80_demand(model, [key, para_value])

        else:
            print_cmd('\nERROR: Could not recognize preset parameter entry {0}\n '.format(key))


# adjust the demand data
# para is one entry of paras: [key, flow]
# TODO: this is a quick hack to make it work.
def adjust_I80_demand(model, para):

    para_name = para[0]

    print_cmd('Setting flows: {0}'.format(para))

    # convert to 30 states
    flow_ratio = []
    for i in range(0, len(para[1])):
        flow_ratio.append(para[1][i])
        flow_ratio.append(para[1][i])

    if len(flow_ratio) != 30:
        print_cmd('\nERROR: adjust_I80_demand now only support 15 flow variables.\n')

    state_type = model.getType("GKTrafficState")

    if para_name.split('_')[0] == 'main':

        main_entrance = findSection(model, 21216)
        truck_ratio = 0.27

        main_count_from_EB3 = [159, 181, 183, 176, 181,
                           180, 175, 206, 187, 195,
                           174, 176, 185, 173, 160,
                           167, 159, 123, 174, 196,
                           168, 132, 149, 173, 143,
                           147, 157, 128, 130, 161]
        # change the input ratio (0.9~1.1) to flow veh/hr
        main_flow = []
        for i in range(0, len(flow_ratio)):
            main_flow.append( flow_ratio[i]*12.0*main_count_from_EB3[i] )

        for state in model.getCatalog().getObjectsByType(state_type).itervalues():

            # set the car state and the order must be correct
            if state.getVehicle().getId() == 53:
                # find which state is this
                state_name = str(state.getName())
                # in this format: cong_state#_car/truck
                middle_name = state_name.split('_')[1]
                # print_cmd('middle_name: {0}'.format(middle_name))
                state_id = int( middle_name[5:] )
                # print_cmd('Set the car state {0} with with new flow {1}'.format(state_name, main_flow[state_id-1]*(1-truck_ratio)))
                state.setEntranceFlow(main_entrance, None, float(main_flow[state_id-1]*(1-truck_ratio)))

            # set the truck state and the order must be correct
            elif state.getVehicle().getId() == 56:
                # find which state is this
                state_name = str(state.getName())
                # in this format: cong_state#_car/truck
                middle_name = state_name.split('_')[1]
                # print_cmd('middle_name: {0}'.format(middle_name))
                state_id = int( middle_name[5:] )
                # print_cmd('Set the truck state {0} with new flow {1}'.format(state_name, main_flow[state_id-1]*truck_ratio))
                state.setEntranceFlow(main_entrance, None, float(main_flow[state_id-1]*truck_ratio))



    elif para_name.split('_')[0] == 'onramp':

        onramp_4 = findSection(model, 1357)








# get the averaged data from the simulation
# The output should be the averaged count and speed data (across 10 replications) of each detector.
# This should be the same format as the loaded validation data.
# output: a dict, avg_data[detector_name] = [[speed],[count]]
def get_averaged_detector_data(model, avg_result, plugin):

    # compute the result
    calculate_status = plugin.calculateResult( avg_result )

    if calculate_status == GGetramModule.eOKCalculateResult:
        print_cmd('Retrieving average data finished.')
    elif calculate_status == GGetramModule.eFailCalculateResult:
        # at 5th iteration failed.
        print_cmd('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        print_cmd('$$$$$$$$$$$ ERROR: Retrieving average data failed.')
        print_cmd('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    elif calculate_status == GGetramModule.eSimRequiredCalculateResult:
        print_cmd('$$$$$$$$$$$ ERROR: Retrieving average data failed. Simulation required.')

    # Not sure if this line is needed
    plugin.readResult(avg_result)

    # read the detector data out
    avg_data = read_detector_data(model, [avg_result])

    return avg_data


# Read the detector data
# This function tests extracting data from the detector from a specific replication
# input: data_origin, a list of replications: [a replication or the average_result(which is a subtype of GKReplication)]
# output: a dict, avg_data[detector_name] = [[speed mph],[count/detection interval]]
def read_detector_data(model, data_origin):

    avg_data = {}

    det_type = model.getType("GKDetector")
    # read data for each replication and then the average
    for replication in data_origin:

        print_cmd('\nReading Replication data: {0}'.format(replication.getName()))

        # get the column id
        speedColumn = det_type.getColumn(GK.BuildContents(GKColumnIds.eSpeed, replication, None))
        countColumn = det_type.getColumn(GK.BuildContents(GKColumnIds.eCount, replication, None))

        # read each detector
        for det in model.getCatalog().getObjectsByType(det_type).itervalues():

            det_name = str(det.getName())
            # print_cmd('----Reading Detector {0}...'.format(det_name))

            # add to dict
            # flow, speed
            avg_data[det_name] = [[],[]]

            speedData = det.getDataValueTS(speedColumn)
            countData = det.getDataValueTS(countColumn)

            if countData.size() == 0 or speedData.size() == 0 or countData.size() != speedData.size():
                print_cmd('ERROR: Detector {0} has no data available'.format(det_name))
            else:
                # print_cmd('----size of data is: {0}'.format(countData.size()))
                # TODO: DONE: the speed data returned from AIMSUN is in km/h; 1 km/h = 0.62137 mph
                for interval in range(countData.size()):
                    avg_data[det_name][0].append(speedData.getValue(GKTimeSerieIndex(interval))[0]*KMH2MPH)
                    avg_data[det_name][1].append(countData.getValue(GKTimeSerieIndex(interval))[0])

                    if _show_detector_data:
                        # print_cmd('--------interval {0}: speed {1}; count {2}'.format(interval,
                        #                                           avg_data[det_name][0][-1],
                        #                                           avg_data[det_name][1][-1]))
                        pass
                # print_cmd('----Detector {0} data:{1}'.format(det.getName(), avg_data[det.getName()])

    return avg_data


# evaluate the objective function
# this function compares the averaged simulation results with the validation data
# output one single valued measure, here we use RMSE
# the avg_data must have the same structure as the valid_dict: valid[det_name] = [[speed],[count]]
def evaluate_obj_val(avg_data, valid_dict, det_weight):

    rms_speed = 0
    rms_count = 0

    print_cmd('\nEvaluating objective value...')
    # print_cmd('avg_data keys: {0}'.format(avg_data.keys())
    # print_cmd('valid_data keys: {0}'.format(valid_dict.keys())

    # adds up the error for speed and count for each detector
    # Note the detector in the main entrance (EB3) is not used
    for key in valid_dict:
        valid_speed = np.array(valid_dict[key][0])
        valid_count = np.array(valid_dict[key][1])
        avg_speed = np.array(avg_data[key][0])
        avg_count = np.array(avg_data[key][1])

        # the following is added to deal with np.nan values
        # print 'before: {0}'.format(valid_speed-avg_speed)
        tmp_array = np.power(valid_speed-avg_speed, 2)
        # print 'after: {0}'.format(tmp_array)

        rms_speed += det_weight[key]*np.sqrt(np.nansum( tmp_array )/len(valid_speed))

        tmp_array = np.power(valid_count-avg_count, 2)
        rms_count += np.sqrt(np.nansum( tmp_array )/len(valid_count))

    print_cmd('Evaluated objective: (RMS_speed, RMS_flow): ({0}, {1})'.format(rms_speed, rms_count))

    return (rms_speed, rms_count)


# NOTE: This function can be used in two ways:
# ---- 1. generate validation data when valid_data is None
# ---- 2. Given any parameters, simulate and output objective function value, and data
def simulate_rep_from_paras(model, experiment, paras,
                             simulator, avg_result, plugin,
                            valid_data, det_weight,
                            name):

    if valid_data is not None:
        print_cmd('\n------------Simulating with {0} paras-------------------'.format(name))
    else:
        print_cmd('\n------------Generating validation data------------------')

    # set new parameters
    set_new_paras(model, experiment, paras)

    simulate_experiment(simulator, avg_result)

    sim_data = get_averaged_detector_data(model, avg_result, plugin)

    if valid_data is not None:
        obj_value = evaluate_obj_val(sim_data, valid_data, det_weight)
    else:
        obj_value = [0, 0]

    return [obj_value, sim_data]




#=================================================================================================
# The following functions are used for interfacing the python script with OptQuest Solver, and deal with
# preset parameter files.
# this function reads the parameters
# return: paras a dict; paras['car'] = (speed_accept, dev_speed_accept) or (reaction_time, max_accl)
# paras is a dict, with keys 'car','truck','ramp3','ramp4', each item has different length following the off->on order
def read_optquest_paras(parafile):

    paras = {}

    wait_time = 0
    timeout = 0
    while not exists(parafile):
        time.sleep(0.1)
        wait_time += 1
        timeout += 1
        if wait_time >= 10: # sleep 1 second
            print_cmd('Waiting for paras...')
            wait_time = 0

        if timeout >= 200: # 20 s
            # assume finished optimization
            return None

    if exists(parafile):

        paras = paras_reader(parafile)

    # delete the file once read
    os.remove(parafile)

    # print_cmd('Have read paras:\n')
    # for key in paras.keys():
    #     print_cmd('---- {0}: {1}'.format(key, paras[key]))

    return paras


# This function read optquest solutions
# OptQuest will stop when its max iteration hits
# solution is also a dict, has the same structure as paras
def read_optquest_solution(solutionfile):

    # Have not finished optimization
    if not exists(solutionfile):
        return None

    # a solution has been converged
    else:

        solution = paras_reader(solutionfile)

        # delete the file once read
        # os.remove(solutionfile)

        print_cmd('Have read solution:\n')
        for key in solution.keys():
            print_cmd('---- {0}: {1}'.format(key, solution[key]))

        return solution


# This function read preset_paras, which will be used to set the experiment
def read_preset_paras(presetparasfile):

    if not exists(presetparasfile):
        print_cmd('\nWARNING: No preset paras set \n ---- could not find file {0}\n'.format(presetparasfile))
        return None

    else:
        preset_paras = paras_reader(presetparasfile)

        print_cmd('Reading the preset paras...\n')
        for key in preset_paras.keys():
            pass
            # print_cmd('---- {0}: {1}'.format(key, preset_paras[key]))

        return preset_paras


# a universal paras reader, for preset_paras, paras, and solutions
# it assumes the file exists
# All types of paras are in the format:
# paras['car/truck/section_paraname'] = [list]
# Gerneral rule is to use the same name as in the GK manual
# input: each line: key, the rest is a list of values
def paras_reader(parafile):

    paras = {}

    f = open(parafile, 'r')

    for line in f:
        line = line.strip()
        items = line.split(',')

        # first item is the key
        paras[items[0]] = []

        for i in range(1, len(items)):
            paras[items[0]].append(float(items[i]))

    f.close()

    return paras


# this function writes the objective function value to file simval.txt
def write_simval(simval, simvalfile):

    f = open(simvalfile, 'w')
    f.write(str.format("{0:.16f}", simval))
    f.close()
    print_cmd('Wrote objective value: {0}'.format(simval))

    return 0



#=================================================================================================
# The following functions loads the previous results: opt_log.csv file to avoid redo the simulations
def load_previous_opt_solutions(opt_step_file):

    if opt_step_file is None:
        return None

    f = open(opt_step_file, 'r')

    paras_list = []
    obj_value_list = []

    for line in f:

        paras = {}

        line = line.strip()
        items = line.split(',')

        # parse each line. first item is the counter, skip
        for i in range(1,len(items)):

            # a key
            if not is_number(items[i]):
                if items[i] == 'RMS':
                    # get the objective value and continue
                    obj_val = float(items[i+1])
                    break
                else:
                    # register key
                    key_name = items[i]
                    paras[key_name] = []

            else:
                # a number, append to the latest key
                paras[key_name].append(float(items[i]))

        # print_cmd('loaded paras: {0}'.format(paras))

        paras_list.append(paras)
        obj_value_list.append(obj_val)

    return [paras_list, obj_value_list]



# the following function tests if the string is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False



# this function tests if a para has been simulated; if so, return obj, otherwise return None
def try_get_obj_from_previous_solutions( solution_list, new_para ):

    # if not solution list
    if solution_list is None:
        return None

    for i in range(0, len(solution_list[0])):

        if new_para == solution_list[0][i]:
            return solution_list[1][i]

    # if none found
    return None




#=================================================================================================
# The following functions handles the logging data.
# this function logs the optimization step
# iteration_counter, car_para1, car_para2..., truck_para1, truck_para2, RMS_speed, RMS_count
def start_opt_logger(path_str, start_time):

    global opt_logger, opt_logger_file

    file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_opt_log'
    opt_logger_file = open(path_str + file_name + '.csv', 'wb')
    opt_logger = csv.writer(opt_logger_file)


def log_opt_step(opt_solution, iteration_counter, paras, RMS):

    list = []
    list.append(iteration_counter)

    #first append key, then parameters
    for key in paras.keys():
        list.append(key)
        for item in paras[key]:
            list.append(item)

    # Last two are reserved for the RMS result
    list.append('RMS')
    list.append(RMS[0])
    list.append(RMS[1])

    opt_solution.append(list)

    # write in file
    if opt_logger is not None:
        opt_logger.writerow(list)


# save logging of the optimization data
def stop_opt_log():

    opt_logger_file.close()


# save the optimization result; and the data generated
# save logging of the optimization data
def save_solution_data(solution, data, path_name, start_time, name):

    file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_sol_'
    f = open(path_name + file_name + name + '.csv', 'wb')
    writer = csv.writer(f)
    # the solution: key1, value, key2, value
    if solution is not None:
        list = []
        for key in solution.keys():
            list.append(key)
            for item in solution[key]:
                list.append(str(item))
        writer.writerow(list)
    else:
        # just to save the validation data. true parameters are unknown
        writer.writerow(['this saves {0} data with unknown parameters'.format(name)])

    for key in data:
        tmp_line = []
        tmp_line.append(key)
        # write speed
        for item in data[key][0]:
            tmp_line.append(item)
        writer.writerow(tmp_line)

        tmp_line = []
        tmp_line.append(key)
        # write count
        for item in data[key][1]:
            tmp_line.append(item)
        writer.writerow(tmp_line)
        # writer.writerow([str(key), str(data[key][0]), str(data[key][1])])

    f.close()



# start logger for the cmd output
def start_cmd_logger(path_str, start_time):

    global cmd_logger

    file_name = start_time.strftime("%Y%m%d_%H%M%S") + '_cmd_log'
    cmd_logger = open(path_str + file_name + '.txt', 'wb')

    # return cmd_logger


# command logger header
# para_list =[ (name, paras) ]
def cmd_logger_header(description,
                      g_num_iter, g_num_rep, seed_list,
                      obj_func,
                      det_for_validation,
                      main_entrance_id, main_entrance_discount,
                      paras_list):

    print_cmd('Calibration experiment description:\n ---- {0}\n'.format(description))
    print_cmd('Calibration Configuration:')
    print_cmd('---- Objective function is:           Minimize    {0}xRMS_speed + {1}xRMS_count'.format(obj_func[0], obj_func[1]))
    print_cmd('---- Detectors used for validation:   {0}'.format(det_for_validation))
    print_cmd('---- Number of iterations:            {0}'.format(g_num_iter))
    print_cmd('---- Main Entrance is:                {0}'.format(main_entrance_id))
    print_cmd('---- Main Entrance flow = VerMac EB3 x{0}'.format(main_entrance_discount))
    print_cmd('---- Number of replications {0}, with seeds: {1}\n'.format(g_num_rep, seed_list))

    print_cmd('\nParameters:')
    for para_tup in paras_list:
        print_cmd('-- {0} paras:'.format((para_tup[0])))
        for key in para_tup[1].keys():
            print_cmd('---- {0}: {1}'.format(key, para_tup[1][key]))


# print out on the cmd, and save the cmd output to a file
# If log_file is None, will only print
def print_cmd(line_str):

    if cmd_logger is None:
        print line_str
    else:
        # save every time
        # Hence even if the file is not correctly closed, those lines are still saved.
        cmd_logger.write(line_str + '\n')
        print line_str

# stop logger for the cmd output
def stop_cmd_logger():

    cmd_logger.close()


# print out result
def print_results(paras_list, obj_val_list):
    print_cmd('\n\n===========================================================================')
    print_cmd('====================Calibration Finished===================================')
    # Forget about beautiful printout. Just log information
    print_cmd('Parameters: \n')
    for para in paras_list:
        print_cmd('---- {0}_paras:'.format(para[0]))
        for key in para[1].keys():
            print_cmd('-------- {0}:    {1}'.format(key, para[1][key]))

    print_cmd('\nObjective values: \n')
    for obj_val in obj_val_list:
        print_cmd(  '{0} objective value:\n---- RMS_speed: {0}      \n---- RMS_count: {1}'.format(obj_val[0],
                                                                                                obj_val[1][0],
                                                                                                obj_val[1][1]))


# This function plots the optimization steps.
# input: the opt_logger as we update in every iteration,
#        the obj_ratio (1,0) (speed, count), and the obj_val using default parameters used as a baseline
# TODO: this function somehow fails due to matplotlib issue in windows
def plot_opt_steps(opt_logger, obj_ratio, obj_default_val):

    global fig_handle

    # compute the objective values to be plotted
    steps = []
    obj_values = []
    for paras in opt_logger:
        steps.append(float(paras[0]))
        obj_values.append(float(paras[-2])*obj_ratio[0] + float(paras[-1]*obj_ratio[1]))

    default_baseline = float(obj_default_val[0])*obj_ratio[0] + float(obj_default_val[1])*obj_ratio[1]

    if fig_handle is None:
        fig_handle, = plt.plot(obj_values)
        fig_handle.set_xlabel('Iteration step')
        fig_handle.set_ylabel('Objective value')
        fig_handle.set_title('Optimization progress (default val {0})'.format(default_baseline))
    else:
        # update and plot the figure
        # fig_handle.set_xdata(steps)
        fig_handle.set_ydata(obj_values)

    plt.draw()
    plt.show()


# just print out the objective value
# obj_list is a list of [(name, value)] would like to print out and compare
def print_opt_steps(opt_logger, obj_ratio, obj_list):

    # compute the objective values to be plotted
    steps = []
    obj_values = []
    for paras in opt_logger:
        steps.append(float(paras[0]))
        obj_values.append(float(paras[-2])*obj_ratio[0] + float(paras[-1]*obj_ratio[1]))

    print_cmd('\nOpt Steps:   {0}'.format(obj_values))
    print_cmd('Opt optimal: {0}'.format(np.min(np.array(obj_values))))
    for item in obj_list:
        compare_baseline = float(item[1][0])*obj_ratio[0] + float(item[1][1])*obj_ratio[1]
        print_cmd('Opt {0}: {1}'.format(item[0], compare_baseline))



# adjust the demand data
# paras['ramp_on1_car']
# no longer used
def adjust_ramps_I80_EB_full(model, paras):

    state_type = model.getType("GKTrafficState")

    # print 'paras: {0}'.format(paras)

    # originally 500 veh/hr
    onramp_3_1 = findSection(model, 21192)    # [0, 500] veh/hr
    onramp_3_2 = findSection(model, 21201)

    # first off ramp at junction 3
    diverge_3_1_main = findSection(model, 3412)
    diverge_3_1_to = findSection(model, 3399)
    diverge_3_1_off = findSection(model, 343) # [0,2]

    # second off ramp at junction 3
    diverge_3_2_main = findSection(model, 3400)
    diverge_3_2_to = findSection(model, 3401)
    diverge_3_2_off = findSection(model, 1039)

    # onramp of 4
    onramp_4 = findSection(model, 1357)

    # originally 3%
    diverge_4_main = findSection(model, 40962)
    diverge_4_to = findSection(model, 3248)
    diverge_4_off = findSection(model, 1501)    # [0, 20]

    for state in model.getCatalog().getObjectsByType(state_type).itervalues():
        # print_cmd('state.getVehicle(): {0}'.format(state.getVehicle().getId() ))
        if state.getVehicle().getId() == 53:

            # ramp3, first off and on
            state.setTurningPercentage(diverge_3_1_main, diverge_3_1_off, None, float(paras['ramp3'][0]))
            state.setTurningPercentage(diverge_3_1_main, diverge_3_1_to, None, 100 - float(paras['ramp3'][0]))
            state.setEntranceFlow(onramp_3_1, None, float(paras['ramp3'][1]))

            # ramp3, second off
            state.setTurningPercentage(diverge_3_2_main, diverge_3_2_off, None, float(paras['ramp3'][2]))
            state.setTurningPercentage(diverge_3_2_main, diverge_3_2_to, None, 100 - float(paras['ramp3'][2]))

            # ramp4, off and on
            state.setTurningPercentage( diverge_4_main, diverge_4_off, None, float(paras['ramp4'][0]) )
            state.setTurningPercentage( diverge_4_main, diverge_4_to, None, 100 - float(paras['ramp4'][0]))
            state.setEntranceFlow(onramp_4, None, float(paras['ramp4'][1]))

            # main entrance flow
            # state.setEntranceFlow






# add a noise model for the validation data.
# noise_mode is speed_noise_model[det] = [bias, normal_distr_dev]
def add_noise_to_data(data, speed_noise_model, count_noise_model):

    new_data = {}

    for det in data.keys():

        speed = data[det][0]
        count = data[det][1]

        new_speed = []
        new_count = []

        for i in range(0, len(speed)):
            new_value = speed[i] + speed_noise_model[det][0] + np.random.normal(0, speed_noise_model[det][1], 1)
            new_speed.append( np.max([0, new_value]) )
            new_value = count[i] + count_noise_model[det][0] + np.random.normal(0, count_noise_model[det][1], 1)
            new_count.append( np.max([0, new_value])  )

        new_data[det] = [ new_speed, new_count ]

    return new_data












