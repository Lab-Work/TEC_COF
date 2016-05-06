__author__ = 'Yanning Li'
# Sep 20, 2015
# This code is to automate the Fundamental Diagram calibration for the freeway and the onramp.
# Key parameters: vf, qmax, rho_c, rho_m
# The calibration is done on a merge with a one-lane freeway and a one-lane ramp merging to a one-lane freeway.
# Calibrating the freeway parameters:
#   - slowly increasing the freeway inflow (100~2500) while keep onramp flow all 0. For calibrating the freeflow side
#   - Keep the upstream freeway maximal inflow; gradually increase the onramp flow to cause congestion. This is to
#     calibrate the congested side.
# Calibrating the onramp parameters:
#   - Keep upstream freeway inflow 0, and slowly increase onramp in flow from (0,2500). Freeflow side.
#   - Keep upstream max inflow, increase the upstream freeway inflow to create congestion on the onramp.

# Features:
# 1. Generate desired demand using a new function generateDemandFD(), in which we can define the demand
# 2. Run 10 replication. Save the speed and flow data for each replication from two detectors located
#    in the middle of the upstream freeway and the downstream freeway.
# 3. A function plot the fundamental diagram scatter plot for the upstream freeway and the onramp link.
# 4. Print the linearly fitted parameters for the two fundamental diagrams.


from AIMSUNFUNCTIONS_V3 import *
from datetime import timedelta


# number of replications
g_num_rep = 10
seed_list = [3503, 23798, 28860, 12358, 6370, 14452, 21488, 7893, 586, 30447]

# Configuration===========================================================
parsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_traffic_control_simulation/FD_calibration/paras.txt'
flowsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_traffic_control_simulation/FD_calibration/flows.txt'
turnsFilePath = 'c:/Users/TrafficControl/Google Drive/AIMSUN_traffic_control_simulation/FD_calibration/turns.txt'
logger_path = 'c:/Users/TrafficControl/Google Drive/AIMSUN_traffic_control_simulation/FD_calibration/Logs/'
data_path = 'c:/Users/TrafficControl/Google Drive/AIMSUN_traffic_control_simulation/FD_calibration/det_data/'

us_freeway_id = 329
ds_freeway_id = 330
onramp_id = 390

calibrated_link = 'freeway'
traffic_state = 'congflow'
all_det_strs =   ['freeway_FD']

# detector used for validation
det_used = all_det_strs
det_interval = 60


# Put one line of description of the simulation
description = 'Calibrating the fundamental diagram of the AIMSUN links.\n'

_debug_FD = True


# Load a network using ANGConsole
def main(argv):

    # keep track of start time:
    start_time = datetime.now()

    # Start the cmd output logger
    start_cmd_logger(logger_path, start_time)
    start_opt_logger(logger_path, start_time)

    # print python version
    # print '\nPython interpreter version: {0} \n'.format(sys.version)

    print_cmd('\n\n========================================================================')
    print_cmd('=========================AIMSUN FD calibration=========================')
    print_cmd('========================================================================')
    print_cmd('FD calibration started at {0}\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    # Auto write the header to the log file
    print_cmd('Calibrating the {0} side of the fundamental diagram of one lane {1} link\n'.format(traffic_state, calibrated_link))



    if (len(argv) < 2):
        print_cmd('Usage: aconsole.exe -script SCRIPT ANG_FILE')
        return -1
    else:
        # Start the Aimsun console
        console = ANGConsole()
        if console.open(argv[1]):

            print_cmd('\nAimsun opening {0} ...\n'.format(argv[1]))

            #========================================================
            # Read demand
            # Get the current Aimsun  model
            model = GKSystem.getSystem().getActiveModel()

            demand = addDemandFD(model, traffic_state+'_demand')

            print_cmd('Finished setting Demand. \n')

            #========================================================
            # Setup model
            # create scenario
            # if exists, then just load
            scenario = setup_scenario(model, traffic_state, demand)

            # create experiment
            experiment = model.getCatalog().findByName( QString(traffic_state+'_exp') )
            if experiment is None or not experiment.isA(QString("GKExperiment")):
                experiment = GKSystem.getSystem().newObject("GKExperiment", model, -1, True)
                print_cmd('ERROR: No traffic experiment named {0}. Creating a new one...\n'.format(traffic_state+'_exp'))

                # attach the new experiment to folder
                folder = getScenariosFolder(model)
                folder.append(experiment)

            experiment.setScenario(scenario)

            # create replications
            avg_result = setup_replication(model, experiment, g_num_rep, seed_list)


            # save the ang file
            console.save("{0}_{1}.ang".format(calibrated_link, traffic_state))

            # create simulator
            simulator = create_simulator(model)
            # plugin is a module which can compute the average for the GKExperimentResult object
            plugin = GKSystem.getSystem().getPlugin( "GGetram" )


            simulate_experiment(simulator, avg_result)

            calculate_status = plugin.calculateResult( avg_result )
            plugin.readResult(avg_result)

            sim_data = {}
            for det in det_used:
                sim_data[det] = read_detector_data_FD(model, experiment, det)

                # save data
                saveDetectionData(calibrated_link, traffic_state, det, sim_data[det])

            # plot the scatter figure

            # plotFD(calibrated_link,  sim_data)



            print_cmd(  'Fundamental Diagram Calibration ended')
            print_cmd(  '===========================================================================')


            print_cmd(  '\nAimsun is now closing...\n')

            # stop logger and save into file
            stop_cmd_logger()

            console.close()

        else:
            console.getLog().addError("Could not open")


# generate the demand
def addDemandFD(model, demand_name):

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

    # Generate the demand as needed
    pars_dict, flows_dict ,turns_dict = generateDemandFD()

    for state_name in pars_dict.keys():

        state = createState(model, state_name, pars_dict, flows_dict, turns_dict, None, 1)

        if _debug_FD:
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

    return demand



# generate demand pars, flows, and turns dict
def generateDemandFD():

    car_string = '53 Car'
    start_time = datetime(2015,9,20,0,0,0)
    duration = datetime(10,10,10,0,30,0)

    pars_dict = {}
    flows_dict = {}
    turns_dict = {}

    flow = 100

    for i in range(0,24):

        state = 'state_{0}'.format(i)

        # set the paras
        pars_dict[state] = list()
        pars_dict[state].append((car_string, '{0}'.format(start_time.time()),'{0}'.format(duration.time())))
        start_time = start_time + timedelta(0,60*30)

        # set the flow
        flows_dict[state] = list()

        # For freeway freeflow
        if calibrated_link == 'freeway' and traffic_state == 'freeflow':
            flows_dict[state].append(('{0}'.format(us_freeway_id), '{0}'.format(flow)))
            flows_dict[state].append(('{0}'.format(onramp_id), '0'))
            flow += 100

        # For freeway congflow
        if calibrated_link == 'freeway' and traffic_state == 'congflow':
            flows_dict[state].append(('{0}'.format(us_freeway_id), '2400'))
            flows_dict[state].append(('{0}'.format(onramp_id), '{0}'.format(flow)))
            flow += 40

        # For onramp freeflow
        if calibrated_link == 'onramp' and traffic_state == 'freeflow':
            flows_dict[state].append(('{0}'.format(us_freeway_id), '0'))
            flows_dict[state].append(('{0}'.format(onramp_id), '{0}'.format(flow)))
            flow += 90

        # For onramp congflow
        if calibrated_link == 'onramp' and traffic_state == 'congflow':
            flows_dict[state].append(('{0}'.format(us_freeway_id), '{0}'.format(flow)))
            flows_dict[state].append(('{0}'.format(onramp_id), '2200'))
            flow += 100


    return pars_dict, flows_dict, turns_dict




# each time can only read one detector
def read_detector_data_FD(model, experiment, detector_name):

    avg_data = {}

    det_type = model.getType("GKDetector")
    # read data for each replication and then the average
    for replication in experiment.getReplications():

        if str(replication.getName()) != 'average_result':

            print_cmd('Reading Replication data: ID {0}'.format(replication.getId()))

            # get the column id
            speedColumn = det_type.getColumn(GK.BuildContents(GKColumnIds.eSpeed, replication, None))
            countColumn = det_type.getColumn(GK.BuildContents(GKColumnIds.eCount, replication, None))

            # read each detector
            for det in model.getCatalog().getObjectsByType(det_type).itervalues():

                det_name = str(det.getName())
                # print_cmd('----Reading Detector {0}...'.format(det_name))

                if det_name == detector_name:

                    # add to dict
                    # flow, speed
                    avg_data[replication] = [[],[]]

                    speedData = det.getDataValueTS(speedColumn)
                    countData = det.getDataValueTS(countColumn)

                    if countData.size() == 0 or speedData.size() == 0 or countData.size() != speedData.size():
                        print_cmd('ERROR: Detector {0} has no data available'.format(det_name))
                    else:
                        # print_cmd('----size of data is: {0}'.format(countData.size()))
                        # TODO: DONE: the speed data returned from AIMSUN is in km/h; 1 km/h = 0.62137 mph
                        for interval in range(countData.size()):
                            avg_data[replication][0].append(speedData.getValue(GKTimeSerieIndex(interval))[0]*KMH2MPH)
                            avg_data[replication][1].append(countData.getValue(GKTimeSerieIndex(interval))[0])

                            if _debug_FD:
                                # print_cmd('--------interval {0}: speed {1}; count {2}'.format(interval,
                                #                                           avg_data[replication][0][-1],
                                #                                           avg_data[replication][1][-1]))
                                pass
                        # print_cmd('----Detector {0} data:{1}'.format(det.getName(), avg_data[det.getName()])

    return avg_data




# This function saves the detection data for each detector and state
# input:
#       calibrated_link: the current link that is calibrated
#       traffic_state: the traffic state of the FD
#       det_name: the detector name
#       data: the avg data, a dict, ['replication'] = [speed (mph), count (/interval)]
# output:
#       each row: [speed (miles/hr), flow (veh/hr)]
def saveDetectionData(link, state, det_name, data):
    # data is a list of string
    data_to_write = []

    for rep in data.keys():

        count = np.array(data[rep][1])
        for row in range(0, len(data[rep][0])):
            line = '{0},{1}\n'.format(data[rep][0][row], count[row]*3600/det_interval)
            data_to_write.append(line)

    # write the data into file
    f = open(data_path  + link + '_' + det_name + '_' + state + '.txt','wb')
    for line in data_to_write:
        f.write(line)
    f.close()



# This function plots the scatter plot of the Fundamental diagram.
# input:
#       link: the link to be calibrated
#       data: A dict of data
#            data['det_name']['rep_name'] = [speed (mph), count (/interval)]
def plotFD(link,  data_all):

    speed = []
    flow = []

    for key in data_all.keys():

        det_name = key

        # stack data
        for rep in data_all[key].keys():

            count = data_all[key][rep][1]

            for row in range(0, len(count)):
                speed.append(data_all[key][rep][0][row])
                flow.append(count[row]*3600/det_interval)

        # comptue the density veh/mile
        density = np.array(flow)/np.array(speed)


        # speed
        fig_window = plt.figure(figsize=(15,10))
        fig = fig_window.add_subplot(111)

        plt.plot(density, flow, '.b')
        # plt.plot(SB1_cong_density, flow_cong, '*r')
        plt.title('FD for {0} from det {1}'.format(link, det_name), fontsize=24)
        plt.xlabel('Traffic density (veh/mile)', fontsize=24)
        plt.ylabel('Traffic flow (veh/hr)', fontsize=24)
        # plt.xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) , np.max(x) + 0.1*(np.max(x)-np.min(x))])
        plt.xlim([0, 30])
        plt.ylim([0, 2000])
        plt.grid(True)
        fig.tick_params(axis='both', which='major', labelsize=20)

        plt.draw()

    plt.show()





if __name__ == "__main__":
    sys.exit(main(sys.argv))
