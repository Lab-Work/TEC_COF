__author__ = 'Yanning Li'

# This scripts assumes the follows:
# - The network contains only one merge: a two-lane freeway and an one-lane onramp merges to the two-lane downstream freeway
# - The onramp is single lane, and freeways are 2 lanes
# - A work zone exists in one segment of the downstream freeway which reduces the number of lanes to 1
# - Each detector across all lanes at its location
# - Only one meter exist which is at the exit of the onramp
# - The meter is preconfigured as a flow meter, we only need to update is flow as needed.
# - We will use the upstream most sensors on the upstream freeway and the onramp, the sensor right before work zone, and
#   the sensors at the junction to run the MPC.



from AAPI import *

_debug = True
MILE2METER = 1609.34

# This function is called when this module is correctly loaded
def AAPILoad():
    AKIPrintString( "AAPILoad" )
    return 0

# This function is called to perform any initialization
def AAPIInit():
    AKIPrintString( "AAPIInit" )

    # export the network to MATLAB
    exportNetwork()

    return 0


# This function is called before every time step
# 1. It should first check if a new signal file is written by MATLAB
# 2. If written, update is signal array.
# 3. Set the meter according to the signal array.
def AAPIManage(time, timeSta, timeTrans, acycle):
    # AKIPrintString( "AAPIManage" )

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
    return 0

def AAPIUnLoad():
    # AKIPrintString( "AAPIUnLoad" )
    return 0




# This function read the network configuration information and write to a file for MATLAB
# 1. It read the network topology, including the links and juncs
# 2. It read the detector information and save them in the global variable
# 3. It read the meter information
# 4. All read information will be exported to a txt file, which can be read by MATLAB.
def exportNetwork():

    # find out what unit the network is using:
    # 1-Metric; 0-English
    unit = AKIInfNetGetUnits()

    # net is a list of strings
    net = []

    # read link section information
    num_section = AKIInfNetNbSectionsANG()
    if _debug:
        AKIPrintString( 'The number of sections:{0}\n'.format(num_section))

    us_fwy_id = None
    ds_fwy_id = None
    onramp_id = None

    for i in range(0, num_section):
        link_id = AKIInfNetGetSectionANGId(i)

        AKIPrintString('Link Id: {0}'.format(link_id))

        # read link information
        link_info = AKIInfNetGetSectionANGInf(link_id)

        if _debug:
            AKIPrintString( '--- Link {0}: speedLimit, {1} kph; capacity {2} veh/h; length {3} meters; num_lanes {4}'.format(
                link_info.id, link_info.speedLimit, link_info.capacity, link_info.length, link_info.nbCentralLanes
            ) )

            AKIPrintString( '--- Link {0}: num_turns {1}'.format(
                link_info.id, link_info.nbTurnings
            ) )

        # if link_info.nbTurnings != 0:
        #
        #     for j in range(0, link_info.nbTurnings ):
        #
        #         AKIPrintString('------ {0} turns to {1}'.format(link_info.nbTurnings,
        #                                                     AKIInfNetGetIdSectionANGDestinationofTurning(link_id, j)))

        link_type = None
        lengthKM = None
        vf= None
        w= None
        kc_pl= None
        km_pl= None
        qmax_pl= None
        v_min= None
        v_max = None
        # Identify the link based on the assumptions we made, save them into file
        if link_info.nbCentralLanes > 1 and link_info.nbTurnings != 0:
            # the upstream freeway
            us_fwy_id = link_id
            link_type = 'freeway'

            # set corresponding parameters which need to be calibrated
            # vf, w, kc_pl, km_pl, qmax_pl, v_min, v_max

        if link_info.nbCentralLanes == 1:
            # the onramp link
            onramp_id = link_id
            link_type = 'onramp'

            # set corresponding parameters which need to be calibrated
            # vf, w, kc_pl, km_pl, qmax_pl, v_min, v_max

        if link_info.nbCentralLanes > 1 and link_info.nbTurnings == 0:
            # the downstream freeway
            ds_fwy_id = link_id
            link_type = 'freeway'

            # set corresponding parameters which need to be calibrated
            # vf, w, kc_pl, km_pl, qmax_pl, v_min, v_max

        if unit == 1:
            lengthKM = link_info.length/1000.0
        elif unit == 0:
            lengthKM = link_info.length*MILE2METER/1000.0


        line = 'link;link_id,{0};num_lanes,{1};lengthKM,{2};link_type,{3};vf,{4};' + \
                'w,{5};kc_pl,{6};km_pl,{7};qmax_pl,{8};v_min,{9};v_max,{10}\n'.format(link_id,
                link_info.nbCentralLanes, lengthKM, link_type, vf, w, kc_pl, km_pl, qmax_pl, v_min, v_max )

        net.append(line)


    # collect the junction info (only one junction)
    # junc_id = AKIInfNetGetJunctionId(0)
    ratio = [1,1]
    line = 'junc;junc_id,{0};inlink,{1},{2};outlink:{3};assignment_ratio,{4},{5};junc_type,onrampjunc'.format(
        1, us_fwy_id, onramp_id, ds_fwy_id, ratio[0], ratio[1]
    )

    net.append(line)



    # # read junc section information
    # num_junc = AKIInfNetNbJunctions()
    # if _debug:
    #     AKIPrintString('The number of junctions:{0}'.format(num_junc))
    # for i in range(0, num_junc):
    #     junc_id = AKIInfNetGetJunctionId(i)
    #     num_turns = AKIInfNetGetNbTurnsInNode(junc_id)
    #     if _debug:
    #         AKIPrintString('--- junction {0}: num_turns: {1}'.format(junc_id, num_turns))
    #
    #     for j in range(0, num_turns):
    #
    #         turn_info = AKIInfNetGetTurnInfo(junc_id, j)
    #         if _debug:
    #             AKIPrintString('------ Turn {0}: from {1} to {2}'.format(turn_info.id,
    #                                                                      turn_info.originSectionId,
    #                                                                      turn_info.destinationSectionId))
    #             AKIPrintString('------ Turn2 {0}: from {1} to {2}'.format(turn_info.id,
    #                                                                      AKIInfNetGetOriginSectionInTurn(junc_id, j),
    #                                                                      AKIInfNetGetDestinationSectionInTurn(junc_id, j)   ))



    # read detector information
    # num_det = AKIDetGetNumberDetectors()
    # if _debug:
    #     AKIPrintString( 'Number of detectors:{0}'.format(num_det))
    # # for each detector
    # for i in range(0, num_det):
    #     det_id = AKIDetGetIdDetector(i)
    #
    #     det_info = AKIDetGetPropertiesDetectorById(det_id)
    #
    #     if _debug:
    #         AKIPrintString('The detector {0}: '.format(det_id))



    # read meter information
    num_meter = ECIGetNumberMeterings()
    if _debug:
        AKIPrintString( 'Number of meters: {0}'.format(num_meter))

    # for each meter
    for i in range(0, num_meter):

        meter_id = ECIGetMeteringIdByPosition(i)
        meter_name = ECIGetMeteringNameById(meter_id)

        # 1-green; 2-flow; 3-delay; 4-ALINEA; 5-greenByLane
        meter_type = ECIGetTypeMeteringById(meter_id)

        if _debug:
            AKIPrintString('Meter {0} has name {1} and type {2}'.format(meter_id, meter_name, meter_type))

        flow_max = doublep()
        flow_min = doublep()
        flow = doublep()

        if meter_type == 2:
            # try to set the flow as 0
            current_config = ECIGetParametersFlowMeteringById(meter_id, 0, flow_max, flow, flow_min)


    # export the net to txt file
    # f = open(files['net'],'wb')
    # for line in net:
    #     f.write(line)
    #
    # f.close()
    #
    # # flag write done
    # f = open(files['net_write_done'],'wb')
    # f.write('net write done\n')
    # f.close()







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
