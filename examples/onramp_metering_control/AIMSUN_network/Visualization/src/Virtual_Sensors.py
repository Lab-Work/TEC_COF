# Author: Juan Carlos Martinez
# 3/21/2016


import numpy as np
import sys
import os
from os.path import exists
import time
import math
import csv
from collections import OrderedDict
import random
import warnings
import matplotlib.pyplot as plt
from default_parameters import default_parameters
from scipy.interpolate import InterpolatedUnivariateSpline
warnings.filterwarnings("error")

"""
    class virtual_sensors
    
        PUBLIC
        
            __init__(self, work_zone, time_step, directory, sections, main_grid_order):
                loads the default simulation and virtual sensor parameters
                loads the trajetories data
                organizes the data
                
            generate_virtual_sensors_data()
                parses the information of the sensors specified in the to_generate file
                builds a crosses matrix that holds the crossing information for a sensor
                generates the virtual_sensor based on that crossing matrix
                save the results in the appropriate files
                
            generate_true_states_data(grid_resolution, plot)
                generate maps for speed and density
                generate true state for travel time # PENDING
                save the results in the appropriate files and plot if specified
                
            
        PRIVATE
        
            __load_data
                loads the trajectories data
                
            __organize_data
                organizes the data by the following dictionary structure:
                self.__organized_data[replication][vehicle][timestamp]
                the value for each key is an array with the following entries:
                [ent, section, section_distance, travelled_distance, lane_idx, speed, acceleration]

            __build_volume_sensor
                build the sensor crosses and readings for each simulation replication
                save the generated data in the appropriate files

            __build_sensor_crosses
                build the sensor crosses matrix
                use binary search to find the crossing time for the upstream bound of the reading area
                iterate after the found timestamp to find the crossing of the downstream bound

            __add_lane_to_sensor_location
                compute the locations of the upstream and downstream visibility boundaries of a sensor

            __build_volume_sensor_readings
                obtain noise parameters, check for occlusion and build sensor readings

            __occlusion_check
                check if occlusion occurs for a given vehicle cross on a sensor visibility region

            __build_travel_time_sensor
                build sensor crossings at the upstream and downstream sensor
                build the travel time pair readings
                save the generated data in the appropriate file

            __build_travel_time_sensor_readings
               obtain travel time noise parameters
               filter vehicles recorded based on penetration rate
               add the travel times of individual vehicles and compute means for each interval

            __build_true_speed_density_states
                set sensor locations (cells) along the main road
                call helper function to build true state at each cell
                plot if specified by the user

            __build_true_speed_density_states_readings
                organize the sensor locations
                generate no-noise sensor readings for each location along the road
                stack the data into speed and density stacks (each row is a time interval)   
    
            
        """


class Virtual_Sensors:

    def __init__(self, work_zone, log_dir, data_dir, sections, main_grid_order, replications,
                        aimsun_start_dur_step=[55800,9000, 1.0], space_start_end=None):
        """

        :param work_zone: str, 'I57', 'I80'
        :param time_step:
        :param log_dir: the top level directory of the logs.
        :param data_dir: the top level directory of the actual data, default E:\\Workzone\\
        :param sections:
        :param main_grid_order:
        :param aimsun_start_dur_step: [start, dur, step ] the start timestamp, duration and step length in seconds
        :param space_start_end: the absolute start and end location of the space domain under investigation
        :return:
        """

        time0 = time.time()
        print('\n')
        print('Constructing virtual_sensors class...')
        
        self.__work_zone = work_zone
        self.__time_step = aimsun_start_dur_step[2]
        self.__sections = sections
        self.__main_grid_order = main_grid_order
        self.__space_start = space_start_end[0]
        self.__space_end = space_start_end[1]
        self.replications = replications

        self.__start_timestamp = aimsun_start_dur_step[0]
        self.__end_timestamp = self.__start_timestamp +  aimsun_start_dur_step[1]
        
        # default parameters for the replication and the sensor types
        from default_parameters import default_parameters
        self.__default_parameters = default_parameters

        # Make the directories a dict with keys as the replication id
        self.__trajectories_file = OrderedDict()
        self.__virtual_sensors_data_folder = OrderedDict()
        self.__true_state_data_folder = OrderedDict()
        self.__virtual_sensors_to_generate_file = OrderedDict()
        self.__virtual_sensors_generated_file =  OrderedDict()

        if sys.platform == 'win32':
            for rep in self.replications:
                self.__trajectories_file[rep] = data_dir + 'Trajectory_data\\' + '{0}_rep{1}.csv'.format(self.__work_zone, rep)
                self.__virtual_sensors_data_folder[rep] = data_dir + 'Virtual_sensor_data\\' + self.__work_zone + '\\' + 'rep{0}'.format(rep) + '\\'
                self.__true_state_data_folder[rep] = data_dir + 'Result\\' + self.__work_zone + '\\' + 'rep{0}'.format(rep) + '\\'
                self.__virtual_sensors_to_generate_file[rep] = log_dir + 'Virtual_sensor_data\\' +  self.__work_zone + '_' + 'rep{0}'.format(rep)  + '_to_generate.txt'
                self.__virtual_sensors_generated_file[rep] =  log_dir + 'Virtual_sensor_data\\' +  self.__work_zone + '_' + 'rep{0}'.format(rep)  + '_generated.txt'
        elif sys.platform == 'darwin':
            for rep in self.replications:
                self.__trajectories_file[rep] = data_dir + 'Trajectory_data/' + '{0}_rep{1}.csv'.format(self.__work_zone, rep)
                self.__virtual_sensors_data_folder[rep] = data_dir + 'Virtual_sensor_data/' + self.__work_zone + '/' + 'rep{0}'.format(rep) + '/'
                self.__true_state_data_folder[rep] = data_dir + 'Result/' + self.__work_zone + '/' + 'rep{0}'.format(rep) + '/'
                self.__virtual_sensors_to_generate_file[rep] = log_dir + 'Virtual_sensor_data/' +  self.__work_zone +'_'+ 'rep{0}'.format(rep)  + '_to_generate.txt'
                self.__virtual_sensors_generated_file[rep] =  log_dir + 'Virtual_sensor_data/' +  self.__work_zone + '_' + 'rep{0}'.format(rep)  + '_generated.txt'
        else:
            raise Exception('Error: run "import sys; sys.platform" to check the platform. ')

        # Do Not load data until generating virtual sensor data and the true states.
        self.__raw_data = None


    def generate_virtual_sensors_data(self):

        time0 = time.time()
        print('Status: Generating virtual sensors data...')

        for rep in self.replications:

            if not exists(self.__virtual_sensors_to_generate_file[rep]):
                print('Status: no additional sensor data need to be generated for Replication {0}.'.format(rep))
                continue

            print('Status: generating sensor data for replication {0}...'.format(rep))

            if self.__raw_data is None:
                time0 = time.time()
                self.__load_data()
                time1 = time.time()
                print('Status: Loaded trajectory data. Took {0:.2f} s. Organizing...'.format(time1-time0))
                self.__organize_data()

                time1 = time.time()
                elapsed_time = '{0:.2f}'.format(time1 - time0)
                print('Elapsed time: ' + elapsed_time + ' s')


            with open(self.__virtual_sensors_to_generate_file[rep]) as sensors_configuration:
                for sensor in sensors_configuration:

                    sensor_parameters = self.__parse_sensor_line(sensor)
                    time_s = time.time()

                    # call appropriate helper functions
                    if (
                        sensor_parameters['type'] == 'icone' or
                        sensor_parameters['type'] == 'radar' or
                        sensor_parameters['type'] == 'rtms' or
                        sensor_parameters['type'] == 'vol_gt'
                        ):
                        self.__build_volume_sensor(sensor_parameters)
                    elif (
                          sensor_parameters['type'] == 'bluetooth' or
                          sensor_parameters['type'] == 'tt_gt'
                          ):
                        self.__build_travel_time_sensor(sensor_parameters)

                    time_e = time.time()
                    print('Status: -- generated {0} data. Took {1:.2f}.'.format(sensor_parameters['type'], time_e - time_s))
                    
        time1 = time.time()
        elapsed_time = '{0:.2f}'.format(time1 - time0)
        print('Elapsed time: ' + elapsed_time + ' s')


    def generate_true_states_data(self, grid_resolution, rep, queue_threshold, plot=False):

        if self.__raw_data is None:
            time0 = time.time()
            print('Status: loading data...')
            self.__load_data()
            time1 = time.time()
            print('Status: Loaded trajectory data. Took {0:.2f} s. Organizing...'.format(time1-time0))
            self.__organize_data()

            time1 = time.time()
            elapsed_time = '{0:.2f}'.format(time1 - time0)
            print('Elapsed time: ' + elapsed_time + ' s')


        time0 = time.time()
        print('\n')
        print('Generating true states data...')
        self.__build_true_flow_density_speed_states(grid_resolution, rep, queue_threshold, plot)
        time1 = time.time()
        elapsed_time = '{0:.2f}'.format(time1 - time0)
        print('Elapsed time: ' + elapsed_time + ' s')


    def __parse_sensor_line(self, sensor):
        """
        This function parses the sensor file whose data need to be generated
        :param sensor: the sensor line
        :return:
        """
        # Parse, load and update sensor parameters
        sensor_parameters = OrderedDict()
        sensor_parameters['custom'] = OrderedDict()
        for couple in sensor.split(';'):
            category = couple.split(':', 1)[0].rstrip()
            value = couple.split(':', 1)[1].rstrip()
            if value[0] == '[':
                value =  [float(i) for i in value.strip('[]').split(',')]
            elif value.find('.') != -1:
                value = float(value)
            elif value.isdigit():
                value = int(value)
            elif value == 'True':   # the occlusion
                value = True
            elif value == 'False':
                value = False

            if ( category != 'type' and category != 'id' and category != 'section' and category != 'distance' and
                 category != 'section_up' and category != 'section_down' and category != 'distance_up' and category != 'distance_down'
                ):
                sensor_parameters['custom'][category] = value
            else:
                sensor_parameters[category] = value
        sensor_parameters['default'] = self.__default_parameters[sensor_parameters['type']].copy()
        if sensor_parameters['custom']:
            for category in sensor_parameters['custom']:
                sensor_parameters['default'][category] = sensor_parameters['custom'][category]

        print('Status: -- generating data for sensor {0}...'.format(sensor_parameters['id']))

        return sensor_parameters


    def __load_data(self):

        self.__raw_data = np.zeros((0,11))
        for rep in self.replications:
            self.__raw_data = np.concatenate( [self.__raw_data,
                                               np.genfromtxt(self.__trajectories_file[rep], dtype='str', delimiter=',')])

        # timestamps = [float(i) for i in self.__raw_data[:,7]]
        # self.__start_timestamp = min(timestamps)
        # self.__end_timestamp = max(timestamps)


    def __organize_data(self):
        
        self.__organized_data = OrderedDict()
        
        # variable that keeps track how far in the network a vehicle has traversed
        length_prev_sections_and_tapers = 0.0
        
        for idx, line in enumerate(self.__raw_data):
            
            replication = int(line[0])
            vehicle = int(line[1])
            ent = int(line[2])
            section = int(line[3])

            # we only input the few sections that covers the road segment under estimation.
            if section not in self.__sections.keys():
                continue

            lane_idx = int(line[4])
            timestamp = float(line[7])
            speed = float(line[8])
            travelled_distance = float(line[9])
            acceleration = float(line[10])
            
            if replication not in self.__organized_data:
                self.__organized_data[replication] = OrderedDict()
            if vehicle not in self.__organized_data[replication]:
                self.__organized_data[replication][vehicle] = OrderedDict()

            
            if idx: 
                prev_line = self.__raw_data[idx-1,:]
                prev_replication = int(prev_line[0])
                prev_vehicle = int(prev_line[1])
                prev_section = int(prev_line[3])
                prev_timestamp = float(prev_line[7])
                prev_travelled_distance = float(prev_line[9])
                
                if replication != prev_replication or vehicle != prev_vehicle:
                    length_prev_sections_and_tapers = 0.0
                    
                elif section != prev_section:

                    # reconstruct the trajectories on the taper
                    length_prev_sections_and_tapers = length_prev_sections_and_tapers + self.__sections[prev_section]['length']
                    prev_section_remaining_length = length_prev_sections_and_tapers - prev_travelled_distance
                    
                    for node_section in self.__sections:
                        # find the right taper based on the sections upstream and downstream
                        # zero acceleration at the tapers is assumed (constant velocity)
                        if 'connections' in self.__sections[node_section]:    
                            if (
                                prev_section in self.__sections[node_section]['connections']['upstream'] and
                                section in self.__sections[node_section]['connections']['downstream']
                                ):      
                                node_lane_idx = lane_idx
                                node_speed_m_s =(travelled_distance - prev_travelled_distance)/(timestamp - prev_timestamp)
                                node_speed_km_h = float('{0:.2f}'.format(node_speed_m_s*3600/1000))
                                node_acceleration_m_s2 = 0.0
                                node_timestamps = np.arange(prev_timestamp + self.__default_parameters['time_step'], timestamp, self.__default_parameters['time_step'])       
                                for node_timestamp in node_timestamps:
                                    node_timestamp = float('{0:.2f}'.format(node_timestamp))
                                    # avoid floating point rounding error
                                    if abs(timestamp - node_timestamp) < self.__default_parameters['time_step']/2:
                                        pass
                                    else:
                                        node_section_distance = float('{0:.2f}'.format(node_speed_m_s*(node_timestamp - prev_timestamp) - prev_section_remaining_length))
                                        node_travelled_distance = float('{0:.2f}'.format(prev_travelled_distance + prev_section_remaining_length + node_section_distance))
                                        # if the vehicle is still on the previous section
                                        if node_section_distance < 0:
                                            selected_section = prev_section
                                            selected_section_distance = float('{0:.2f}'.format(self.__sections[prev_section]['length'] + node_section_distance))   
                                        else:
                                            selected_section = node_section
                                            selected_section_distance = node_section_distance
                                        self.__organized_data[replication][vehicle][node_timestamp] = [-1, selected_section, selected_section_distance, node_travelled_distance,
                                                                                                       node_lane_idx, node_speed_km_h, node_acceleration_m_s2]
                                length_prev_sections_and_tapers = length_prev_sections_and_tapers + self.__sections[node_section]['length']
                                # break and resume with next section
                                break
                            
            section_distance = float('{0:.2f}'.format(travelled_distance - length_prev_sections_and_tapers))
            # set trajectory points that have section_distance exceeding the specified section distance to have the maximum
            # allowable value
            # this takes care of possible inconcistencies in the trajectory data
            if section_distance > self.__sections[section]['length']:
                section_distance = self.__sections[section]['length']
            
            self.__organized_data[replication][vehicle][timestamp] = [ent, section, section_distance, travelled_distance, lane_idx, speed, acceleration]


    def __build_volume_sensor(self, sensor_parameters):

        for replication in self.__organized_data:
            
            # build crosses and readings
            sensor_crosses = self.__build_sensor_crosses(sensor_parameters, replication)
            sensor_readings = self.__build_volume_sensor_readings(sensor_crosses, sensor_parameters)
            
            sensor_file = self.__virtual_sensors_data_folder[replication] + str(sensor_parameters['id']) + '.txt'

            if not os.path.exists( os.path.dirname(sensor_file) ):
                os.makedirs(os.path.dirname(sensor_file))
            # os.makedirs(os.path.dirname(sensor_file), exist_ok=True)

            with open(sensor_file,'w+') as file:
                for interval in sensor_readings:
                    line = str([interval, float('%.2f' % sensor_readings[interval]['sms']),
                                sensor_readings[interval]['count']]).strip('[]').replace(' ','')
                    file.write(line + '\n')
            with open(self.__virtual_sensors_generated_file[replication], 'a+') as file:
                line = self.__sensordict2string( str(sensor_parameters['id']), sensor_parameters )
                file.write(line + '\n')


    def __build_sensor_crosses(self, sensor_parameters, replication):

        #v_hat = []
        #v_tilda = []

        # ordered section list of the main road segment (dictionary is built in constructor)
        sections_list = self.__main_grid_order
        
        alpha, theta = np.radians(sensor_parameters['default']['alpha']), np.radians(sensor_parameters['default']['theta'])
        s, L = sensor_parameters['default']['s'], sensor_parameters['default']['L']
        sensor_section, sensor_section_distance = sensor_parameters['section'], sensor_parameters['distance']
        
        sensor_location = OrderedDict()
        sensor_crosses = np.atleast_2d(np.empty((0,5)))

        #vehicles = []
        
        for vehicle in self.__organized_data[replication]:

            # binary search over the timestamps of a vehicle
            try:
                timestamps = list(self.__organized_data[replication][vehicle].keys())
                min_idx, max_idx = 0, len(timestamps) - 1
                sensor_crossed = False
                
                while min_idx <= max_idx and not sensor_crossed:
                    
                    mid_idx = int(round((max_idx - min_idx)/2) + min_idx)
                    
                    timestamp, next_timestamp = timestamps[mid_idx], timestamps[mid_idx + 1]
                    line = self.__organized_data[replication][vehicle][timestamp]
                    next_line = self.__organized_data[replication][vehicle][next_timestamp]
                    section, section_distance, travelled_distance, lane_idx, speed = line[1], line[2], line[3], line[4], line[5]

                    next_section, next_section_distance, next_travelled_distance, next_speed = next_line[1], next_line[2], next_line[3], next_line[5]
                    
                    if lane_idx not in sensor_location:
                        sensor_location = self.__add_lane_to_sensor_location(sensor_location, lane_idx, s, L, alpha, theta, sensor_section, sensor_section_distance, sections_list)

                    # check if upstream boundary is crossed
                    upstream_cross_before, upstream_cross_after = False, False
                    # take care of vehicles that cross sections that are not on the main line

                    if section in sections_list:
                        upstream_cross_before = ( sections_list.index(section) < sections_list.index(sensor_location[lane_idx]['upstream']['section']) or
                                                  ( section == sensor_location[lane_idx]['upstream']['section'] and
                                                    section_distance + self.__default_parameters['veh_length']/2 < sensor_location[lane_idx]['upstream']['distance']
                                                  )
                                                )  
                        upstream_cross_after = ( sections_list.index(sensor_location[lane_idx]['upstream']['section']) < sections_list.index(next_section) or
                                                 ( next_section == sensor_location[lane_idx]['upstream']['section']  and
                                                   sensor_location[lane_idx]['upstream']['distance'] <= next_section_distance + self.__default_parameters['veh_length']/2
                                                 )
                                                )                         
                    upstream_cross = upstream_cross_before and upstream_cross_after


                    if upstream_cross:
                        distance_to_upstream = sensor_location[lane_idx]['upstream']['distance'] - (section_distance + self.__default_parameters['veh_length']/2)
                        upstream_section_idx = sections_list.index(sensor_location[lane_idx]['upstream']['section'])
                        while distance_to_upstream < 0:
                            upstream_section_idx = upstream_section_idx - 1
                            distance_to_upstream = distance_to_upstream + self.__sections[sections_list[upstream_section_idx]]['length']
                        cross_speed_m_s = (next_travelled_distance - travelled_distance)/(next_timestamp - timestamp)
                        enter_timestamp = float('{0:.2f}'.format(timestamp + distance_to_upstream/cross_speed_m_s))
                        enter_speed_km_h = speed

                        # iterate over the timestamps following the upstream cross timestamp to find the timestamp at
                        # the crossing of the downstream boundary
                        for idx in range (mid_idx, len(timestamps)):
                            
                            timestamp, next_timestamp = timestamps[idx], timestamps[idx + 1]
                            line = self.__organized_data[replication][vehicle][timestamp]
                            next_line = self.__organized_data[replication][vehicle][next_timestamp]
                            section, section_distance, travelled_distance, lane_idx, speed = line[1], line[2], line[3], line[4], line[5]
                            next_section, next_section_distance, next_travelled_distance = next_line[1], next_line[2], next_line[3]

                            if lane_idx not in sensor_location:
                                sensor_location = self.__add_lane_to_sensor_location(sensor_location, lane_idx, s, L, alpha, theta, sensor_section, sensor_section_distance, sections_list)

                            # check if downstream boundary is crossed
                            downstream_cross_before, downstream_cross_after = False, False
                            if section in sections_list:
                                downstream_cross_before = ( sections_list.index(section) < sections_list.index(sensor_location[lane_idx]['downstream']['section']) or
                                                            ( section == sensor_location[lane_idx]['downstream']['section'] and
                                                              section_distance + self.__default_parameters['veh_length']/2 < sensor_location[lane_idx]['downstream']['distance']
                                                            )
                                                          )  
                                downstream_cross_after = ( sections_list.index(sensor_location[lane_idx]['downstream']['section']) < sections_list.index(next_section) or
                                                           ( next_section == sensor_location[lane_idx]['downstream']['section'] and
                                                             sensor_location[lane_idx]['downstream']['distance'] < next_section_distance + self.__default_parameters['veh_length']/2
                                                       )
                                                      )                         
                            downstream_cross = downstream_cross_before and downstream_cross_after
                        
                            if downstream_cross:
                                distance_to_downstream = sensor_location[lane_idx]['downstream']['distance'] - (section_distance + self.__default_parameters['veh_length']/2)
                                downstream_section_idx = sections_list.index(sensor_location[lane_idx]['downstream']['section'])
                                while distance_to_downstream < 0:
                                    downstream_section_idx = downstream_section_idx - 1
                                    distance_to_downstream = distance_to_downstream + self.__sections[sections_list[downstream_section_idx]]['length']
                                cross_speed_m_s = (next_travelled_distance - travelled_distance)/(next_timestamp - timestamp)
                                exit_timestamp = float('{0:.2f}'.format(timestamp + distance_to_upstream/cross_speed_m_s))
                                exit_speed_km_h = speed
                                
                                cross_lane_idx = lane_idx
                                recorded_speed_km_h = cross_speed_m_s*3600/1000
                                #recorded_speed_km_h = (enter_speed_km_h + exit_speed_km_h)/2
                                #v_hat.append(recorded_speed_km_h)
                                #v_tilda.append(cross_speed_m_s*3600/1000)
                                #vehicles.append(vehicle)

                                cross_data = [enter_timestamp, exit_timestamp, vehicle, cross_lane_idx, recorded_speed_km_h]
                                sensor_crosses = np.vstack((sensor_crosses, cross_data))
                                sensor_crossed = True
                                break
                                                                   
                    elif min_idx == max_idx:
                        break
                    elif upstream_cross_before:
                        min_idx = mid_idx + 1
                    else:
                        max_idx = mid_idx - 1
                            
            except IndexError:
                pass

        return sensor_crosses
    

    def __add_lane_to_sensor_location(self, sensor_location, lane_idx, s, L, alpha, theta,
                                      sensor_section, sensor_section_distance, sections_list):
     
        sensor_location[lane_idx] = OrderedDict()
        d = (s + lane_idx*L - L/2)*np.tan(alpha)
        b = (s + lane_idx*L - L/2)*(np.tan(alpha + theta) - np.tan(alpha))
        a = (s + lane_idx*L - L/2)*(np.tan(alpha) - np.tan(alpha - theta))

        # add upstream crossing location
        sensor_location[lane_idx]['upstream'] = OrderedDict()
        upstream_section = sensor_section
        upstream_section_distance = sensor_section_distance - d - b
        while upstream_section_distance < 0:
            upstream_section_idx = sections_list.index(upstream_section)
            upstream_section = sections_list[upstream_section_idx - 1]
            upstream_section_distance = upstream_section_distance + self.__sections[upstream_section]['length']
        sensor_location[lane_idx]['upstream']['section'] = upstream_section
        sensor_location[lane_idx]['upstream']['distance'] = upstream_section_distance

        # add downstream crossing location
        sensor_location[lane_idx]['downstream'] = OrderedDict()
        downstream_section = sensor_section
        downstream_section_distance = sensor_section_distance - d + a
        while upstream_section_distance < 0:
            downstream_section_idx = sections_list.index(downstream_section)
            downstream_section = sections_list[downstream_section_idx - 1]
            downstream_section_distance = downstream_section_distance + self.__sections[downstream_section]['length']
        downstream_distance = sensor_section_distance - d + a
        sensor_location[lane_idx]['downstream']['section'] = downstream_section
        sensor_location[lane_idx]['downstream']['distance'] = downstream_section_distance

        return sensor_location
    

    def __build_volume_sensor_readings(self, crosses, sensor_parameters):

        # obtain noise parameters
        random.seed()
        noise_type, p_occlusion_accept, occlusion_config = \
            sensor_parameters['default']['noise_type'], \
            sensor_parameters['default']['p_occlusion_accept'], \
            sensor_parameters['default']['occlusion']

        aggregation_sec, awake_sec = sensor_parameters['default']['aggregation_sec'], sensor_parameters['default']['awake_sec']
        v_range, v_threshold = sensor_parameters['default']['v_range'], sensor_parameters['default']['v_threshold']
        p_missing_ff, p_missing_cf = sensor_parameters['default']['p_missing_ff'], sensor_parameters['default']['p_missing_cf']
        if noise_type == 'relative':
            v_bias_ff, v_bias_cf = sensor_parameters['default']['v_bias_ff'], sensor_parameters['default']['v_bias_cf']
            v_accuracy_p_ff, v_accuracy_p_cf = sensor_parameters['default']['v_accuracy_p_ff'], sensor_parameters['default']['v_accuracy_p_cf']
        elif noise_type == 'absolute':
            v_noise_sigma_ff, v_noise_sigma_cf = sensor_parameters['default']['v_noise_sigma_ff'], sensor_parameters['default']['v_noise_sigma_cf']
            v_noise_mu_ff, v_noise_mu_cf = sensor_parameters['default']['v_noise_mu_ff'], sensor_parameters['default']['v_noise_mu_cf']


        curr_r_occlusion_accept = random.uniform(p_occlusion_accept[0], p_occlusion_accept[1])/100
        curr_r_missing_ff, curr_r_missing_cf = random.uniform(p_missing_ff[0], p_missing_ff[1])/100, random.uniform(p_missing_cf[0], p_missing_cf[1])/100

        if noise_type == 'relative':
            curr_v_bias_ff, curr_v_bias_cf = random.uniform(v_bias_ff[0], v_bias_ff[1]), random.uniform(v_bias_cf[0], v_bias_cf[1])
            curr_v_accuracy_r_ff, curr_v_accuracy_r_cf = random.uniform(v_accuracy_p_ff[0], v_accuracy_p_ff[1])/100, random.uniform(v_accuracy_p_cf[0], v_accuracy_p_cf[1])/100
        if noise_type == 'absolute':
            curr_v_noise_sigma_ff, curr_v_noise_sigma_cf = random.uniform(v_noise_sigma_ff[0], v_noise_sigma_ff[1]), random.uniform(v_noise_sigma_cf[0], v_noise_sigma_cf[1])
            curr_v_noise_mu_ff, curr_v_noise_mu_cf = random.uniform(v_noise_mu_ff[0], v_noise_mu_ff[1]), random.uniform(v_noise_mu_cf[0], v_noise_mu_cf[1])

        sensor = OrderedDict()
        intervals = aggregation_sec*(np.arange(math.floor(self.__start_timestamp/aggregation_sec),math.floor(self.__end_timestamp/aggregation_sec + 1)) -
                                     math.floor(self.__start_timestamp/aggregation_sec))

        for interval in intervals:
            sensor[interval] = OrderedDict()
            sensor[interval]['crossed?'] = False
            sensor[interval]['inverse_speed_sum'] = 0
            sensor[interval]['count'] = 0
            sensor[interval]['flow'] = 0
            sensor[interval]['sms'] = float('inf')
            sensor[interval]['density'] = 0

        if crosses.any():
            max_lane_idx = int(np.max(crosses[:,3]))
            # iterate over every lane_idx available
            for lane_idx in range (1, max_lane_idx + 1):
                for veh_cross in crosses[crosses[:,3] == float(lane_idx), :]:
                    # the interval is set based on the enter_timestamp
                    enter_timestamp = veh_cross[0]
                    if (enter_timestamp % aggregation_sec) <= awake_sec:
                        # check for occlusion
                        clear_vision = True
                        if lane_idx != 1 and curr_r_occlusion_accept < 1:
                            clear_vision = self.__occlusion_check(crosses, veh_cross, curr_r_occlusion_accept, occlusion_config)
                        if clear_vision:
                            interval = aggregation_sec*(math.floor(enter_timestamp/aggregation_sec) - math.floor(self.__start_timestamp/aggregation_sec))
                            if not sensor[interval]['crossed?']:
                                sensor[interval]['crossed?'] = True
                            speed = veh_cross[4]
                            # add noise
                            if noise_type == 'relative':
                                if speed >= v_threshold:
                                    r_missing, v_bias, v_accuracy_r = curr_r_missing_ff, curr_v_bias_ff, curr_v_accuracy_r_ff
                                else:
                                    r_missing, v_bias, v_accuracy_r = curr_r_missing_cf, curr_v_bias_cf, curr_v_accuracy_r_cf
                                noisy_speed = abs(random.gauss(speed + v_bias,(1 - v_accuracy_r)*speed))
                            elif noise_type == 'absolute':
                                if speed >= v_threshold:
                                    r_missing, v_noise_mu, v_noise_sigma = curr_r_missing_ff, curr_v_noise_mu_ff, curr_v_noise_sigma_ff
                                else:
                                    r_missing, v_noise_mu, v_noise_sigma = curr_r_missing_cf, curr_v_noise_mu_cf, curr_v_noise_sigma_cf
                                noisy_speed = abs(speed + random.gauss(v_noise_mu,v_noise_sigma))
                            retain_data_check = random.uniform(0, 1)

                            # update data
                            if v_range[0] <= noisy_speed and noisy_speed <= v_range[1] and retain_data_check >= r_missing:
                                if noisy_speed:
                                    sensor[interval]['inverse_speed_sum'] += 1/noisy_speed
                                sensor[interval]['count'] += 1
                                sensor[interval]['flow'] = int(sensor[interval]['count']*3600/aggregation_sec)
                                if sensor[interval]['inverse_speed_sum']:
                                    # print( sensor[interval]['inverse_speed_sum'] )
                                    sensor[interval]['sms'] = float('{0:.10f}'.format(sensor[interval]['count']/sensor[interval]['inverse_speed_sum']))
                                    print( sensor[interval]['sms'] )
                                    sensor[interval]['density'] = float('{0:.2f}'.format(sensor[interval]['flow']/sensor[interval]['sms']))


        # tell difference between zero count and no veh read because of drop rate
        for interval in sensor:
            if sensor[interval]['crossed?'] and not sensor[interval]['count']:
                sensor[interval]['count'] = -1
                sensor[interval]['flow'] = -1
                sensor[interval]['density'] = -1
                sensor[interval]['sms'] = -1

        # drop the data based on the data missing parameter
        for interval in sensor:
            # if freeflow
            if sensor[interval]['sms'] != -1 and sensor[interval]['sms'] >= v_threshold:
                r_missing = curr_r_missing_ff
            # if congested flow
            elif sensor[interval]['sms'] !=-1 and sensor[interval]['sms'] < v_threshold:
                r_missing = curr_r_missing_cf

            retain_data_check = random.uniform(0,1)
            if retain_data_check <= r_missing:
                # throw away data
                sensor[interval]['sms'] = -1
                sensor[interval]['count'] = -1
                sensor[interval]['flow'] = -1
                sensor[interval]['density'] = -1

        return sensor


    def __occlusion_check(self, crosses, veh_cross, curr_r_occlusion_accept, occlusion_config):

        clear_vision = True

        if occlusion_config is True:

            lane_idx = int(veh_cross[3])
            veh_cross_time = veh_cross[1] - veh_cross[0]
            occlusion_acceptance = curr_r_occlusion_accept*veh_cross_time
            # Check if occlusion occurs on any of the inner lanes
            for inner_lane_idx in range(1, lane_idx):
                # Update the clear_vision boolean by checking that there is no overlap on any of the lanes. This happens when:
                # current vehicle exit_time < enter_time of any vehicle with a lower lane_idx (e2 < s1) OR
                # current vehicle enter_time > exit_time of any vehicle with a lower lane_index (e1 < s2)
                clear_vision = clear_vision and np.logical_or(veh_cross[1] < crosses[crosses[:,3] == float(inner_lane_idx),0] - occlusion_acceptance,
                                                              veh_cross[0] > crosses[crosses[:,3] == float(inner_lane_idx),1] + occlusion_acceptance
                                                              ).all()
        return clear_vision


    def __build_travel_time_sensor(self, sensor_parameters, save2file=True):

        up_sensor_parameters, down_sensor_parameters = sensor_parameters.copy(), sensor_parameters.copy()
        up_sensor_parameters['section'], down_sensor_parameters['section'] = sensor_parameters['section_up'], sensor_parameters['section_down']
        up_sensor_parameters['distance'], down_sensor_parameters['distance'] = sensor_parameters['distance_up'], sensor_parameters['distance_down']

        all_pair_readings = OrderedDict()

        for replication in self.__organized_data:

            # build crosses at upstream and downstream sensors
            sensor_up_crosses = self.__build_sensor_crosses(up_sensor_parameters, replication)
            sensor_down_crosses = self.__build_sensor_crosses(down_sensor_parameters, replication)

            # build the sensor readings based on the sensor crossings
            pair_readings = self.__build_travel_time_sensor_readings(sensor_up_crosses, sensor_down_crosses, sensor_parameters, replication)

            if save2file is True:

                # sensor_file = self.__virtual_sensors_data_folder + str(replication) + '\\' + str(sensor_parameters['id']) + '.txt'
                sensor_file = self.__virtual_sensors_data_folder[replication] + str(sensor_parameters['id']) + '.txt'

                if not os.path.exists( os.path.dirname(sensor_file) ):
                    os.makedirs(os.path.dirname(sensor_file))
                #os.makedirs(os.path.dirname(sensor_file), exist_ok=True)


                with open(sensor_file,'w+') as file:
                    for interval in pair_readings:
                        # if travel time measurement exist
                        if 'travel_time' in pair_readings[interval].keys() and not np.isnan(pair_readings[interval]['travel_time']):
                            value = pair_readings[interval]['travel_time']
                        else:
                            value = np.nan
                        line = ','.join( [ str(i) for i in [ interval, value ] ])
                        file.write(line + '\n')

                with open(self.__virtual_sensors_generated_file[replication], 'a+') as file:
                    line = self.__sensordict2string( str( sensor_parameters['id'] ), sensor_parameters )
                    file.write(line + '\n')

            else:
                travel_time = []
                latest_reading = np.nan

                for interval in range(0, 9000, 5):

                    if interval in pair_readings.keys() and 'travel_time' in pair_readings[interval].keys() and not np.isnan(pair_readings[interval]['travel_time']):

                        travel_time.append( pair_readings[interval]['travel_time'] )
                        # update the latest reading
                        latest_reading = pair_readings[interval]['travel_time']
                    else:
                        travel_time.append( latest_reading )


                all_pair_readings[replication] = travel_time

        return all_pair_readings


    def __build_travel_time_sensor_readings(self, sensor_up_crosses, sensor_down_crosses, sensor_parameters, replication):

        # obtain noise parameters
        random.seed()
        aggregation_sec = sensor_parameters['default']['aggregation_sec']
        p_penetration = sensor_parameters['default']['p_penetration']
        t_noise_sigma, t_noise_mu = sensor_parameters['default']['t_noise_mu'], sensor_parameters['default']['t_noise_mu']
        t_sigma, t_mu = random.uniform(t_noise_sigma[0],t_noise_sigma[1]), random.uniform(t_noise_mu[0],t_noise_mu[1])
        curr_r_penetration = random.uniform(p_penetration[0],p_penetration[1])/100
        
        pair = OrderedDict()
        # add the travel time of each penetrated vehicle
        for vehicle in self.__organized_data[replication]:
            if random.uniform(0,1) <= curr_r_penetration:
                try:
                    up_timestamp = sensor_up_crosses[sensor_up_crosses[:,2] == vehicle,0]
                    down_timestamp = sensor_down_crosses[sensor_down_crosses[:,2] == vehicle,0]
                    noisy_up_timestamp = abs(up_timestamp + random.gauss(t_mu, t_sigma))
                    noisy_down_timestamp = abs(down_timestamp + random.gauss(t_mu, t_sigma))
                    interval = aggregation_sec*(math.floor(up_timestamp/aggregation_sec) - math.floor(self.__start_timestamp/aggregation_sec))
                    if interval not in pair:
                        pair[interval] = OrderedDict()
                        pair[interval]['up_times_list'] = []
                        pair[interval]['down_times_list'] = []
                    if up_timestamp and down_timestamp:
                        pair[interval]['up_times_list'].append(noisy_up_timestamp)
                        pair[interval]['down_times_list'].append(noisy_down_timestamp)
                except:
                    pass
        # compute the means of vehicle travel times for each interval
        for interval in pair:
            # pair[interval]['up_time'] = float('%.2f' % np.mean(pair[interval]['up_times_list']))
            # pair[interval]['down_time'] = float('%.2f' % np.mean(pair[interval]['down_times_list']))

            # return the mean travel time
            if len(pair[interval]['up_times_list']) != 0:
                pair[interval]['travel_time'] = float('%.2f' % (np.mean( np.array(pair[interval]['down_times_list']) -
                                                                         np.array(pair[interval]['up_times_list']))) )
            else:
                # otherwise not travel time recoreded
                pair[interval]['travel_time'] = np.nan
            
        return pair


    def __build_true_flow_density_speed_states(self, resolution, replication, queue_threshold, plot=False):

        """
        This function generates matrices that hold the true traffic state on
        a road network based on Edie's definitions
        :param resolution: Tuple [aggregation in seconds, aggregation distance in meters]
        :param replication: Replication number
        :param plot: Boolean (True if plots are requested, False otherwise)
        :return:
        """

        # extract resolution and construct cell space and time boundaries
        agg_sec = resolution[0]
        agg_dist = resolution[1]
        cell_area = agg_sec*agg_dist
        main_len = 0.0
        section_bounds = [0.0]
        for section in self.__main_grid_order:
            main_len = main_len + self.__sections[section]['length']
            section_bounds.append(float('{0:.2f}'.format(main_len)))

        # get the space cell grids only for the section investigation
        num_pt = round( (self.__space_end - self.__space_start)/agg_dist ) + 1
        space_cell_bounds = np.linspace( self.__space_start, self.__space_end, num_pt )

        # get the time cell grids only for the time domain under investigation
        num_pt = round( (self.__end_timestamp - self.__start_timestamp)/agg_sec ) + 1
        time_cell_bounds = np.linspace( self.__start_timestamp, self.__end_timestamp, num_pt ) - self.__start_timestamp

        # initialize an empty matrix for the space sum and time sum on each cell
        # for n cell boundaries on an axis, there are n-1 cells on that axis
        spaces_sum_mat = np.zeros((len(space_cell_bounds) - 1, (len(time_cell_bounds) - 1)))
        times_sum_mat = np.zeros((len(space_cell_bounds) - 1, (len(time_cell_bounds) - 1)))

        real_space_cell_bounds = space_cell_bounds
        # space_cell_bounds = np.unique(list(space_cell_bounds) + list(section_bounds))

        # iterate over every vehicle
        for vehicle in self.__organized_data[replication]:

            # get timestamps, sections and distances in arrays
            timestamps = np.array(list(self.__organized_data[replication][vehicle].keys()))
            sections = np.array([self.__organized_data[replication][vehicle][timestamp][1] for
                                 timestamp in self.__organized_data[replication][vehicle]])
            distances = np.array([self.__organized_data[replication][vehicle][timestamp][2] for
                                  timestamp in self.__organized_data[replication][vehicle]])

            # get idxs of the sections that match the sections on the main grid
            # adjust section distances to be absolute to the beginning of main grid
            # extract main absolute distances and timestamps
            main_idxs = []
            for section in self.__main_grid_order:
                idxs = np.array(np.where(sections == section)[0]).tolist()
                up_sect_len = self.__get_up_sect_len(section)
                distances[idxs] = distances[idxs] + up_sect_len
                timestamps[idxs] = timestamps[idxs] - self.__start_timestamp
                main_idxs.extend(idxs)
            main_distances = distances[[main_idxs]]
            main_timestamps = timestamps[[main_idxs]]

            # get rid of duplicate distances and timestamps
            dup_idxs = [idx for idx, item in enumerate(main_distances) if item in main_distances[:idx]]
            main_distances = np.delete(main_distances,dup_idxs)
            main_timestamps = np.delete(main_timestamps,dup_idxs)

            print('checking vehicle: {0}'.format(vehicle))
            print('length of vehicle traj {0}'.format( len(distances) ))
            print('length of main_distances: {0}'.format( len(main_distances) ))

            # Make sure the vehicle was on the main freeway
            if len(main_distances) > 1:
                # obtain cell bounds that are relevant to current vehicle
                veh_space_cell_bounds = space_cell_bounds[space_cell_bounds>=main_distances[0]-agg_dist]
                veh_space_cell_bounds = veh_space_cell_bounds[veh_space_cell_bounds<=main_distances[-1]+agg_dist]
                veh_time_cell_bounds = time_cell_bounds[time_cell_bounds>=main_timestamps[0]-agg_sec]
                veh_time_cell_bounds = veh_time_cell_bounds[veh_time_cell_bounds<=main_timestamps[-1]+agg_sec]

                # set a 1st order interpolation/extrapolation for the absoulte distances and timestamps
                order = 1
                # interpolate/extrapolate on space given time
                space_spl = InterpolatedUnivariateSpline(main_timestamps, main_distances, k=order)
                # interpolate/extrapolate on time given space
                time_spl = InterpolatedUnivariateSpline(main_distances, main_timestamps, k=order)
                time_at_space_bounds = self.__interpolate_space_crosses(time_spl,veh_space_cell_bounds)
                space_at_time_bounds = self.__interpolate_time_crosses(space_spl,veh_time_cell_bounds)
                crosses = (np.sort(time_at_space_bounds + space_at_time_bounds, axis=0)).tolist()
                crosses = [x for i, x in enumerate(crosses) if not i or x!=crosses[i-1]]

                # update values on spaces_sum_mat and times_sum_mat
                for i in range(1,len(crosses)-2):
                    ti, xi = crosses[i][0], crosses[i][1]
                    tf, xf = crosses[i+1][0], crosses[i+1][1]
                    if ti >= time_cell_bounds[0] and xi >= real_space_cell_bounds[0] and ti < time_cell_bounds[-1] and xi < real_space_cell_bounds[-1]:
                        delta_t = tf - ti
                        delta_x = xf - xi
                        t_idx, d_idx = int((ti-time_cell_bounds[0])/agg_sec), int((xi-real_space_cell_bounds[0])/agg_dist)
                        times_sum_mat[d_idx,t_idx] = times_sum_mat[d_idx,t_idx] + delta_t
                        spaces_sum_mat[d_idx,t_idx] = spaces_sum_mat[d_idx,t_idx] + delta_x

        # compute flow, density and speed based on edie's definitions
        spaces_sum_mat = np.matrix(spaces_sum_mat)
        times_sum_mat = np.matrix(times_sum_mat)
        q_mat = np.true_divide(spaces_sum_mat, cell_area)   # veh/s
        k_mat = np.true_divide(times_sum_mat, cell_area)    # veh/m

        with np.errstate(divide='ignore', invalid='ignore'):
            v_mat = np.true_divide(spaces_sum_mat, times_sum_mat)   # m/s
            v_mat[v_mat == np.inf] = np.nan


        # save the true states on the appropriate files
        true_density_file = self.__true_state_data_folder[replication] + 'truestate_{0}s{1}m_density.txt'.format(agg_sec, agg_dist)
        true_speed_file = self.__true_state_data_folder[replication] + 'truestate_{0}s{1}m_speed.txt'.format(agg_sec, agg_dist)
        true_queue_file = self.__true_state_data_folder[replication] + 'truestate_{0}s{1}m_queue.txt'.format(agg_sec, agg_dist)
        true_traveltime_file = self.__true_state_data_folder[replication] + 'truestate_{0}s{1}m_traveltime.txt'.format(agg_sec, agg_dist)

        # fill in the nan values, using no update
        dim_space, dim_time = v_mat.shape
        for x in range(0, dim_space):
            for t in range(0, dim_time):
                if np.isnan(v_mat[x,t]):
                    # replace the nan value by not updating
                    j = t-1
                    while j!= -1 and np.isnan(v_mat[x,j]):
                        j -= 1

                    if j != -1:
                        # if found a value that is not nan
                        v_mat[x,t] = v_mat[x,j]


        # get the end of the queue
        queue_mat = np.zeros((1, dim_time))

        for step in range(0, dim_time):
            index = v_mat[:, step] <= queue_threshold

            if sum(index) <= 2: # use 2 to suppress false alarms
                # if less or equal then 2 cells are in congestion, it may be caused by noise.
                queue_mat[0, step] = 0
            else:
                queue_mat[0, step] = \
                    agg_dist*( dim_space - np.argmax(index) )


        # get the travel time, use a perfect sensor
        # ----------- perfect BT ----------- #
        # TODO: this is a quick hack for generating the true travel time for I57
        # sensor_paras = OrderedDict()
        #
        # sensor_paras['id'] = 'BTperfect'
        # sensor_paras['type'] = 'bluetooth'
        # sensor_paras['section_up'] = 51413
        # sensor_paras['distance_up'] = self.__space_start
        # sensor_paras['section_down'] = 51413
        # sensor_paras['distance_down'] = self.__space_end
        #
        # sensor_paras['default'] = {}
        # sensor_paras['default']['alpha'] = 0
        # sensor_paras['default']['theta'] = 0
        # sensor_paras['default']['s'] = 1
        # sensor_paras['default']['L'] = 3.5
        #
        # sensor_paras['default']['aggregation_sec'] = 5
        # sensor_paras['default']['p_penetration'] = [100,100]
        # sensor_paras['default']['t_noise_sigma'] = [0.0,0.0]
        # sensor_paras['default']['t_noise_mu'] = [0.0,0.0]
        #
        # tt_mat_with_time = self.__build_travel_time_sensor(sensor_paras, save2file=False)
        # tt_mat = np.matrix(tt_mat_with_time[51595])

        tt_mat = queue_mat*0.0

        np.savetxt( true_speed_file, v_mat.T, delimiter=',')
        np.savetxt( true_density_file, k_mat.T, delimiter=',')
        np.savetxt( true_queue_file, queue_mat.T, delimiter=',')
        np.savetxt( true_traveltime_file, tt_mat.T, delimiter=',')

        # plot figures
        if plot:

            extend_array = [0, self.__end_timestamp - self.__start_timestamp, 0, main_len]
            
            plt.figure()
            plt.imshow(spaces_sum_mat, cmap='jet', origin='lower', aspect='auto', interpolation='none', extent=extend_array)
            plt.title(self.__work_zone + ' Cell Travelled Distance Map')
            plt.xlabel('Time [s]')
            plt.ylabel('Distance [m]')
            cbar = plt.colorbar()
            cbar.set_label('Travelled Distance [veh*m]')
            
            plt.figure()
            plt.imshow(times_sum_mat, cmap='jet', origin='lower', aspect='auto', interpolation='none', extent=extend_array)
            plt.title(self.__work_zone + ' Cell Travelled Time Map')
            plt.xlabel('Time [s]')
            plt.ylabel('Distance [m]')
            cbar = plt.colorbar()
            cbar.set_label('Travelled Time [veh*hr]')

            plt.figure()
            qplt = plt.imshow(q_mat, cmap='jet', origin='lower', aspect='auto', interpolation='none', extent=extend_array)
            plt.title(self.__work_zone + ' Flow Map')
            plt.xlabel('Time [s]')
            plt.ylabel('Distance [m]')
            cbar = plt.colorbar()
            cbar.set_label('Flow [veh/hr]')

            plt.figure()
            kplt = plt.imshow(k_mat, cmap='jet', origin='lower', aspect='auto', interpolation='none', extent=extend_array)
            plt.title(self.__work_zone + ' Density Map')
            plt.xlabel('Time [s]')
            plt.ylabel('Distance [m]')
            cbar = plt.colorbar()
            cbar.set_label('Density [veh/km]')

            plt.figure()
            vplt = plt.imshow(v_mat, cmap='jet_r', origin='lower', aspect='auto', interpolation='none', extent=extend_array)
            plt.title(self.__work_zone + ' Speed Map')
            plt.xlabel('Time [s]')
            plt.ylabel('Distance [m]')
            cbar = plt.colorbar()
            cbar.set_label('Speed [km/h]')
            
            plt.show()


    def __interpolate_space_crosses(self, time_spl, space_cell_bounds):

        """
        This function interpolates on space given a time spline
        :param time_spl: time spline that takes space
        :param space_cell_bounds: space cell bounds
        :return vector of tuples time_space_cross
        """

        time_space_cross = []
        for space_bound in space_cell_bounds:
            time_space_cross.append([float('{0:.2f}'.format(float(time_spl(space_bound)))),float('{0:.2f}'.format(space_bound))])
        return time_space_cross

    def __interpolate_time_crosses(self, space_spl, time_cell_bounds):

        """
        This function interpolates on time given a space spline
        :param space_spl: space spline that takes time
        :param time_cell_bounds: time cell bounds
        :return vector of tuples time_space_cross
        """

        time_space_cross = []
        for time_bound in time_cell_bounds:
            time_space_cross.append([float('{0:.2f}'.format(time_bound)), float('{0:.2f}'.format(float(space_spl(time_bound))))])
        return time_space_cross


        
    def __get_up_sect_len(self, section):

        """
        This function returns the cumulative lengths of the main grid
        sections that are upstream of the section given
        :param section: section id on the main grid
        :return: float up_sect_len
        """

        up_sect_len = 0.0
        for main_sect in self.__main_grid_order:
            if section == main_sect:
                return up_sect_len
            else:
                up_sect_len = up_sect_len + self.__sections[main_sect]['length']

        raise Exception('Error: The main grid does not contain the section {0}'.format(section))

    @staticmethod
    def __sensordict2string( sensor_id, sensor_att):
        """
        This function converts the sensor key-value store to a string in the following format.
        Entries are separated by ; The entry name and value are separated by :, first entry is id:
        :param sensor_id: the sensor id
        :param sensor_att: key-value store of the sensors.
                general_keys: id, type, section, distance...
                default: Dict, first copied from default parameters and then overwritten by custom values. Use this dict for parameters
                custom: Dict, the set of parameters that have been written to default.
        :return: a string in the above defined format
        """

        line = []
        # append the first id entry
        entry = ':'.join( ['id', sensor_id ] )
        line.append(entry)

        # append other attribute entries
        for key in sensor_att.keys():

            if key != 'id':
                # Do Not repeat id entry

                if key == 'custom':
                    # skip the custom entries since they have been written to default sets.
                    continue
                elif key == 'default':
                    # add all parameters in default dict
                    for para in sensor_att['default']:
                        entry = ':'.join( [ para, str(sensor_att['default'][para]) ] )
                        line.append(entry)
                else:
                    # other keys, such as type, distance...
                    entry = ':'.join( [ key, str(sensor_att[key]) ] )
                    line.append(entry)

        # join all entries and return
        return ';'.join( line )
