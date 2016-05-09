# Juan Carlos Martinez
# Last Updated: 2/22/2016
# About the file:
#   - This file sets the 'default' parameters for the virtual sensors
#     Ultimately, a dictionary called default_parameters is constructed
# Notes:
#   - This file is meant to be modified by the user to adjust to his/her needs
#   - The parameters of a single sensor can be deviated from the 'default'
#     paramters by adjusting the sensor configuration file

# Code begins here

from collections import OrderedDict
default_parameters = OrderedDict()

default_parameters['veh_length'] = 4.5
default_parameters['time_step'] = 0.8


# -- Volume Sensor Configuration -- #
# The following comments specify and explain the default parameters that need
# to be set for construction a 'volume' sensor

## Deployment Parameters ##
# alpha                      : Angle (CW) between the sensor beam and the road cross-section [deg]
# theta                      : Angle (CW) of visibility on each side of sensor beam [deg]
# s                          : Road shoulder [m]
# L                          : Lane width [m]

## General Noise Parameters ##
# noise_type                 : Noise type ['relative' or 'absolute']
# occlusion                  : Occlusion [True or False]
# p_occlusion_accept         : Percent range of a vehicle's crossing time that may overlap with another vehicle's
#                              crossing time without causing occlusion
# aggregation_sec            : Aggregation period [s]
# awake_sec                  : Detection phase [s]
# v_range                    : Speed range [kph]
# v_threshold                : Speed threshhold for free flow vs. congested flow [kph]
# p_missing_ff               : Percent missing range in free flow []
# p_missing_cf               : Percent missing range in congested flow []

## Relative Noise Parameters ##
# v_bias_ff                  : Speed bias range in free flow [kph]
# v_bias_cf                  : Speed bias range in congested flow [kph]
# v_accuracy_ff              : Speed accuracy percent range in free flow []
# v_accuracy_cf              : Speed accuracy percent range in congested flow []

## Absolute Noise Parameters ##
# v_noise_sigma_ff           : Speed noise std. range in free flow [kph]
# v_noise_sigma_cf           : Speed noise std. range in congested flow [kph]
# v_noise_mu_ff              : Speed noise mean range in free flow [kph]
# v_noise_mu_cf              : Speed noise mean range in congested flow [kph]


# ------ Volume Ground Truth ------ #
default_parameters['vol_gt'] = OrderedDict()

default_parameters['vol_gt']['alpha'] = 0
default_parameters['vol_gt']['theta'] = 0
default_parameters['vol_gt']['s'] = 1  
default_parameters['vol_gt']['L'] = 3.5        

default_parameters['vol_gt']['noise_type'] = 'absolute'
default_parameters['vol_gt']['p_occlusion_accept'] = [100,100]
default_parameters['vol_gt']['aggregation_sec'] = 30
default_parameters['vol_gt']['awake_sec'] = 30
default_parameters['vol_gt']['v_range'] = [0,105]
default_parameters['vol_gt']['v_threshold'] = 0
default_parameters['vol_gt']['p_missing_ff'] = [0.0,0.0]
default_parameters['vol_gt']['p_missing_cf'] = [0.0,0.0]

default_parameters['vol_gt']['v_bias_ff']  = [0,0]
default_parameters['vol_gt']['v_bias_cf'] = [0,0]
default_parameters['vol_gt']['v_accuracy_p_ff'] = [0,0]
default_parameters['vol_gt']['v_accuracy_p_cf'] = [0,0]

default_parameters['vol_gt']['v_noise_sigma_ff'] = [0,0]
default_parameters['vol_gt']['v_noise_sigma_cf'] = [0,0]
default_parameters['vol_gt']['v_noise_mu_ff'] = [0,0]
default_parameters['vol_gt']['v_noise_mu_cf'] = [0,0]
# --------------------------------- #


# ------------- icone ------------- #
default_parameters['icone'] = OrderedDict()

default_parameters['icone']['alpha'] = 30
default_parameters['icone']['theta'] = 5
default_parameters['icone']['s'] = 1
default_parameters['icone']['L'] = 3.5

default_parameters['icone']['noise_type'] = 'relative'
default_parameters['icone']['p_occlusion_accept'] = [30,30]
default_parameters['icone']['aggregation_sec'] = 30
default_parameters['icone']['awake_sec'] = 30
default_parameters['icone']['v_range'] = [0,105]
default_parameters['icone']['v_threshold'] = 10
default_parameters['icone']['p_missing_ff'] = [0.0,0.0]
default_parameters['icone']['p_missing_cf'] = [0.0,0.0]

default_parameters['icone']['v_bias_ff']  = [1,5]
default_parameters['icone']['v_bias_cf'] = [1,5]
default_parameters['icone']['v_accuracy_p_ff'] = [10,10]
default_parameters['icone']['v_accuracy_p_cf'] = [20,20]

default_parameters['icone']['v_noise_sigma_ff'] = [0,50]
default_parameters['icone']['v_noise_sigma_cf'] = [2,3.6]
default_parameters['icone']['v_noise_mu_ff'] = [0,0]
default_parameters['icone']['v_noise_mu_cf'] = [2,3]
# --------------------------------- #


# ------------- radar ------------- #
default_parameters['radar'] = OrderedDict()

default_parameters['radar']['alpha'] = 30
default_parameters['radar']['theta'] = 5
default_parameters['radar']['s']= 1
default_parameters['radar']['L'] = 4

default_parameters['radar']['noise_type'] = 'relative'
default_parameters['radar']['p_occlusion_accept'] = [30,30]
default_parameters['radar']['aggregation_sec'] = 30
default_parameters['radar']['awake_sec'] = 30
default_parameters['radar']['v_range'] = [2,105]
default_parameters['radar']['v_threshold'] = 40
default_parameters['radar']['p_missing_ff'] = [0.2,1.0]
default_parameters['radar']['p_missing_cf'] = [1.0,20]

default_parameters['radar']['v_bias_ff']  = [1,5]
default_parameters['radar']['v_bias_cf'] = [1,5]
default_parameters['radar']['v_accuracy_p_ff'] = [10,10]
default_parameters['radar']['v_accuracy_p_cf'] = [20,20]

default_parameters['radar']['v_noise_sigma_ff'] = [1,1]
default_parameters['radar']['v_noise_sigma_cf'] = [1,1]
default_parameters['radar']['v_noise_mu_ff'] = [0,0]
default_parameters['radar']['v_noise_mu_cf'] = [0,0]
# --------------------------------- #


# ------------- rtms -------------- #
default_parameters['rtms'] = OrderedDict()

default_parameters['rtms']['alpha'] = 30
default_parameters['rtms']['theta'] = 5
default_parameters['rtms']['s']= 1
default_parameters['rtms']['L'] = 4

default_parameters['rtms']['noise_type'] = 'relative'
default_parameters['rtms']['p_occlusion_accept'] = [30,30]
default_parameters['rtms']['aggregation_sec'] = 30
default_parameters['rtms']['awake_sec'] = 30
default_parameters['rtms']['v_range'] = [0,100]
default_parameters['rtms']['v_threshold'] = 40
default_parameters['rtms']['p_missing_ff'] = [0.2,1.0]
default_parameters['rtms']['p_missing_cf'] = [1.0,20]

default_parameters['rtms']['v_bias_ff']  = [1,5]
default_parameters['rtms']['v_bias_cf'] = [1,5]
default_parameters['rtms']['v_accuracy_p_ff'] = [10,10]
default_parameters['rtms']['v_accuracy_p_cf'] = [20,20]

default_parameters['rtms']['v_noise_sigma_ff'] = [1,1]
default_parameters['rtms']['v_noise_sigma_cf'] = [1,1]
default_parameters['rtms']['v_noise_mu_ff'] = [0,0]
default_parameters['rtms']['v_noise_mu_cf'] = [0,0]
# --------------------------------- #



#  Travel Time Sensor Configuration #
# The following comments specify and explain the default parameters that need
# to be set for construction a 'travel time' sensor

## Deployment Parameters ##
# alpha                      : Angle (CW) between the sensor beam and the road cross-section [deg] (Set to 0 [Zero])
# theta                      : Angle (CW) of visibility on each side of sensor beam [deg] (Set to 0 [Zero])
# s                          : Road shoulder [m]
# L                          : Lane width [m]

## General Noise Parameters ##
# aggregation_sec            : Aggregation period [s]
# veh_length                 : Average vehicle length [m]
# p_penetration              : Penetration percentage range []
# t_noise_sigma              : Time recording noise sigma range [s]
# t_noise_mu                 : Time recording noise mean range [s]


# --- Travel Time Ground Truth ---- #
default_parameters['tt_gt'] = OrderedDict()

default_parameters['tt_gt']['alpha'] = 0
default_parameters['tt_gt']['theta'] = 0
default_parameters['tt_gt']['s'] = 1  
default_parameters['tt_gt']['L'] = 3.5        

default_parameters['tt_gt']['aggregation_sec'] = 30
default_parameters['tt_gt']['p_penetration'] = [100,100]
default_parameters['tt_gt']['t_noise_sigma'] = [0,0]
default_parameters['tt_gt']['t_noise_mu'] = [0,0]
# --------------------------------- #

# ----------- bluetooth ----------- #
default_parameters['bluetooth'] = OrderedDict()

default_parameters['bluetooth']['alpha'] = 0
default_parameters['bluetooth']['theta'] = 0
default_parameters['bluetooth']['s'] = 1  
default_parameters['bluetooth']['L'] = 3.5 

default_parameters['bluetooth']['aggregation_sec'] = 30
default_parameters['bluetooth']['p_penetration'] = [2,20]
default_parameters['bluetooth']['t_noise_sigma'] = [1,1]
default_parameters['bluetooth']['t_noise_mu'] = [0,0]
# --------------------------------- #
