__author__ = 'Yanning Li'
"""
This is the script for visualizing the simulated traffic states using trajectory data generated from AIMSUN.
"""

from Virtual_Sensors import Virtual_Sensors
from collections import OrderedDict
import matplotlib.pyplot as plt
from os.path import exists
import numpy as np
import sys

# =====================================================================================================
# Configure the folders and directories
# Uncontrolled: I01 (freeway), I02 (onramp)
# optimal control: I11 (freeway), I12 (onramp)

workzone = 'I02'

# 'C:\\Users\\Carlos\\Documents\\GitHub\\IDOT_Code\\I57_configurations_input.txt'
log_dir = '/Users/Yanning/Dropbox/TEC_COF/examples/onramp_metering_control/AIMSUN_network/Visualization/'
data_dir = '/Users/Yanning/Dropbox/TEC_COF/examples/onramp_metering_control/AIMSUN_network/Visualization/Data/'

queue_threshold = 17.88 # speed in m/s for specifying the queue.
grid_res = (5,10)

unit = 'metric' # metric to be used for visualization; m/s, or mph
limit = [0,35]  # limit used for visualization,

# =====================================================================================================
# No modification after this.
# =====================================================================================================

def load_workzone(topo_file, grid_res):
    """
    This function loads the topology of the network
    :param topo_file: the file path
    :param grid_res: (seconds, meters) tuple for the resolution of the visualization
    :return: t_grid, x_grid, topology of the network including the location of ramps, freeway sec order and replications
    """
    start_dur_step = None
    t_grid = None
    x_grid = None

    workzone_topo = { 'sections':OrderedDict(),
                      'loc_onramp':None,
                      'loc_offramp':None,
                      'fwy_sec_order':None,
                      'replications':None}

    f_topo = open(topo_file,'r')

    for line in f_topo:

        if line[0] == '#':
            continue

        if 'sec_id' in line:
            items = line.strip().split(';')
            sec_id = int( items[0].split(':')[1] )
            workzone_topo['sections'][sec_id] = {}

            # add length
            workzone_topo['sections'][sec_id][ items[1].split(':')[0] ] = float( items[1].split(':')[1] )

            if 'upstream' in line:
                workzone_topo['sections'][sec_id]['connections'] = {}
                for i in range(2,len(items)):
                    entry = items[i].split(':')
                    if entry[0] == 'upstream':
                        workzone_topo['sections'][sec_id]['connections']['upstream'] = [int(j) for j in entry[1].split(',')]
                    elif entry[0] == 'downstream':
                        workzone_topo['sections'][sec_id]['connections']['downstream'] = [int(j) for j in entry[1].split(',')]

        elif 'fwy_sec_order' in line:
            items = line.strip().split(':')
            workzone_topo['fwy_sec_order'] = [int(i) for i in items[1].split(',')]

        elif 'loc_start_end' in line:
            items = line.strip().split(':')
            start = float( items[1].split(',')[0] )
            end = float( items[1].split(',')[1] )

            # get x_grid
            x_grid = np.arange(end, start, - grid_res[1])
            if x_grid[-1] != start:
                x_grid = np.concatenate( [ x_grid, np.array([start]) ] )
            x_grid = x_grid[::-1]

            # round every digit to two decimals
            for i in range(0, len(x_grid)):
                x_grid[i] = round(x_grid[i], 2)

        elif 'time_start_dur_step' in line:
            start_dur_step = [float(i) for i in line.split(':')[1].split(',')]
            t_grid = np.arange(0, start_dur_step[1], grid_res[0] )
            if t_grid[-1] != start_dur_step[1]:
                t_grid = np.concatenate( [t_grid, np.array([start_dur_step[1]])] )

            # round every digit to two decimals
            for i in range(0, len(t_grid)):
                t_grid[i] = round(t_grid[i], 2)

        elif 'replications' in line:
            items = line.strip().split(':')
            workzone_topo['replications'] = [ int(i) for i in items[1].split(',') ]

    f_topo.close()

    # the following entries must be set in the topology file
    if t_grid is None or x_grid is None or workzone_topo['replications'] is None:
        raise Exception('Error: loc_start_end, simulation_duration, replications must be set in workzone topology.')
    else:
        t_grid = t_grid.tolist()
        x_grid = x_grid.tolist()

    return t_grid, x_grid, workzone_topo, start_dur_step


def plot_true_speed_for_rep(grid_res, rep, unit='metric', limit=[0,40]):
    """
    This function plots the true speed profile in the specified unit
    :param unit: 'metric', 'imperial'; respectively 'm, s, m/s', and 'mile, hour, mph'
    :param limit: The limit of the colorbar in above units
    :return: A figure profile with x-axis being the time, and y-axis being the space. (flow direction upwards)
    """

    # set the result dir to read from depending on the platform
    if sys.platform == 'darwin':
        result_dir = data_dir + '/Result/{0}/rep{1}/'.format(workzone, rep)
    elif sys.platform == 'win32':
        result_dir = data_dir + 'Result\\{0}\\rep{1}\\'.format(workzone, rep)
    else:
        raise Exception('Check platform.')

    # read the result from file
    true_speed_file = result_dir + 'truestate_{0}s{1}m_speed.txt'.format(grid_res[0], grid_res[1])

    speed_data = np.genfromtxt(true_speed_file, delimiter=',')
    speed_data = np.matrix(speed_data).T

    if unit == 'metric':
        # all internal values are in metric, so plot directly
        speed = np.flipud( speed_data )
        unit_str = 'm/s'
    elif unit == 'imperial':
        speed = np.flipud( __metric2imperial(speed_data, 'speed') )
        # limit = self.__metric2imperial(np.array(limit), 'speed')
        unit_str = 'mph'
    else:
        raise Exception('Error: Unrecognized unit for plotting speed.')

    fig = plt.figure( figsize=(18,8), dpi=100 )
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    im = ax.imshow(speed, cmap=plt.get_cmap('jet_r'),
                    interpolation='nearest',
                    aspect='auto',
                    vmin=limit[0], vmax=limit[1])
    ax.autoscale(False)


    ax.set_title('True speed'.format(rep, grid_res[0], grid_res[1]))
    plt.xlabel('Time (min)')
    plt.ylabel('Space, traffic direction $\mathbf{\Rightarrow}$')
    cax = fig.add_axes([0.95, 0.25, 0.01, 0.5])
    fig.colorbar(im, cax=cax, orientation='vertical')

    plt.draw()


def __metric2imperial(value = np.zeros((1,1)) , option='speed'):
    """
    A utility function which converts the metric (m, s, m/s) to imperial (mile, hour, m/h)
    :param value: float, np.array, or np.matrix. the to be converted value
    :param option: 'speed', 'density'
    :return: converted value
    """

    if type(value) is float or type(value) is np.float64 \
            or type(value) is np.ndarray or type(value) is np.matrix:
        if option == 'speed':
            return value*3600.0/1609.34
        elif option == 'density':
            return value*1609.34
        elif option == 'distance':
            return value/1609.34
        else:
            raise Exception('Error: Unrecognized unit conversion option.')
    else:
        print(type(value))
        raise Exception('Error: Unrecognized value type for unit conversion.')


def main(argv):

    # load the work zone configuration.
    time_grid, space_grid, workzone_topo, aimsun_start_dur_step = \
        load_workzone( log_dir + workzone + '_topology.txt', grid_res)

    # Generate a Virtual_Sensor class for processing trajectory
    data_generator = Virtual_Sensors(workzone, log_dir, data_dir,
                                     workzone_topo['sections'],
                                     workzone_topo['fwy_sec_order'],
                                     workzone_topo['replications'],
                                     aimsun_start_dur_step,
                                     [space_grid[0], space_grid[-1]])

    for rep in workzone_topo['replications']:

        # find the result dir to save to, depending on the platform
        if sys.platform == 'darwin':
            result_dir = data_dir + '/Result/{0}/rep{1}/'.format(workzone, rep)
        elif sys.platform == 'win32':
            result_dir = data_dir + 'Result\\{0}\\rep{1}\\'.format(workzone, rep)
        else:
            raise Exception('Check platform.')
        true_speed_file = result_dir + 'truestate_{0}s{1}m_speed.txt'.format(grid_res[0], grid_res[1])

        # if the true state data has been generated before, then skip
        if not exists( true_speed_file ):
            # generate the true state data
            print('Status: Generating true states for Replication {0} at grid {1}s x {2}m...'.format(rep,
                                                                                                     grid_res[0],
                                                                                                     grid_res[1]))
            data_generator.generate_true_states_data(grid_res, rep, queue_threshold)
        else:
            print('Status: True state for Replication {0} at grid {1}s x {2}m has been previously generated.'.format(rep, grid_res[0], grid_res[1]))

        # plot the generated true speed
        plot_true_speed_for_rep(grid_res, rep, unit=unit, limit=limit)

    plt.show()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

