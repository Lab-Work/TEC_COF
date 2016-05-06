__author__ = 'Yanning'


import numpy as np
import matplotlib.pyplot as plt

# link_str = 'freeway'
link_str = 'onramp'

data_path_freeflow = '/Users/Yanning/GoogleDrive/AIMSUN_traffic_control_simulation/FD_calibration/det_data/'+ link_str +'_' + link_str + '_FD_freeflow.txt'
data_path_congflow = '/Users/Yanning/GoogleDrive/AIMSUN_traffic_control_simulation/FD_calibration/det_data/'+ link_str +'_' + link_str + '_FD_congflow.txt'

f_free = open(data_path_freeflow)
f_cong = open(data_path_congflow)
# next(f)

speed_free = []
flow_free = []

for line in f_free:

    line = line.strip()
    items = line.split(',')

    # print line

    if len(items[0]) == 0 or len(items[1]) == 0:
        continue

    else:
        speed_free.append(float(items[0]))
        flow_free.append(float(items[1]))

f_free.close()


speed_cong = []
flow_cong = []

for line in f_cong:

    line = line.strip()
    items = line.split(',')

    # print line

    if len(items[0]) == 0 or len(items[1]) == 0:
        continue

    else:
        speed_cong.append(float(items[0]))
        flow_cong.append(float(items[1]))

f_cong.close()

speed_free = np.array(speed_free)
flow_free = np.array(flow_free)
density_free = flow_free/speed_free

speed_cong = np.array(speed_cong)
flow_cong = np.array(flow_cong)
density_cong = flow_cong/speed_cong

# fig = plt.figure(figsize=(15,10))
# plt.plot(speed)

# ======================================
# linear fitting of the freeflow
# Our model is y = a * x, so things are quite simple, in this case...
# x needs to be a column vector instead of a 1D vector for this, however.
density_free = density_free[:,np.newaxis]
vf, _, _, _ = np.linalg.lstsq(density_free, flow_free)


# ======================================
# The cong side fitting
# force the line to go through a point
q_max = 2100.0
rho_c = q_max/vf[0]

# print q_max, rho_c

# remove the few free flow points to avoid bias
cong_index = (density_cong >= 50)

# coe = np.polyfit(density_cong[cong_index], flow_cong[cong_index], 1)

density_shifted = density_cong[cong_index] - rho_c
flow_shifted = flow_cong[cong_index] - q_max

density_shifted = density_shifted[:, np.newaxis]
w, _, _, _ = np.linalg.lstsq( density_shifted, flow_shifted  )

# ======================================
# compute the parameters
# w = coe[0]
# rho_c = coe[1]/(vf-w)
# q_max = vf*rho_c
# rho_m = -coe[1]/w
rho_m = rho_c - q_max/w[0]

# print 'The congested line is q = {0}*rho + {1}'.format(coe[0], coe[1])


# visualize the result
fig_window = plt.figure(figsize=(15,10))
fig = fig_window.add_subplot(111)
# free flow fitting
# plt.scatter(density_free, flow_free)
plt.plot(density_free, vf*density_free, 'r-', linewidth=2.0)

# congested fitting
plt.scatter(density_cong[cong_index], flow_cong[cong_index], color='b', marker='o')
plt.scatter(density_cong[~cong_index], flow_cong[~cong_index], color='r', marker='*')

#
dens = np.linspace(0, 400, 100)
plt.plot(dens, w[0]*(dens - rho_m), 'r-', linewidth=2.0)

plt.title('Fundamental diagram for single lane ' + link_str, fontsize=24)
plt.xlabel('Traffic density (veh/mile)', fontsize=24)
plt.ylabel('Traffic flow (veh/hr)', fontsize=24)
# plt.xlim([ np.min(x) - 0.1*(np.max(x)-np.min(x)) , np.max(x) + 0.1*(np.max(x)-np.min(x))])
plt.xlim([0, 300])
plt.ylim([0, 2800])
plt.grid(True)
fig.tick_params(axis='both', which='major', labelsize=20)



fig.text(200, 2000, 'vf:   {0} mph\n rho_c:  {1} veh/mile\n rho_m:  {2} veh/mile;\n q_max:  {3} veh/hr\n w:  {4} mph'.
         format(np.round(vf[0], 2),
                np.round(rho_c, 2),
                np.round(rho_m, 2),
                np.round(q_max, 2),
                np.round(w[0],2)),
         style='italic',
         bbox={'facecolor':'red', 'alpha':0.5, 'pad':10},
         fontsize=16)

fig_window = plt.figure(figsize=(15, 10))
fig = fig_window.add_subplot(111)

plt.scatter( density_cong, np.array(speed_cong) )


plt.draw()





plt.show()









