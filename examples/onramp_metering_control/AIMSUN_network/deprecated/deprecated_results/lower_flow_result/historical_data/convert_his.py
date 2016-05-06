# This file is to convert the historical data to the data format we defined
# from ISO time, flow means (flow dev) to 
# link_id, upstream, time, flow, speed


us_link = 329
ds_link = 330
on_link = 390

time_interval = 60

data_to_write = []

file_list = ['his_us_in', 'his_ds_out', 'his_on_in']


for file in file_list:

	time_stamp = 0

	name = file.split('_')

	if name[1] == 'us':
		link = 329
	elif name[1] == 'ds':
		link = 330
	elif name[1] == 'on':
		link = 390

	if name[2] == 'in':
		bound = 'upstream'
	elif name[2] == 'out':
		bound = 'downstream'

	f = open(file + '.txt', 'r')

	next(f)

	for line in f:
		items = line.split('\t')
		flow = float(items[1].split(' ')[0]) # veh/hr

		time_stamp = time_stamp + time_interval
		row = '{0},{1},{2},{3},NaN\n'.format(link, bound, time_stamp, flow)
		data_to_write.append(row)

	f.close()


f = open('historical_data.txt','wb')
for row in data_to_write:
	f.write(row)
f.close()












