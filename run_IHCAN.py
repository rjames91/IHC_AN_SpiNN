import numpy
import model_launch_framework
import pylab as plt
from vision.spike_tools.vis.vis_tools import plot_output_spikes

data = []
n_ihcan=352

for i in range(n_ihcan):
    data.append([])
    drnl_data=numpy.fromfile("./c_models/IHC_AN_float_MAP_14Copy/load_files/load" + str(i+1) +"_1",dtype='float32')
    #insure audio data can be divided evenly into 100 sample segements
    data[i].append(drnl_data[0:int(numpy.floor(len(drnl_data)/100)*100)])


#create framework of connected model vertices and run
samples = model_launch_framework.run_model(
    data, n_chips=numpy.ceil(n_ihcan/10.0),n_ihcan=10,fs=44100,resample_factor=1)


#convert to spike train
spike_trains=list()
spike_index=0
ihc_index=0

#obtain list of IHCAN outputs

ihc_output = [samples[x:x+2*int(numpy.floor(len(data[0][0])/100))*100] for x in xrange(0, len(samples), 2*int(numpy.floor(len(data[0][0])/100))*100)]

for ihc in ihc_output:
    #obtain fibre response
    HSR=[ihc[x:x+100] for x in xrange(0,len(ihc),200)]
    HSR = [item for sublist in HSR for item in sublist]
    #find non zero indicies for fibre and record in spike trains list
    #spike_trains.append((spike_index,numpy.nonzero(HSR)))
    idxs = numpy.nonzero(HSR)
    if len(idxs[0])==0:
        print "spike_index {} from ihc{} did not record in full, re-simulate this particular channel".format(spike_index,ihc_index)
    for i in idxs[0]:
        spike_trains.append((spike_index, i))

    #increment spike_index
    spike_index=spike_index+1

    #obtain fibre response
    LSR=[ihc[x:x+100] for x in xrange(100,len(ihc),200)]
    LSR = [item for sublist in LSR for item in sublist]
    # find non zero indicies for fibre and record in spike trains list
    #spike_trains.append((spike_index, numpy.nonzero(LSR)))
    idxs = numpy.nonzero(LSR)
    if len(idxs[0])==0:
        print "spike_index {} from ihc{} did not record in full, re-simulate this particular channel".format(spike_index,ihc_index)
    for i in idxs[0]:
        spike_trains.append((spike_index, i))

    # increment spike_index
    spike_index = spike_index + 1

    # increment ihc_index
    ihc_index=ihc_index + 1


#spike_trains=numpy.load("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains.npy")
#plot results
#plt.figure()
#plot_output_spikes(spike_trains,plotter=plt,markersize=1)
#plt.show()

# Save the results
#numpy.save("/home/rjames/Dropbox (The University of Manchester)/EarProject/spike_trains_v2.npy", spike_trains)
numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/results.csv", samples, fmt="%f", delimiter=",")
#numpy.savetxt("/home/rjames/Dropbox (The University of Manchester)/EarProject/complete.txt",[1],fmt="%f")

