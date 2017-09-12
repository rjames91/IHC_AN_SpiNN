import spinnaker_graph_front_end as g

from IHCAN_vertex import IHCANVertex
import model_binaries

from pacman.model.constraints.placer_constraints\
    .chip_and_core_constraint import ChipAndCoreConstraint
from pacman.model.graphs.machine import MachineEdge

from spinnman.model.enums.cpu_state import CPUState

from spinn_front_end_common.utilities import globals_variables
from spinn_machine.utilities.progress_bar import ProgressBar


import numpy
import logging
import time

logger = logging.getLogger(__name__)


def run_model(
        data, n_chips=None,n_ihcan=0,fs=44100,resample_factor=1):

    # Set up the simulation
    g.setup(n_chips_required=n_chips, model_binary_module=model_binaries)

    # Get the number of cores available for use
    n_cores = 0
    machine = g.machine()

    # Create a OME for each chip
    boards = dict()

    #changed to lists to ensure data is read back in the same order that verticies are instantiated
    ihcans=list()

    cf_index=0
    count=0
    for chip in machine.chips:
        if count >= n_chips:
            break
        else:
            boards[chip.x, chip.y] = chip.ip_address

            for j in range(n_ihcan):
                ihcan=IHCANVertex(data[j][:],fs,resample_factor)
                g.add_machine_vertex_instance(ihcan)
                # constrain placement to local chip
                ihcan.add_constraint(ChipAndCoreConstraint(chip.x, chip.y))
                #ihcans[chip.x, chip.y,j] = ihcan
                ihcans.append(ihcan)

            count=count+1

# Run the simulation
    g.run(None)

    # Wait for the application to finish
    txrx = g.transceiver()
    app_id = globals_variables.get_simulator()._app_id
    #logger.info("Running {} worker cores".format(n_workers))
    logger.info("Waiting for application to finish...")
    running = txrx.get_core_state_count(app_id, CPUState.RUNNING)
    while running > 0:
        time.sleep(0.5)
        error = txrx.get_core_state_count(app_id, CPUState.RUN_TIME_EXCEPTION)
        watchdog = txrx.get_core_state_count(app_id, CPUState.WATCHDOG)
        if error > 0 or watchdog > 0:
            error_msg = "Some cores have failed ({} RTE, {} WDOG)".format(
                error, watchdog)
            raise Exception(error_msg)
        running = txrx.get_core_state_count(app_id, CPUState.RUNNING)

    # Get the data back
    samples = list()
    progress = ProgressBar(len(ihcans), "Reading results")

    for ihcan in ihcans:
        samples.append(ihcan.read_samples(g.buffer_manager()))
        progress.update()
    progress.end()
    samples = numpy.hstack(samples)

    # Close the machine
    g.stop()

    print "channels running: ",len(ihcans)/5.0
    print "output data: {} fibres with length {}".format(len(ihcans)*2,len(samples))
    if(len(samples) != len(ihcans)*2*numpy.floor(len(data[0][0])/100)*100*(1.0/resample_factor)):
        print "samples length {} isn't expected size {}".format(len(samples),len(ihcans)*2*numpy.floor(len(data[0][0])/100)*100*(1.0/resample_factor))

    return samples
