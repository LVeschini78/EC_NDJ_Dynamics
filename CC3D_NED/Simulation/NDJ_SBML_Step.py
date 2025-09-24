from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
from random import uniform
from NDJ_V9 import *

class NDJ_SBML_Step(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
            
        # Antimony model string stored in NDJ_V9                       
        options = {'relative': 1e-10, 'absolute': 1e-12, 'step':1}
        self.set_sbml_global_options(options)
        step_size = 10
        
        #****************************************************
        # Keep only few cells to evaluate isolated cells
        # for cell in self.cell_list_by_type(self.EC):
            # if cell.id not in [20,50,51,52,47,65,70]: 
                # self.delete_cell(cell)
        #****************************************************
        
        # Apply model string to cells        
        for cell in self.cell_list_by_type(self.EC):

            self.add_antimony_to_cell(model_string=model_string, model_name='NDS', cell=cell, step_size=step_size)

            cell.sbml.NDS['D4'] = uniform(0,3000)
            cell.sbml.NDS['N1'] = uniform(0,3000)
            cell.sbml.NDS['J1'] = uniform(0,3000)
            # cell.sbml.NDS['N4'] = uniform(0,3000)
            cell.sbml.NDS['tmDN'] = uniform(100,3000)
            # cell.sbml.NDS['tmJN'] = uniform(100,3000)

            #PS
            cell.sbml.NDS['ps_bd4'] = ps_bd4
            cell.sbml.NDS['ps_bn1'] = ps_bn1
            # cell.sbml.NDS['ps_bn4'] = ps_bn4
            cell.sbml.NDS['ps_ind'] = ps_i
            cell.sbml.NDS['ps_ja'] = ps_ja
            cell.sbml.NDS['ps_Kdni'] = ps_Kdni
            cell.sbml.NDS['ps_ni'] = ps_ni
            cell.sbml.NDS['ps_kconv'] = ps_kconv
            cell.sbml.NDS['cr'] = cr
            cell.sbml.NDS['h1_auto '] = auto
            
    def step(self, mcs):
        if mcs > 20:
            self.timestep_sbml()
            
# Simulate treatment with DAPT introduced at mcs 100        
        # if mcs == 100:
            # for cell in self.cell_list_by_type(self.EC):
                # cell.sbml.NDS['ps_ni'] = 2100
 