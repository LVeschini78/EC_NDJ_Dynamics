from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *

class NDJ_Interactions(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
      # For each cell evaluate total amount of imput D4 and J1 (iD4/iJ1)
      # as function of common surface with neighbours.   

      if mcs> 20:
        
        for cell in self.cell_list_by_type(self.EC):
            
            cell.sbml.NDS['tmDN'] = 0
            # cell.sbml.NDS['tmJN'] = 0
            cell.sbml.NDS['dD4'] = 0
            cell.sbml.NDS['dJ1'] = 0
            cell_N1 = cell.sbml.NDS['N1']
            # cell_N4 = cell.sbml.NDS['N4']
            tot_C_area=0
            
            nb =self.get_cell_neighbor_data_list(cell)
            tot_C_area = nb.common_surface_area_by_type()[1]
                           
            for neighbor_tuple in nb:    
                nb_cell = neighbor_tuple[0]
                common_area = neighbor_tuple[1]                
                C_area_fr = common_area / tot_C_area
                                               
                if nb_cell:
                    nb_D4 = nb_cell.sbml.NDS['D4']                    
                    nb_J1 = nb_cell.sbml.NDS['J1']
                    d4_contrib = nb_D4 * C_area_fr
                    j1_contrib = nb_J1 * C_area_fr
                    n1_contrib = cell_N1 * C_area_fr
                    # n4_contrib = cell_N4 * C_area_fr
                    
                    
                    kJC= 1 
                    kpDJ1 = (d4_contrib/(1+d4_contrib+(kJC*j1_contrib))) # Competition of J1 for N1
                    kpDJ2 = 1-kpDJ1
                    nb_cell.sbml.NDS['dD4'] += min((kpDJ1*d4_contrib),(kpDJ1*n1_contrib)) #D4 donated in trans
                    nb_cell.sbml.NDS['dJ1'] += min((kpDJ2*d4_contrib),(kpDJ2*n1_contrib)) #J1 donated in trans                  
                    cell.sbml.NDS['tmDN'] += min((kpDJ1*d4_contrib),(kpDJ1*n1_contrib))                    
                   
                    # cell.sbml.NDS['tmJN'] += min((kpDJ2*j1_contrib),(kpDJ2*n4_contrib))
                    # nb_cell.sbml.NDS['dJ1'] += min((kpDJ2*j1_contrib),(kpDJ2*n4_contrib))
            #print(cell.sbml.NDS['tmDN'])            
            #print(cell.sbml.NDS['tmJN'])            
        #print("XXXXXXXXXX")                