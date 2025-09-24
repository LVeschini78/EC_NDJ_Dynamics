from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *

class NDJ_CellMorph(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):        
        
        for cell in self.cell_list_by_type(self.EC):
                    
            surface = cell.surface
            if surface < 120:                
                cell.targetSurface = 120 
                cell.lambdaSurface = 1
            elif surface > 300:
                cell.targetSurface = surface*0.9 
                cell.lambdaSurface = 1
            else:
                cell.targetSurface = surface 
                cell.lambdaSurface = 1
                
                           
        # Rtrieve ECPT data and assign to cells
        # NOTCH_Dict = {}
        # with open('Path/data_file.csv') as f:             
            # for line in f:               
                # line= line.split()[0]
                # line=line.split(',')
                # NOTCH_Dict[int(line[2])]={"cy":float(line[0]),"nu":float(line[1])}

        # for cell in self.cell_list:
            # # if cell.id in NOTCH_Dict:
                # # # cell.sbml.NDS["mN1"]= 1+NOTCH_Dict[cell.id]["nu"]
                # # # cell.sbml.NDS["mDL4"]= 1-NOTCH_Dict[cell.id]["nu"]
                # # cell.sbml.NDS["NICD"]= NOTCH_Dict[cell.id]["cy"]
                # # print(str(cell.id) + " found")
            # # else:
                # # # cell.sbml.NDS["mN1"]= 1
                # # # cell.sbml.NDS["mDL4"]= 1
                # # cell.sbml.NDS["NICD"]= 0

    def step(self,mcs):
        pass                     
            
    def finish(self):
        
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):        
        return

