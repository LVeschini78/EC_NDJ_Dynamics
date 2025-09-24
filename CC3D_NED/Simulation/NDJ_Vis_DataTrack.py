from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
from datetime import datetime
from NDJ_V9 import *

import pickle
# now = datetime.now()
# dt = now.strftime("%d%m%Y_%H_%M")

# Export data file path

tracking_filename = '/Users/lory/Desktop/1_Chesnais_et_al_2022_NOTCH Modelling/EC_Connect_V10_222022/EC_Connect_V10/Simulation/Data'
dat = {'step': [], 'cell': [], 'nicd': [], 'n4icd': [], 'hes': [], 'hey': []}

class NDJ_Vis_DataTrack(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.create_scalar_field_cell_level_py("D4")
        self.create_scalar_field_cell_level_py("J1")
        self.create_scalar_field_cell_level_py("N1")
        # self.create_scalar_field_cell_level_py("N4")
        self.create_scalar_field_cell_level_py("NICD")
        # self.create_scalar_field_cell_level_py("N4ICD")
        self.create_scalar_field_cell_level_py("H1")
        # self.create_scalar_field_cell_level_py("Hy")
        self.create_scalar_field_cell_level_py("tmDN")
        # self.create_scalar_field_cell_level_py("tmJN")

    def start(self):
        # self.trk_file = tracking_filename+ 'QP_' + str(ps_bd4) + '_' + str(ps_bn1) + '_' + str(replicates) + '_' + str(ps_ja) + '_' + str(cr)+ '.csv'


        # with open (self.trk_file,"w") as trk:
            # trk.write('bd4,bn1,jag,cr,rep,step,cell,nicd,n4icd,hes,hey')
            # trk.write('\n')
            
        # initialize setting for Histogram
        # self.plot_winh = self.add_new_plot_window(title='Histogram of HES', x_axis_title='Number of Cells',
                                                 # y_axis_title='Counts')
        # _alpha is transparency 0 is transparent, 255 is opaque
        # self.plot_winh.add_histogram_plot(plot_name='Hist 1', color='green', alpha=100)
        

        self.plot_win = self.add_new_plot_window(title='HES1 Neighborood',
                                                 x_axis_title='HES1 - Center Cell',
                                                 y_axis_title='HES1 - Neighbour Cell', 
                                                 x_scale_type='linear', 
                                                 y_scale_type='linear',
                                                 grid=False)
        self.plot_win2 = self.add_new_plot_window(title='HES1 One Cell',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='NICD', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.plot_all = self.add_new_plot_window(title='HES1 all Cells',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='HES1', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)                                                
                                                 
        self.plot_win.add_plot("Cell", style='Lines', color='red', size=8)      
        self.plot_win.add_plot("Cell1", style='Dots', color='red', size=1)
        self.plot_win.add_plot("Cell2", style='Dots', color='yellow', size=1)
        self.plot_win.add_plot("Cell3", style='Dots', color='green', size=1)
        self.plot_win.add_plot("Cell4", style='Lines', color='orange', size=4)
        self.plot_win.add_plot("Cell5", style='Lines', color='white', size=4)
        self.plot_win.add_plot("Cell6", style='Lines', color='pink', size=4)
        # self.plot_win.add_plot("Cell8", style='Dots', color='cadmiumyellow', size=3)
        # self.plot_win.add_plot("Cell9", style='Dots', color='deeppink1', size=3)
        self.plot_win2.add_plot("NICD", style='Lines', color='green', size=3)
        self.plot_all.add_plot("HES1", style='Dots', color='green', size=3)
        # print(dir(self.plot_win.plot_data))
        # print(self.plot_win.plot_data)
        
    def step(self, mcs):
        
        

        if mcs >10 and mcs % 1 == 0:

            d4 = self.field.D4
            n1 = self.field.N1
            j1 = self.field.J1
            nICD = self.field.NICD
            hES = self.field.H1
            tmdn = self.field.tmDN


            for cell in self.cell_list_by_type(self.EC):
                d4[cell] = cell.sbml.NDS['D4']/5000
                n1[cell] = cell.sbml.NDS['N1']/5000
                j1[cell] = cell.sbml.NDS['J1']/5000
                nICD[cell] = cell.sbml.NDS['NICD']/100
                hES[cell] = cell.sbml.NDS['H1']/100
                tmdn[cell] = cell.sbml.NDS['tmDN']/5000
                
                if mcs >= 200:
                    # self.plot_all.add_data_point("HES1", mcs, hES[cell])
                    
                    if cell.id == 75:
                        
                        # self.plot_win.add_data_point("Cell", mcs, hES[cell])
                        # neigh_n = len(self.get_cell_neighbor_data_list(cell))
                        nbs = self.get_cell_neighbor_data_list(cell)
                        tr_neigh, sa = nbs[1]
                        tr_neigh2, sa2 = nbs[2]
                        tr_neigh3, sa3 = nbs[3]
                        self.plot_win.add_data_point("Cell1", hES[cell], hES[tr_neigh])
                        self.plot_win.add_data_point("Cell2", hES[cell], hES[tr_neigh2])
                        self.plot_win.add_data_point("Cell3", hES[cell], hES[tr_neigh3])
                        
                        # self.plot_win2.add_data_point("NICD", mcs, hES[cell])
                if mcs >= 4850:
                    self.plot_all.add_data_point("HES1", mcs, hES[cell])
                    
            if mcs == 5000:
                
                # print(self.plot_win.plot_data)
                nb_filename =f'/Users/lory/Desktop/1_Chesnais_et_al_2022_NOTCH Modelling/NOTCH_Dynamics_V3/new_fig3_data/{data_file_head}_nb.pkl'
                with open(nb_filename, 'wb') as f:
                     pickle.dump(self.plot_win.plot_data, f)
                     
                all_filename =f'/Users/lory/Desktop/1_Chesnais_et_al_2022_NOTCH Modelling/NOTCH_Dynamics_V3/new_fig3_data/{data_file_head}_all.pkl'
                with open(all_filename, 'wb') as f:
                     pickle.dump(self.plot_all.plot_data, f)
                # self.plot_win.savePlotAsData(f'/Users/lory/Desktop/1_Chesnais_et_al_2022_NOTCH Modelling/NOTCH_Dynamics_V3/new_fig3_data/{data_file_head}_nb.csv',1)
                # self.plot_all.savePlotAsData(f'/Users/lory/Desktop/1_Chesnais_et_al_2022_NOTCH Modelling/NOTCH_Dynamics_V3/new_fig3_data/{data_file_head}_all.csv',1)
                
                    
                    
                    
