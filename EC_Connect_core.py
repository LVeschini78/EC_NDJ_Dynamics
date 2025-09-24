from cc3d.core.simservice.CC3DSimService import CC3DSimService
from cc3d.core.PySteppables import SteppableBasePy
from cc3d.core.PyCoreSpecs import (PottsCore, 
                                   Metadata, 
                                   CellTypePlugin,
                                   NeighborTrackerPlugin,
                                   PIFInitializer)

import numpy as np
from random import uniform
from scipy.stats import ks_2samp as ks
import json
import os

class ECConnectCore:
    """
    A class to manage and run EC_Connect simulations using CC3D with programmatic control.
    """

    def __init__(self, store: str, exp_name: str, params: dict) -> None:
        self.store = store
        self.exp_name = exp_name
        self.params = params
        self.initialise_core()

    # -------------------------------------------------------
    # SIMULATION CORE
    # -------------------------------------------------------
    def initialise_core(self):
        """
        Initialize the simulation core.
        Users must define the simulation specifications and register them here.
        """
        self.specs = []  # Initialize the specs list

        # Define Potts Core Specs
        simCore = PottsCore(dim_x=270, 
                            dim_y=270, 
                            dim_z=1, 
                            steps=2000, 
                            fluctuation_amplitude=5.0,
                            neighbor_order=2,
                            boundary_x='Periodic',
                            boundary_y='Periodic')
        self.specs.append(simCore)

        # Define Metadata Specs
        simMeta = Metadata(num_processors=1, debug_output_frequency=0)
        self.specs.append(simMeta)

        # Define Cell Type Plugin Specs
        simCells = CellTypePlugin()
        simCells.cell_type_append("EC", type_id=1, frozen=True)
        self.specs.append(simCells)

        simNeighborTracker = NeighborTrackerPlugin()
        self.specs.append(simNeighborTracker)
        # Define PIF Initializer Specs
        # Options used in the paper are HU_CT11.piff (medium VEGF HUVEC),
        # HU_VE2.piff (high VEGF HUVEC), QuadPattern.piff (Regular square cells)
        curr_wd = os.getcwd()
        piff_path = os.path.join(curr_wd, "CC3D_NED","Simulation", "HU_VE2.piff")

        pif_initializer = PIFInitializer(pif_name=piff_path)
        self.specs.append(pif_initializer)

        # Define output directory for simulation results.
        # It does not really matter as we store the data in a zarr file
        # But it is a requirement of CC3D
        out_dir = os.path.join(curr_wd, "Output_CC3D")

        # Instantiate the core simulation
        self.core_sim = CC3DSimService(output_dir=out_dir)
        self.core_sim.register_specs(self.specs)

        # Register Steppables (defined below)
        self.core_sim.register_steppable(ECConnectStep(store=self.store,
                                                    exp_name=self.exp_name,
                                                    exp_params=self.params,
                                                    field_names=['H1', 'NICD'],
                                                    frequency=1))

        # Run the simulation
        self.core_sim.run()
        self.core_sim.init()
        self.core_sim.start()
        

    # -------------------------------------------------------
    # DEFINE DEDICATED METHODS FOR PROGRAMMATIC CONTROL OF THE SIMULATION STEPS
    # -------------------------------------------------------
    def step(self) -> None:
        """
        Run the simulation for a single step without writing output.
        """
        # Set simulation input and run a single step
        self.core_sim.sim_input = ['skip_write']
        self.core_sim.step()
        # return self.core_sim.sim_output


    def step_write(self) -> object:
        """
        Run the simulation for a single step and write output.
        """
        #Set simulation input and run a single step with output
        self.core_sim.sim_input = ['write']
        self.core_sim.step()
        return self.core_sim.sim_output
        
# -------------------------------------------------------
# DEFINE A STEPPABLE CLASS TO ENCODE SIMULATION BEHAVIOUR
# -------------------------------------------------------
class ECConnectStep(SteppableBasePy):
    """
    A steppable to manage EC_Connect simulation steps.
    This class handles the interactions between the spatial context and the subcellular
    Notch-Delta-Jagged signaling model.
    Additionally, it writes simulation data to a specified store.
    """

    def __init__(self, store: str, exp_name: str, exp_params: dict, field_names: list, frequency: int = 1) -> None:
        super().__init__(frequency)
        # We need to scale down the parameters to feed to the SBML model 
        scaled_params = {k: round(v / 100, 4) for k, v in exp_params.items()}
        self.exp_params = scaled_params
        self.store = store
        self.exp_name = exp_name
        self.field_names = field_names
        # Set up the model string for the SBML simulation
        self.set_model_string()
        self.cumulative_data = []

        # Set up the store manager to handle data storage
        # Defined below
        self.store_manager = StoreManager(store) 
        if not os.path.exists(self.store):
            self.store_manager.create_store()
        self.store_manager.setup_exp(expname=exp_name, exparams=scaled_params)

        # Create scalar fields for each field name to facilitate access and fast write
        # Only works with CC3D version >= 4.7
        for field in self.field_names:
            self.create_scalar_field_py(field)

    def start(self):
            
        # Set global SBML options                      
        options = {'relative': 1e-10, 'absolute': 1e-12, 'step':1}
        self.set_sbml_global_options(options)
        # Each MCS the SBML will be run for 10 steps of 72s each (12 min total)
        step_size = 10
        self.timeseries = {}
        self.tracked_nb_id = 75 # Track a specific cell and its neighbors
        self.time_tracked = 200 # Time to track in MCS

        # Apply model string to cells
        for cell in self.cell_list_by_type(self.EC):

            self.add_antimony_to_cell(model_string=self.NED_model, model_name='NDS', cell=cell, step_size=step_size)

            cell.sbml.NDS['D4'] = uniform(0,3000)
            cell.sbml.NDS['N1'] = uniform(0,3000)
            cell.sbml.NDS['J1'] = uniform(0,3000)
            cell.sbml.NDS['tmDN'] = uniform(100,3000)

            for param in self.exp_params.keys():
                cell.sbml.NDS[param] = self.exp_params[param]
    # A method to retrieve fields in a loop for writing to store
    # It is the only method that works with CC3D < 4.7
    def retrieve_fields_loop(self):
        shelve_vals = np.zeros((self.dim.x, self.dim.y, self.dim.z, 2 + len(self.field_names)))  
        for x, y, z in self.every_pixel():
            cell = self.cell_field[x,y,z]
            if cell:
                shelve_vals[x,y,z,0] = cell.type
                shelve_vals[x,y,z,1] = cell.id
            
                shelve_vals[x,y,z,2] = cell.sbml.NDS['H1']
                shelve_vals[x,y,z,3] = cell.sbml.NDS['NICD']
        
        return shelve_vals
    
    def write_to_store(self, ks_stat: float , p_value: float, dict_data:dict=None):
        """
        Write simulation data to the store.
        """
        shelve_vals = self.retrieve_fields_loop()

        # Write data to the store
        self.store_manager.write_data(shelve_vals, ks_stat, p_value, dict_data=dict_data)
        
    def step(self, mcs: int):
        tracked_nb_id = self.tracked_nb_id
        time_tracked = self.time_tracked
        # Apply a warmup phase of 20 MCS
        if mcs > 20:
            # Step the SBML model for all cells
            self.timestep_sbml()
            
            # Calculate ligand dynamics, and track specific cells for reporting
            for cell in self.cell_list_by_type(self.EC):           
                cell.sbml.NDS['tmDN'] = 0
                cell.sbml.NDS['dD4'] = 0
                cell.sbml.NDS['dJ1'] = 0
                cell_N1 = cell.sbml.NDS['N1']
                tot_C_area=0
                # Calculate total common surface area for the cell
                nb = self.get_cell_neighbor_data_list(cell)
                tot_C_area = nb.common_surface_area_by_type()[1]
                if cell.id == tracked_nb_id and mcs > 1001 - time_tracked:
                    if 'center' in self.timeseries:
                        self.timeseries['center'].append(cell.sbml.NDS['H1'])
                    else:
                        self.timeseries['center'] = [cell.sbml.NDS['H1']]

                for neighbor_tuple in nb:       
                    nb_cell = neighbor_tuple[0]
                    common_area = neighbor_tuple[1]                
                    C_area_fr = common_area / tot_C_area
                                                
                    if nb_cell:
                        # Calculate the contribution of each ligand for the neighbour
                        nb_D4 = nb_cell.sbml.NDS['D4']                    
                        nb_J1 = nb_cell.sbml.NDS['J1']
                        d4_contrib = nb_D4 * C_area_fr
                        j1_contrib = nb_J1 * C_area_fr
                        n1_contrib = cell_N1 * C_area_fr
                                            
                        kJC= 1 
                        kpDJ1 = (d4_contrib/(1+d4_contrib+(kJC*j1_contrib))) # Competition of J1 for N1
                        kpDJ2 = 1-kpDJ1
                        nb_cell.sbml.NDS['dD4'] += min((kpDJ1*d4_contrib),(kpDJ1*n1_contrib)) #D4 donated in trans
                        nb_cell.sbml.NDS['dJ1'] += min((kpDJ2*d4_contrib),(kpDJ2*n1_contrib)) #J1 donated in trans                  
                        cell.sbml.NDS['tmDN'] += min((kpDJ1*d4_contrib),(kpDJ1*n1_contrib))

                        if cell.id == tracked_nb_id and mcs > 1001 - time_tracked:
                            if f'nb_{nb_cell.id}' in self.timeseries:
                                self.timeseries[f'nb_{nb_cell.id}'].append(nb_cell.sbml.NDS['H1'])
                            else:
                                self.timeseries[f'nb_{nb_cell.id}'] = [nb_cell.sbml.NDS['H1']]

        # At the end of the simulation (MCS 1000) or during the last 150 MCS collect data for KS statistic
        if mcs in range(850,1001):
            task = self.external_input[0]
            # print(f"Task: {task}")
            if task == 'write':
                # Get the ref data from the json file huvec_data.json
                with open(f'huvec_data.json', 'r') as f:
                    ref_data = json.load(f)
                
                ref_data = ref_data['huvec_data'] 
                c_data = self.cumulative_data
                ac_data = np.array(c_data).flatten()
                ar_data = np.array(ref_data).flatten()
                # print("shapes", ac_data.shape, ar_data.shape)
                ks_stat, p_value = ks(ac_data, ar_data)
                # print("KS Stat:", ks_stat, "p_value:", p_value)
                series = {'sim_data': ac_data.tolist(), 
                          'ref_data': ar_data.tolist(),
                          'nb_track_data': self.timeseries}
                self.write_to_store(ks_stat=ks_stat, p_value=p_value, dict_data=series)
                self.external_output = [ks_stat, p_value]

            else:
                curr_data =[cell.sbml.NDS['H1']/100 for cell in self.cell_list_by_type(self.EC)]
                self.cumulative_data.append(curr_data)

    # A method to set the NED model string
    def set_model_string(self):
        self.NED_model = f"""
        ## Utility functions

        function Hp(Vm,S,hM,h) 
            Vm*((pow(S/hM,h))/(1+pow(S/hM,h))) 
        end
        
        function Hc(Vm,A,R,Q,Ka,ha,Kr,hr,Kq,hq) 
            Vm*((pow(A/Ka,ha))/(1+pow(A/Ka,ha)+pow(R/Kr,hr)+pow(Q/Kq,hq))) 
        end
        
        function Hn(bp,S,hM,h)
            bp/(1+pow(S/hM,h))
        end
        
        function Ma(S,Kd)
            S*Kd
        end

        function MM(S,vmax,Km)
            vmax*S/(Km+S)
        end
        
        model NDS()
        
        # Timescale 1 tick= 1', in CC3D each MCS runs the SBML for 10 tiks= 10'               
        
        # To PS:
        h1_auto = 0.0; # 0.0 for no  H1 autoregulation, 1.0 for H1 autoregulation
        ps_bd4 = 0.0;
        ps_bn4 = 0.0;
        ps_bn1 = 0.0;
        ps_bj1 = 0.0;
        ps_ind = 0.0;
        ps_ja = 0.0;
        ps_ni = 2000;
        ps_Kdni = 0.0; 
        ps_kconv = 0.0;
        cr = 0;
        ps_K_nui = 0.0;
        
        ## NICD Production Degradation
        Vmni= 100; hMni:= ps_ni; Kdni:= ps_Kdni; tmDN=0; NICD = 0;
        NI_prod: => NICD; Hp(Vmni,(tmDN+cmDN),hMni,4); 
        NI_deg: NICD => ; Ma(NICD,Kdni);
        
        ## Ligands/receptors production and degradation
        Amp=10 # amplification mRNA=>Protein     

        # NOTCH1 expression, can be modifed to include HES1 inducibility
        N1=0; VmN1= Amp*1; hMn1=50; Kdn1=0.05; bN1:= Amp*ps_bn1; sc=25; #bN * scaling 25/min

        # N1 inducibility by HES1
        N1_prod: => N1; sc*(bN1);  # + Hp(VmN1,H1,hMn1,4)); 
        N1_deg: N1=> ; Ma(N1,Kdn1)+(tmDN*Hp(Kdn1,(tmDN+cmDN),hMni,4));;
        
        # DLL4 expression repressed by HES1
        D4=0; VmD4= Amp*1; hMd4=50; Kdd4=0.05; bD4:= Amp*ps_bd4; dD4=0;
        
        D4_prod: =>D4; sc*Hn(bD4,H1,hMd4,4); 
        D4_deg: D4=> ; Ma(D4+(dD4*D4/(1+D4)),Kdd4); 
          
        # Jagged1 expression basal and inducible by HES1
        J1=0; VmJ1:= Amp*ps_ind; hMj1:= ps_ja; Kdj1= 0.05; dJ1= 0; bJ1:= Amp*ps_bj1;
        
        J1_prod: =>J1; sc*(bJ1+Hp(VmJ1,H1,hMj1,4)); 
        J1_deg: J1=> ; Ma(J1+(J1*dJ1/(1+J1)),Kdj1);
        
        ## Transcription Factors expression
        # h1_auto 0.0 for no autoregulation, 1.0 for autoregulation
        
        # HES 1 Autoregulation with delay 
        vm_prna = 0.12; hh1 = 8; K = 0.5;
        pmRNA_prod: => pmRNA; h1_auto * Hc(vm_prna, NICD*0.01, pmRNA, HES, K, 1, K, 1, K, hh1); # Competitive hill function for inactive mRNA transcription
        
        kmrp = 0.1; kmrd = 0.05;
        mRNA_prod: pmRNA => mRNA; h1_auto * MM(pmRNA,kmrp,K); # Michaelis-Menten reaction for mRNA processing
        mRNA_deg: mRNA => ; h1_auto * MM(mRNA,kmrd,K); # Michaelis-Menten reaction for mRNA degradation
        
        # Michaelis-Menten reactions for protein translation, processing, and degradation
        kphp = 0.1; khp = 0.1; khd = 0.075; 
        pHES_prod:  => pHES; h1_auto * MM(mRNA,kphp,K); 
        HES_prod: pHES => HES;  h1_auto * MM(pHES,khp,K); 
        HES_deg: HES => ; h1_auto * MM(HES,khd,K); 

        # Hes scaling
        # H1 = 0
        H1_scaling: => H1; h1_auto * (HES * 100 - H1);
        
        # HES 1 Without autoregulatory feedback
        hMh1=35; 
        Kdh1 = 0.09; bH1 = 0.0; VmH1 = 10;
        H1_prod: =>H1; (1.0 - h1_auto) * (bH1 + Hp(VmH1,NICD,hMh1,2)); 
        H1_deg: H1=> ; (1.0 - h1_auto) * Ma(H1,Kdh1);

        ## Cis ligand-receptor interactions
        cDNi:= piecewise(D4, D4<N1, N1); # Chose the minimum of D4 and N1    
        cJNi:= piecewise(J1, J1<N1, N1); # Chose the minimum of J1 and N1
        
        # Relative proportion of cis interactions after trans interactions scaled by cr
        crd:=cr*(1-(tmDN/(300+tmDN))); 
        kcd= 0.10; # Degradation rate of cis complexes
        kconv:= ps_kconv; # Conversion rate of cDN complex to transcriptionally active cmDN

        cDN_complex_prod: D4+N1 => cDN;  Ma(cDNi,crd); 
        cmDN_prod_deg: cDN => cmDN; Ma(cDN,kconv) - 0.05 * cmDN;
        cDN_deg: cDN => ; Ma(cDN,kcd);
        
        cJN_complex_prod: J1+N1 => cJN; Ma(cJNi,crd);  
        cJN_deg: cJN => ; Ma(cJN,kcd);

        end
        """

# -------------------------------------------------------
# DEFINE A STORE MANAGER CLASS TO HANDLE DATA STORAGE
# -------------------------------------------------------
import zarr
class StoreManager():
    """ 
    A class to manage data storage for experiments.
    This class handles the creation of a store, setting up experiments,
    retrieving experiment parameters, and writing data to the store.
    It uses the zarr library for efficient storage and retrieval of data.
    """
    def __init__(self, store:str):
        self.store = store
        self.expname = None
        self.exparams = None
        
    def create_store(self):
        with zarr.open(self.store, mode='w') as root:
            data_group = root.create_group('experiments')
    
    def get_experiments(self):
        with zarr.open(self.store) as root:
            experiments = root['experiments']
            return list(experiments.keys())

    def setup_exp(self, expname:str, exparams:dict):
        with zarr.open(self.store) as root:
            if expname in self.get_experiments():
                return
            else:
                self.expname = expname
                self.exparams = exparams
                exp_group = root['experiments'].require_group(expname)
                exp_group.attrs['params'] = exparams

    def set_expname(self, expname:str):
        with zarr.open(self.store) as root:
            if expname in self.get_experiments():
                self.expname = expname
                self.exparams = root['experiments'][expname].attrs['params']
            else:
                print(f'Experiment {expname} not found')
                self.expname = None
                
    def get_exp_params(self, expname:str):
        with zarr.open(self.store) as root:
            if expname in self.get_experiments():
                exp_group = root['experiments'][expname]
                return exp_group.attrs['params']
            else:
                print(f'Experiment {expname} not found')
                return None

    def write_data(self, data:np.ndarray, ks_stat: float, p_value: float, dict_data:dict=None):     
        if self.expname:
            with zarr.open(self.store) as root:
                exp_group = root['experiments']
                exp_group = exp_group[self.expname]
                exp_group.require_dataset('data', shape=data.shape, dtype=data.dtype, data=data, compression='blosc', overwrite=True)
                exp_group.attrs['ks'] = [ks_stat, p_value]
                if dict_data is not None:
                    for key, value in dict_data.items():
                        exp_group.attrs[key] = value                
            
            print("_________________________________")
            print(f"Experiment: {self.expname}")
            print(f"Params: {self.exparams}") 
            print(f"Ks Stat: {ks_stat}, p_value: {p_value}")
            print("_________________________________") 
        else:
            print('No experiment set???')