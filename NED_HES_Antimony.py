import multiprocessing as mp
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te

class NED_runner():
    def __init__(self, autoregulation:bool=False, preset:bool=True):
        # Complete Antimony Model definition for the NDS model
        # NED_string = f
        
        self.utility_functions = """
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
        
        """

        self.model_specs = """
        model NDS()
        
        # Timescale 1 tick= 1', in CC3D each MCS runs the SBML for 10 tiks= 10'               
        # TODO: remove To PS params
        # To PS:
        ps_bd4 = 0.0;
        ps_bn1 = 0.0;
        ps_bJ1 = 0.0;
        ps_ind = 0.0;
        ps_ja = 0.0;
        ps_ni = 2000;
        ps_Kdni = 0.0; 
        ps_kconv = 0.0;
        cr = 0;
        """

        self.NICD_reg = """
        ## NICD Production Degradation
        Vmni= 100; hMni:= ps_ni; Kdni:= ps_Kdni; tmDN=0; NICD = 0;
        NI_prod: => NICD; Hp(Vmni,(tmDN+cmDN),hMni,4); 
        NI_deg: NICD => ; Ma(NICD,Kdni);
        """

        self.ligands_receptors_reg = """
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
        J1=0; VmJ1:= Amp*ps_ind; hMj1:= ps_ja; Kdj1= 0.05; dJ1= 0; bJ1:= Amp*ps_bJ1;
        
        J1_prod: =>J1; sc*(bJ1+Hp(VmJ1,H1,hMj1,4)); 
        J1_deg: J1=> ; Ma(J1+(J1*dJ1/(1+J1)),Kdj1);
        
        """
        self.HES_auto_orig = """
        #HES1 With Autoregulatory feedback
        iH1=0; H1=0; VmH1= Amp*1; hMh1=50; Kdh1=0.09; bH1= 0.0; K_nuI:= ps_K_nui;
        =>iH1; bH1+Hc(VmH1,NICD,0.00000005*H1,hMh1,2); #=>iH1; bH1+Hp(VmH1,(NICD/(1+(H1/200))),hMh1,2);                
        iH1=>H1; Ma(iH1,K_nuI);
        H1=> ; Ma(H1,Kdh1);
        """

        self.HES_auto = """
        # HES 1 Autoregulation with delay 
        Kni = 0.5; Khe =0.5; 
        kdel = 0.55; kmm = 0.55; hh1 = 8; Mdh = 0.05;
        vm_prna = 1.4; kdh = 0.35;
        
        pmRNA_prod:  => pmRNA; Hc(vm_prna,(NICD*0.01),HES,pmRNA,Kni,1,Khe,hh1,Kni,1); # Competitive hill function for inactive mRNA production
        mRNA_prod: pmRNA => mRNA; MM(pmRNA,kdel,kmm); # Michaelis-Menten for active mRNA production
        pH1_prod: mRNA => pHES; MM(mRNA,kdel,kmm); # Michaelis-Menten for inactive protein production
        HES_prod: pHES => HES; MM(pHES,kdel,kmm); # Michaelis-Menten for active protein production
        HES_deg: HES => ; MM(HES,kdh,Mdh); # Michaelis-Menten for protein degradation

        # Hes scaling
        # H1 = 0
        H1_scaling: => H1; (HES * 100 - H1);
        """
        # CORRECT ONE
        self.NED_string_compact = """
        vm_prna = 0.1; hh1 = 8; K = 0.5; # Parameters for the competitive hill function
        pmRNA_prod: => pmRNA; Hc(vm_prna, NICD, pmRNA, HES, K, 1, K, 1, K, hh1); # Competitive hill function for inactive mRNA transcription
        kmrp = 0.1; kmrd = 0.05; # Parameters for mRNA processing and degradation
        mRNA_prod: pmRNA => mRNA; MM(pmRNA,kmrp,K); # Michaelis-Menten reaction for mRNA processing
        mRNA_deg: mRNA => ; MM(mRNA,kmrd,K); # Michaelis-Menten reaction for mRNA degradation
        kphp = 0.1; khp = 0.1; khd = 0.075; # Parameters for inactive and active protein translation and degradation
        pHES_prod:  => pHES; MM(mRNA,kphp,K); # Michaelis-Menten reaction for inactive protein translation
        HES_prod: pHES => HES;  MM(pHES,khp,K); # Michaelis-Menten reaction for active protein translation
        HES_deg: HES => ; MM(HES,khd,K); # Michaelis-Menten reaction for active protein degradation
        """
        self.HES_no_auto = """
        # HES 1 Without autoregulatory feedback
        hMh1=35; 
        Kdh1 = 0.09; bH1 = 0.0; VmH1 = 10;
        H1_prod: =>H1; (bH1 + Hp(VmH1,NICD,hMh1,2)); 
        H1_deg: H1=> ; Ma(H1,Kdh1);
        """

        self.cis_interactions = """
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
        if autoregulation is True:
            NED_string = self.utility_functions + self.model_specs + self.NICD_reg + self.ligands_receptors_reg + self.HES_auto + self.cis_interactions
        else:
            NED_string = self.utility_functions + self.model_specs + self.NICD_reg + self.ligands_receptors_reg + self.HES_no_auto + self.cis_interactions

        self.NED = te.loada(NED_string)
        if preset is True:            
            self.NED['tmDN'] = 0
            self.NED['ps_bd4'] = 1.0
            self.NED['ps_bn1'] = 1.0
            self.NED['ps_bJ1'] = 0.0
            self.NED['ps_ind'] = 0.7
            self.NED['ps_ja'] = 50
            self.NED['ps_ni'] = 2000
            self.NED['ps_Kdni'] = 1
            self.NED['ps_kconv'] = 0.0
            self.NED['cr'] = 0.0

        self.__results = None
        self.__values = None

    def display_component(self, component:str):
        print(getattr(self, component))

    def list_model_components(self):
        # List all attributes of the class instance
        attributes = dir(self)
        # Filter out methods and built-in attributes
        model_components = [attr for attr in attributes if not callable(getattr(self, attr)) and not attr.startswith("__") and not attr.startswith("_")]
        print("Model components:")
        for component in model_components:
            print(component)
        print("-------------------")
        
    def run(self, st:int, params:list[str] = ['time','NICD']):
            result = self.NED.simulate(0, st, st, params)
            return result
    
    def scan(self, 
             st:int, 
             par2sc:str, 
             rng:tuple, 
             points:int, 
             observe:str,
             time2upp:bool = False,
             legend:bool = False):
        results = []
        stabilities = []
        self.__values = np.linspace(rng[0], rng[1], points)
        for val in tqdm(self.__values):
            self.NED.reset()
            self.NED[par2sc] = val
            result = self.run(st, ['time', observe])
            results.append(result)
            if legend is True:
                # Find steady state
                self.NED.steadyState()
                eigen = self.NED.getFullEigenValues()
                stable = all(np.real(eigen) < 0)
                stabilities.append("Stable" if stable else "Unstable")
            if time2upp is True:
                half_upper = np.max(result[:,1]) / 2
                t2upp = np.argmax(result[:,1] > half_upper)
                print(f'val: {val.round(2)}, time to 1/2 max: {t2upp}, max_value: {np.max(result[:,1])}')
        self.__results = np.array(results)
        
        colormap = plt.get_cmap('jet')
        plt.figure(figsize=(4,2))
        for i,result in enumerate(self.__results):
            col_num = i/len(self.__results)
            colour = colormap(col_num)
            plt.plot(result[:,0], 
                     result[:,1], 
                     color=colour,
                     label=f"{self.__values[i].round(2)}, {stabilities[i]}" if legend is True else None)
            
        plt.xlabel('Time')
        plt.ylabel(observe)
        if legend is True:
            legend = plt.legend(loc='upper left', bbox_to_anchor=(1, 1.1),prop={'size': 8})
            legend.set_title(par2sc)
        plt.show()
    
    def scan_bu(self,
                par2sc:str, 
                rng:tuple, 
                observe:str,
                intercept:float = 0):
        results = []
        self.__values = np.linspace(rng[0], rng[1], 500)
        for val in tqdm(self.__values):
            self.NED.reset()
            self.NED[par2sc] = val
            result = self.run(1000, ['time', observe])
            bottom = np.min(result[:,1][500:999])
            upper = np.max(result[:,1][500:999])
            results.append([val, bottom, upper])

        self.__results = np.array(results)
        # plot the results
        colourmap = plt.get_cmap('jet')
        plt.figure(figsize=(4,2))
        plt.plot(self.__results[:,0], self.__results[:,1], color=colourmap(0.5), label='Bottom')
        plt.plot(self.__results[:,0], self.__results[:,2], color=colourmap(0.95), label='Upper')
        if intercept != 0.0:
            # Find the intercept
            dist_list = np.abs(self.__results[:,1] - intercept)
            print(f'Minimum distance to intercept: {np.min(dist_list)}')
            print(f'max value: {np.max(self.__results[:,2])}')
            if np.min(dist_list) > 0.01 * np.max(self.__results[:,0]):
                print('No intercept found')
            else:
                idx = np.argmin(dist_list)
                plt.plot([self.__results[idx,0], self.__results[idx,0]], [0, intercept], color='green', lw=0.5)
                plt.plot([0, self.__results[idx,0]], [intercept, intercept], color='green', lw=0.5)
                # Set a text label under the axis with the intercept value
                plt.text(self.__results[idx,0], 10, f'{int(self.__results[idx,0])}', color='green', fontsize=10)
        
        plt.xlim(rng)
        plt.ylim(0, np.max(self.__results[:,2] * 1.1))
        plt.xlabel(par2sc)
        plt.ylabel(observe)
        # plt.legend()
        plt.show()

class HES_runner:
    def __init__(self):
        self.NED_string_compact = """
        vm_prna = 0.1; hh1 = 8; K = 0.5; # Parameters for the competitive hill function
        pmRNA_prod: => pmRNA; Hc(vm_prna, NICD, pmRNA, HES, K, 1, K, 1, K, hh1); # Competitive hill function for inactive mRNA transcription
        kmrp = 0.1; kmrd = 0.05; # Parameters for mRNA processing and degradation
        mRNA_prod: pmRNA => mRNA; MM(pmRNA,kmrp,K); # Michaelis-Menten reaction for mRNA processing
        mRNA_deg: mRNA => ; MM(mRNA,kmrd,K); # Michaelis-Menten reaction for mRNA degradation
        kphp = 0.1; khp = 0.1; khd = 0.075; # Parameters for inactive and active protein translation and degradation
        pHES_prod:  => pHES; MM(mRNA,kphp,K); # Michaelis-Menten reaction for inactive protein translation
        HES_prod: pHES => HES;  MM(pHES,khp,K); # Michaelis-Menten reaction for active protein translation
        HES_deg: HES => ; MM(HES,khd,K); # Michaelis-Menten reaction for active protein degradation
        """
        self.NED_string = """
        NICD = 1.0;

        # NOTE: We replace scaled (NICD * 0.01) with NICD
        # Competitive hill function for inactive mRNA transcription        
        vm_prna = 0.1; hh1 = 8; K = 0.5;
        pmRNA_prod:  => pmRNA; vm_prna * NICD/K / (1 + NICD/K + pmRNA/K + HES^hh1/K^hh1);
        
        # Michaelis-Menten reactions for mRNA processing and degradation
        kmrp = 0.1; kmrd = 0.05;
        mRNA_prod: pmRNA => mRNA; kmrp * pmRNA / (K + pmRNA);
        mRNA_deg: mRNA => ; kmrd * mRNA / (K + mRNA); 

        # Michaelis-Menten reactions for protein translation, processing, and degradation
        kphp = 0.1; khp = 0.1; khd = 0.075;
        pHES_prod:  => pHES; kphp * mRNA / (K + mRNA);
        HES_prod: pHES => HES;  khp * pHES / (K + pHES); 
        HES_deg: HES => ; (khd * HES) / (K + HES); 
        """
        self.NED = te.loada(self.NED_string)
        
    def display_modelstring(self):
        print(self.NED_string)

    def compute_task(self, args):
        val1, val2, par1, par2, observe = args
        NED = te.loada(self.NED_string)  # Create NED object within the worker process
        NED.reset()
        NED[par1] = val1
        NED[par2] = val2
        if par1 != 'NICD' and par2 != 'NICD':
            NED['NICD'] = 1.0
        try:
            result = NED.simulate(0, 1000,1000, ['time', observe])
        except Exception as e:
            print(f"Error in compute_task with val1={val1}, val2={val2}: {e}")
            return np.array([val1, val2, -1.0, -1.0, -1.0, -1.0])
        res_slice = result[200:999]        
        topvalue = np.max(res_slice[:,1])
        # print("top", topvalue)
        std = NED.steadyState()
        # print("steady", std)
        eigen = NED.getFullEigenValues()
        stable = all(np.real(eigen) < 0)
        stbl = 1.0 if stable else 0.0
        amplitude = (np.max(res_slice[:,1]) - np.min(res_slice[:,1])).round(2) if stable is False else 0
        min_idx = np.where((res_slice[1:-1] < res_slice[:-2]) & (res_slice[1:-1] < res_slice[2:]))[0]
        if len(min_idx) > 2:
            period = int(res_slice[min_idx[-1],0] - res_slice[min_idx[-2],0])
        else:
            period = -1
        # print(np.array([val1, val2, topvalue, stbl]))
        return np.array([val1, val2, topvalue, stbl, amplitude, period])
        # except Exception as e:
        #     print(f"Error in compute_task with val1={val1}, val2={val2}: {e}")
        #     return None

    def parallel_compute(self, values1, values2, par1, par2, observe, num_cpus=1):

        tasks = [(val1, val2, par1, par2, observe) for val1 in values1 for val2 in values2]
        
        with mp.Pool(processes=num_cpus) as pool:
            results = list(tqdm(pool.imap(self.compute_task, tasks), total=len(tasks)))
        
        return results
    def add_colourbar(self, ax, vmin, vmax, cmap):
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        cbar = plt.colorbar(sm, ax=ax)
        # Define custom ticks
        ticks = [vmin]
        if vmax > 1.2 and ax == 'ax1':
            ticks.append(1.0)
            ticks.append(0.5)
            ticks.append(vmax)
        elif ax == 'ax3' or ax == 'ax4':
            ticks.append(vmin)
            ticks.append(vmax)
        return cbar
        # Set the ticks and labels
        # cbar.set_ticks(ticks)
        # cbar.set_ticklabels([f'{tick:.2f}' for tick in ticks])
        # return cbar

    def phase_portrait_parallel(self, par1, par2, rng1, rng2, observe, num_cpus=1, points=20, target=1.0):
        values1 = np.linspace(rng1[0], rng1[1], points).round(3)
        values2 = np.linspace(rng2[0], rng2[1], points).round(3)
        results = self.parallel_compute(values1, values2, par1, par2, observe, num_cpus)
        self.results = np.array(results)
        # plot the results
        colourmap = plt.get_cmap('jet')
        plt.figure(figsize=(4,4))
        fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        for ax in [ax1, ax2, ax3, ax4]:
            ax.set_box_aspect(1)

        for res in self.results:
            if res[2] == -1.0:
                continue
            # ax1 = plt.subplot(121)
            colour = 'green' if res[3] == 1 else 'red'
            if np.max(self.results[:,2]) == 0:
                scaled_topvalue = 0
            else:
                scaled_topvalue = res[2] / np.max(self.results[:,2])
            
            topcolour = colourmap(scaled_topvalue)
             
            # Plot maxima
            ax1.scatter(res[0], res[1], 
                color=topcolour, 
                edgecolor='none', 
                s=22)
            # overplot a small black dot for topvalues = 1.0±0.01
            # if abs(res[2]-target) < 0.05: 
            #     ax1.scatter(res[0], res[1], 
            #         color='black', 
            #         edgecolor='none', 
            #         s=5)
            # Highlight points of interest with a black circle (no fill)
            if abs(res[2]-target) < 0.05: 
                for ax in [ax1, ax2, ax3, ax4]:
                    ax.scatter(res[0], res[1], 
                        facecolors='none', 
                        edgecolors='black', 
                        linewidths=0.5,
                        s=22, zorder=10)
                    
            ax1.set_title('Maxima')
            ax1.set_xlabel(par1)
            ax1.set_ylabel(par2)

            # Add colourbar
            # a1_bar = self.add_colourbar(ax1, np.min(self.results[:,2]), np.max(self.results[:,2]), colourmap)

            # Plot stability
            ax2.scatter(res[0], res[1],
                color=colour, 
                edgecolor='none', 
                s=22)
            ax2.set_title('Stability')
            ax2.set_xlabel(par1)
            ax2.set_ylabel(par2)
            
            # Plot amplitude
            amplitude = res[4]
            max_amplitude = np.max(self.results[:,4])
            scaled_amplitude = amplitude / max_amplitude if max_amplitude > 0 else 0
            ampcolour = colourmap(scaled_amplitude)
            ax3.scatter(res[0], res[1],
                color=ampcolour, 
                edgecolor='none', 
                s=22)
            ax3.set_title('Amplitude')
            ax3.set_xlabel(par1)
            ax3.set_ylabel(par2)
            # a3_bar = self.add_colourbar(ax3, np.min(self.results[:,4]), np.max(self.results[:,4]), colourmap)
            # plot period
            period = res[5]
            min_period = np.min(self.results[:,5])
            max_period = np.max(self.results[:,5])
            scaled_period = (period - min_period) / (max_period - min_period) if max_period > min_period else 0
            periodcolour = colourmap(scaled_period)
            ax4.scatter(res[0], res[1],
                color=periodcolour, 
                edgecolor='none', 
                s=22)

            ax4.set_title('Period')
            ax4.set_xlabel(par1)
            ax4.set_ylabel(par2)
        
        a1_bar = self.add_colourbar(ax1, np.min(self.results[:,2]), np.max(self.results[:,2]), colourmap)
        a3_bar = self.add_colourbar(ax3, np.min(self.results[:,4]), np.max(self.results[:,4]), colourmap)
        a4_bar = self.add_colourbar(ax4, np.min(self.results[:,5]), np.max(self.results[:,5]), colourmap)
        # vmin=np.min(self.results[:,2])
        # vmax=np.max(self.results[:,2])
        # norm = plt.Normalize(vmin=vmin, vmax=vmax)
        # # Add colourbar
        # sm = plt.cm.ScalarMappable(cmap=colourmap, norm=norm)
        # sm._A = []
        # cbar = plt.colorbar(sm, ax=ax1)
        # # Define custom ticks
        # ticks = [vmin]
        # if vmax > 1.2:
        #     ticks.append(1.0)
        # ticks.append(0.5)
        # ticks.append(vmax)
        
        # # Set the ticks and labels
        # cbar.set_ticks(ticks)
        # cbar.set_ticklabels([f'{tick:.2f}' for tick in ticks])
        plt.tight_layout()
        plt.show()
    
    def run(self, st:int, params:list[str] = ['time','NICD']):
            result = self.NED.simulate(0, st, st*2, params)
            return result
    
    def scan(self, 
             st2p:int, 
             par2sc:str, 
             rng:tuple, 
             points:int, 
             observe:str,
             time2upp:bool = False,
             legend:bool = False):
        results = []
        stabilities = []
        frequencies = []
        amplitudes = []        
        st = 3000
        self.__values = np.linspace(rng[0], rng[1], points)
        for val in tqdm(self.__values):
            self.NED.reset()
            self.NED[par2sc] = val
            result = self.run(st, ['time', observe])
            results.append(result)
            if legend is True:
                self.NED.steadyState()
                eigen = self.NED.getFullEigenValues()
                stable = all(np.real(eigen) < 0)
                s_str = "S"
                stabilities.append(f"{s_str:>6}" if stable else "NS")

                if stable is False:         
                    res_slice = result[4000:5999]
                    # find the local maxima indexes i.e. where x-1 < x > x+1
                    max_idx = np.where((res_slice[1:-1] > res_slice[:-2]) & (res_slice[1:-1] > res_slice[2:]))[0]
                    min_idx = np.where((res_slice[1:-1] < res_slice[:-2]) & (res_slice[1:-1] < res_slice[2:]))[0]
                    amplitude = (np.max(res_slice[:,1]) - np.min(res_slice[:,1])).round(2)
                    # print(max_idx)
                    # find the period
                    if len(min_idx) > 1:
                        period = int(res_slice[min_idx[-1],0] - res_slice[min_idx[-2],0])
                    else:
                        period = '>600'
                else:
                    period = 'NA'
                    period = f"{period:>5}"
                    amplitude = 0
                frequencies.append(period)
                amplitudes.append(amplitude)

            if time2upp is True:
                half_upper = np.max(result[:,1]) / 2
                t2upp = np.argmax(result[:,1] > half_upper)
                print(f'val: {val.round(2)}, time to 1/2 max: {t2upp}, max_value: {np.max(result[:,1])}')
        self.__results = np.array(results)
        
        colormap = plt.get_cmap('jet')
        plt.figure(figsize=(4,2))
        for i,result in enumerate(self.__results):
            col_num = i/len(self.__results)
            colour = colormap(col_num)
            lbl = lbl = f"{self.__values[i]:>5.2f} {stabilities[i]:>5} {amplitudes[i]:>5.2f} {frequencies[i]:>5}" if legend is True else None
            plt.plot(result[:,0], 
                     result[:,1], 
                     color=colour,
                     label=lbl)
        plt.xlim(0, st2p)
        plt.xlabel('Time')
        plt.ylabel(observe)
        if legend is True:
            ctitles = [par2sc , 'Stb', 'Amp', 'Freq']
            title_str = f"{ctitles[0]:>13} {ctitles[1]} {ctitles[2]} {ctitles[3]}"
            legend = plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.1),prop={'size': 8})
            legend.set_title(title_str)
        plt.show()
    
    def scan_ap(self,
                par2sc:str, 
                rng:tuple, 
                observe:str,
                feat:str='amp'):
        results = []
        self.__values = np.linspace(rng[0], rng[1], 100)
        for val in tqdm(self.__values):
            self.NED.reset()
            self.NED[par2sc] = val
            result = self.run(3000, ['time', observe])    
            res_slice = result[4000:5999]
            # find the local maxima/minima indexes i.e. where x-1 < x > x+1 or x-1 > x < x+1
            max_idx = np.where((res_slice[1:-1] > res_slice[:-2]) & (res_slice[1:-1] > res_slice[2:]))[0]
            min_idx = np.where((res_slice[1:-1] < res_slice[:-2]) & (res_slice[1:-1] < res_slice[2:]))[0]
            amplitude = (np.max(res_slice[:,1]) - np.min(res_slice[:,1])).round(2)
            # find the period
            if abs(amplitude) < 0.01:
                period = 0
            else:
                if len(min_idx) > 1:
                    period = int(res_slice[min_idx[-1],0] - res_slice[min_idx[-2],0])
                else:
                    period = 510
            results.append([val, amplitude, period])
        self.results = np.array(results)
        # plot the results
        colourmap = plt.get_cmap('jet')
        plt.figure(figsize=(4,2))
        if feat == 'amp':
            plt.plot(self.results[:,0], self.results[:,1], color=colourmap(0.5), label='Amplitude')
        else:
            plt.plot(self.results[:,0], self.results[:,2], color=colourmap(0.5), label='Period')
     
        plt.xlim(rng)
        maxy= np.max(self.results[:,1]) if feat == 'amp' else np.max(self.results[:,2])
        plt.ylim(0, maxy * 1.1)
        plt.xlabel(par2sc)
        plt.ylabel('Amplitude' if feat == 'amp' else 'Period')
        # plt.legend()
        plt.show()
    
    def scan_bu(self,
                par2sc:str, 
                rng:tuple, 
                observe:str,
                intercept:float = 0):
        results = []
        self.__values = np.linspace(rng[0], rng[1], 500)
        for val in tqdm(self.__values):
            self.NED.reset()
            self.NED[par2sc] = val
            result = self.run(5000, ['time', observe])
            bottom = np.min(result[:,1][8000:9999])
            upper = np.max(result[:,1][8000:9999])
            results.append([val, bottom, upper])

        self.__results = np.array(results)
        # plot the results
        colourmap = plt.get_cmap('jet')
        plt.figure(figsize=(4,2))
        plt.plot(self.__results[:,0], self.__results[:,1], color=colourmap(0.5), label='Bottom')
        plt.plot(self.__results[:,0], self.__results[:,2], color=colourmap(0.95), label='Upper')
        if intercept != 0.0:
            # Find the intercept
            dist_list = np.abs(self.__results[:,1] - intercept)
            print(f'Minimum distance to intercept: {np.min(dist_list)}')
            print(f'max value: {np.max(self.__results[:,2])}')
            if np.min(dist_list) > 0.01 * np.max(self.__results[:,0]):
                print('No intercept found')
            else:
                idx = np.argmin(dist_list)
                plt.plot([self.__results[idx,0], self.__results[idx,0]], [0, intercept], color='green', lw=0.5)
                plt.plot([0, self.__results[idx,0]], [intercept, intercept], color='green', lw=0.5)
                # Set a text label under the axis with the intercept value
                plt.text(self.__results[idx,0], 10, f'{int(self.__results[idx,0])}', color='green', fontsize=10)
        
        plt.xlim(rng)
        plt.ylim(0, np.max(self.__results[:,2] * 1.1))
        plt.xlabel(par2sc)
        plt.ylabel(observe)
        # plt.legend()
        plt.show()


if __name__ == '__main__':
    mp.set_start_method('spawn')
    runner = HES_runner()
    runner.phase_portrait_parallel('NICD', 'hh1', (0, 1), (0.3, 0.8), 'HES', num_cpus=4)  # Set the number of CPUs to use