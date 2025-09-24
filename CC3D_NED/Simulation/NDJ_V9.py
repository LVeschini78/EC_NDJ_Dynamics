
# Parameters to scan = Default value 
# To add more parameters to the scan add below in model AND in NDJ_SBML_Step     
sq = 1.1 # Track wether using square cells or real cell shapes
cr =  0.1 #0.05 #0.06 Relative proportion of cis interaction
ps_bd4 = 1.0 # Base Dll4 production rate. Increase ps_bd4 (e.g. 1.2)to simulate effect of high VEGF
ps_bn1 = 1.0# Base Notch1 production rate
ps_bn4 = 0.0 # Base Notch4 Jagged 1 production rate (if 0 J1 & N4 can only be produced by H1 induction)
ps_bj1 =0.0
ps_i = 0.0 # J1/N4 production 0 OFF, 1 ON (1= Vmax same as D4&N1, can have values different than 1 for weaker stronger ind)
ps_ja = 50 # K of Hill function for Jagged 1 (inducible) production. N1 induction & Dll4 repression have K =50
ps_Kdni = 1 # NICD degradation, default 1 (0.002/min).
ps_ni = 2000 #K of Hill function for NICD production. Increase to model effect of DAPT (default 2000, ki50=2100)
ps_kconv = 0.08 # 0.078 # Extent of productive cis interactions (fraction of cDN converted to cmDN)
auto = 0.0 # HES1 autoregulation Active/inactive 1/0
replicates = 1.0

data_file_head = f'cr{cr}_ji{ps_i}_aut{auto}_lyt{sq}'

model_string = f"""
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
        ps_bJ1 = 0.0;
        ps_ind = 0.0;
        ps_ja = 0.0;
        ps_ni = 2000;
        ps_Kdni = 0.0; 
        ps_kconv = 0.0;
        cr = 0;
        
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
        J1=0; VmJ1:= Amp*ps_ind; hMj1:= ps_ja; Kdj1= 0.05; dJ1= 0; bJ1:= Amp*ps_bJ1;
        
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