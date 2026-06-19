model_string = """
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
        
        # To Scanned Parameters:
        auto = 0;   # If 1, use the auto-regulatory feedback of HES1, if 0, use the standard model without feedback
        ps_bd4 = 0.0;
        ps_bn1 = 0.0;
        ps_bJ1 = 0.0;
        ps_J1ind = 0.0;
        ps_ja = 0.0;
        ps_ni = 2000;
        ps_Kdni = 0.0; 
        ps_kconv = 0.0;
        cr = 0;

        ## NICD Production Degradation
        Vmni= 100/4; hMni:= ps_ni; Kdni:= ps_Kdni/4; tmDN=0; NICD = 0; 
        NI_prod: => NICD; Hp(Vmni,(tmDN+cmDN),hMni,4);                 # Eq (1)
        NI_deg: NICD => ; Ma(NICD,Kdni);                               # Eq (2)
        
        # Ligands/receptors production and degradation
        # NOTCH1 expression, can be modifed to include HES1 inducibility
        N1=0; VmN1= 250; hMn1=50; Kdn1=0.05; bN1:= 250;               

                                   # N1 inducibility by HES1
        N1_prod: => N1; Ma(1,bN1); # + Hp(VmN1,H1,hMn1,4));            # Eq (3)
        N1_deg: N1=> ; Ma(N1,Kdn1)+(tmDN*Hp(Kdn1,(tmDN+cmDN),hMni,4)); # Eq (4)
        
        # DLL4 expression repressed by HES1
        D4=0; hMd4=50; Kdd4=0.05; bD4:= 250; dD4=0;                    
        
        D4_prod: =>D4; Hn(bD4,H1,hMd4,4);                              # Eq (5)
        D4_deg: D4=> ; Ma(D4 + (dD4*D4/(1+D4)), Kdd4);                 # Eq (6)
          
        # Jagged1 expression basal and inducible by HES1
        J1=0; VmJ1:= amp * ps_J1ind; hMj1:= ps_ja; Kdj1= 0.05; dJ1= 0; bJ1:= 0;
        
        J1_prod: =>J1; Ma(1,bJ1) + Hp(VmJ1,H1,hMj1,4);                  # Eq (7) 
        J1_deg: J1=> ; Ma(J1 + (J1*dJ1/(1+J1)), Kdj1);                  # Eq (8)

        # HES 1 Without autoregulatory feedback
        H1 = 0
        hMh1=35; 
        Kdh1 = 0.075; bH1 = 0.0; VmH1 = 7.5;
        H1_prod: =>H1; auto * (bH1 + Hp(VmH1,NICD,hMh1,2));             # Eq (9)
        H1_deg: H1=> ; Ma(H1,Kdh1);                                     # Eq (10)


        # HES 1 WITH autoregulatory feedback
        vm_prna = 0.1; hh1 = 8; K = 0.5;
        pmRNA_prod: => pmRNA; (1-auto) * Hc(vm_prna, NICD, pmRNA, HES, K, 1, K, 1, K, hh1); # Eq (11)
        
        kmrp = 0.1; kmrd = 0.05;                                        
        mRNA_prod: pmRNA => mRNA; MM(pmRNA,kmrp,K);                     # Eq (12)
        mRNA_deg: mRNA => ; MM(mRNA,kmrd,K);                            # Eq (13)
        
        kphp = 0.1; khp = 0.1; khd = 0.075; 
        pHES_prod:  => pHES; MM(mRNA,kphp,K);                           # Eq (14)
        HES_prod: pHES => HES;  MM(pHES,khp,K);                         # Eq (15)
        HES_deg: HES => ; MM(HES,khd,K);                                # Eq (16)

        # Hes scaling
        # H1 = 0
        H1_scaling: => H1; (1-auto) * (HES * 100 - H1);                 # Eq (17)

        ## Cis ligand-receptor interactions
        cDNi:= piecewise(D4, D4<N1, N1); # Chose the minimum of D4 and N1    
        cJNi:= piecewise(J1, J1<N1, N1); # Chose the minimum of J1 and N1
        
        # Relative proportion of cis interactions after trans interactions scaled by cr
        crd:=cr*(1-(tmDN/(300+tmDN)));
        kcd= 0.10; # Degradation rate of cis complexes
        kconv:= ps_kconv; # Conversion rate of cDN complex to transcriptionally active cmDN

        cDN_complex_prod: D4+N1 => cDN;  Ma(cDNi,crd); 
        cmDN_prod_deg: cDN => cmDN; Ma(cDN,kconv) - 0.05 * cmDN;                 # Eq (18)
        cDN_deg: cDN => ; Ma(cDN,kcd);                                           # Eq (19)
        
        cJN_complex_prod: J1+N1 => cJN; Ma(cJNi,crd);                            # Eq (20)
        cJN_deg: cJN => ; Ma(cJN,kcd);                                           # Eq (21)

        end

        """
        
