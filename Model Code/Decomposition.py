#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Package Import#
import numpy as np
import pandas as pd
import math


# In[2]:


##Pool for Plant Residue/Coarse Woody Debris(Litter/Lignin)##
class CWD_pool():   #Class of litter fall, three attribute: C, N, Lignin#
    def __init__(self,obser_data,sample_site,name):   #Initialization#
        if name == "litter":
            self.C = obser_data.loc[sample_site,"litter_C"]   #C Content in Plant Residue/Coarse Woody Debris#
            self.N = obser_data.loc[sample_site,"litter_N"]   #N Content in Plant Residue/Coarse Woody Debris#
            self.lignin = obser_data.loc[sample_site,"litter_lignin"]   #Lignin percent in Plant Residue/Coarse Woody Debris#

    def AddCarbon(self,addc):    #Add Carbon#
        self.C = self.C+addc   
    def SubCarbon(self,subc):    #Subtract Carbon#
        if self.C > subc:
            self.C = self.C-subc
        else:
            self.C = 0    
    def AddNitrogen(self,addn):   #Add Nitrogen#
        self.N = self.N + addn
    def SubNitrogen(self,subn):   #Subtract Nitrogen#
        if self.N > 0:
            self.N = self.N-subn
        else:
            self.N = 0          
    def CNRatio(self):            #CN Ratio#
        if self.N == 0:
            print ('NA')
        else:
            CN = self.C/self.N
        return CN


# In[3]:


##Pools for SOM##
class SOM_pool(): #Class:Different soil carbon pool, two commom attribute: C, N, two special attribute for Microbe: CU, NU #
    def __init__(self,obser_data,sample_site,name):   #Pool Size Initialization#
        self.C = obser_data.loc[sample_site,name+'_C']   #C Content in Pool#
        self.N = obser_data.loc[sample_site,name+'_N']   #N Content in Pool#
        self.CU = 0   #Usable Carbon for Microbial Pools#
        self.NU = 0   #Usable Nitrogen for Microbial Pools#
        
    def AddCarbon(self,addc):      #Add Carbon#
        self.C = self.C+addc
    def SubCarbon(self,subc):      #Subtract Carbon#
        if self.C > subc:
            self.C = self.C-subc
        else:
            self.C = 0
    def AddNewCarbon(self,addcu):   #Carbon Substrate for Microbial Pools#
        self.CU = self.CU+addcu      
    def SubNewCarbon(self,subcu):   #Empty Carbon Substrate for Microbial Pools#
        if self.CU > subcu:
            self.CU = self.CU-subcu
        else:
            self.CU = 0  
    def AddNitrogen(self,addn):     #Add Nitrogen#
        self.N = self.N+addn   
    def SubNitrogen(self,subn):     #Subtract Nitrogen#
        if self.N > 0:
            self.N = self.N-subn
        else:
            self.N = 0        
    def AddNewNitrogen(self,addnu): #Nitrogen Substrate for Microbial Pools#
        self.NU = self.NU+addnu  
    def SubNewNitrogen(self,subnu): #Empty Nitrogen Substrate for Microbial Pools#
        if self.NU > 0:
            self.NU = self.NU-subnu
        else:
            self.NU = 0
    def CNRatio(self):   #CN Ratio#
        CN = self.C/self.N
        return CN


# In[4]:


##Available Nitrogen pool##
class av_Nitro():
    def __init__(self,obser_data,sample_site):   #Initialization Available Nitrogen Pool, two attributes: ammonium, nitrate#
        self.beta = 0.05
        self.NH4 = obser_data.loc[sample_site,"NH4"] #Soil NH4 Content#
        self.NO3 = obser_data.loc[sample_site,"NO3"] #Soil NO3 Content#
    def avn(self):   #Available Nitrogen Content#
        avn_nitro = self.NH4+self.NO3 #Available Nitrogen Content = Soil NH4 + Soil NO3#
        return avn_nitro
    def imm(self,Nmin):   #Immobilize Available Nitrogen to Supplement deficient N in Decomposition, Nmin: N lacking for decomposition#
        avnitro = self.NH4+self.NO3
        if avnitro > 0:   #Calculating NH4/NO4 content#
            f_INH4 = (self.NH4/avnitro)/(self.beta+((1-self.beta)*self.NH4/avnitro)) # fraction of immobilization
            INH4 = min(f_INH4*Nmin,self.NH4)   #INH4 : NH4 Immobilized#
            INO3 = min(Nmin-INH4,self.NO3)     #INO3 : NO3 Immobilized#
        else:
            INH4 = 0
            INO3 = 0
        self.NH4 = self.NH4-INH4 #Updating NH4#
        self.NO3 = self.NO3-INO3 #Updaitng NO3#
    def Addmine(self,addmine):   #Mineralizing Organic Nitrogen into Available Nitrogen#
        self.NH4 = self.NH4+addmine   #Mineralization: All N mineralize to NH4# 
    def nitri(self,nitri_rate):       #Nitrification: NH4 tranfer into NO3#
        nitrification = self.NH4*nitri_rate  #Nitrification: NH4 Nitrify into NO3#
        self.NH4 = self.NH4-nitrification    #Updating NH4#
        self.NO3 = self.NO3+nitrification    #Updating NO3#


# In[5]:


##Organic Nitrogen pool##
class Org_Nitro():       #Organic Nitrogen Pool (This class is not used in this model currently)#
    def __init__(self):
        self.beta = 0.05
        self.OrgN = 0  
        
    def AddOrgNitrogen(self,addOrgN):   #Add Organic Nitrogen#
        self.OrgN = self.OrgN+addOrgN
        return self.OrgN 
    def mine(self,mine_rate):   #Organic Nitrogen Mineralization#
        mineralization = self.OrgN*mine_rate
        return mineralization


# In[6]:


##Pool/Flow in litter/Soil Substrate##  #Pool for MB and Substrate for litter and soil#
class initial_content():      #Carbon Flow from the litter/soil pool: determine whether exist C flow from CWD and Soil#
    def __init__(self):
        self.litterzero = 0   #Carbon flow out from CWD pool(1st day): Updata once#
        self.littersub  = 0   #Carbon flow out from CWD: Update daily# ###Q: we assigned litter.C*kbase_litter as the initial value of littersub
        #self.littersub  = litter.C*kfrag_litter
        self.soilzero   = 0   #Carbon flow out from soil(1st day): Update once#
        self.soilsub    = 0   #Carbon flow out from soil: Update daily# ###Q: we assign NOM.C*kmax_NOM to soilsub in cal.py file
        #self.soilsub   = NOM.C*kNOM   
        self.soilzeropot = 0  #Potential Carbon flow out from soil(1st day)#
        self.soilsubpot = 0   #Potential Carbon flow out fron soil: Update daily#
        self.microSMB1Zero  = 0   #Initial C Pool size of SMB1#  ###Q:change to initmicroSMB1 or not?# NO, we assigned SMB1 and SMB2 initial value in cal.py file 
        self.microSMB2Zero  = 0   #Initial C Pool size of SMB2#
    def LitterZero(self,substrate):   #(1st)Carbon flow out form CWD pool#
        self.litterzero = substrate
    def LitterSub(self,substrate):    #Daily Updating Carbon flow out from CWD#
        self.littersub = substrate   
    def SoilZero(self,substrate):     #(1st)Carbon flow out from Soil#
        self.soilzero = substrate
    def SoilSub(self,substrate):      #Daily Updating Carbon flow out from Soil#
        self.soilsub = substrate   
    def SoilZeroPot(self,substrate):  #(1st)Potential Carbon flow out from soil pool# ###Q:what does these potential stand for#
        self.soilzeropot = substrate
    def SoilSubPot(self,substrate):   #Daily Updating potential Carbon flow out from soil pool#
        self.soilsubpot = substrate  
    def MicroSMB1Zero(self,poolsize):     #Initial pool size of SMB1#
        self.microSMB1Zero = poolsize   
    def MicroSMB2Zero(self,poolsize):     #Initial pool size of SMB2#
        self.microSMB2Zero = poolsize


# In[7]:


#Class for litter-lignin-soil CO2 seperation#
class lls_CO2():                                      # litligsoil_CO2 in cal.py
    def __init__(self):                               # This is class attribute
        # self.labelLignin_pct = 0.004/(0.1*0.19)
        self.labelLignin_pct = 0.031                  ###QQQ: 20210608: Lab incibation(label lignin)
        self.litter_CO2  = 0
        self.lignin_CO2  = 0
        self.soil_CO2    = 0
        self.litterC_MB1 = 0
        self.ligninC_MB1 = 0
        self.soilC_MB1   = 0
        self.litterC_MB2 = 0
        self.ligninC_MB2 = 0
        self.soilC_MB2   = 0
        
    def litterCO2(self,CO2):                         # This is class method
        self.litter_CO2 = self.litter_CO2 + CO2
    def ligninCO2(self,CO2):             
        self.lignin_CO2 = self.lignin_CO2 + CO2
    def soilCO2(self,CO2):               
        self.soil_CO2 = self.soil_CO2 + CO2
    
    def litterToMB1(self,carbon):            
        self.litterC_MB1 = self.litterC_MB1 + carbon
    def ligninToMB1(self,carbon):            
        self.ligninC_MB1 = self.ligninC_MB1 + carbon
    def soilToMB1(self,carbon):              
        self.soilC_MB1 = self.soilC_MB1 + carbon    
    def litterToMB2(self,carbon):            
        self.litterC_MB2 = self.litterC_MB2 + carbon
    def ligninToMB2(self,carbon):            
        self.ligninC_MB2 = self.ligninC_MB2 + carbon
    def soilToMB2(self,carbon):              
        self.soilC_MB2 = self.soilC_MB2 + carbon    


# In[8]:


### Corse Woody Debris fragmentation and allocation
def CWDFrag_OrgInAllocLab (sample_site,obser_data,litter,AOM1,AOM2,AOM3,CN_AOM1,CN_AOM2,CN_AOM3):
    
    ###QQQ: For lab version: CWD allocate to AOM pool directly because the added material was not whole plant, and no N involved in this processes#     # 
    ##Environmental Factors Scalar##
    #EnVirFactor 1: temp# 
    temp = obser_data.loc[sample_site,"temp"]   #Unit: celsius degree# 
    if temp < (-10):  #Celsius < 10,t_scale=0
        t_scale = 0
    else:
        t_scale = np.exp(308.56*((1.0/71.02)-(1.0/((temp+273.15)-227.13))))   #t_scale = 0.8817
    #EnVirFactor 2: water# 
    water_cont = obser_data.loc[sample_site,"water_cont"]   #Volumetric Water Content: (Unit:%)#
    water_sat  = obser_data.loc[sample_site,"water_sat" ]
    wfp = water_cont/water_sat
    w_scale = (1-np.exp(-wfp))/0.632
    w_scale = min(1,max(0,w_scale))   #w_scale=0.895#
    
    #kfrag_litter = (1/kbase_litter)/365*t_scale*w_scale    ###QQQ: 20210601: Only for lab incubation: CWD allocation at one time. Decomposition Rate of Plant Residue(Loss Ratio Every Day)#
    litter_loseC_lab = litter.C                ###QQQ: 20210601: Only for lab incubation: CWD allocation at one time. Total Carbon Loss from Plant Residue(CWD) every day# 
    litter_loseN_lab = litter.N                ###QQQ: 20210601: Only for lab incubation: CWD allocation at one time.Total Nitrogen Loss from Plant Residue(CWD) every day# 
    
    ###QQQ: 20210601: For lab incubation version
    #First: Plant residues Allocating to AOM3(Lignin) Pool#
    litter_loseC_AOM3_lab  = litter_loseC_lab * (litter.lignin/100)      ###QQQ:20210608: for lab (lignin,Not changed currently)
    litter_loseN_AOM3_lab  = litter_loseN_lab * 0 
    litter_loseC_AOM12_lab = litter_loseC_lab-litter_loseC_AOM3_lab
    litter_loseN_AOM12_lab = litter_loseN_lab-litter_loseN_AOM3_lab
    #Second: Plant Residues Allocating to AOM1/AOM2(Litter) Pools#
    #One:When C/N loss from Plant Residue#
    if litter_loseC_AOM12_lab >0 and litter_loseN_AOM12_lab >0:
        if litter_loseC_AOM12_lab/litter_loseN_AOM12_lab <= CN_AOM1:  #CN Ratio of Plant Residue < AOM1(Liable Litter pool)
            litter_loseC_AOM1_lab = litter_loseC_AOM12_lab      #All Carbon allocating to AOM1#
            litter_loseN_AOM1_lab = litter_loseN_AOM12_lab      #All Nitrogen allocating to AOM2#
            litter_loseC_AOM2_lab = 0
            litter_loseN_AOM2_lab = 0    
        elif CN_AOM1 < litter_loseC_AOM12_lab/litter_loseN_AOM12_lab <= CN_AOM2: #CN Ratio of Plant Residue ∈ [AOM1,AOM2]#
            alloc_AOM1_lab = (CN_AOM2*litter_loseN_AOM12_lab-litter_loseC_AOM12_lab)/(CN_AOM2-CN_AOM1)
            alloc_AOM2_lab = (litter_loseC_AOM12_lab-CN_AOM1*litter_loseN_AOM12_lab)/(CN_AOM2-CN_AOM1)
            allocRatio_AOM1_lab = alloc_AOM1_lab/(alloc_AOM1_lab+alloc_AOM2_lab) #CN Proportion to AOM1#
            allocRatio_AOM2_lab = alloc_AOM2_lab/(alloc_AOM1_lab+alloc_AOM2_lab) #CN Proportion to AOM2#
            litter_loseC_AOM1_lab = litter_loseC_AOM12_lab*(allocRatio_AOM1_lab)/(allocRatio_AOM1_lab+allocRatio_AOM2_lab)
            litter_loseC_AOM2_lab = litter_loseC_AOM12_lab*(allocRatio_AOM2_lab)/(allocRatio_AOM1_lab+allocRatio_AOM2_lab)
            litter_loseN_AOM2_lab = litter_loseC_AOM2_lab/150                  ###QQQ: 20210602: For lab version: CN_AOM2=150, CN_AOM1=15~25
            litter_loseN_AOM1_lab = litter_loseN_AOM12_lab - litter_loseN_AOM2_lab
        elif litter_loseC_AOM12_lab/litter_loseN_AOM12_lab > CN_AOM2: #CN Ratio of Plant Residue > AOM2#
            litter_loseC_AOM1_lab = 0
            litter_loseN_AOM1_lab = 0
            litter_loseC_AOM2_lab = litter_loseC_AOM12_lab #All Carbon allocating to AOM2#
            litter_loseN_AOM2_lab = litter_loseN_AOM12_lab #All Nitrogen allocating to AOm2#
    #Two:When no C/N loss from Plant Residue#
    if litter_loseC_AOM12_lab==0 and litter_loseN_AOM12_lab==0: #No Plant Residue input#
        litter_loseC_AOM1_lab = 0
        litter_loseN_AOM1_lab = 0
        litter_loseC_AOM2_lab = 0
        litter_loseN_AOM2_lab = 0
    
    ###QQQ:20210601: For lab incubation###
    AOM1.AddCarbon(litter_loseC_AOM1_lab)                   ###QQQ:20210601: Updating Carbon content in AOM1#
    AOM1.AddNitrogen(litter_loseN_AOM1_lab)                 ###QQQ:20210601:Updating Nitrogen content in AOM1#
    AOM2.AddCarbon(litter_loseC_AOM2_lab)                   ###QQQ:20210601:Updating Carbon Content in AOM2#
    AOM2.AddNitrogen(litter_loseN_AOM2_lab)                 ###QQQ:20210601:Updating Nitrogen Content in AOM2#
    #AOM3.AddCarbon(litter_loseC_AOM3_lab*0.031)            ###QQQ:20210608:Label Lignin cal(Hide 97% lignin carbon)  ###QQQ:20210601:Updating Carbon content in AOM3#
    #AOM3.AddNitrogen(litter_loseN_AOM3_lab*0.031)          ###QQQ:20210608:label Lignin cal(Hide 97% lignin carbon)  ###QQQ:20210601:Updating Nitrogen content in AOM3#
    AOM3.AddCarbon(litter_loseC_AOM3_lab)                   ###QQQ:20210908:Use all lignin to cal, because of lignin protection effect
    AOM3.AddNitrogen(litter_loseN_AOM3_lab)                 ###QQQ:20210908:Use all lignin to cal, because of lignin protection effect
    
    litter.SubCarbon(litter_loseC_lab)
    litter.SubNitrogen(litter_loseN_lab)
    
    # return(alloc_AOM1_lab,alloc_AOM2_lab,allocRatio_AOM1_lab,allocRatio_AOM2_lab)


# In[9]:


### Decomposition - Tentative decomposition
##Function for Tentative Decomposition Procedure##
def SoilDecom_pot(sample_site,t,obser_data,potimm,flag_SMB1,flag_SMB2,lig_p1,lig_p2,             #Scalar(6)#   ###QQQ: 20210908: lignin decay log para #
                  AOM1,AOM2,AOM3,SMB1,SMB2,SMR,NOM,MOM,initCont,                   #Class(9)#
                  CN_SMB1,CN_SMB2,CN_SMR,CN_NOM,CN_MOM,                            #CN Ratio(5)#
                  kmax_AOM1,kmax_AOM2,kmax_AOM3,kmax_SMR,kmax_NOM,kmax_MOM,        #Decay Rate(6)#
                  fAOM1_SMB1,fAOM2_SMB1,fAOM_MOM,fAOM3_SMB1,fNOM_SMB1,fNOM_MOM,    #Fraction(6)#
                  fE_AOM1,fE_AOM2,fE_AOM3,fE_NOM,fA_SMB1,fA_SMB2,fM_SMB1,fM_SMB2): #Microbe(8)#
    
    #Environmental Factors Scalar#
    #EnVirFactor 1: temp#
    temp    = obser_data.loc[sample_site,"temp"]                 #temp(unit:℃)# 
    t_scale = 4.89*np.exp(-3.432+0.1*temp*(1-0.5*temp/36.9))     #t_scale=0.7698#
    #EnVirFactor 2: water#
    theta     = obser_data.loc[sample_site,"water_cont"]         #Volumetric water content(unit:%)#
    theta_fc  = obser_data.loc[sample_site,"water_fc"]           #Field water content#
    theta_sat = obser_data.loc[sample_site,"water_sat"]          #Saturated water content#
    if theta <= theta_fc:
        w_scale = (1-np.exp(-theta/theta_sat))/(1-np.exp(-theta_fc/theta_sat))
    elif theta > theta_fc:
        w_scale = 1.0044-(0.0044/(np.exp(-5*(((theta/theta_sat)-(theta_fc/theta_sat))/(1-theta_fc/theta_sat))))) #w_scale=1.0#
    #EnVirFactor 3: clay#
    Pclay = obser_data.loc[sample_site, "Pclay"]                #Clay Proportion:Unit(percent without %)#
    clay_scale = 1-(0.75*Pclay/100)                             #clay_scale = 0.722#

    ###Step1: ASOM(Active Soil OM Pool) Allocating to PSOM,SMB1,SMB2###
    kNOM = 1/kmax_NOM/365*t_scale*w_scale    #Decay rate of NOM pool#
    #Decomposition_Tentative:CN Loss & Transfer#
    C_leaveNOM_pot    = NOM.C*kNOM                                #Potential Carbon Leave from NOM Pool#
    C_NOM_MOM_pot     = C_leaveNOM_pot*fNOM_MOM                   #Carbon transfer from NOM to MOM#
    C_NOM_SMB1_pot    = C_leaveNOM_pot*(1-fNOM_MOM)*(fNOM_SMB1)   #Carbon transfer from NOM to SMB1#
    C_NOM_SMB2_pot    = C_leaveNOM_pot*(1-fNOM_MOM)*(1-fNOM_SMB1) #Carbon transfer from NOM to SMB2#
    C_NOM_SMB1_E_pot  = C_NOM_SMB1_pot*fE_NOM                     #Carbon Used by SMB1#
    C_NOM_SMB2_E_pot  = C_NOM_SMB2_pot*fE_NOM                     #Carbon Used by SMB2#
    C_NOM_SMB1_rh_pot = C_NOM_SMB1_pot*(1-fE_NOM)                 #Carbon directly loss as CO2#
    C_NOM_SMB2_rh_pot = C_NOM_SMB2_pot*(1-fE_NOM)                 #Carbon directly loss as CO2#
    N_leaveNOM_pot    = NOM.N*kNOM                                #Potential Nitrogen leave from NOM pool#
    N_NOM_MOM_pot     = N_leaveNOM_pot*fNOM_MOM                   #Nitrogen transfer from NOM to MOM#
    N_NOM_SMB1_pot    = N_leaveNOM_pot*(1-fNOM_MOM)*fNOM_SMB1     #Nitrogen transfer from NOM to AOM1#
    N_NOM_SMB2_pot    = N_leaveNOM_pot*(1-fNOM_MOM)*(1-fNOM_SMB1) #Nitrogen transfer from NOM to AOM2#
    #Substrate content(When No litter Imput): CN loss from Soil#
    initCont.SoilSubPot(C_leaveNOM_pot)
    if t==1: #Initial Substrate#
        initCont.SoilZeroPot(C_leaveNOM_pot)
    #N Coupling:Immobilization/Mineralization#
    #Determing N Immobilization Potential in NOM Decomposition--bimm_NOM#
    N_need = (C_NOM_MOM_pot/CN_MOM)+(C_NOM_SMB1_E_pot/CN_SMB1)+(C_NOM_SMB2_E_pot/CN_SMB2) #Actual Nitrogen demand from MOM to NOM#
    if N_need >= N_leaveNOM_pot: 
        bimm_NOM = 1 #Immobilization#
    else:
        bimm_NOM = 0 #Mineralization#
    #Updating Potebtial Immobilization#
    potimm = potimm+(N_need-N_leaveNOM_pot) #Potential N Immobilization Content along Decomposition#
    
    ###Step2: PSOM allocating to ASOM###
    kMOM = 1/kmax_MOM/365*t_scale*w_scale*clay_scale #Decay rate of MOM pool#
    #Decomposition_Tentative:CN Loss & Transfer#
    C_leaveMOM_pot = MOM.C*kMOM 
    C_MOM_NOM_pot  = C_leaveMOM_pot
    N_leaveMOM_pot = MOM.N*kMOM
    N_MOM_NOM_pot  = N_leaveMOM_pot
    #N Coupling:Immobilization/Mineralization#
    #Determing N immobilization potential in MOM Allocation Overall--bimm_MOM#
    N_need = C_MOM_NOM_pot/CN_NOM   #Actual Nitrogen demand from MOM to NOM#
    if N_need > N_leaveMOM_pot:                          
        bimm_MOM = 1 #Immobilization#
    else:
        bimm_MOM = 0 #Mineralization#
    #Updating Potebtial Immobilization#
    potimm = potimm+(N_need-N_leaveMOM_pot) #Potential N Immobilization Content along Decomposition#
   
    ###Step3: AOM1 allocating to SMB1,SMB2###
    kAOM1 = 1/kmax_AOM1/365*t_scale*w_scale #Decay rate of AOM1 pool#
    #Decomposition_Tentative:CN Loss & Transfer#
    C_leaveAOM1_pot    = AOM1.C*kAOM1
    C_AOM1_SMB1_pot    = C_leaveAOM1_pot*fAOM1_SMB1
    C_AOM1_SMB2_pot    = C_leaveAOM1_pot*(1-fAOM1_SMB1)
    C_AOM1_SMB1_E_pot  = C_AOM1_SMB1_pot*fE_AOM1
    C_AOM1_SMB2_E_pot  = C_AOM1_SMB2_pot*fE_AOM1
    C_AOM1_SMB1_rh_pot = C_AOM1_SMB1_pot*(1-fE_AOM1)
    C_AOM1_SMB2_rh_pot = C_AOM1_SMB2_pot*(1-fE_AOM1)
    N_leaveAOM1_pot    = AOM1.N*kAOM1
    N_AOM1_SMB1_pot    = N_leaveAOM1_pot*fAOM1_SMB1
    N_AOM1_SMB2_pot    = N_leaveAOM1_pot*(1-fAOM1_SMB1)
    #N Coupling:Immobilization/Mineralization#
    #Determing N immobilization potential in AOM1 Allocation Overall--bimm_AOM1#
    N_need = (C_AOM1_SMB1_E_pot/CN_SMB1)+(C_AOM1_SMB2_E_pot/CN_SMB2) #Actual Nitrogen demand from AOM1 to SMB#
    if N_need > N_leaveAOM1_pot: #Actual Nitrogen demand from MOM to NOM#
        bimm_AOM1 = 1   # 1: immobilization
    else:                        
        bimm_AOM1 = 0   # 0: mineralization
    #Updating Potebtial Immobilization#
    potimm = potimm+(N_need-N_leaveAOM1_pot) #Potential N Immobilization Content along Decomposition#

    ###Step4: AOM2 Allocating to SMB1,SMB2###
    kAOM2 = 1/kmax_AOM2/365*t_scale*w_scale * math.exp(-3*AOM3.C/(AOM2.C+AOM3.C))            ###QQQ: 20210908: add lignin protection mechanism #Decay rate of AOM2 pool#
    #Decomposition_Tentative:CN Loss & Transfer#
    C_leaveAOM2_pot    = AOM2.C*kAOM2
    C_AOM2_SMB1_pot    = C_leaveAOM2_pot*fAOM2_SMB1
    C_AOM2_SMB2_pot    = C_leaveAOM2_pot*(1-fAOM2_SMB1)
    C_AOM2_SMB1_E_pot  = C_AOM2_SMB1_pot*fE_AOM2
    C_AOM2_SMB2_E_pot  = C_AOM2_SMB2_pot*fE_AOM2
    C_AOM2_SMB1_rh_pot = C_AOM2_SMB1_pot*(1-fE_AOM2)
    C_AOM2_SMB2_rh_pot = C_AOM2_SMB2_pot*(1-fE_AOM2)
    N_leaveAOM2_pot    = AOM2.N*kAOM2
    N_AOM2_SMB1_pot    = N_leaveAOM2_pot*fAOM2_SMB1
    N_AOM2_SMB2_pot    = N_leaveAOM2_pot*(1-fAOM2_SMB1)
    #N Coupling:Immobilization/Mineralization#
    #Determing N immobilization potential in AOM2 Allocation Overall--bimm_AOM2#
    N_need = (C_AOM2_SMB1_E_pot/CN_SMB1)+(C_AOM2_SMB2_E_pot/CN_SMB2) #Actual Nitrogen demand from AOM2 to SMB#
    if N_need > N_leaveAOM2_pot:
        bimm_AOM2 = 1   # 1: immobilization
    else:                        
        bimm_AOM2 = 0   # 0: mineralization
    #Updating Potebtial Immobilization#
    potimm = potimm+(N_need-N_leaveAOM2_pot) #Potential N Immobilization Content along Decomposition#
    
    ###Step5: AOM3 Allocating to SMB1,SMB2###        
    kAOM3 = 1/kmax_AOM3/365*t_scale*w_scale                  #Decay rate of AOM3 pool#
    #Decomposition_Tentative:CN Loss & Transfer#
    C_leaveAOM3_pot    = AOM3.C*kAOM3
    ###QQQ: 20210606: For Lab Version 
    C_AOM3_MOM_pot     = C_leaveAOM3_pot*(fAOM_MOM)           ###QQQ: 20210606: For Lab Version        
    C_AOM3_SMB1_pot    = C_leaveAOM3_pot*(1-fAOM_MOM)*fAOM3_SMB1
    C_AOM3_SMB2_pot    = C_leaveAOM3_pot*(1-fAOM_MOM)*(1-fAOM3_SMB1)
    C_AOM3_SMB1_E_pot  = C_AOM3_SMB1_pot*fE_AOM3
    C_AOM3_SMB2_E_pot  = C_AOM3_SMB2_pot*fE_AOM3
    C_AOM3_SMB1_rh_pot = C_AOM3_SMB1_pot*(1-fE_AOM3)
    C_AOM3_SMB2_rh_pot = C_AOM3_SMB2_pot*(1-fE_AOM3)
    N_leaveAOM3_pot    = AOM3.N*kAOM3
    N_AOM3_MOM_pot     = N_leaveAOM3_pot*(fAOM_MOM)          ###QQQ: 20210606: For Lab Version 
    N_AOM3_SMB1_pot    = N_leaveAOM3_pot*(1-fAOM_MOM)*fAOM3_SMB1
    N_AOM3_SMB2_pot    = N_leaveAOM3_pot*(1-fAOM_MOM)*(1-fAOM3_SMB1)
    #N Coupling:Immobilization/Mineralization#
    #Determing N immobilization potential in AOM3 Allocation Overall--bimm_AOM3#
    N_need = N_AOM3_MOM_pot/CN_MOM+(C_AOM3_SMB1_E_pot)/CN_SMB1+(C_AOM3_SMB2_E_pot)/CN_SMB2  ###QQQ: 20210606: For Lab Version  #Actual Nitrogen demand from AOM3 to SMB#
    if N_need > N_leaveAOM3_pot:
        bimm_AOM3 = 1   # 1: immobilization
    else:                        
        bimm_AOM3 = 0   # 0: mineralization
    #Updating Potebtial Immobilization# 
    potimm = potimm+(N_need-N_leaveAOM3_pot) #Potential N Immobilization Content along Decomposition#
    
    #Carbon loss content from Litter Pool#            ###QQQ: 20210601: For lab incubation#
    initCont.LitterSub(C_leaveAOM1_pot+C_leaveAOM2_pot+C_leaveAOM3_pot)
    if t==1: #Initial Carbon Loss Content from Plant Residue(1st)#
        initCont.LitterZero(C_leaveAOM1_pot+C_leaveAOM2_pot+C_leaveAOM3_pot)   #Litter Substrate in the First Day#
        
    ###Step6: SMB2 Allocating to SMR###
    ##Situation 1: Soil Decomposition##
    if initCont.litterzero == 0:
        micro_CP   = 0.5
        death_rate = 0.1 
        #Determing Decomposition Period# 
        if flag_SMB2 == 0: #Microbial anundance Increase Period#
            fM_SMB2 = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
            fM_SMB2 = min(0.8,fM_SMB2)
            death_rate = death_rate
        if flag_SMB2 == 1: #Microbial anundance Decline Period#
            fM_SMB2 = fM_SMB2*(SMB2.C/(micro_CP*initCont.microSMB1Zero))
            fM_SMB2 = min(0.8,fM_SMB2)
            death_rate = (1+death_rate)
        ##Microbial Activities & CN Transfer in Decomposition##    
        C_SMB2_A_pot  = SMB2.CU*(fA_SMB2)                   #Assimilation
        C_SMB2_E_pot  = SMB2.CU*(1-fA_SMB2)                 #Egestion
        C_SMB2_rh_pot = SMB2.CU*(fA_SMB2*fM_SMB2)           #Respiration
        C_SMB2_P_pot  = SMB2.CU*(fA_SMB2*(1-fM_SMB2))       #Procreation
        C_SMB2_D_pot  = C_SMB2_P_pot*death_rate             #Death
        N_SMB2_A_pot  = SMB2.NU*(fA_SMB2)
        N_SMB2_E_pot  = SMB2.NU*((1-fA_SMB2)+(fA_SMB2*fM_SMB2))
        N_SMB2_rh_pot = SMB2.NU*0
        N_SMB2_P_pot  = SMB2.NU*(fA_SMB2*(1-fM_SMB2))
        N_SMB2_D_pot  = N_SMB2_P_pot*death_rate
    ##Situation 2: Litter Decomposition##
    if initCont.litterzero != 0:
        substrate_CP = 0.75            #Substrate_CP(hyper-parameter)
        micro_CP     = 0.5             #Pool Size_CP(hyper-parameter)
        death_rate   = 0.1            #Death rate(hyper-parameter)
        #Determing Decomposition Period#
        if flag_SMB2 == 0: #Microbial anundance Increase Period#
            #if initCont.littersub >= (substrate_CP*initCont.litterzero):
                #fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
                #death_rate   = (death_rate/10)
            #if initCont.littersub < (substrate_CP*initCont.litterzero):
                #fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
                #death_rate   = death_rate
            fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
            death_rate   = (death_rate)
        if flag_SMB2 == 1: #Microbial anundance Decline Period#
            if SMB2.C >= initCont.microSMB2Zero:
                fM_SMB2_flag = fM_SMB2
                death_rate   = (1+death_rate)
            if SMB2.C < initCont.microSMB2Zero:
                fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
                death_rate   = (1+death_rate/10)             
        ##Microbial Activity & CN Transfer in Decomposition##
        C_SMB2_A_pot  = SMB2.CU*(fA_SMB2)                     #Assimilation
        C_SMB2_E_pot  = SMB2.CU*(1-fA_SMB2)                   #Egestion
        C_SMB2_rh_pot = SMB2.CU*(fA_SMB2*fM_SMB2_flag)        #Respiration
        C_SMB2_P_pot  = SMB2.CU*(fA_SMB2*(1-fM_SMB2_flag))    #Procreation
        C_SMB2_D_pot  = C_SMB2_P_pot*death_rate               #Death
        N_SMB2_A_pot  = SMB2.NU*(fA_SMB2)
        N_SMB2_E_pot  = SMB2.NU*((1-fA_SMB2)+(fA_SMB2*fM_SMB2_flag))
        N_SMB2_rh_pot = SMB2.NU*0
        N_SMB2_P_pot  = SMB2.NU*(fA_SMB2*(1-fM_SMB2_flag))
        N_SMB2_D_pot  = N_SMB2_P_pot*death_rate
    #Determing N immobilization potential in SMB2 Allocation Overall--bimm_SMB2#
    #N Coupling:Immobilization/Mineralization#
    C_SMB2_SMR_pot = (C_SMB2_E_pot+C_SMB2_D_pot)
    N_SMB2_SMR_pot = (N_SMB2_E_pot+N_SMB2_D_pot)
    N_need = C_SMB2_SMR_pot/CN_SMR  #Actual Nitrogen demand from SMB2 to SMR#
    if N_need > N_SMB2_SMR_pot:
        bimm_SMB2 = 1     # 1: immobilization
    else:
        bimm_SMB2 = 0     # 0: mineralization
    #Updating Potebtial Immobilization# 
    potimm = potimm+(N_need-N_SMB2_SMR_pot)  #Potential N Immobilization Content along Decomposition#
            
    ###Step7: SMB1 Allocating to SMR###
    ##Situation 1: Soil Decomposition##
    if initCont.litterzero == 0:
        micro_CP = 1.1
        death_rate = 0.1
        #Determing Decomposition Period# 
        if flag_SMB1 == 0:  #Microbial anundance Increase Period#
            fM_SMB1 = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
            fM_SMB1 = min(0.8,fM_SMB1)
            death_rate = (death_rate)
        if flag_SMB1 == 1:  #Microbial anundance Decline Period#
            fM_SMB1 = fM_SMB1*(SMB1.C/(micro_CP*initCont.microSMB1Zero))
            fM_SMB1 = min(0.8,fM_SMB1)
            death_rate = (1+death_rate)
        ##Microbial Activity in Decomposition##
        C_SMB1_A_pot  = SMB1.CU*(fA_SMB1)                     #Assimilation
        C_SMB1_E_pot  = SMB1.CU*(1-fA_SMB1)                   #Egestion
        C_SMB1_rh_pot = SMB1.CU*(fA_SMB1*fM_SMB1)             #Respiration
        C_SMB1_P_pot  = SMB1.CU*(fA_SMB1*(1-fM_SMB1))         #Procreation
        C_SMB1_D_pot  = C_SMB1_P_pot*death_rate               #Death
        N_SMB1_A_pot  = SMB1.NU*(fA_SMB1)
        N_SMB1_E_pot  = SMB1.NU*((1-fA_SMB1)+(fA_SMB1*fM_SMB1))
        N_SMB1_rh_pot = SMB1.NU*0
        N_SMB1_P_pot  = SMB1.NU*(fA_SMB1*(1-fM_SMB1))       
        N_SMB1_D_pot  = N_SMB1_P_pot*death_rate
    ##Situation 2: Litter Decomposition##    
    if initCont.litterzero != 0:
        micro_CP = 1.25   #Pool Size#
        death_rate = 0.1  
        #Determing Decomposition Period# 
        if flag_SMB1 == 0: #Microbial anundance Increase Period#
            fM_SMB1_flag = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
            death_rate = death_rate
        if flag_SMB1 == 1: #Microbial anundance Decline Period#
            fM_SMB1_flag = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
            death_rate = (1+death_rate)
        ##Microbial Activity & CN Transfer in Decomposition##
        C_SMB1_A_pot  = SMB1.CU*(fA_SMB1)                     #Assimilation
        C_SMB1_E_pot  = SMB1.CU*(1-fA_SMB1)                   #Egestion     
        C_SMB1_rh_pot = SMB1.CU*(fA_SMB1*fM_SMB1_flag)        #Respiration
        C_SMB1_P_pot  = SMB1.CU*(fA_SMB1*(1-fM_SMB1_flag))    #Procreation
        C_SMB1_D_pot  = C_SMB1_P_pot*death_rate               #Death
        N_SMB1_A_pot  = SMB1.NU*(fA_SMB1)
        N_SMB1_E_pot  = SMB1.NU*((1-fA_SMB1)+(fA_SMB1*fM_SMB1_flag))
        N_SMB1_rh_pot = SMB1.NU*0
        N_SMB1_P_pot  = SMB1.NU*(fA_SMB1*(1-fM_SMB1_flag))
        N_SMB1_D_pot  = N_SMB1_P_pot*death_rate
    #Determing N immobilization potential in SMB1--bimm_SMB1#
    #N Coupling:Immobilization/Mineralization#
    C_SMB1_SMR_pot = (C_SMB1_E_pot+C_SMB1_D_pot)
    N_SMB1_SMR_pot = (N_SMB1_E_pot+N_SMB1_D_pot)
    N_need = C_SMB1_SMR_pot/CN_SMR  #Actual Nitrogen demand from SMB1 to SMR#
    if N_need > N_SMB1_SMR_pot:                          
        bimm_SMB1 = 1  # 1: immobilization
    else:
        bimm_SMB1 = 0  # 0: mineralization
    #Updating Potebtial Immobilization# 
    potimm = potimm+(N_need-N_SMB1_SMR_pot)   #Potential N Immobilization Content along Decomposition#
          
    ###Step8: SMR Allocating to ASOM###
    #Decomposition_Tentative:CN Loss & Transfer#
    kSMR = 1/kmax_SMR/365*t_scale*w_scale   #Decay Rate#
    C_leaveSMR_pot = SMR.C*kSMR
    C_SMR_NOM_pot  = C_leaveSMR_pot
    N_leaveSMR_pot = SMR.N*kSMR
    N_SMR_NOM_pot  = N_leaveSMR_pot
    #Determing N immobilization potential in SMR--bimm_SMR#
    #N coupling: Immobilization/Mineralization#
    N_need = C_SMR_NOM_pot/CN_NOM   #Actucal N demand from SMR to NOM#
    if N_need > N_leaveSMR_pot:                               
        bimm_SMR = 1   #Immobilization#
    else:                        
        bimm_SMR = 0   #Mineralization#
    #Updating Potebtial Immobilization#
    potimm = potimm + (N_need-N_leaveSMR_pot)    #Potential N Immobilization Content along Decomposition#
    
    return (bimm_AOM1,bimm_AOM2,bimm_AOM3,bimm_SMR,bimm_NOM,bimm_MOM,potimm)


# In[10]:


### Decomposition - Actual decomposition
def SoilDecom_real (sample_site,t,obser_data,Nmin,potimm,flag_SMB1,flag_SMB2,Rh,lig_p1,lig_p2,     #Scalar(7)#  ###QQQ: 20210908: add 2 logistic lignin para to control lignin decomposition ###
                    AOM1,AOM2,AOM3,SMB1,SMB2,SMR,NOM,MOM,avN,initCont,litligsoil_CO2,#Class(10)#
                    bimm_AOM1,bimm_AOM2,bimm_AOM3,bimm_SMR,bimm_NOM,bimm_MOM,        #N scalar(6)#
                    CN_SMB1,CN_SMB2,CN_SMR,CN_NOM,CN_MOM,                            #CN Ratio#
                    kmax_AOM1,kmax_AOM2,kmax_AOM3,kmax_SMR,kmax_NOM,kmax_MOM,        #Decay Rate(6)#
                    fAOM1_SMB1,fAOM2_SMB1,fAOM_MOM,fAOM3_SMB1,fNOM_SMB1,fNOM_MOM,    #Fraction(6)#
                    fE_AOM1,fE_AOM2,fE_AOM3,fE_NOM,fA_SMB1,fA_SMB2,fM_SMB1,fM_SMB2): #Microbe(8)#
    
    ##Environmental Factors##
    #EnVirFactor 1: temp# 
    temp    = obser_data.loc[sample_site,"temp"]   #temp(unit celsius degree )#
    t_scale = 4.89*np.exp(-3.432+0.1*temp*(1-0.5*temp/36.9))
    #EnVirFactor 2: water# 
    theta    = obser_data.loc[sample_site,"water_cont"]   #volumetric water content(unit:%)#
    theta_fc  = obser_data.loc[sample_site,"water_fc"]    #Field water content#
    theta_sat = obser_data.loc[sample_site,"water_sat"]   #Saturated water content#
    if theta <= theta_fc:
        w_scale = (1-np.exp(-theta/theta_sat))/(1-np.exp(-theta_fc/theta_sat))
    elif theta > theta_fc:
        w_scale = 1.0044-(0.0044/(np.exp(-5*(((theta/theta_sat)-(theta_fc/theta_sat))/(1-theta_fc/theta_sat)))))
    #EnVirFactor 3: clay# 
    Pclay = obser_data.loc[sample_site,"Pclay"]   #Clay Proportion:Unit(percent without %)#
    clay_scale = 1-(0.75*Pclay/100)
    #EnVirFactor 4: Nitrogen Scalar-Modify the k value#
    if initCont.litterzero != 0:
        avn_opt = (initCont.littersub/25 + initCont.soilsub/25)    #Optimal Nitrogen = 25*litter_loseC (Unit:g N)# #avn_opt = 2.0#
    if initCont.litterzero == 0:
        avn_opt = (initCont.soilsub/25)   #Optimal Nitrogen = 25*litter_loseC (Unit:g N)# #avn_opt = 2.0#
    #Nitrogen Scalar1-Immobilization scaler#
    if (potimm>0):   #Immobilization in Decomposition#
        ni_scale = avN.avn()/potimm   #f(N)=avn/Nimm# (potimm=Nimm)
    else:                         
        ni_scale = 1
    ni_scale = min(1,max(0.5,ni_scale))          ###Q: Change 0.1 to 0.5###
    #ni_scale = fmin(1., fmax(0.1, ni_scale))#   ###Q
    #Nitrogen Scalar2-Mineralization scaler#
    if avN.avn() >= 1.5*avn_opt:   #Mineralization Occured when Nitrogen is over Optimal Content#
        nm_scale = 1-(avN.avn()-avn_opt)/avn_opt   #avn,avN.NH4,avN.NO3:(Unit:g N)#(gN/m2/day)
    elif avn_opt/2 <= avN.avn() < 1.5*avn_opt:   #Nitrogen is in the Optimal Interval#
        nm_scale = 1
    elif avN.avn() < avn_opt/2:   #Mineralization Occured when Nitrogen is not enough#
        nm_scale = 1+(0.5*avn_opt-avN.avn())/avn_opt
    nm_scale = min(1.5,max(0.8,nm_scale))       ###Q
    #if(avn > avn_opt) nm_scale = 1. - (avn - avn_opt)/avn_opt;
    #else if(avn <= avn_opt && avn >=0.5*avn_opt) nm_scale = 1.;
    #else nm_scale = 1. + (0.5*avn_opt - avn)/avn_opt;
    #nm_scale = fmin(1.5, fmax(0.8, nm_scale));
    #Nitrification Process#
    avN.nitri(0.1)   #1% NH4 is nitrified into NO3 Every Day#
    
    ###STEP1: ASOM(Soil Available Carbon Pool) Allocating to PSOM/SMB1/SMB2###
    #Decay Rate#
    if bimm_NOM == 1:
        kNOM = 1/kmax_NOM/365*t_scale*w_scale*ni_scale
    else:
        kNOM = 1/kmax_NOM/365*t_scale*w_scale*nm_scale
    ##Decomposition: C&N Loss & Transfer##
    C_leaveNOM    = NOM.C*kNOM                               #Carbon Leave from NOM Pool#
    C_NOM_MOM     = C_leaveNOM*fNOM_MOM                      #Carbon transfer from NOM to MOM#
    C_NOM_SMB1    = C_leaveNOM*(1-fNOM_MOM)*fNOM_SMB1        #Carbon transfer from NOM to SMB1#
    C_NOM_SMB2    = C_leaveNOM*(1-fNOM_MOM)*(1-fNOM_SMB1)    #Carbon transfer from NOM to SMB2#
    C_NOM_SMB1_E  = C_NOM_SMB1*fE_NOM                        #Carbon Used by SMB1#
    C_NOM_SMB2_E  = C_NOM_SMB2*fE_NOM                        #Carbon Used by SMB2#
    C_NOM_SMB1_rh = C_NOM_SMB1*(1-fE_NOM)                    #Carbon directly loss as CO2#
    C_NOM_SMB2_rh = C_NOM_SMB2*(1-fE_NOM)                    #Carbon directly loss as CO2#
    N_leaveNOM    = NOM.N*kNOM                               #Nitrogen leave from NOM pool#
    N_NOM_MOM     = N_leaveNOM*fNOM_MOM                      #Nitrogen transfer from NOM to MOM#
    N_NOM_SMB1    = N_leaveNOM*(1-fNOM_MOM)*fNOM_SMB1        #Nitrogen transfer from NOM to AOM1#
    N_NOM_SMB2    = N_leaveNOM*(1-fNOM_MOM)*(1-fNOM_SMB1)    #Nitrogen transfer from NOM to AOM2#
    #Substrate content(Without litter Input): CN loss from Soil#
    initCont.SoilSub(C_leaveNOM)
    if t==1:
        initCont.SoilZero(C_leaveNOM)   #Soil Substrate in the First Day#
    ##Updating Upstream(ASOM) Pool##
    NOM.SubCarbon(C_leaveNOM)
    NOM.SubNitrogen(N_leaveNOM)
    Rh.append(C_NOM_SMB1_rh)
    Rh.append(C_NOM_SMB2_rh)
    litligsoil_CO2.soilCO2(C_NOM_SMB1_rh + C_NOM_SMB2_rh)  ###Q: litter-lignin-soil seperated CO2
    litligsoil_CO2.soilToMB1(C_NOM_SMB1_E)   ###Q
    litligsoil_CO2.soilToMB2(C_NOM_SMB2_E)   ###Q
    
    #N Coupling:Immobilization/Mineralization#
    #From ASOM to PSOM#
    if N_NOM_MOM < C_NOM_MOM/CN_MOM:  #Immobilization#
        Nimmobilize = C_NOM_MOM/CN_MOM-N_NOM_MOM
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream(MOM) Pool#
        MOM.AddCarbon(C_NOM_MOM)
        MOM.AddNitrogen(N_NOM_MOM+Nimmobilize)
    elif N_NOM_MOM >= C_NOM_MOM/CN_MOM: #Mineralization#
        Nmineralize = N_NOM_MOM-C_NOM_MOM/CN_MOM
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream(MOM) Pool#
        MOM.AddCarbon(C_NOM_MOM)
        MOM.AddNitrogen(C_NOM_MOM/CN_MOM)
    #From ASOM to SMB1#
    if N_NOM_SMB1 < C_NOM_SMB1_E/CN_SMB1:   #Immobilization#
        Nimmobilize = C_NOM_SMB1_E/CN_SMB1-N_NOM_SMB1
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream(SMB1) Pool#
        SMB1.AddNewCarbon(C_NOM_SMB1_E)
        SMB1.AddNewNitrogen(N_NOM_SMB1+Nimmobilize)
    elif N_NOM_SMB1 >= C_NOM_SMB1_E/CN_SMB1:   #Mineralization#
        Nmineralize = N_NOM_SMB1-C_NOM_SMB1_E/CN_SMB1
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream(SMB1) Pool#
        SMB1.AddNewCarbon(C_NOM_SMB1_E)
        SMB1.AddNewNitrogen(C_NOM_SMB1_E/CN_SMB1)
    #From ASOM to SMB2#
    if N_NOM_SMB2 < C_NOM_SMB2_E/CN_SMB2:   #Immobilization#
        Nimmobilize = C_NOM_SMB2_E/CN_SMB2-N_NOM_SMB2
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream(SMB2) Pool#
        SMB2.AddNewCarbon(C_NOM_SMB2_E)
        SMB2.AddNewNitrogen(N_NOM_SMB2+Nimmobilize)
    elif N_NOM_SMB2 >= C_NOM_SMB2_E/CN_SMB2:   #Mineralization#
        Nmineralize = N_NOM_SMB2-(C_NOM_SMB2_E/CN_SMB2)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream(SMB2) Pool#
        SMB2.AddNewCarbon(C_NOM_SMB2_E)
        SMB2.AddNewNitrogen(C_NOM_SMB2_E/CN_SMB2)
        
    ###STEP2: PSOM Allocating to ASOM###
    #Decay Rate#
    if bimm_MOM == 1:
        kMOM = 1/kmax_MOM/365*t_scale*w_scale*ni_scale*clay_scale
    else:
        kMOM = 1/kmax_MOM/365*t_scale*w_scale*nm_scale*clay_scale
    ##Decomposition: CN Loss & Transfer##
    C_leaveMOM = MOM.C*kMOM   #Carbon leave from MOM Pool#
    C_MOM_NOM  = C_leaveMOM   #Carbon transfer from MOM to NOM#
    N_leaveMOM = MOM.N*kMOM   #Nitrogen Leave from MOM Pool#
    N_MOM_NOM  = N_leaveMOM   #Nitrogen leave from MOM Pool#
    #Updating Upstream(MOM) Pool#
    MOM.SubCarbon(C_leaveMOM)
    MOM.SubNitrogen(N_leaveMOM)
    #N Coupling:Immobilization/Mineralization#
    #From PSOM to ASOM#
    if N_MOM_NOM < C_MOM_NOM/CN_NOM: #Immobilization#
        Nimmobilize = (C_MOM_NOM/CN_NOM)-N_MOM_NOM
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        NOM.AddCarbon(C_MOM_NOM)
        NOM.AddNitrogen(N_MOM_NOM+Nimmobilize)
    elif N_MOM_NOM >= C_MOM_NOM/CN_NOM: #Mineralization#
        Nmineralize = N_MOM_NOM-(C_MOM_NOM/CN_NOM)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        NOM.AddCarbon(C_MOM_NOM)
        NOM.AddNitrogen(C_MOM_NOM/CN_NOM)
    
    ###STEP3: AOM1 Allocating to SMB1,SMB2###
    #Decay Rate#
    if bimm_AOM1 == 1:
        kAOM1 = 1/kmax_AOM1/365*t_scale*w_scale*ni_scale
    else:
        kAOM1 = 1/kmax_AOM1/365*t_scale*w_scale*nm_scale
    ##Decomposition: CN Loss & Transfer##
    C_leaveAOM1    = AOM1.C*kAOM1                #Carbon leave from AOM1#
    C_AOM1_SMB1    = C_leaveAOM1*fAOM1_SMB1      #Carbon transfer from AOM1 to SMB1#
    C_AOM1_SMB2    = C_leaveAOM1*(1-fAOM1_SMB1)  #Carbon transfer from AOM1 to SMB2#
    C_AOM1_SMB1_E  = C_AOM1_SMB1*fE_AOM1         #Carbon used by SMB1#
    C_AOM1_SMB2_E  = C_AOM1_SMB2*fE_AOM1         #Carbon used by SMB2#
    C_AOM1_SMB1_rh = C_AOM1_SMB1*(1-fE_AOM1)     #Carbon directly loss as CO2#
    C_AOM1_SMB2_rh = C_AOM1_SMB2*(1-fE_AOM1)     #Carbon directly loss as CO2#
    N_leaveAOM1    = AOM1.N*kAOM1                #Nitrogen leave from AOM1#
    N_AOM1_SMB1    = N_leaveAOM1*fAOM1_SMB1      #Nitrogen transfer into SMB1#
    N_AOM1_SMB2    = N_leaveAOM1*(1-fAOM1_SMB1)  #Nitrogen transfer into SMB2#
    #Updating Upstream(MOM) Pool#
    AOM1.SubCarbon(C_leaveAOM1)
    AOM1.SubNitrogen(N_leaveAOM1)
    Rh.append(C_AOM1_SMB1_rh)
    Rh.append(C_AOM1_SMB2_rh)
    litligsoil_CO2.litterCO2(C_AOM1_SMB1_rh + C_AOM1_SMB2_rh)    ###Q: litter-lignin-soil seperation 
    litligsoil_CO2.litterToMB1(C_AOM1_SMB1_E)   ###Q
    litligsoil_CO2.litterToMB2(C_AOM1_SMB2_E)   ###Q
    
    #N Coupling:Immobilization/Mineralization#
    #From AOM1 to SMB1#
    if N_AOM1_SMB1 < C_AOM1_SMB1_E/CN_SMB1:  #Immobilization#
        Nimmobilize = (C_AOM1_SMB1_E/CN_SMB1)-N_AOM1_SMB1
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        SMB1.AddNewCarbon(C_AOM1_SMB1_E)
        SMB1.AddNewNitrogen(N_AOM1_SMB1+Nimmobilize)
    elif N_AOM1_SMB1 >= C_AOM1_SMB1_E/CN_SMB1:  #Mineralization#
        Nmineralize = N_AOM1_SMB1-(C_AOM1_SMB1_E/CN_SMB1)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        SMB1.AddNewCarbon(C_AOM1_SMB1_E)
        SMB1.AddNewNitrogen(C_AOM1_SMB1_E/CN_SMB1)
    #From AOM1 to SMB2#
    if N_AOM1_SMB2 < C_AOM1_SMB2_E/CN_SMB2:  #Immobilization#
        Nimmobilize = (C_AOM1_SMB2_E/CN_SMB2)-N_AOM1_SMB2
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        SMB2.AddNewCarbon(C_AOM1_SMB2_E)
        SMB2.AddNewNitrogen(N_AOM1_SMB2+Nimmobilize)
    elif N_AOM1_SMB2 >= C_AOM1_SMB2_E/CN_SMB2:  #Mineralization#
        Nmineralize = N_AOM1_SMB2-(C_AOM1_SMB2_E/CN_SMB2)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        SMB2.AddNewCarbon(C_AOM1_SMB2_E)
        SMB2.AddNewNitrogen(C_AOM1_SMB2_E/CN_SMB2)

    ###STEP4: AOM2 Allocating to SMB1,SMB2###
    #Decay Rate#
    if bimm_AOM2 == 1:
        kAOM2 = 1/kmax_AOM2/365*t_scale*w_scale*ni_scale * math.exp(-3 * (AOM3.C/(AOM2.C+AOM3.C+AOM1.C)))   ###QQQ: 20210908: lignin protection effect mechanism
        #kAOM2 = 1/kmax_AOM2/365*t_scale*w_scale*ni_scale     ###QQQ: 20210908: lignin protection ettect mechanism
    else:
        kAOM2 = 1/kmax_AOM2/365*t_scale*w_scale*nm_scale * math.exp(-3* (AOM3.C/(AOM2.C+AOM3.C+AOM1.C)))    ###QQQ: 20210908: lignin protection effect mechanism
        # kAOM2 = 1/kmax_AOM2/365*t_scale*w_scale*ni_scale     ###QQQ: 20210908: lignin protection effect mechanism
    ##Decomposition: CN Loss & Transfer##
    C_leaveAOM2    = AOM2.C*kAOM2                    #Carbon leave from AOM2#
    C_AOM2_SMB1    = C_leaveAOM2*fAOM2_SMB1          #Carbon transfer from AOM2 into SMB1#
    C_AOM2_SMB2    = C_leaveAOM2*(1-fAOM2_SMB1)      #Carbon transfer from AOM2 into SMB2#
    C_AOM2_SMB1_E  = C_AOM2_SMB1*fE_AOM2             #Carbon used by SMB1#
    C_AOM2_SMB2_E  = C_AOM2_SMB2*fE_AOM2             #Carbon used by SMB2#
    C_AOM2_SMB1_rh = C_AOM2_SMB1*(1-fE_AOM2)         #Carbon directly loss as CO2#
    C_AOM2_SMB2_rh = C_AOM2_SMB2*(1-fE_AOM2)         #Carbon directly loss as CO2#
    N_leaveAOM2    = AOM2.N*kAOM2                    #Nitrogen leave from AOM2#
    N_AOM2_SMB1    = N_leaveAOM2*fAOM2_SMB1          #Nitrogen transfer from AOM2 into SMB1#
    N_AOM2_SMB2    = N_leaveAOM2*(1-fAOM2_SMB1)      #Nitrogen transfer from AOM2 into SMB2#
    ##Updating Upstream/AOM2 Pool##
    AOM2.SubCarbon(C_leaveAOM2)
    AOM2.SubNitrogen(N_leaveAOM2)
    Rh.append(C_AOM2_SMB1_rh)
    Rh.append(C_AOM2_SMB2_rh)
    litligsoil_CO2.litterCO2(C_AOM2_SMB1_rh + C_AOM2_SMB2_rh)   ###Q: litter-lignin-soil: CO2 seperation
    litligsoil_CO2.litterToMB1(C_AOM2_SMB1_E)   ###Q
    litligsoil_CO2.litterToMB2(C_AOM2_SMB2_E)   ###Q
    
    #N Coupling:Immobilization/Mineralization#
    #From AOM2 to SMB1#
    if N_AOM2_SMB1 < C_AOM2_SMB1_E/CN_SMB1:  #Immobilization#
        Nimmobilize = (C_AOM2_SMB1_E/CN_SMB1)-N_AOM2_SMB1
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin - Nimmobilize
        #Updating Downstream Pool#
        SMB1.AddNewCarbon(C_AOM2_SMB1_E)
        SMB1.AddNewNitrogen(N_AOM2_SMB1+Nimmobilize)
    elif N_AOM2_SMB1 >= C_AOM2_SMB1_E/CN_SMB1:  #Mineralization#
        Nmineralize = N_AOM2_SMB1-(C_AOM2_SMB1_E/CN_SMB1)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        SMB1.AddNewCarbon(C_AOM2_SMB1_E)
        SMB1.AddNewNitrogen(C_AOM2_SMB1_E/CN_SMB1)
    #From AOM2 to SMB2#
    if N_AOM2_SMB2 < C_AOM2_SMB2_E/CN_SMB2:  #Immobilization#
        Nimmobilize = (C_AOM2_SMB2_E/CN_SMB2)-N_AOM2_SMB2
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        SMB2.AddNewCarbon(C_AOM2_SMB2_E)
        SMB2.AddNewNitrogen(N_AOM2_SMB2+Nimmobilize)
    elif N_AOM2_SMB2 >= C_AOM2_SMB2_E/CN_SMB2:  #Mineralization#
        Nmineralize = N_AOM2_SMB2-(C_AOM2_SMB2_E/CN_SMB2)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        SMB2.AddNewCarbon(C_AOM2_SMB2_E)
        SMB2.AddNewNitrogen(C_AOM2_SMB2_E/CN_SMB2)

    ###STEP5: AOM3 Allocating to SMB1,SMB2###
    #Decay Rate#
    #lig_p1 = 4                # determine the maximun lignin   (rate 0.08 = 4) to (20)
    #lig_p2 = 0.9              # determine the slop of lignin    (0.9: vertical line - 0.99: linear line)
    if bimm_AOM3 == 1:                         
        kAOM3 = 1/kmax_AOM3/365*t_scale*w_scale*ni_scale * (lig_p1/(1+pow(lig_p2,-t))+0.75)         ###QQQ:20210908: lignin decay logistic function
    else:
        kAOM3 = 1/kmax_AOM3/365*t_scale*w_scale*nm_scale * (lig_p1/(1+pow(lig_p2,-t))+0.75)         ###QQQ:20210908: lignin decay logistic function
    ##Decomposition: CN Loss & Transfer##
    C_leaveAOM3    = AOM3.C*kAOM3                #Carbon leave from AOM3#
    ###QQQ: 20201607: For Lab Incubation
    C_AOM3_MOM     = C_leaveAOM3*fAOM_MOM        ###QQQ: 20210607: For lab incubation: Carbon transfer from AOM3 into MOM#
    C_AOM3_SMB1    = C_leaveAOM3*(1-fAOM_MOM)*fAOM3_SMB1      #Carbon transfer from AOM3 into SMB1#
    C_AOM3_SMB2    = C_leaveAOM3*(1-fAOM_MOM)*(1-fAOM3_SMB1)  #Carbon transfer from AOM3 into SMB2#
    C_AOM3_SMB1_E  = C_AOM3_SMB1*fE_AOM3         #Carbon used by SMB1#
    C_AOM3_SMB2_E  = C_AOM3_SMB2*fE_AOM3         #Carbon used by SMB2# 
    C_AOM3_rh      = C_leaveAOM3*(1-fAOM_MOM)*(1-fE_AOM3)     ###Q: 20210517: Simplify CO2 emission from AOM3 pool
    C_AOM3_SMB1_rh = C_AOM3_SMB1*(1-fE_AOM3)     #Carbon directly loss as CO2#
    C_AOM3_SMB2_rh = C_AOM3_SMB2*(1-fE_AOM3)     #Carbon directly loss as CO2#
    N_leaveAOM3    = AOM3.N*kAOM3                #Nitrogen loss from AOM3#
    N_AOM3_MOM     = AOM3.N*(1-fAOM_MOM)         ###QQQ: 20210607: For lab incubation
    N_AOM3_SMB1    = N_leaveAOM3*(1-fAOM_MOM)*fAOM3_SMB1      #Nitrogen transfer from AOM3 into SMB1#
    N_AOM3_SMB2    = N_leaveAOM3*(1-fAOM_MOM)*(1-fAOM3_SMB1)  #Nitrogen transfer from AOM3 into SMB2#
    ##Updating Upstream/AOM3 Pool##
    AOM3.SubCarbon(C_leaveAOM3)
    AOM3.SubNitrogen(N_leaveAOM3)
    Rh.append(C_AOM3_SMB1_rh)
    Rh.append(C_AOM3_SMB2_rh)
    #litligsoil_CO2.ligninCO2(C_AOM3_SMB1_rh+ C_AOM3_SMB2_rh)   ###Q: litter-lignin-soil CO2 seperation
    litligsoil_CO2.ligninCO2(C_AOM3_rh)                         ###Q: 20210517: Simplify CO2 emission from AOM3 pool
    litligsoil_CO2.ligninToMB1(C_AOM3_SMB1_E)   ###Q
    litligsoil_CO2.ligninToMB2(C_AOM3_SMB2_E)   ###Q
    
    #N Coupling:Immobilization/Mineralization#
    #From AOM3 to MOM# ###QQQ: 20210607: For lab incubation
    if N_AOM3_MOM < C_AOM3_MOM/CN_MOM:  #Immobilization#      ###QQQ: 20210607: For lab incubation
        Nimmobilize = (C_AOM3_MOM/CN_MOM)-N_AOM3_MOM
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        MOM.AddNewCarbon(C_AOM3_MOM)
        MOM.AddNewNitrogen(N_AOM3_MOM+Nimmobilize)
    elif N_AOM3_MOM >= C_AOM3_MOM/CN_MOM:  #Mineralization#
        Nmineralize = N_AOM3_MOM-(C_AOM3_MOM/CN_MOM)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        MOM.AddNewCarbon(C_AOM3_MOM)
        MOM.AddNewNitrogen(C_AOM3_MOM/CN_MOM)                ###QQQ: 20210607: For lab incubation
    #From AOM3 to SMB1#
    if N_AOM3_SMB1 < C_AOM3_SMB1_E/CN_SMB1:  #Immobilization#
        Nimmobilize = (C_AOM3_SMB1_E/CN_SMB1)-N_AOM3_SMB1
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        SMB1.AddNewCarbon(C_AOM3_SMB1_E)
        SMB1.AddNewNitrogen(N_AOM3_SMB1+Nimmobilize)
    elif N_AOM3_SMB1 >= C_AOM3_SMB1_E/CN_SMB1:  #Mineralization#
        Nmineralize = N_AOM3_SMB1-(C_AOM3_SMB1_E/CN_SMB1)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        SMB1.AddNewCarbon(C_AOM3_SMB1_E)
        SMB1.AddNewNitrogen(C_AOM3_SMB1_E/CN_SMB1)
    #From AOM3 to SMB2#
    if N_AOM3_SMB2 < C_AOM3_SMB2_E/CN_SMB2:  #Immobilization#
        Nimmobilize = (C_AOM3_SMB2_E/CN_SMB2)-N_AOM3_SMB2
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream Pool#
        SMB2.AddNewCarbon(C_AOM3_SMB2_E)
        SMB2.AddNewNitrogen(N_AOM3_SMB2+Nimmobilize)
    elif N_AOM3_SMB2 >= C_AOM3_SMB2_E/CN_SMB2:  #Mineralization#
        Nmineralize = N_AOM3_SMB2-(C_AOM3_SMB2_E/CN_SMB2)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream Pool#
        SMB2.AddNewCarbon(C_AOM3_SMB2_E)
        SMB2.AddNewNitrogen(C_AOM3_SMB2_E/CN_SMB2)
    
    #Carbon loss content from Litter Pool#            ###QQQ: 20210601: For lab incubation#
    initCont.LitterSub(C_leaveAOM1+C_leaveAOM2+C_leaveAOM3)
    if t==1: #Initial Carbon Loss Content from Plant Residue(1st)#
        initCont.LitterZero(C_leaveAOM1+C_leaveAOM2+C_leaveAOM3)   #Litter Substrate in the First Day#
    
    
    ###STEP6: SMB2 Allocating to SMR###
    ##Situation 1: Soil Decomposition##
    if initCont.litterzero == 0:
        micro_CP   = 0.5
        death_rate = 0.1
        #Determing Decomposition Period#   ###Q: Can we simplify this?
        if SMB2.C >= (micro_CP*initCont.microSMB1Zero):  #initCont.microSMB2 is a fixed value, representing SMB2 pool size at 0 day#
            flag_SMB2 = 1
        if flag_SMB2 == 0:  #Microbial anundance Increase Period#
            fM_SMB2 = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
            fM_SMB2 = min(0.75,fM_SMB2)
            death_rate = death_rate
        if flag_SMB2 == 1:  #Microbial anundance Decline Period#
            fM_SMB2 = fM_SMB2*(SMB2.C/(micro_CP*initCont.microSMB1Zero))
            fM_SMB2 = max(0.25,min(0.8,fM_SMB2))
            death_rate = (1+death_rate)
        ##Microbial Activities & CN Transfer in Decomposition##
        C_SMB2_A  = SMB2.CU*(fA_SMB2)                   #Assimilation
        C_SMB2_E  = SMB2.CU*(1-fA_SMB2)                 #Egestion
        C_SMB2_rh = SMB2.CU*(fA_SMB2*fM_SMB2)           #Respiration
        C_SMB2_P  = SMB2.CU*(fA_SMB2*(1-fM_SMB2))       #Procreation
        C_SMB2_D  = C_SMB2_P*death_rate                 #Death
        N_SMB2_A  = SMB2.NU*(fA_SMB2)
        N_SMB2_E  = SMB2.NU*((1-fA_SMB2)+(fA_SMB2*fM_SMB2))
        N_SMB2_rh = SMB2.NU*0
        N_SMB2_P  = SMB2.NU*(fA_SMB2*(1-fM_SMB2))
        N_SMB2_D  = N_SMB2_P*death_rate
    ##Situation 2: Litter Decomposition##
    if initCont.litterzero != 0:
        substrate_CP = 0.75            #Substrate_CP
        micro_CP     = 0.5             #Pool Size_CP
        death_rate   = 0.1
        #Determing Decomposition Period#
        if SMB2.C >= (micro_CP*initCont.microSMB1Zero):         #Critical point#
            flag_SMB2 = 1                                   #Decline Period#
        if flag_SMB2 == 0:  #Microbial anundance Increase Period#
            #if initCont.littersub >= (substrate_CP*initCont.litterzero):
                #fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
                #death_rate   = (death_rate/10)
            #if initCont.littersub < (substrate_CP*initCont.litterzero):
                #fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
                #death_rate   = death_rate
            fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
            death_rate   = death_rate
        if flag_SMB2 == 1:  #Microbial anundance Decline Period#
            if SMB2.C >= initCont.microSMB2Zero:  #Decline Period 1#
                fM_SMB2_flag = fM_SMB2
                death_rate   = (1+death_rate)
            if SMB2.C < initCont.microSMB2Zero:   #Decline Period 2#
                fM_SMB2_flag = fM_SMB2*(SMB2.C/initCont.microSMB2Zero)
                death_rate   = (1+death_rate/10)             
        ##Microbial Activity & CN Transfer in Decomposition##
        C_SMB2_A  = SMB2.CU*(fA_SMB2)                     #Assimilation
        C_SMB2_E  = SMB2.CU*(1-fA_SMB2)                   #Egestion
        C_SMB2_rh = SMB2.CU*(fA_SMB2*fM_SMB2_flag)        #Respiration
        C_SMB2_P  = SMB2.CU*(fA_SMB2*(1-fM_SMB2_flag))    #Procreation
        C_SMB2_D  = C_SMB2_P*death_rate                   #Death
        N_SMB2_A  = SMB2.NU*(fA_SMB2)
        N_SMB2_E  = SMB2.NU*((1-fA_SMB2)+(fA_SMB2*fM_SMB2_flag))
        N_SMB2_rh = SMB2.NU*0
        N_SMB2_P  = SMB2.NU*(fA_SMB2*(1-fM_SMB2_flag))
        N_SMB2_D  = N_SMB2_P*death_rate
    #Updating Upstream/SMB2 Pool#
    SMB2.SubNewCarbon(C_SMB2_E+C_SMB2_rh+C_SMB2_P)
    SMB2.SubNewNitrogen(N_SMB2_E+N_SMB2_rh+N_SMB2_P)
    SMB2.AddCarbon(C_SMB2_P)
    SMB2.AddNitrogen(N_SMB2_P)
    SMB2.SubCarbon(C_SMB2_D)
    SMB2.SubNitrogen(N_SMB2_D)
    Rh.append(C_SMB2_rh)
    litterToMB2_pct = litligsoil_CO2.litterC_MB2/(litligsoil_CO2.litterC_MB2+litligsoil_CO2.ligninC_MB2+litligsoil_CO2.soilC_MB2) ###Q
    ligninToMB2_pct = litligsoil_CO2.ligninC_MB2/(litligsoil_CO2.litterC_MB2+litligsoil_CO2.ligninC_MB2+litligsoil_CO2.soilC_MB2) ###Q
    soilToMB2_pct   = litligsoil_CO2.soilC_MB2/(litligsoil_CO2.litterC_MB2+litligsoil_CO2.ligninC_MB2+litligsoil_CO2.soilC_MB2)   ###Q
    litligsoil_CO2.litterCO2(C_SMB2_rh*litterToMB2_pct)    ###Q
    litligsoil_CO2.ligninCO2(C_SMB2_rh*ligninToMB2_pct)    ###Q:20210517 Change back to: litligsoil_CO2.ligninCO2(C_SMB2_rh*ligninToMB2_pct)
    litligsoil_CO2.soilCO2(C_SMB2_rh*soilToMB2_pct)        ###Q

    #N Coupling:Immobilization/Mineralization#
    #From SMB2 to SMR#
    C_SMB2_SMR = (C_SMB2_E+C_SMB2_D) #Carbon Transfer into SMB2#
    N_SMB2_SMR = (N_SMB2_E+N_SMB2_D) #Nitrogen Transfer into SMB2#
    if N_SMB2_SMR < C_SMB2_SMR/CN_SMR: #Immobilization#
        Nimmobilize = (C_SMB2_SMR/CN_SMR)-N_SMB2_SMR
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream/SMR Pool#
        SMR.AddCarbon(C_SMB2_SMR)
        SMR.AddNitrogen(N_SMB2_SMR+Nimmobilize)
    elif N_SMB2_SMR >= C_SMB2_SMR/CN_SMR:  #Mineralization#
        Nmineralize = N_SMB2_SMR-(C_SMB2_SMR/CN_SMR)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream/SMR Pool#
        SMR.AddCarbon(C_SMB2_SMR)
        SMR.AddNitrogen(C_SMB2_SMR/CN_SMR)

    ###STEP7: SMB1 Allocating to SMR ###
    ##Situation 1: Soil Decomposition##
    if initCont.litterzero == 0:
        micro_CP = 1.1
        death_rate = 0.1
        #Determing Decomposition Period#
        if SMB1.C >= (micro_CP*initCont.microSMB1Zero):     #Critical point#
            flag_SMB1 = 1                               #Decline Period#
        if flag_SMB1 == 0:  #Microbial anundance Increase Period#
            fM_SMB1 = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
            fM_SMB1 = min(0.75,fM_SMB1)
            death_rate = (death_rate)
        if flag_SMB1 == 1:  #Microbial anundance Decline Period#
            fM_SMB1 = fM_SMB1*(SMB1.C/(micro_CP*initCont.microSMB1Zero))
            fM_SMB1 = max(0.25,min(0.8,fM_SMB1))
            death_rate = (1+death_rate)
        ##Microbial Activity in Decomposition##
        C_SMB1_A  = SMB1.CU*(fA_SMB1)                     #Assimilation
        C_SMB1_E  = SMB1.CU*(1-fA_SMB1)                   #Egestion
        C_SMB1_rh = SMB1.CU*(fA_SMB1*fM_SMB1)             #Respiration
        C_SMB1_P  = SMB1.CU*(fA_SMB1*(1-fM_SMB1))         #Procreation
        C_SMB1_D  = C_SMB1_P*death_rate                   #Death
        N_SMB1_A  = SMB1.NU*(fA_SMB1)
        N_SMB1_E  = SMB1.NU*((1-fA_SMB1)+(fA_SMB1*fM_SMB1))
        N_SMB1_rh = SMB1.NU*0
        N_SMB1_P  = SMB1.NU*(fA_SMB1*(1-fM_SMB1))       
        N_SMB1_D  = N_SMB1_P*death_rate
    ##Situation 2: Litter Decomposition##
    if initCont.litterzero != 0:
        micro_CP = 1.25      #Pool Size#
        death_rate = 0.1
        #Determing Decomposition Period# 
        if SMB1.C >= (micro_CP*initCont.microSMB1Zero):           # Critical point
            flag_SMB1 = 1                                     # Start to decrease
        if flag_SMB1 == 0:   #Microbial anundance Increase Period#
            fM_SMB1_flag = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
            death_rate = death_rate
        if flag_SMB1 == 1:   #Microbial anundance Decline Period#
            #fM_SMB1_flag = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
            #death_rate = (1+death_rate)
            if SMB1.C >= initCont.microSMB1Zero:  #Decline Period 1#  ###Q: 20210421: SMB1 pool size and flux become 0 in old version
                fM_SMB1_flag = fM_SMB1
                death_rate   = (1+death_rate)
            if SMB1.C < initCont.microSMB1Zero:  #Decline Period 2#
                fM_SMB1_flag = fM_SMB1*(SMB1.C/initCont.microSMB1Zero)
                death_rate = (1+death_rate/10)
        ##Microbial Activity & CN Transfer in Decomposition##
        C_SMB1_A  = SMB1.CU*(fA_SMB1)                     #Assimilation
        C_SMB1_E  = SMB1.CU*(1-fA_SMB1)                   #Egestion     
        C_SMB1_rh = SMB1.CU*(fA_SMB1*fM_SMB1_flag)        #Respiration
        C_SMB1_P  = SMB1.CU*(fA_SMB1*(1-fM_SMB1_flag))    #Procreation
        C_SMB1_D  = C_SMB1_P*death_rate                   #Death
        N_SMB1_A  = SMB1.NU*(fA_SMB1)
        N_SMB1_E  = SMB1.NU*((1-fA_SMB1)+(fA_SMB1*fM_SMB1_flag))
        N_SMB1_rh = SMB1.NU*0
        N_SMB1_P  = SMB1.NU*(fA_SMB1*(1-fM_SMB1_flag))
        N_SMB1_D  = N_SMB1_P*death_rate
        litterToMB1_pct = litligsoil_CO2.litterC_MB1/(litligsoil_CO2.litterC_MB1+litligsoil_CO2.ligninC_MB1+litligsoil_CO2.soilC_MB1) ###Q
        ligninToMB1_pct = litligsoil_CO2.ligninC_MB1/(litligsoil_CO2.litterC_MB1+litligsoil_CO2.ligninC_MB1+litligsoil_CO2.soilC_MB1) ###Q
        soilToMB1_pct   = litligsoil_CO2.soilC_MB1/(litligsoil_CO2.litterC_MB1+litligsoil_CO2.ligninC_MB1+litligsoil_CO2.soilC_MB1)   ###Q
        litligsoil_CO2.litterCO2(C_SMB1_rh*litterToMB1_pct)    ###Q
        litligsoil_CO2.ligninCO2(C_SMB1_rh*ligninToMB1_pct)    ###Q: 20210517: change back to: litligsoil_CO2.ligninCO2(C_SMB1_rh*ligninToMB1_pct)
        litligsoil_CO2.soilCO2(C_SMB1_rh*soilToMB1_pct)        ###Q
        
    ##Updating Upstream Pool/SMB1 pool##
    SMB1.SubNewCarbon(C_SMB1_E+C_SMB1_rh+C_SMB1_P)
    SMB1.SubNewNitrogen(N_SMB1_E+N_SMB1_rh+N_SMB1_P)
    SMB1.AddCarbon(C_SMB1_P)
    SMB1.AddNitrogen(N_SMB1_P)
    SMB1.SubCarbon(C_SMB1_D)
    SMB1.SubNitrogen(N_SMB1_D)
    Rh.append(C_SMB1_rh)
    #N Coupling:Immobilization/Mineralization#
    #From SMB1 to SMR#
    C_SMB1_SMR = (C_SMB1_E+C_SMB1_D)
    N_SMB1_SMR = (N_SMB1_E+N_SMB1_D)
    if N_SMB1_SMR < C_SMB1_SMR/CN_SMR: #Immobilization#
        Nimmobilize = (C_SMB1_SMR/CN_SMR)-N_SMB1_SMR
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream/SMR Pool#
        SMR.AddCarbon(C_SMB1_SMR)
        SMR.AddNitrogen(N_SMB1_SMR+Nimmobilize)
    elif N_SMB1_SMR >= C_SMB1_SMR/CN_SMR: #Mineralization#
        Nmineralize = N_SMB1_SMR-(C_SMB1_SMR/CN_SMR)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream/SMR Pool#
        SMR.AddCarbon(C_SMB1_SMR)
        SMR.AddNitrogen(C_SMB1_SMR/CN_SMR)
     
    ###STEP8: SMR Allocating to NOM###
    #Decay rate#
    if bimm_SMR == 1:
        kSMR = 1/kmax_SMR/365*t_scale*w_scale*ni_scale
    else:
        kSMR = 1/kmax_SMR/365*t_scale*w_scale*nm_scale
    #Decomposition: CN Loss & Transfer#
    C_leaveSMR = SMR.C*kSMR   #Carbon loss from SMR#
    C_SMR_NOM  = C_leaveSMR   #Carbon transfer from SMR to NOM#
    N_leaveSMR = SMR.N*kSMR   #Nitrogen loss from SMR#
    N_SMR_NOM  = N_leaveSMR   #Nitrogen transfer from SMR to NOM#
    #Updating Upstream/SMR Pool#
    SMR.SubCarbon(C_leaveSMR)
    SMR.SubNitrogen(N_leaveSMR)
    #N coupling: Immobilization/Mineralization#
    #From SMR to NOM#
    if N_SMR_NOM < C_SMR_NOM/CN_NOM:  #Immobilization#
        Nimmobilize = (C_SMR_NOM/CN_NOM)-N_SMR_NOM
        Nimmobilize = min(Nimmobilize,avN.avn())
        avN.imm(Nimmobilize)
        Nmin = Nmin-Nimmobilize
        #Updating Downstream/NOM Pool#
        NOM.AddCarbon(C_SMR_NOM)
        NOM.AddNitrogen(N_SMR_NOM+Nimmobilize)
    elif N_SMR_NOM >= C_SMR_NOM/CN_NOM:  #Mineralization#
        Nmineralize = N_SMR_NOM-(C_SMR_NOM/CN_NOM)
        avN.NH4 = avN.NH4+Nmineralize
        Nmin = Nmin+Nmineralize
        #Updating Downstream/NOM Pool#
        NOM.AddCarbon(C_SMR_NOM)
        NOM.AddNitrogen(C_SMR_NOM/CN_NOM)
    
##########################################################################################################################
    #Add up all CO2#
    rh = C_AOM1_SMB1_rh+C_AOM1_SMB2_rh+C_AOM2_SMB1_rh+C_AOM2_SMB2_rh+C_AOM3_SMB1_rh+C_AOM3_SMB2_rh+C_SMB2_rh+C_SMB1_rh+C_NOM_SMB1_rh+C_NOM_SMB2_rh
    
    Rh.append(rh)
    Rh.append(litligsoil_CO2.litter_CO2+(litligsoil_CO2.lignin_CO2*0.69))  ###QQQ: 20210908: lignin cal: total lignin*0.69 (69% unlabeled litter-lignin )
    Rh.append(litligsoil_CO2.lignin_CO2*0.31*0.1*1000)                     ###QQQ: 20210908: lignin cal: total lignin*0.31 (31% synthesized lignin) * 1/10 (labeled lignin C) * 1000 (convert to ug/g/d) 
    Rh.append(litligsoil_CO2.soil_CO2)
    
    #Rh.append(litligsoil_CO2.litter_CO2+(litligsoil_CO2.lignin_CO2*10*2.2))  ###QQQ: 20210608: (label lignin) AOM1+AOM2+AOM3-litter lignin
    #Rh.append(litligsoil_CO2.lignin_CO2*1000)                                ###QQQ: 20210608: (label lignin) Only add labeled lignin-C in model (rate:ug/g/d)
    #Rh.append(litligsoil_CO2.soil_CO2)
    
    #Rh.append((litligsoil_CO2.litter_CO2 + litligsoil_CO2.lignin_CO2*(1-litligsoil_CO2.labelLignin_pct)))  ###Q: 03292021: Non-lignin portion rate
    #Rh.append((litligsoil_CO2.lignin_CO2 * litligsoil_CO2.labelLignin_pct*(0.1)*1000))                     ###Q: 03292021: synthesized lignin rate (ug/g/d)
    #Rh.append((litligsoil_CO2.soil_CO2))
    
    #Rh.append(rh)
    #Rh.append((litligsoil_CO2.litter_CO2))                                                 ###Q: 05172021: litter CO2_Simple version
    #Rh.append((litligsoil_CO2.lignin_CO2*litligsoil_CO2.labelLignin_pct*(0.1)*1000))       ###Q: 05172021: lignin CO2_Simple version
    #Rh.append((litligsoil_CO2.soil_CO2))
    
    
    return (Rh,flag_SMB1,flag_SMB2,fM_SMB1,fM_SMB2)


# In[ ]:




