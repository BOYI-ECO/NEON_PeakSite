#!/usr/bin/env python
# coding: utf-8

# # Main Function: Simulation Lab incubation

# In[1]:


##Sample site & Lab incubation days##
site = 2
day	 = 600

#Package import#
import numpy as np
import pandas as pd
from Decomposition_Peak import CWD_pool
from Decomposition_Peak import SOM_pool
from Decomposition_Peak import av_Nitro
from Decomposition_Peak import CWDFrag_OrgInAllocLab		       ###QQQ:20210602
from Decomposition_Peak import SoilDecom_pot
from Decomposition_Peak import SoilDecom_real
from Decomposition_Peak import initial_content
from Decomposition_Peak import lls_CO2

## Date for Lab measurement ##
day_selected = np.int32(np.loadtxt('day_obser.txt'))
day_selected = day_selected
# Read (Optimization) parameters getting from Matlab #
decom_para_1 = np.loadtxt('para.txt')					   # Paras need to be optimized #
decom_para = pd.read_csv("decom_para.csv",index_col=0)	   # Paras #
 
#Looping for Each Site#
for sample_site in range(1,site):
	#Input data#
	obser_data = pd.read_csv("obser_data.csv",index_col=0)
	#Parameters-MCMC#
	#Parameters#
	CN_AOM1	   = decom_para_1[0]
	CN_AOM2	   = decom_para_1[1]
	CN_AOM3	   = decom_para_1[2]
	CN_SMB1	   = decom_para_1[3]
	CN_SMB2	   = decom_para_1[4]
	CN_SMR	   = decom_para_1[5]
	CN_NOM	   = decom_para_1[6]
	CN_MOM	   = decom_para_1[7]
	kmax_AOM1  = decom_para_1[8]
	kmax_AOM2  = decom_para_1[9]
	kmax_AOM3  = decom_para_1[10]
	kmax_SMR   = decom_para_1[11]
	kmax_NOM   = decom_para_1[12]
	kmax_MOM   = decom_para_1[13]
	fAOM1_SMB1 = decom_para_1[14]
	fAOM2_SMB1 = decom_para_1[15]
	fAOM_MOM   = decom_para_1[16]
	fAOM3_SMB1 = decom_para_1[17]
	fNOM_SMB1  = decom_para_1[18]
	fNOM_MOM   = decom_para_1[19]
	fE_AOM1	   = decom_para_1[20]
	fE_AOM2	   = decom_para_1[21]
	fE_AOM3	   = decom_para_1[22]
	fE_NOM	   = decom_para_1[23]
	fA_SMB1	   = decom_para_1[24]
	fM_SMB1	   = decom_para_1[25]
	fA_SMB2	   = decom_para_1[26]
	fM_SMB2	   = decom_para_1[27]
	lig_p1	   = decom_para_1[28]
	lig_p2	   = decom_para_1[29]
	u	       = decom_para_1[30]				###QQQ: 20220310(1): day,u,theta,ymax, for peak sites
	theta	   = decom_para_1[31]
	ymax	   = decom_para_1[32]
	#Class initialization#
	litter = CWD_pool(obser_data,sample_site,"litter")
	AOM1   = SOM_pool(obser_data,sample_site,'AOM1')
	AOM2   = SOM_pool(obser_data,sample_site,'AOM2')
	AOM3   = SOM_pool(obser_data,sample_site,'AOM3')
	SMB1   = SOM_pool(obser_data,sample_site,'SMB1')
	SMB2   = SOM_pool(obser_data,sample_site,'SMB2')
	SMR	   = SOM_pool(obser_data,sample_site,'SMR')
	NOM	   = SOM_pool(obser_data,sample_site,'NOM')
	MOM	   = SOM_pool(obser_data,sample_site,'MOM')
	avN	   = av_Nitro(obser_data,sample_site)
	initCont   = initial_content()
	initCont.MicroSMB1Zero(SMB1.C)
	initCont.MicroSMB2Zero(SMB2.C)
	#initCont.LitterSub(litter.C*kbase_litter)
	#initCont.SoilSub(NOM.C*kmax_NOM)
	#Data Collection#
	Rh=[]
	flag_SMB1 = 0
	flag_SMB2 = 0


	CWDFrag_OrgInAllocLab(sample_site,obser_data,litter,AOM1,AOM2,AOM3,CN_AOM1,CN_AOM2,CN_AOM3)	               ###QQQ:20210602
	#Looping for Each Inbubation Day#
	for t in range(1,day):
		potimm = 0
		Nmin = 0
		litligsoil_CO2 = lls_CO2()

		


		bimm_AOM1,bimm_AOM2,bimm_AOM3,bimm_SMR,bimm_NOM,bimm_MOM,potimm = SoilDecom_pot(sample_site,t,obser_data,potimm,flag_SMB1,flag_SMB2,lig_p1,lig_p2,	   ###QQQ: 20210908
																						AOM1,AOM2,AOM3,SMB1,SMB2,SMR,NOM,MOM,initCont, 
																						CN_SMB1,CN_SMB2,CN_SMR,CN_NOM,CN_MOM, 
																						kmax_AOM1,kmax_AOM2,kmax_AOM3,kmax_SMR,kmax_NOM,kmax_MOM, 
																						fAOM1_SMB1,fAOM2_SMB1,fAOM_MOM,fAOM3_SMB1,fNOM_SMB1,fNOM_MOM, 
																						fE_AOM1,fE_AOM2,fE_AOM3,fE_NOM,fA_SMB1,fA_SMB2,fM_SMB1,fM_SMB2)

		Rh,flag_SMB1,flag_SMB2,fM_SMB1,fM_SMB2 = SoilDecom_real (u,theta,ymax,sample_site,t,obser_data,Nmin,potimm,flag_SMB1,flag_SMB2,Rh,lig_p1,lig_p2,	   ###QQQ: 20220310(2): day,u,theta,ymax, for peak sites
																 AOM1,AOM2,AOM3,SMB1,SMB2,SMR,NOM,MOM,avN,initCont,litligsoil_CO2, 
																 bimm_AOM1,bimm_AOM2,bimm_AOM3,bimm_SMR,bimm_NOM,bimm_MOM, 
																 CN_SMB1,CN_SMB2,CN_SMR,CN_NOM,CN_MOM, 
																 kmax_AOM1,kmax_AOM2,kmax_AOM3,kmax_SMR,kmax_NOM,kmax_MOM, 
																 fAOM1_SMB1,fAOM2_SMB1,fAOM_MOM,fAOM3_SMB1,fNOM_SMB1,fNOM_MOM, 
																 fE_AOM1,fE_AOM2,fE_AOM3,fE_NOM,fA_SMB1,fA_SMB2,fM_SMB1,fM_SMB2)


	##Data Output##
	Rh_array = np.array(Rh)
	Rh_Day = Rh_array.reshape(-1,14)
	Rh_Day_selected = Rh_Day[day_selected,:]
	Rh_Day_selected_frame = pd.DataFrame(Rh_Day_selected[:,-3:])
	Rh_Day_selected_frame_lss = pd.melt(Rh_Day_selected_frame).iloc[:,-1]
	np.savetxt('Rh_sum.txt',Rh_Day_selected_frame_lss)




