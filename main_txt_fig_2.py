'''
Created on Aug 30, 2022

@author: garmin
'''
'''
Using one irradiance to match Liefer data. Using liefer to validate macromolecular allocation in CFM

@author: Keisuke, gabby
'''

from pylab import *
from af001_energy_calculation import *
from numpy import *
from matplotlib import *
from matplotlib import pyplot
from sklearn.metrics import r2_score
from matplotlib.pyplot import figure, show, xticks, yticks, xlim, ylim, plot, xlabel, ylabel, title, stackplot, legend, scatter, bar, margins
import matplotlib.patches as mpat
from cmath import nan


#Liefer experimental notes: 85 umol/m2s, 12:12 light:dark cycle, 
#producing a bar graph to compare macromolecular allocation data and model output

What_is_limiting=1    #0: P-limiting  1:N-limiting

DPI = 300

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#Function beging here
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
def kkI(healey85data,What_is_limiting):   #this function calculate for the same irradiance
    #=============================
    # Parameter sets
    #=============================

    m=3.79146798299876E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
    Pmax=0.00320513285659728
    OT=0.00863364097132997
    Ynphoto_chl=3.56099164557551          #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
    Cnbiosynth=4.34728279914354E-10        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Nconst_protein=4.45336898828389E-15     #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable=6212.59249917364        #(s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other=5.44534485638617E-17               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max=25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
    Cessential=1.51786753491048E-15             #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
    
#==============================

    I = 85                              #healey85data.Lightintensity
    #plotcolor=healey85data.plotcolor

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Parameters
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    if What_is_limiting==0:
    #For P-limiting case
        Pin=0.002  #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.2    #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=1.00*10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)
    
    elif What_is_limiting==1:
    #For N-limiting case
        Pin=0.02  #(mol/m3) Phosphorus concentration in the incoming medium (Healey 1985)
        Nin=0.05    #(mol/m3) Nitrate concentration in the incoming medium (Healey 1985)
        Qc=10**(-12)/12      #(molC/cell) biomass C per cell (196-18)(from Healey 1985)

    E3=evalue()
    E=E3.E
    C=6
    
    Ddmax=1.6
    #Dmax=4
    Dstep=0.001
    Dd=arange(Dstep,Ddmax+Dstep,Dstep)       #(h-1) growth rate
    #print(Dd)
    U=arange(0,Ddmax/Dstep,1).astype(int)
    D=Dd/(3600*24)
    
    #------------------------------
    #Photosynthesis
    #------------------------------

    Pchl=Pmax*(1-exp(-OT*I)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
    Pchl=Pchl/2
    ls=D*Qc                    #(molC s-1) Biomass synthesis rate (193-25) ABIO
    #------------------------------
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Key parameters: parameterization ideas -> Kei 193-28
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    Molar_mass_DNA_AT_average=308.47        #(g mol-1) Molar mass AT average (from "05 Review nucleic acid composition.xlsx")
    Molar_mass_DNA_CG_average=308.97        #(g mol-1) Molar mass CG average (from "05 Review nucleic acid composition.xlsx")
    
    Molar_mass_RNA_AT_average=317.47        #(g mol-1) Molar mass AT average (from "05 Review nucleic acid composition.xlsx")
    Molar_mass_RNA_CG_average=324.97        #(g mol-1) Molar mass CG average (from "05 Review nucleic acid composition.xlsx")

    #================================
    #E coli
    #================================
    CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)

    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)

    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]

#New part++++++++++++++++++++++++++++++++++++++++++++++++++
    AT=1-CG
    AU=1-CG
    #Values per P in mol/mol from "05 Review nucleic acid composition.xlsx"
    #RNA  
    C_CG_RNA = 19/2
    N_CG_RNA = 8/2
    
    C_AU_RNA = 19/2
    N_AU_RNA = 7/2
  
    #DNA  
    C_CG_DNA = 19/2
    N_CG_DNA = 8/2
    
    C_AT_DNA = 20/2
    N_AT_DNA = 7/2

    
    YnucacidN_P = N_CG_RNA*CG + N_AU_RNA*AU #same for RNA and DNA
    YnucacidP_N = 1/YnucacidN_P
    
    YdnaC_N = (C_CG_DNA*CG + C_AT_DNA*AT)/(N_CG_DNA*CG + N_AT_DNA*AT)
    YrnaC_N = (C_CG_RNA*CG + C_AU_RNA*AU)/(N_CG_RNA*CG + N_AU_RNA*AU)
    #print(YdnaC_N, YrnaC_N,YnucacidN_P,YdnaC_N*YnucacidN_P,YrnaC_N*YnucacidN_P)
#++++++++++++++++++++++++++++++++++++++++++++++++++
    #print(YnucacidN_P,YdnaC_N,YrnaC_N)
    
    DNAmb=2.1269                   #(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    Avogadro=6.022*10**23           #(molecules mol-1) Avogadro constant
    Pdna_const=DNAmb*2*10**6/Avogadro                #(molP cell-1) Constant part of DNA in phosphorus 
    Prna_const=Pdna_const*RNA_DNA_molar_ratio       #(molP cell-1) Constant part of RNA in phosphorus
    #* Make sure to multiply by 2 as they are base PAIRs"
    # print(YnucacidP_N,'b')
    # YnucacidP_N = 1/3.75
    Ndna_const=Pdna_const/YnucacidP_N      #(molN cell-1) Constant part of DNA in nitrogen
    Nrna_const=Ndna_const*RNA_DNA_molar_ratio   #(molN cell-1) Constatn part of RNA in phosphorus
    #print(Nrna_const,'a')
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Calculation
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    Chl=((1+E)*ls+m)/Pchl       #(molC chl cell-1) Chlrophyll concentration (193-25) 
    Nphoto=Chl*Ynphoto_chl  #(molN cell-1) Photosynthesis related protein nitrogen (193-25)
    Nbiosynth=D*Cnbiosynth             #(molN cell-1) various part of biosynthesis related protein in N (193-37)
    Nprotein=Nphoto+Nconst_protein+Nbiosynth    #(molN cell-1) All the proteins in N (193-26)
    Nrna_variable=Nprotein*D*Cnrna_variable        #(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
    Ndna_variable=Ndna_const*Dd/1.2*(18.3-7.6)/7.6        #(molN cell-1) variable part of nitrogen in DNA (193-26) Increasing ratio based on Bremmer 1996
    Ndna_variable=0*Dd                     #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend. 
   
    #print(Dd/1.2*(18.3-7.6)/7.6)
    #Ndnarna=Nconst_dnarna+Nnucacid_variable
    Nchl=Chl*4/55                           #(molN cell-1) Chlorophyll nitrogen (actually almost negligiable)
    Pthylakoid=Chl*Ypthylakoid_chl          #(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
    Prna_variable=Nrna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in RNA (193-26)
    Pdna_variable=Ndna_variable*YnucacidP_N     #(molP cell-1) variable part of phosphorus in DNA (193-26)
    
    #=================================
    #Total calculation
    #=================================
    Qn_max=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl+Nstore_max      #(molN cell-1)  nitrogen content in the cell (193-26)                                 #(molN cell-1) total phosphorus content in the cell (193-26)                                                                             #(molP cell-1) total phosphorus content in the cell (193-26)
    Qn_min=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl            #(molN cell-1) total nitrogen in the cell without storage
    Qp_min=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const      #(molP cell-1) total phosphorus in the cell without storage
    
    #=================================
    #Vector preparation
    #=================================
    Nstore=zeros(size(Dd))
    X=zeros(size(Dd))
    Qn_test=zeros(size(Dd))
    Qp_test=copy(X)
    Qp=copy(X)
    Qn=copy(X)
    Pstore=copy(X)
    Limitation=copy(X)
    #=================================
    #Population calculation
    #=================================
    Xn_max=Nin/Qn_min
    Xp_max=Pin/Qp_min
    for i in U:
        if Xn_max[i]>Xp_max[i]:
            X[i]=Xp_max[i]
            Qp[i]=Qp_min[i]
            Qn_test[i]=Nin/X[i]
            Limitation[i]=0
            if Qn_test[i]<Qn_max[i]:
                Qn[i]=Qn_test[i]
                Nstore[i]=Qn_test[i]-Nprotein[i]-Nrna_variable[i]-Nrna_const-Ndna_variable[i]-Ndna_const-Nchl[i]  #(molN cell-1) Nitrogen storage in the cell
            else:
                Qn[i]=Qn_max[i]
                Nstore[i]=Nstore_max
        else:
            X[i]=Xn_max[i]
            Qn[i]=Qn_min[i]
            Qp_test[i]=Pin/X[i]
            Limitation[i]=1
            if Qp_test[i]<Qp_max:
                Qp[i]=Qp_test[i]
            else:
                Qp[i]=Qp_max
                Qp[i]=Qp_min[i]+Qp_max
            Pstore[i]=Qp[i]-Pconst_other-Pthylakoid[i]-Prna_variable[i]-Prna_const-Pdna_variable[i]-Pdna_const   #(molP cell-1) Stored phosphorus in the cell
            if Pstore[i]<0:
                Pstore[i]=0
            
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc   #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc*14*10**6/(12*10**3)    #(ug N / mg C) biomass N to C ratio (164-20)
    PtoCplot=Qp/Qc*30.97*10**6/(12*10**3)    #(ug P / mg C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp*14*10**6/(30.97*10**6)        #(ug N /ug P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc         #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6     #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 3 (unit adjustment)
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Nunit=1/Qc#*14*10**6/(12*10**3)          #((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
    Punit=1/Qc#*30.97*10**6/(12*10**3)       #((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
    Numbertoarray=ones(size(Dd))            #(dimensionless) Number to array converter
    
    NunitData = 14*10**6/(12*10**3)
    PunitData = 30.97*10**6/(12*10**3)
    
    
    #=======================================
    #Calculation of carbon usage (195-16)
    #=======================================
    CNprotein=3.820219393   #(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
    
    #-------------------
    #C protein
    #-------------------
    Cphoto=Nphoto*CNprotein     #(molC cell-1) carbon in photosystem protein (195-16)
    Cbiosynth=Nbiosynth*CNprotein   #(molC cell-1) carbon in biosynthesis protein (195-16)
    Cconst_protein=Nconst_protein*CNprotein  #(molC cell-1) carbon in other protein assumed constant (195-16)
 
    #----------------------
    #C chlorophyll
    #----------------------
    Cchl=Chl                    #(molC cell-1) carbon in chlorophyll (195-16)
    
    #----------------------
    #C DNA RNA
    #----------------------
    Crna_const=Nrna_const*YrnaC_N       #(molC cell-1) carbon in variable part of RNA (195-16)
    Crna_variable=Nrna_variable*YrnaC_N     #(molC cell-1) carbon in variable part of RNA (195-16)
    
    Cdna_const=Ndna_const*YdnaC_N       #(molC cell-1) carbon in constant part of DNA (195-16)
    Cdna_variable=Ndna_variable*YdnaC_N     #(molC cell-1) carbon in variable part of DNA (195-16)
    
    Cnstore=Nstore*YcyanoC_N        #(molC cell-1) carbon in nitrogen storage (cyanophycin)
    CthylakoidPG=Pthylakoid*YpgC_P           #(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
    Pstore_plot=Pstore*Punit 
    Prna_variable_plot=Prna_variable*Punit      #(ug P/mgC) Phosphorus in variable part of RNA (193-37)
    Pdna_variable_plot=Pdna_variable*Punit      #(ug P/mgC) Phosphorus in variable part of DNA (193-37)
    
    Pthylakoid_plot=Pthylakoid*Punit     #(ug P/ mgC) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)(193-33)
    Pconst_other_plot=Pconst_other*Punit*Numbertoarray       #(ug P/ mgC) PHosphorus in other parts (ex. phospholipid in outer membrane, ATP, ADP, etc. (assuming constant) (193-33)
    
    #PtoCplot2=Pdna_plot+Prna_plot+Pthylakoid_plot+Pconst_other_plot   #For calculation check
    Pdna_const_plot=Pdna_const*Punit*Numbertoarray    #(ug P/ mgC) Phosphorus in constant part of DNA
    Prna_const_plot=Prna_const*Punit*Numbertoarray    #(ug P/ mgC) Phosphorus in constant part of RNA
    
    
    #print(Prna_variable_plot)
    PtoCplot=Pconst_other_plot+Pdna_const_plot+Pdna_variable_plot+Pthylakoid_plot+Prna_variable_plot+Prna_const_plot+Pstore_plot
    
    #print(Pstore_plot)
    #=======================================
    #C for plot in
    #=======================================
    percentorratio=100       #100: percent, 1:ratio
    Cphoto_plot=Cphoto/Qc*percentorratio           
    Cbiosynth_plot=Cbiosynth/Qc*percentorratio
    Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numbertoarray
    Cchl_plot=Cchl/Qc*percentorratio
    Crna_const_plot=Crna_const/Qc*percentorratio*Numbertoarray
    Crna_variable_plot=Crna_variable/Qc*percentorratio
    Cdna_const_plot=Cdna_const/Qc*percentorratio*Numbertoarray
    Cdna_variable_plot=Cdna_variable/Qc*percentorratio
    Cother_plot=1*percentorratio-Cphoto_plot-Cbiosynth_plot-Cconst_protein_plot-Cchl_plot\
             -Crna_const_plot-Crna_variable_plot-Cdna_const_plot-Cdna_variable_plot
    #Cother_plot=Cother_plot/Qc*percentorratio
    Cessential_plot=Cessential/Qc*percentorratio*Numbertoarray
    Cnstore_plot=Cnstore/Qc*percentorratio
    CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio
    
    #print(Cother_plot)
    #=======================================
    #N for plot ***(Mainly from 193-33)***
    #=======================================
    Nphoto_plot=Nphoto*Nunit    #(ug N/ mgC) Photosynthesis related protein nitrogen (193-25)(193-33)
    #Nproteinsynth_plot=Nproteinsynth*Nunit      #(ug N/ mgC) protein synthesis related protein in N (193-25)(193-33)
    Nbiosynth_plot=Nbiosynth*Nunit      #(ug N/ mgC) biosynthesis related protein in N (193-37)
    Nconst_protein_plot=Nconst_protein*Nunit*Numbertoarray    #(ug N/ mgC) constant protein pool in nitrogen (193-25)(193-33)
    Nchl_plot=Nchl*Nunit        #(ug N/ mgC) Chlorophyll nitrogen (actually almost negligiable) (193-33)
#    Nconst_other_plot=Nconst_other*Nunit*Numbertoarray        #(ug N/ mgC) Other kinds of nitrogen assuming constant (193-33)
    
    Nrna_variable_plot=Nrna_variable*Nunit      #(ug N/ mgC) Nitrogen in Variable part of nucleic acid (193-37)
    Ndna_variable_plot=Ndna_variable*Nunit      #(ug N/ mgC) Nitrogen in Variable part of nucleic acid (193-37)
    
    Ndna_const_plot=Ndna_const*Nunit*Numbertoarray    #(ug N/ mgC) Nitrogen in constant part of DNA
    Nrna_const_plot=Nrna_const*Nunit*Numbertoarray    #(ug N/ mgC) Nitrogen in constant part of RNA
    
    Nstore_plot=Nstore*Nunit      #(ug N/ mgC) Nitrogen in storage
    
    
    Pstore_plot=Pstore*Punit 
    Prna_variable_plot=Prna_variable*Punit      #(ug P/mgC) Phosphorus in variable part of RNA (193-37)
    Pdna_variable_plot=Pdna_variable*Punit      #(ug P/mgC) Phosphorus in variable part of DNA (193-37)
    
    Pthylakoid_plot=Pthylakoid*Punit     #(ug P/ mgC) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)(193-33)
    Pconst_other_plot=Pconst_other*Punit*Numbertoarray       #(ug P/ mgC) PHosphorus in other parts (ex. phospholipid in outer membrane, ATP, ADP, etc. (assuming constant) (193-33)
    
    #PtoCplot2=Pdna_plot+Prna_plot+Pthylakoid_plot+Pconst_other_plot   #For calculation check
    Pdna_const_plot=Pdna_const*Punit*Numbertoarray    #(ug P/ mgC) Phosphorus in constant part of DNA
    Prna_const_plot=Prna_const*Punit*Numbertoarray    #(ug P/ mgC) Phosphorus in constant part of RNA
    
    
    #print(Prna_variable_plot)
    PtoCplot=Pconst_other_plot+Pdna_const_plot+Pdna_variable_plot+Pthylakoid_plot+Prna_variable_plot+Prna_const_plot+Pstore_plot
    
    #print(Pstore_plot)
    #---------------------------------------------------
    #C other: Here revised to include Nstore reduction
    #---------------------------------------------------

    #print(Qc)
    Cother_without_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-CthylakoidPG
    
    Cother_with_full_Nstore=Qc-Cphoto-Cbiosynth-Cconst_protein-Cchl\
            -Crna_const-Crna_variable-Cdna_const-Cdna_variable\
            -Cessential-Cnstore-CthylakoidPG
            
    Cother=Cother_with_full_Nstore            
    
    Nstore_reduce=logical_and(Cother_without_Nstore>0, Cother_with_full_Nstore<0)
    
    Cnstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]
    Nstore0=copy(Nstore)
    Nstore[Nstore_reduce]=Cother_without_Nstore[Nstore_reduce]/YcyanoC_N
    Cother[Nstore_reduce]=0
    Qn[Nstore_reduce]=Qn[Nstore_reduce]+Nstore[Nstore_reduce]-Nstore0[Nstore_reduce]
    
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For plotting 1
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    BiomassC=12*X*Qc   #(mg C L-1) Biomass concentration
    NtoCplot=Qn/Qc#*14*10**6/(12*10**3)    #(mol N / mol C) biomass N to C ratio (164-20)
    
    #PtoCplot=Qp/Qc#*30.97*10**6/(12*10**3)    #(mol P / mol C) biomass P to C ratio (164-20)
    NtoPplot=Qn/Qp#*14*10**6/(30.97*10**6)        #(mol N /mol P) biomass N to P ratio (164-20)
    ChltoC0=Chl/Qc         #(mol C chl mol C -1) Chlorophyll to carbon ratio
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
    ChltoCplot=ChltoC0/12/1000*Mchl/55*10**6     #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
    

    Dd[Cother<0]=nan
    
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Plot free parameters #Here units adjusted for paper
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    from decimal import Decimal
    def sci(Numb):
        Numb='%.2E' % Decimal(Numb)
        return Numb

    #import ebenezeer data
    data= genfromtxt('Chloropicon.csv', delimiter=',')
    growth_rate= data[1:46,13]
    gr_mari=growth_rate[0:15]
    gr_maure=growth_rate[15:30]
    gr_rosco=growth_rate[30:]
    
    NP= data[1:46,9]
    NP_mari=NP[0:15]
    NP_maure=NP[15:30]
    NP_rosco=NP[30:]
     
     
    C= data[1:46,16]
    C_mari=C[0:15]
    C_maure=C[15:30]
    C_rosco=C[30:]
    
      #N
    N= data[1:46,17]
    N_mari=N[0:15]
    N_maure=N[15:30]
    N_rosco=N[30:]
    
    #N:C
    NC_mari=N_mari/C_mari
    NC_maure=N_maure/C_maure
    NC_rosco=N_rosco/C_rosco
    
    #P
    P= data[1:46,18]
    P_mari=P[0:15]
    P_maure=P[15:30]
    P_rosco=P[30:]
    
    #PC
    PC_mari=P_mari/C_mari
    PC_maure=P_maure/C_maure
    PC_rosco=P_rosco/C_rosco
    
    NP_mari=N_mari/P_mari
    NP_maure=N_maure/P_maure
    NP_rosco=N_rosco/P_rosco
    
   
    protein=data[1:46,20]
    carb=data[1:46,21]
    lipid=data[1:46,22]
    chlorophyll=data[1:46,23]
    RNA=data[1:46,24]
    DNA=data[1:46,25]
    Pinlipid=data[1:46,26]
    polyphosphate=data[1:46,27]
    
    
        ##### find carbon allocation elements ###### MARI    N_mari=N[0:15]

    RNA_mari=RNA[0:15]
    DNA_mari=DNA[0:15]
    Pro_mari=protein[0:15]
    Pig_mari=chlorophyll[0:15]
    Carb_mari=carb[0:15]
    Lip_mari=lipid[0:15]
    Pinlipid_mari=Pinlipid[0:15]
    poly_mari=polyphosphate[0:15]
    
    #C fractions
    RNA_C=(RNA_mari/C_mari)
    DNA_C=(DNA_mari/C_mari)
    Pro_C=(Pro_mari/C_mari)
    Pig_C=(Pig_mari/C_mari)
    Carb_c=(Carb_mari/C_mari)
    Lip_c=(Lip_mari/C_mari)

    
    
    RNA_Cmari=RNA_C*0.34
    DNA_Cmari=DNA_C*0.36
    Pro_Cmari=Pro_C*0.53
    Pig_Cmari=Pig_C*0.74
    Carb_Cmari=Carb_c*0.40
    Lip_Cmari=Lip_c*0.25
    total=RNA_Cmari+DNA_Cmari+Pro_Cmari+Pig_Cmari+Carb_Cmari+Lip_Cmari


    RNA_Cmari=RNA_Cmari/total
    DNA_Cmari=DNA_Cmari/total
    Pro_Cmari=Pro_Cmari/total
    Pig_Cmari=Pig_Cmari/total
    Carb_Cmari=Carb_Cmari/total
    Lip_Cmari=Lip_Cmari/total
    
    ##### find nitrogen allocation elements #####
    RNA_Nmari=RNA_mari/N_mari
    DNA_Nmari=DNA_mari/N_mari
    Pro_Nmari=Pro_mari/N_mari
    Pig_Nmari=Pig_mari/N_mari
    
    RNA_Nmari=RNA_Nmari*0.1612
    DNA_Nmari=DNA_Nmari*0.1684
    Pro_Nmari=Pro_Nmari*0.16
    Pig_Nmari=Pig_Nmari*0.063
       
    total=RNA_Nmari+DNA_Nmari+Pro_Nmari+Pig_Nmari
    
    RNA_Nmari=RNA_Nmari/total
    DNA_Nmari=DNA_Nmari/total
    Pro_Nmari=Pro_Nmari/total
    Pig_Nmari=Pig_Nmari/total
    


    Growt=size(gr_mari)  
    mari_growth_plot = zeros(Growt)
    mari_NC= zeros(Growt)
    mari_PC= zeros(Growt)
    mari_NP= zeros(Growt)
    
      #carbon mari
    mari_RNA_plot = zeros(Growt)
    mari_DNA_plot = zeros(Growt)
    mari_Pro_plot = zeros(Growt)
    mari_Pig_plot = zeros(Growt)
    mari_Carb_plot = zeros(Growt)
    mari_Lip_plot = zeros(Growt)
    mari_resid_plot= zeros(Growt)
    
    #nitrogen
    mari_RNA_Nplot = zeros(Growt)
    mari_DNA_Nplot = zeros(Growt)
    mari_Pro_Nplot = zeros(Growt)
    mari_Pig_Nplot = zeros(Growt)
    
    #phosphorus
    mari_RNA_Pplot = zeros(Growt)
    mari_DNA_Pplot = zeros(Growt)
    mari_Lip_Pplot = zeros(Growt)
    mari_Res_Pplot = zeros(Growt)

    
    for i in arange(Growt):
        j=where(gr_mari==gr_mari.min())[0][0]
        mari_growth_plot[i] = gr_mari[j]
        mari_NC[i]=NC_mari[j]
        mari_PC[i]=PC_mari[j]
        mari_NP[i]=NP_mari[j]
        #####################
        mari_RNA_plot[i] = RNA_Cmari[j]
        mari_DNA_plot[i]=DNA_Cmari[j]
        mari_Pro_plot[i]=Pro_Cmari[j]
        mari_Pig_plot[i]=Pig_Cmari[j]
        mari_Carb_plot[i]=Carb_Cmari[j]
        mari_Lip_plot[i]=Lip_Cmari[j]
        #########################
        mari_RNA_Nplot[i] = RNA_Nmari[j]
        mari_DNA_Nplot[i]=DNA_Nmari[j]
        mari_Pro_Nplot[i]=Pro_Nmari[j]
        mari_Pig_Nplot[i]=Pig_Nmari[j]
        #########################]
        gr_mari[j] = 1e30        

   
#import liefer data
    liefer=genfromtxt('data_clean2.csv',delimiter=',')

    #separate into different species
    th_ps=liefer[1:16,2:]
    th_we=liefer[18:33,2:]
    os_ta=liefer[35:50,2:]
    micro=liefer[52:67,2:]


    #find growth rates
    growth_ps=th_ps[:,1]
    growth_we=th_we[:,1]
    growth_os=os_ta[:,1]
    growth_micro=micro[:,1]
 
############################ Elemental Stoichiometry #######################################

    #### N:P values
    NP_ps=th_ps[:,46]
    NP_we=th_we[:,46]
    NP_os=os_ta[:,46]
    NP_micro=micro[:,46]
 
    #### C:N values
    CN_ps=th_ps[:,44]
    CN_we=th_we[:,44]
    CN_os=os_ta[:,44]
    CN_micro=micro[:,44]
    ###convert to N:C
    CN_ps=1/CN_ps
    CN_we=1/CN_we
    CN_os=1/CN_os
    CN_micro=1/CN_micro
 
    #### C:P values
    CP_ps=th_ps[:,45]
    CP_we=th_we[:,45]
    CP_os=os_ta[:,45]
    CP_micro=micro[:,45]

    #convert to P:C values
    CP_ps=1/CP_ps
    CP_we=1/CP_we
    CP_os=1/CP_os
    CP_micro=1/CP_micro
    
    G = size(growth_os)
    ps_growth_plot = zeros(G)
    #stoichiometry
    NC=zeros(G)
    PC=zeros(G)
    NP=zeros(G)
        ##### find carbon allocation elements ######
    RNA_ps=th_ps[:,23]
    DNA_ps=th_ps[:,24]
    Pro_ps=th_ps[:,25]
    Pig_ps=th_ps[:,26]
    Carb_ps=th_ps[:,27]
    Lip_ps=th_ps[:,28]
    
    ##### find nitrogen allocation elements #####
    RNA_Nps=th_ps[:,33]
    DNA_Nps=th_ps[:,34]
    Pro_Nps=th_ps[:,35]
    Pig_Nps=th_ps[:,36]
    
    ##### find phosphorus allocation elements #####
    RNA_Pps=th_ps[:,40]
    DNA_Pps=th_ps[:,41]
    Lip_Pps=th_ps[:,42]
    Res_Pps=th_ps[:,43]
    
    #### C:N values
    CN_ps=th_ps[:,44]
    ###convert to N:C
    CN_ps=1/CN_ps
    #### C:P values
    CP_ps=th_ps[:,45]
    #convert to P:C values
    CP_ps=1/CP_ps
    
    #setting up the arrays
    U = size(growth_ps)
    growth_plot = zeros(U)
    
    #carbon
    ps_RNA_plot = zeros(U)
    ps_DNA_plot = zeros(U)
    ps_Pro_plot = zeros(U)
    ps_Pig_plot = zeros(U)
    ps_Carb_plot = zeros(U)
    ps_Lip_plot = zeros(U)
    
    #nitrogen
    ps_RNA_Nplot = zeros(U)
    ps_DNA_Nplot = zeros(U)
    ps_Pro_Nplot = zeros(U)
    ps_Pig_Nplot = zeros(U)
    
    #phosphorus
    ps_RNA_Pplot = zeros(U)
    ps_DNA_Pplot = zeros(U)
    ps_Lip_Pplot = zeros(U)
    ps_Res_Pplot = zeros(U)
    
    #stoichiometry
    ps_NC=zeros(U)
    ps_PC=zeros(U)
    
    for i in arange(U):
        j = where(growth_ps==growth_ps.min())[0][0]
        ps_growth_plot[i] = growth_ps[j]
        #######################
        ps_RNA_plot[i] = RNA_ps[j]
        ps_DNA_plot[i]=DNA_ps[j]
        ps_Pro_plot[i]=Pro_ps[j]
        ps_Pig_plot[i]=Pig_ps[j]
        ps_Carb_plot[i]=Carb_ps[j]
        ps_Lip_plot[i]=Lip_ps[j]
        #########################
        ps_RNA_Nplot[i] = RNA_Nps[j]
        ps_DNA_Nplot[i]=DNA_Nps[j]
        ps_Pro_Nplot[i]=Pro_Nps[j]
        ps_Pig_Nplot[i]=Pig_Nps[j]
        #########################
        ps_RNA_Pplot[i] = RNA_Pps[j]
        ps_DNA_Pplot[i] = DNA_Pps[j]
        ps_Lip_Pplot[i] = Lip_Pps[j]
        ps_Res_Pplot[i] = Res_Pps[j]
        ########################
        ps_NC[i]=CN_ps[j]
        ps_PC[i]=CP_ps[j]
        ###########################
        growth_ps[j] = 1e30
    
    
    ###Get rid of the data point with a nan value#####
    ps_growth_Cplot=delete(growth_plot,[9])
    ps_Lip_plot=delete(ps_Lip_plot,[9])
    ps_RNA_plot = delete(ps_RNA_plot, [9])
    ps_DNA_plot = delete(ps_DNA_plot, [9])
    ps_Pro_plot = delete(ps_Pro_plot, [9])
    ps_Pig_plot = delete(ps_Pig_plot, [9])
    ps_Carb_plot = delete(ps_Carb_plot, [9])
    
    
    G = size(growth_os)
    we_growth_plot = zeros(G)
    #stoichiometry
    we_NC=zeros(G)
    we_PC=zeros(G)
    we_NP=zeros(G)
    #find carbon allocation elements
    RNA_we=th_we[:,23]
    DNA_we=th_we[:,24]
    Pro_we=th_we[:,25]
    Pig_we=th_we[:,26]
    Carb_we=th_we[:,27]
    Lip_we=th_we[:,28]
    
    ##### find nitrogen allocation elements #####
    RNA_Nwe=th_we[:,33]
    DNA_Nwe=th_we[:,34]
    Pro_Nwe=th_we[:,35]
    Pig_Nwe=th_we[:,36]
    
    ##### find phosphorus allocation elements #####
    RNA_Pwe=th_we[:,40]
    DNA_Pwe=th_we[:,41]
    Lip_Pwe=th_we[:,42]
    Res_Pwe=th_we[:,43]
    
    #### C:N values
    CN_we=th_we[:,44]
    ###convert to N:C
    CN_we=1/CN_we
    #### C:P values
    CP_we=th_we[:,45]
    #convert to P:C values
    CP_we=1/CP_we
    
    ##setting the arrays###
    U = size(growth_we)
    we_growth_plot = zeros(U)
    #carbon
    we_RNA_plot = zeros(U)
    we_DNA_plot = zeros(U)
    we_Pro_plot = zeros(U)
    we_Pig_plot = zeros(U)
    we_Carb_plot = zeros(U)
    we_Lip_plot = zeros(U)
    #nitrogen
    we_RNA_Nplot = zeros(U)
    we_DNA_Nplot = zeros(U)
    we_Pro_Nplot = zeros(U)
    we_Pig_Nplot = zeros(U)
    #phosphorus
    we_RNA_Pplot = zeros(U)
    we_DNA_Pplot = zeros(U)
    we_Lip_Pplot = zeros(U)
    we_Res_Pplot = zeros(U)
    #stoichiometry
    we_NC=zeros(U)
    we_PC=zeros(U)
    
    for i in arange(U):
         j = where(growth_we==growth_we.min())[0][0]
         we_growth_plot[i] = growth_we[j]
         #############################
         we_RNA_plot[i] = RNA_we[j]
         we_DNA_plot[i]=DNA_we[j]
         we_Pro_plot[i]=Pro_we[j]
         we_Pig_plot[i]=Pig_we[j]
         we_Carb_plot[i]=Carb_we[j]
         we_Lip_plot[i]=Lip_we[j]
         ############################
         we_RNA_Nplot[i] = RNA_Nwe[j]
         we_DNA_Nplot[i]=DNA_Nwe[j]
         we_Pro_Nplot[i]=Pro_Nwe[j]
         we_Pig_Nplot[i]=Pig_Nwe[j]
         ##########################
         we_RNA_Pplot[i] = RNA_Pwe[j]
         we_DNA_Pplot[i] = DNA_Pwe[j]
         we_Lip_Pplot[i] = Lip_Pwe[j]
         we_Res_Pplot[i] = Res_Pwe[j]
         ##########################
         we_NC[i]=CN_we[j]
         we_PC[i]=CP_we[j]
        ###########################
         growth_we[j] = 1e30
    
    
    ###Get rid of the data point with a nan value#####
    we_growth_Cplot=delete(growth_plot,[8])
    we_Lip_plot=delete(we_Lip_plot,[8])
    we_RNA_plot = delete(we_RNA_plot, [8])
    we_DNA_plot = delete(we_DNA_plot, [8])
    we_Pro_plot = delete(we_Pro_plot, [8])
    we_Pig_plot = delete(we_Pig_plot, [8])
    we_Carb_plot = delete(we_Carb_plot, [8])
    
    we_growth_Pplot=delete(we_growth_plot,[8])
    we_Lip_Pplot=delete(we_Lip_Pplot,[8])
    we_RNA_Pplot = delete(we_RNA_Pplot, [8])
    we_DNA_Pplot = delete(we_DNA_Pplot, [8])
    we_Res_Pplot = delete(we_Res_Pplot,[8])
    we_PC= delete(we_PC, [8])
    
    G = size(growth_os)
    growth_plot = zeros(G)
    #stoichiometry
    NC=zeros(G)
    PC=zeros(G)
    NP=zeros(G)
     #find carbon allocation elements
    RNA_os=os_ta[:,23]
    DNA_os=os_ta[:,24]
    Pro_os=os_ta[:,25]
    Pig_os=os_ta[:,26]
    Carb_os=os_ta[:,27]
    Lip_os=os_ta[:,28]
    
    ##### find nitrogen allocation elements #####
    RNA_Nos=os_ta[:,33]
    DNA_Nos=os_ta[:,34]
    Pro_Nos=os_ta[:,35]
    Pig_Nos=os_ta[:,36]
    
    ##### find phosphorus allocation elements #####
    RNA_Pos=os_ta[:,40]
    DNA_Pos=os_ta[:,41]
    Lip_Pos=os_ta[:,42]
    Res_Pos=os_ta[:,43]
      
            #setting the arrays up
    U = size(growth_os)
    ot_growth_plot = zeros(U)
    #carbon
    ot_RNA_plot = zeros(U)
    ot_DNA_plot = zeros(U)
    ot_Pro_plot = zeros(U)
    ot_Pig_plot = zeros(U)
    ot_Carb_plot = zeros(U)
    ot_Lip_plot = zeros(U)
    #nitrogen
    ot_RNA_Nplot = zeros(U)
    ot_DNA_Nplot = zeros(U)
    ot_Pro_Nplot = zeros(U)
    ot_Pig_Nplot = zeros(U)
    #phosphorus
    ot_RNA_Pplot = zeros(U)
    ot_DNA_Pplot = zeros(U)
    ot_Lip_Pplot = zeros(U)
    ot_Res_Pplot = zeros(U)
    #stoichiometry
    ot_NC=zeros(U)
    ot_PC=zeros(U)
    
    for i in arange(U):
         j = where(growth_os==growth_os.min())[0][0]
         ot_growth_plot[i] = growth_os[j]
         ###############################
         ot_RNA_plot[i] = RNA_os[j]
         ot_DNA_plot[i]=DNA_os[j]
         ot_Pro_plot[i]=Pro_os[j]
         ot_Pig_plot[i]=Pig_os[j]
         ot_Carb_plot[i]=Carb_os[j]
         ot_Lip_plot[i]=Lip_os[j]
         ################################
         ############################
         ot_RNA_Nplot[i] = RNA_Nos[j]
         ot_DNA_Nplot[i]=DNA_Nos[j]
         ot_Pro_Nplot[i]=Pro_Nos[j]
         ot_Pig_Nplot[i]=Pig_Nos[j]
         ##########################
         ot_RNA_Pplot[i] = RNA_Pos[j]
         ot_DNA_Pplot[i] = DNA_Pos[j]
         ot_Lip_Pplot[i] = Lip_Pos[j]
         ot_Res_Pplot[i] = Res_Pos[j]
         ##########################
         ot_NC[i]=CN_os[j]
         ot_PC[i]=CP_os[j]
         ###########################
         growth_os[j] = 1e30
         
         ##### find carbon allocation elements ######  MAURE
         
    RNA_maure=RNA[15:30]
    DNA_maure=DNA[15:30]
    Pro_maure=protein[15:30]
    Pig_maure=chlorophyll[15:30]
    Carb_maure=carb[15:30]
    Lip_maure=lipid[15:30]
    Pinlipid_maure=Pinlipid[15:30]
    poly_maure=polyphosphate[15:30]

    
    RNA_C=(RNA_maure/C_maure)
    DNA_C=(DNA_maure/C_maure)
    Pro_C=(Pro_maure/C_maure)
    Pig_C=(Pig_maure/C_maure)
    Carb_C=(Carb_maure/C_maure)
    Lip_C=(Lip_maure/C_maure)
    #print(Pro_C)
    
    RNA_Cmaure=RNA_C*0.34
    DNA_Cmaure=DNA_C*0.36
    Pro_Cmaure=Pro_C*0.53
    Pig_Cmaure=Pig_C*0.74
    Carb_Cmaure=Carb_C*0.40
    Lip_Cmaure=Lip_C*0.25
    total=RNA_Cmaure+DNA_Cmaure+Pro_Cmaure+Pig_Cmaure+Carb_Cmaure+Lip_Cmaure

    RNA_Cmaure=RNA_Cmaure/total
    DNA_Cmaure=DNA_Cmaure/total
    Pro_Cmaure=Pro_Cmaure/total
    Pig_Cmaure=Pig_Cmaure/total
    Carb_Cmaure=Carb_Cmaure/total
    Lip_Cmaure=Lip_Cmaure/total
    
    Growth=size(gr_maure)
    maure_growth_plot = zeros(Growth)
    maure_NC= zeros(Growth)
    maure_PC= zeros(Growth)
    maure_NP= zeros(Growth)
    
    ##### find nitrogen allocation elements #####
    RNA_Nmaure=RNA_maure/N_maure
    DNA_Nmaure=DNA_maure/N_maure
    Pro_Nmaure=Pro_maure/N_maure
    Pig_Nmaure=Pig_maure/N_maure
    
    RNA_Nmaure=RNA_Nmaure*0.1612
    DNA_Nmaure=DNA_Nmaure*0.1684
    Pro_Nmaure=Pro_Nmaure*0.16
    Pig_Nmaure=Pig_Nmaure*0.063
       
    total=RNA_Nmaure+DNA_Nmaure+Pro_Nmaure+Pig_Nmaure
    
    RNA_Nmaure=RNA_Nmaure/total
    DNA_Nmaure=DNA_Nmaure/total
    Pro_Nmaure=Pro_Nmaure/total
    Pig_Nmaure=Pig_Nmaure/total
    
    
    ##### find phosphorus allocation elements #####
    RNA_Pmaure=RNA_maure*(1/15.5)
    DNA_Pmaure=DNA_maure*(1/15.9)
    Lip_Pmaure=Pinlipid_maure
    Res_Pmaure=poly_maure
       #carbon maure
    maure_RNA_plot = zeros(Growth)
    maure_DNA_plot = zeros(Growth)
    maure_Pro_plot = zeros(Growth)
    maure_Pig_plot = zeros(Growth)
    maure_Carb_plot = zeros(Growth)
    maure_Lip_plot = zeros(Growth)
    
    #nitrogen
    maure_RNA_Nplot = zeros(Growth)
    maure_DNA_Nplot = zeros(Growth)
    maure_Pro_Nplot = zeros(Growth)
    maure_Pig_Nplot = zeros(Growth)
   
    for i in arange(Growth):
        j=where(gr_maure==gr_maure.min())[0][0]
        maure_growth_plot[i] = gr_maure[j]
        maure_NC[i]=NC_maure[j]
        maure_PC[i]=PC_maure[j]
        maure_NP[i]=NP_maure[j]
        ######################
        maure_RNA_plot[i] = RNA_Cmaure[j]
        maure_DNA_plot[i]=DNA_Cmaure[j]
        maure_Pro_plot[i]=Pro_Cmaure[j]
        maure_Pig_plot[i]=Pig_Cmaure[j]
        maure_Carb_plot[i]=Carb_Cmaure[j]
        maure_Lip_plot[i]=Lip_Cmaure[j]
        #########################
        maure_RNA_Nplot[i] = RNA_Nmaure[j]
        maure_DNA_Nplot[i]=DNA_Nmaure[j]
        maure_Pro_Nplot[i]=Pro_Nmaure[j]
        maure_Pig_Nplot[i]=Pig_Nmaure[j]
        #########################

        gr_maure[j] = 1e30 
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #5.Plot   Here updated version; copy and past to the other ones and use $\mu$ for mu
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    rcParams.update({'font.size': 25,
                     'lines.markersize':12,
                     'lines.markeredgewidth':1})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'xtick.major.pad': 15})
    rcParams.update({'figure.autolayout': True})
    rcParams['figure.figsize']=8,6.5
    rcParams.update({'figure.facecolor':'w'})
    rcParams.update({'lines.linewidth':3})   
    
    
    rcParams.update({'axes.linewidth':1.5})
    rcParams.update({'xtick.major.width':1})
    rcParams.update({'ytick.major.width':1})
    rcParams.update({'mathtext.default': 'regular' })
    
    lowlim=0
    highlim=Ddmax+1e-5
    step=0.2
    AxisFontSize = 30
    
       
    #model_components
    carbs=Cother_plot+Cessential_plot
    lipids=CthylakoidPG_plot
    DNA=Cdna_const_plot+Cdna_variable_plot
    RNA=Crna_const_plot+Crna_variable_plot
    protein=Cnstore_plot+Cconst_protein_plot+Cphoto_plot+Cbiosynth_plot
    pigment=Cchl_plot
    carb_lip=carbs+lipids

    print(ot_growth_plot)
    #observations
    mari_data_carbs=100*mari_Carb_plot
    mari_data_lip=100*mari_Lip_plot
    mari_data_pig=100*mari_Pig_plot
    mari_data_pro=100*mari_Pro_plot
    mari_data_RNA=100*mari_RNA_plot
    mari_data_DNA=100*mari_DNA_plot
    mari_data_res=100*mari_resid_plot
    mari_data_lip_carb=mari_data_carbs+mari_data_lip       #+data_res

    we_data_carbs=100*we_Carb_plot
    we_data_lip=100*we_Lip_plot    
    we_data_pig=100*we_Pig_plot
    we_data_pro=100*we_Pro_plot
    we_data_RNA=100*we_RNA_plot
    we_data_DNA=100*we_DNA_plot
    we_data_lip_carb=we_data_carbs+we_data_lip
    
    ps_data_carbs=100*ps_Carb_plot
    ps_data_lip=100*ps_Lip_plot    
    ps_data_pig=100*ps_Pig_plot
    ps_data_pro=100*ps_Pro_plot
    ps_data_RNA=100*ps_RNA_plot
    ps_data_DNA=100*ps_DNA_plot
    ps_data_lip_carb=ps_data_carbs+ps_data_lip
    
    ot_data_carbs=100*ot_Carb_plot
    ot_data_lip=100*ot_Lip_plot    
    ot_data_pig=100*ot_Pig_plot
    ot_data_pro=100*ot_Pro_plot
    ot_data_RNA=100*ot_RNA_plot
    ot_data_DNA=100*ot_DNA_plot
    ot_data_lip_carb=ot_data_carbs+ot_data_lip
    
    maure_data_carbs=100*maure_Carb_plot
    maure_data_lip=100*maure_Lip_plot
    maure_data_pig=100*maure_Pig_plot
    maure_data_pro=100*maure_Pro_plot
    maure_data_RNA=100*maure_RNA_plot
    maure_data_DNA=100*maure_DNA_plot
    maure_data_lip_carb=maure_data_carbs+maure_data_lip
    
    #print(maure_data_lip_carb)
    
    
    cl_array=array([carb_lip[0],ps_data_lip_carb[5],ot_data_lip_carb[3],we_data_lip_carb[2]])
    pig_array=array([pigment[0],ps_data_pig[5],ot_data_pig[3],we_data_pig[2]])
    pro_array=array([protein[0],ps_data_pro[5],ot_data_pro[3],we_data_pro[2]])
    RNA_array=array([RNA[0],ps_data_RNA[5],ot_data_RNA[3],we_data_RNA[2]])
    DNA_array=array([DNA[0],ps_data_DNA[5],ot_data_DNA[3],we_data_DNA[2]])
    
    total= cl_array[0]+pig_array[0]+pro_array[0]+RNA_array[0]+DNA_array[0]
    species= ['Model','T.pseudonana','O.tauri','T.weissflogii']
    ind = arange(4)#[x for x, _ in enumerate(growth_rate)]
    #print(ind)
    width = 0.45  
    
    figure(1)
    bar(ind+width, cl_array, width=width, label='Other C', color='teal', bottom=pig_array+pro_array+RNA_array+DNA_array)
    bar(ind+width, pig_array, width=width, label='Pigment', color='#DDCC77', bottom=pro_array+RNA_array+DNA_array)
    bar(ind+width, pro_array, width=width, label='Protein', color='#CC6677', bottom=RNA_array+DNA_array)
    bar(ind+width, RNA_array, width=width, label='RNA', color='#33BBEE', bottom=DNA_array)
    bar(ind+width, DNA_array, width=width, label='DNA', color='#4B0082')
    
    xticks(ind + 2*width / 2, species)
    xticks(fontsize=12,ha='center')
    ylabel("C allocation (mol $C_{x}$ mol $C_{Total}$$^{-1}$)")
    xlabel("Species", labelpad=10)
    title('Low growth= 0.001 d$^{-1}$')
    
    
#     print(cl_array)
#     print(pig_array)
#     print(pro_array)
#     print(RNA_array)
#     print(DNA_array)

    N_pigment=Nchl_plot
    N_protein=Nconst_protein_plot+Nbiosynth_plot+Nphoto_plot
    N_RNA=Nrna_variable_plot+Nrna_const_plot
    N_DNA=Ndna_const_plot+Ndna_variable_plot
    
    N_mari_data_pig=mari_Pig_Nplot*mari_NC
    N_mari_data_pro=mari_Pro_Nplot*mari_NC
    N_mari_data_RNA=mari_RNA_Nplot*mari_NC
    N_mari_data_DNA=mari_DNA_Nplot*mari_NC
    
    N_ps_data_pig=ps_Pig_Nplot*ps_NC
    N_ps_data_pro=ps_Pro_Nplot*ps_NC
    N_ps_data_RNA=ps_RNA_Nplot*ps_NC
    N_ps_data_DNA=ps_DNA_Nplot*ps_NC
    
    N_ot_data_pig=ot_Pig_Nplot*ot_NC
    N_ot_data_pro=ot_Pro_Nplot*ot_NC
    N_ot_data_RNA=ot_RNA_Nplot*ot_NC
    N_ot_data_DNA=ot_DNA_Nplot*ot_NC
    
    N_we_data_pig=we_Pig_Nplot*we_NC
    N_we_data_pro=we_Pro_Nplot*we_NC
    N_we_data_RNA=we_RNA_Nplot*we_NC
    N_we_data_DNA=we_DNA_Nplot*we_NC
    
    N_maure_data_pig=maure_Pig_Nplot*maure_NC
    N_maure_data_pro=maure_Pro_Nplot*maure_NC
    N_maure_data_RNA=maure_RNA_Nplot*maure_NC
    N_maure_data_DNA=maure_DNA_Nplot*maure_NC

    N_pig_array=array([N_pigment[0],N_ps_data_pig[5],N_ot_data_pig[3],N_we_data_pig[2]])
    N_pro_array=array([N_protein[0],N_ps_data_pro[5],N_ot_data_pro[3],N_we_data_pro[2]])
    N_RNA_array=array([N_RNA[0],N_ps_data_RNA[5],N_ot_data_RNA[3],N_we_data_RNA[2]])
    N_DNA_array=array([N_DNA[0],N_ps_data_DNA[5],N_ot_data_DNA[3],N_we_data_DNA[2]])
    
    species= ['Model','T.pseudonana','O.tauri','T.weissflogii']
    ind = arange(4)#[x for x, _ in enumerate(growth_rate)]
    #print(ind)
    width = 0.45  
    
    figure(2)
    bar(ind+width, N_pig_array, width=width, label='Pigment', color='#DDCC77', bottom=N_pro_array+N_RNA_array+N_DNA_array)
    bar(ind+width, N_pro_array, width=width, label='Protein', color='#CC6677', bottom=N_RNA_array+N_DNA_array)
    bar(ind+width, N_RNA_array, width=width, label='RNA', color='#33BBEE', bottom=N_DNA_array)
    bar(ind+width, N_DNA_array, width=width, label='DNA', color='#4B0082')
    
    xticks(ind + 2*width / 2, species)
    xticks(fontsize=12,ha='center')
    ylabel("N allocation (mol N mol C$^{-1}$)")
    ylim(0,0.2)
    xlabel("Species", labelpad=10)
    #title('Low growth= 0.001 d$^{-1}$')
    
    

    h_cl_array=array([carb_lip[599],ps_data_lip_carb[9],ot_data_lip_carb[11],we_data_lip_carb[13]])
    h_pig_array=array([pigment[599],ps_data_pig[9],ot_data_pig[11],we_data_pig[13]])
    h_pro_array=array([protein[599],ps_data_pro[9],ot_data_pro[11],we_data_pro[13]])
    h_RNA_array=array([RNA[599],ps_data_RNA[9],ot_data_RNA[11],we_data_RNA[13]])
    h_DNA_array=array([DNA[599],ps_data_DNA[9],ot_data_DNA[11],we_data_DNA[13]])
    
    species_2= ['Model','T.pseudonana','O.tauri','T.weissflogii']
    
    figure (3)
    bar(ind+width, h_cl_array, width=width, label='Other C', color='teal', bottom=h_pig_array+h_pro_array+h_RNA_array+h_DNA_array)
    bar(ind+width, h_pig_array, width=width, label='Pigment', color='#DDCC77', bottom=h_pro_array+h_RNA_array+h_DNA_array)
    bar(ind+width, h_pro_array, width=width, label='Protein', color='#CC6677', bottom=h_RNA_array+h_DNA_array)
    bar(ind+width, h_RNA_array, width=width, label='RNA', color='#33BBEE', bottom=h_DNA_array)
    bar(ind+width, h_DNA_array, width=width, label='DNA', color='#4B0082')
    
    #print(h_cl_array)
    #print(h_pig_array)
    #print(h_pro_array)
    #print(h_RNA_array)
    #print(h_DNA_array)
    
    xticks(ind + 2*width / 2, species_2)
    xticks(fontsize=12,ha='center')
    ylabel("C allocation (mol $C_{x}$ mol $C_{Total}$$^{-1}$)")
    xlabel("Species", labelpad=10)
    title('High growth= ~0.60 d$^{-1}$')
    
    
    h_N_pig_array=array([N_pigment[599],N_ps_data_pig[9],N_ot_data_pig[11],N_we_data_pig[13]])
    h_N_pro_array=array([N_protein[599],N_ps_data_pro[9],N_ot_data_pro[11],N_we_data_pro[13]])
    h_N_RNA_array=array([N_RNA[599],N_ps_data_RNA[9],N_ot_data_RNA[11],N_we_data_RNA[13]])
    h_N_DNA_array=array([N_DNA[599],N_ps_data_DNA[9],N_ot_data_DNA[11],N_we_data_DNA[13]])
    
    figure(4)
    bar(ind+width, h_N_pig_array, width=width, label='Pigment', color='#DDCC77', bottom=h_N_pro_array+h_N_RNA_array+h_N_DNA_array)
    bar(ind+width, h_N_pro_array, width=width, label='Protein', color='#CC6677', bottom=h_N_RNA_array+h_N_DNA_array)
    bar(ind+width, h_N_RNA_array, width=width, label='RNA', color='#33BBEE', bottom=h_N_DNA_array)
    bar(ind+width, h_N_DNA_array, width=width, label='DNA', color='#4B0082')
    
    xticks(ind + 2*width / 2, species_2)
    xticks(fontsize=12,ha='center')
    ylim(0,0.2)
    ylabel("N allocation (mol N mol C$^{-1}$)")
    xlabel("Species", labelpad=10)
    #title('High growth= ~0.60 d$^{-1}$')

    
    
    

##########################
    return

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# Main part
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAA


if What_is_limiting==0:
    #for P limiting data
    from Healey85_data_08_Chl_from_Nlimited_added import healey85C
elif What_is_limiting==1:
    #for N limiting data
    from Healey85_data_09_Nlimited_case import healey85C


healey85data=healey85C()
for a in healey85data:
    kkI(a,What_is_limiting)
    

pyplot.show()
    
    