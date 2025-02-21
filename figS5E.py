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
from matplotlib.pyplot import figure, show, xticks, yticks, xlim, ylim, plot, xlabel, ylabel, title, stackplot, legend, scatter, bar
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
    Cnbiosynth=0.5*4.34728279914354E-10        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Nconst_protein=0.4*4.45336898828389E-15     #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable=19*6212.59249917364        #(s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other=5.44534485638617E-16               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max=1.5*25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
    Cessential=1.51786753491048E-15             #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%

#==============================

    I = 100                               #healey85data.Lightintensity
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
    
    Molar_mass_DNA_AT_average=307.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_DNA_CG_average=307.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    Molar_mass_RNA_AT_average=316.47        #(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    Molar_mass_RNA_CG_average=323.97        #(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
    
    #================================
    #E coli
    #================================
    CG_Ecoli=0.506          #(dimensionless) from [http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016)]
    AT_Ecoli=1-CG_Ecoli     #(dimensionless) 
    
    Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of DNA unit
    Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli     #(g mol-1) Molar mass of RNA unit
    
    RNA_DNA_mass_ratio=20/7.6   #(ug/ug) Bremer and Dennis 1996
    RNA_DNA_mass_ratio=17.844/6.5239  #(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
    
    RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli    #(mol mol-1)
    
    #================================
    #Stoichiometric parameters
    #================================
    YcyanoC_N=2                             #(molC molN) C/N molar ratio of cyanophycin
    YpgC_P=40                           #(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)

    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]

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
    print(YdnaC_N, YrnaC_N,YnucacidN_P,YdnaC_N*YnucacidN_P,YrnaC_N*YnucacidN_P)
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
    print(Nrna_const,'a')
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
    gr_maure=growth_rate[15:30]

    NP= data[1:46,9]
    NP_maure=NP[15:30]

     
     
    C= data[1:46,16]
    C_maure=C[15:30]

      #N
    N= data[1:46,17]
    N_maure=N[15:30]
   
    
    #N:C
    NC_maure=N_maure/C_maure
  
    
    #P
    P= data[1:46,18]
    P_maure=P[15:30]

    
    #PC
    PC_maure=P_maure/C_maure

    NP_maure=N_maure/P_maure
   
   
    protein=data[1:46,20]
    carb=data[1:46,21]
    lipid=data[1:46,22]
    chlorophyll=data[1:46,23]
    RNA=data[1:46,24]
    DNA=data[1:46,25]
    Pinlipid=data[1:46,26]
    polyphosphate=data[1:46,27]
    

    
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
    
     #carbon maure
    maure_RNA_plot = zeros(Growth)
    maure_DNA_plot = zeros(Growth)
    maure_Pro_plot = zeros(Growth)
    maure_Pig_plot = zeros(Growth)
    maure_Carb_plot = zeros(Growth)
    maure_Lip_plot = zeros(Growth)
    
   
   
    
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
        gr_maure[j] = 1e30 
        

    #print(maure_Pro_plot)
    
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
    
    carbs_m=array([carbs[133],carbs[158],carbs[250],carbs[366],carbs[620]])
    lipids_m=array([lipids[133],lipids[158],lipids[250],lipids[366],lipids[620]])
    carb_lip_m=array([carb_lip[133],carb_lip[158],carb_lip[250],carb_lip[366],carb_lip[620]])
    pigment_m=array([pigment[133],pigment[158],pigment[250],pigment[366],pigment[620]])
    protein_m=array([protein[133],protein[158],protein[250],protein[366],protein[620]])
    RNA_m=array([RNA[133],RNA[158],RNA[250],RNA[366],RNA[620]])
    DNA_m=array([DNA[133],DNA[158],DNA[250],DNA[366],DNA[620]])
    
    #observations
    data_carbs=100*maure_Carb_plot
    data_lip=100*maure_Lip_plot
    data_pig=100*maure_Pig_plot
    data_pro=100*maure_Pro_plot
    data_RNA=100*maure_RNA_plot
    data_DNA=100*maure_DNA_plot
    
    
    #print(data_carbs)
    data_lip_carb=data_carbs+data_lip
    
    print(data_lip_carb)
    
    carbs_d=array([data_carbs[3],data_carbs[5],data_carbs[7],data_carbs[9],data_carbs[14]])
    lipids_d=array([data_lip[3],data_lip[5],data_lip[7],data_lip[9],data_lip[14]])
    carbs_lip_d=array([data_lip_carb[3],data_lip_carb[5],data_lip_carb[7],data_lip_carb[9],data_lip_carb[14]])
    pigment_d=array([data_pig[3],data_pig[5],data_pig[7],data_pig[9],data_pig[14]])
    protein_d=array([data_pro[3],data_pro[5],data_pro[7],data_pro[9],data_pro[14]])
    RNA_d=array([data_RNA[3],data_RNA[5],data_RNA[7],data_RNA[9],data_RNA[14]])
    DNA_d=array([data_DNA[3],data_DNA[5],data_DNA[7],data_DNA[9],data_DNA[14]])

    
    growth_rate = ['0.13', '0.16', '0.25', '0.38', '0.62']
    ind = arange(5)#[x for x, _ in enumerate(growth_rate)]
    #print(ind)
    width = 0.35  
    
    bar(ind, carb_lip_m, width=width, label='Other C', color='teal', bottom=pigment_m+protein_m+RNA_m+DNA_m)
    #bar(ind, lipids_m, width=width, label='Lipids', color='#882255', bottom=pigment_m+protein_m+RNA_m+DNA_m)
    bar(ind, pigment_m, width=width, label='Pigment', color='#DDCC77', bottom=protein_m+RNA_m+DNA_m)
    bar(ind, protein_m, width=width, label='Protein', color='#CC6677', bottom=RNA_m+DNA_m)
    bar(ind, RNA_m, width=width, label='RNA', color='#33BBEE', bottom=DNA_m)
    bar(ind, DNA_m, width=width, label='DNA', color='#4B0082')
    
#     bar(ind+width, carbs_d, width=width, label='Other C', color='teal', bottom=lipids_d+pigment_d+protein_d+RNA_d+DNA_d,edgecolor='black')
#     bar(ind+width, lipids_d, width=width, label='Lipids', color='#882255', bottom=pigment_d+protein_d+RNA_d+DNA_d,edgecolor='black')
    bar(ind+width, carbs_lip_d, width=width, label='Other C', color='teal', bottom=pigment_d+protein_d+RNA_d+DNA_d,edgecolor='black')
    bar(ind+width, pigment_d, width=width, label='Pigment', color='#DDCC77', bottom=protein_d+RNA_d+DNA_d,edgecolor='black')
    bar(ind+width, protein_d, width=width, label='Protein', color='#CC6677', bottom=RNA_d+DNA_d,edgecolor='black')
    bar(ind+width, RNA_d, width=width, label='RNA', color='#33BBEE', bottom=DNA_d,edgecolor='black')
    bar(ind+width, DNA_d, width=width, label='DNA', color='#4B0082',edgecolor='black')
    
    xticks(ind + width / 2, growth_rate)
    ylabel("Fraction of C")
    xlabel("Growth rate")
    title('C.maureeniae')
    #legend(loc="upper right")
    #title("Model-Data Comparison")

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
    
    