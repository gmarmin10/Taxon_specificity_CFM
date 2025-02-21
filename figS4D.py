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

    m=2.8*3.79146798299876E-19         #(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
    Pmax=0.00320513285659728
    OT=0.00863364097132997
    Ynphoto_chl=0.5*3.56099164557551          #((molN cell-1)/(molC chl cell-1)) the stoichiometric ratio for cell photosynthetic enzyme (Rubisco etc.) nitrogen to chlorophyll (193-25)
    Cnbiosynth=0.5*4.34728279914354E-10        #(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
    Nconst_protein=0.5*4.45336898828389E-15     #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Nstore_max=2.91679384515998E-15         #(molN cell-1) Constant protein pool in nitrogen (193-25)
    Cnrna_variable=0.1*6212.59249917364        #(s) Constant for Variable part of RNA (193-26)
    Ypthylakoid_chl=0.0281633095303638        #((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
    Pconst_other=5.44534485638617E-20               #(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
    Qp_max=0.7*1.5*25.26/(3.097e16)                                                                              #(molP cell-1) total phosphorus content in the cell (193-26)
    Cessential=1.51786753491048E-15          #(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%

#==============================

    I = 85                               #healey85data.Lightintensity
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
    CG=0.59

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
    #Ndna_variable=0*Dd                     #(molN cell-1) While Bremer and Dennis shows increasing trend, Parrott 1980 shows decreasing trend. 
   
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
    growth_plot = zeros(U)
    #carbon
    RNA_plot = zeros(U)
    DNA_plot = zeros(U)
    Pro_plot = zeros(U)
    Pig_plot = zeros(U)
    Carb_plot = zeros(U)
    Lip_plot = zeros(U)
    #nitrogen
    RNA_Nplot = zeros(U)
    DNA_Nplot = zeros(U)
    Pro_Nplot = zeros(U)
    Pig_Nplot = zeros(U)
    #phosphorus
    RNA_Pplot = zeros(U)
    DNA_Pplot = zeros(U)
    Lip_Pplot = zeros(U)
    Res_Pplot = zeros(U)
    #stoichiometry
    NC=zeros(U)
    PC=zeros(U)
    
    for i in arange(U):
         j = where(growth_os==growth_os.min())[0][0]
         growth_plot[i] = growth_os[j]
         ###############################
         RNA_plot[i] = RNA_os[j]
         DNA_plot[i]=DNA_os[j]
         Pro_plot[i]=Pro_os[j]
         Pig_plot[i]=Pig_os[j]
         Carb_plot[i]=Carb_os[j]
         Lip_plot[i]=Lip_os[j]
         ################################
         ############################
         RNA_Nplot[i] = RNA_Nos[j]
         DNA_Nplot[i]=DNA_Nos[j]
         Pro_Nplot[i]=Pro_Nos[j]
         Pig_Nplot[i]=Pig_Nos[j]
         ##########################
         RNA_Pplot[i] = RNA_Pos[j]
         DNA_Pplot[i] = DNA_Pos[j]
         Lip_Pplot[i] = Lip_Pos[j]
         Res_Pplot[i] = Res_Pos[j]
         ##########################
         NC[i]=CN_os[j]
         PC[i]=CP_os[j]
         NP[i]=NP_os[j]
         ###########################
         growth_os[j] = 1e30

#     
#     growth_plot=growth_plot[5:]
#     NC=NC[5:]
#     PC=PC[5:]
#     NP=NP[5:]
   
    
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
    
 
    figure(2)
    N_pigment=Nchl_plot
    N_protein=Nconst_protein_plot+Nbiosynth_plot+Nphoto_plot
    N_RNA=Nrna_variable_plot+Nrna_const_plot
    N_DNA=Ndna_const_plot+Ndna_variable_plot
    #print(N_protein)
    
    N_pigment_m=array([ N_pigment[0], N_pigment[232], N_pigment[522], N_pigment[749], N_pigment[859]])
    N_protein_m=array([N_protein[0],N_protein[232],N_protein[522],N_protein[749],N_protein[859]])
    N_RNA_m=array([N_RNA[0],N_RNA[232],N_RNA[522],N_RNA[749],N_RNA[859]])
    N_DNA_m=array([N_DNA[0],N_DNA[232],N_DNA[522],N_DNA[749],N_DNA[859]])
     
    #observations
    N_data_pig=NC*Pig_Nplot
    N_data_pro=NC*Pro_Nplot
    N_data_RNA=NC*RNA_Nplot
    N_data_DNA=NC*DNA_Nplot
    #print(N_data_pro)

    N_pigment_d=array([N_data_pig[3],N_data_pig[8],N_data_pig[9],N_data_pig[12],N_data_pig[14]])
    N_protein_d=array([N_data_pro[3],N_data_pro[8],N_data_pro[9],N_data_pro[12],N_data_pro[14]])
    N_RNA_d=array([N_data_RNA[3],N_data_RNA[8],N_data_RNA[9],N_data_RNA[12],N_data_RNA[14]])
    N_DNA_d=array([N_data_DNA[3],N_data_DNA[8],N_data_DNA[9],N_data_DNA[12],N_data_DNA[14]])
 
    growth_rate = ['0.001', '0.23', '0.52', '0.75', '0.86'] 
    ind = arange(5)#[x for x, _ in enumerate(growth_rate)]
    #print(ind)
    width = 0.35  
     
    bar(ind, N_pigment_m, width=width, label='Pigment', color='#DDCC77', bottom=N_protein_m+N_RNA_m+N_DNA_m)
    bar(ind, N_protein_m, width=width, label='Protein', color='#CC6677', bottom=N_RNA_m+N_DNA_m)
    bar(ind, N_RNA_m, width=width, label='RNA', color='#33BBEE', bottom=N_DNA_m)
    bar(ind, N_DNA_m, width=width, label='DNA', color='#4B0082')
    
    
    
    bar(ind+width, N_pigment_d, width=width, label='Pigment', color='#DDCC77', bottom=N_protein_d+N_RNA_d+N_DNA_d,edgecolor='black')
    bar(ind+width, N_protein_d, width=width, label='Protein', color='#CC6677', bottom=N_RNA_d+N_DNA_d,edgecolor='black')
    bar(ind+width, N_RNA_d, width=width, label='RNA', color='#33BBEE', bottom=N_DNA_d,edgecolor='black')
    bar(ind+width, N_DNA_d, width=width, label='DNA', color='#4B0082',edgecolor='black')
    
    #print(N_DNA_d)
    #print(N_DNA_m)
    
    print(N_pigment_m)
    print(N_pigment_d)
    print(N_protein_m)
    print(N_protein_d)
    print(N_RNA_m)
    print(N_RNA_d)
    print(N_DNA_m)
    print(N_DNA_d)
    
    
    xticks(ind + width / 2, growth_rate)
    yticks()
    ylabel("N:C")
    title('O.tauri')
    xlabel("Growth rate")
    #c_leg=mpat.Patch(color='teal',label='Other C',alpha=0.75)
    pig_leg=mpat.Patch(color='#DDCC77',label='Pigment',alpha=0.75)
    Pro_leg=mpat.Patch(color='#CC6677',label='Protein',alpha=0.75)
    RNA_leg=mpat.Patch(color='#33BBEE',label='RNA',alpha=0.75)
    DNA_leg=mpat.Patch(color='#4B0082',label='DNA',alpha=0.75)

    handles=[pig_leg,Pro_leg,RNA_leg,DNA_leg]
    #pyplot.legend(handles=handles,loc='upper center',bbox_to_anchor=(0.65,-0.35), ncol=4,fontsize='x-small',frameon=False)
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
    
    