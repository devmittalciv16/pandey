from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.
import numpy as np
nnn={"image":0}
beamEnergy =  639.5; # in eV
properties = np.array(([1, 0,   0,         0,       0,        0,       28,        0,        0,       0,      0],
                       [2, 0,  31.00,      11.0,    0,      8.31E-24,  312.84,    8.5,      2.22,     0,     0],
                       [3, 1,  15.00,      11,     7.5,    6.55E-24,  224.723,  -30+54,    37.75,   -0.15,  0.05],
                       [4, 1,  26.00,      10,      7.2,    6.55E-24,  224.723,  -32+50,    40.75,   -0.65,  0.46],
                       [5, 1,  23.04,      7.00,    8.2,    6.55E-24,  224.7237, -32+50,    40.75,   -0.25,  0.55], 
                       [6, 1,  20.04,      4.00,    5,      6.55E-24,  224.7237,  -22+38,   35.75,   -0.08,  0.15], 
                       [7, 0,   'Inf',       5.91,    0,      2.33E-24,  182.487,   21.138,   12.10,      0,    0]),
                       dtype=object)

#magicConstant1 = np.array([1,0.9,1.02,1.05,0.9,1,1])  # % for charge scattering factors
v1=1;v2=0.9;v3=1.02;v4=1.05;v5=0.9;v6=1;v7=1;
magicConstant1 = [v1,v2,v3,v4,v5,v6,v7]
v8=1.15;v9=0.30;v10=0.25;v11=1.85;v12=1.55;v13=1.25;v14=1.25; 
magicConstant2 = [v8,v9,v10,v11,v12,v13,v14]   # for magnetic scattering factors%
bob = 'd'; 

def calci():
    
     # 'd' means plots asym ratio, if empty plots plus/minus curve
    farbe = 'k';
    thetaRangeDeg = np.arange(0.0,34.001, 0.2);

    import math
    pi=math.pi
    lambda1 = 12398.42 / beamEnergy; # in Angstrom
    kZero = 2*pi / lambda1; # in Angstrom^-1
    tZeroPlus = 1 /math.sqrt(2) * np.array((1,1j),dtype=complex);
    tZeroMinus = 1 / math.sqrt(2) * np.array((1,-1j),dtype=complex);
    RZERO = 2.81794E-5; # classical electron radius in Angstrom
    NA = 6.022E23; # Avogadro's constant

    # transform thetaRange
    thetaRangeRad = thetaRangeDeg * pi / 180;
    QzRange = 4 * pi / lambda1 * np.sin(thetaRangeRad);

    #read out properties
    nLayers = len(properties);
    thick = properties[:,2];
    chemRough2 = np.square(properties[:,3]);
    magRough2 = np.square(properties[:,4]);
    rho = properties[:,5];
    molW = properties[:,6];
    fc = properties[:,7:9];
    fm = properties[:,9:11];

    #sorted

    # calculate more properties
    rhoAtom = NA * rho/ molW; # in Angstrom^-3
    delta = lambda1**2 / (2 * pi) * RZERO * rhoAtom * fc[:,0];
    delta.astype('complex');
    beta = lambda1**2 / (2 * pi) * RZERO * rhoAtom * fc[:,1];
    beta.astype('complex');
    refrIdx = 1 - delta + 1j * beta;

    A = (-fc[:,0] + 1j*fc[:,1])* magicConstant1 ;
    B = (fm[:,0] + 1j*fm[:,1]) * magicConstant1* magicConstant2 ;

    # find which layers are magnetic
    magIdx = np.where(properties[:,1]==1)[0];
    nonMagIdx = np.where(properties[:,1]==0)[0];

    nonMagIdx.shape

    intPlus= np.zeros((len(thetaRangeRad)));
    intMinus= np.zeros((len(thetaRangeRad))); 
    asymmRatio = np.zeros((len(thetaRangeRad)));
    intdiff = np.zeros((len(thetaRangeRad)));
    iTheta = -1;

    for theta in thetaRangeRad:
      iTheta = iTheta+1;
      #calculate uPlus, uMinus...
      uPlus= np.zeros((nLayers),dtype=complex);
      uMinus= np.zeros((nLayers),dtype=complex);
      uZero = np.zeros((nLayers),dtype=complex);
    
      # ... non-magnetic
      if nonMagIdx.any():
        Kemcho=np.sin(theta)**2 - 2 * delta[nonMagIdx] + 2j* beta[nonMagIdx];
        Kemcho=Kemcho.astype('complex');
        uPlus[nonMagIdx] = np.sqrt(Kemcho);
        uMinus[nonMagIdx] = uPlus[nonMagIdx];
      #... magnetic
      if magIdx.any():
        Kemcho1= np.sin(theta)**2 + 4 * pi / (kZero**2)  * RZERO * rhoAtom[magIdx] * (A[magIdx] + B[magIdx]);
        Kemcho1=Kemcho1.astype('complex');
        uPlus[magIdx] = np.sqrt(Kemcho1);
        Kemcho2= np.sin(theta)**2 + 4 * pi / (kZero**2) * RZERO * rhoAtom[magIdx]*(A[magIdx] - B[magIdx]);
        Kemcho2=Kemcho2.astype('complex');
        uMinus[magIdx] = np.sqrt(Kemcho2);
        
      invFmatrix = np.zeros((2,2,nLayers-1),dtype=complex);
      x1=1j * kZero * np.multiply(uPlus[0:len(uPlus)-1],thick[0:len(thick)-1]);
      x1=x1.astype('complex');
      invFmatrix[0,0,:]=np.exp(x1);
      x2=1j * kZero * np.multiply(uMinus[0:len(uMinus)-1],thick[0:len(thick)-1]);
      x2=x2.astype('complex');
      invFmatrix[1,1,:] = np.exp(x2); 

      Mtt= np.zeros((2,2,nLayers-1),dtype=complex);
      Mtr= np.zeros((2,2,nLayers-1),dtype=complex);
      Mrt= np.zeros((2,2,nLayers-1),dtype=complex);
      Mrr= np.zeros((2,2,nLayers-1),dtype=complex);

      transition = (properties[1:len(properties),1] == 1) + 2*(properties[0:len(properties)-1,1] == 1);
      x3=np.arange(0,nLayers-1)

      for iTop in x3:
        # iTop 1 is between layer 1 and 2 etc.
        # iBot = iTop + 1;

        if (transition[iTop]==0):
          uZeroTop = uPlus[iTop];
          uZeroBot = uPlus[iTop+1];

          Mtt[:,:,iTop] = (2 * uZeroTop / (uZeroTop + uZeroBot) * np.exp( (kZero**2) / 2 * ((uZeroTop - uZeroBot)**2) * chemRough2[iTop+1]) ) * np.identity(2);
          Mrr[:,:,iTop] = (2 * uZeroBot / (uZeroTop + uZeroBot) * np.exp( kZero**2 / 2 * ((uZeroTop - uZeroBot)**2) * chemRough2[iTop+1]) ) * np.identity(2);
          Mtr[:,:,iTop] = ((uZeroBot - uZeroTop) / (uZeroTop + uZeroBot) * np.exp( -2 * (kZero**2) * uZeroTop * uZeroBot * chemRough2[iTop+1]) ) *np.identity(2);
          Mrt[:,:,iTop] = -Mtr[:,:,iTop];


        if (transition[iTop]==1):
          uZero = uPlus[iTop];
          uP = uPlus[iTop+1];
          uM = uMinus[iTop+1];

          factA = (uP**2 + uM**2) / 2 - uZero**2;
          factB = (uP**2 - uM**2) / 2;

          D1p = factA *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * chemRough2[iTop+1]);
          D1m = factA *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * chemRough2[iTop+1]);

          D3p = factA *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * chemRough2[iTop+1]);
          D3m = factA *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * chemRough2[iTop+1]);

          D2p = factB *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * magRough2[iTop+1]);
          D2m = factB *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * magRough2[iTop+1]);

          D4p = factB *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * magRough2[iTop+1]);
          D4m = factB *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * magRough2[iTop+1]);
        
          #CorrectTill Now

          # n -> r
          RzeroUnr= np.zeros((2,2),dtype=complex);
          IVnr= np.zeros((2,2),dtype=complex);
          IVpNr= np.zeros((2,2),dtype=complex);
          TzeroNr = np.zeros((2,2),dtype=complex);

          RzeroUnr[0,0] = (D1p + D2p) / ((uP + uZero)**2) + (D3p - D4p) / ((uM + uZero)**2);
          RzeroUnr[1,1] = RzeroUnr[0,0];
          RzeroUnr[0,1] = -1j * ((D1p + D2p) / (uP + uZero)**2 - (D3p - D4p) / (uM + uZero)**2);
          RzeroUnr[1,0] = -RzeroUnr[0,1];
          RzeroUnr = -0.5 * RzeroUnr;

          IVnr[0,0] = (D1m + D2m) / (uP**2 - uZero**2) + (D3m - D4m) / (uM**2 - uZero**2);
          IVnr[1,1] = IVnr[0,0];
          IVnr[0,1] = -1j * ((D1m + D2m) / (uP**2 - uZero**2) - (D3m - D4m) / (uM**2 - uZero**2));
          IVnr[1,0] = -IVnr[0,1];
          IVnr = 0.5 * IVnr;

          IVpNr[0,0] = (D1m + D2m) / (uP**2 - uZero**2);
          IVpNr[1,1] = (D3m - D4m) / (uM**2 - uZero**2);

          TzeroNr[0,0] = uZero / (uP + uZero);
          TzeroNr[0,1] = -1j * (uZero / (uP + uZero));
          TzeroNr[1,0] = uZero / (uM + uZero);
          TzeroNr[1,1] = 1j * (uZero / (uM + uZero));

          Mrt[:,:,iTop] = np.matmul(np.linalg.inv(IVnr), RzeroUnr);
          Mtt[:,:,iTop] = np.matmul(np.linalg.inv(IVpNr), TzeroNr);
        
          # r -> n

          RzeroUrn= np.zeros((2,2),dtype=complex);
          IVrn = np.zeros((2,2),dtype=complex);
          IVpRn = np.zeros((2,2),dtype=complex);
          TzeroRn= np.zeros((2,2),dtype=complex);

          RzeroUrn[0,0] = (D1p + D2p) / (uP + uZero)**2;
          RzeroUrn[1,1] = (D3p - D4p) / (uM + uZero)**2;

          IVrn = IVpNr;

          IVpRn = IVnr;

          TzeroRn[0,0] = 2 * uP / (uP + uZero);
          TzeroRn[0,1] = 2 * uM / (uM + uZero);
          TzeroRn[1,0] = 1j * 2 * uP / (uP + uZero);
          TzeroRn[1,1] = -1j * 2 * uM / (uM + uZero);

          Mtr[:,:,iTop] = np.matmul(np.linalg.inv(IVrn), RzeroUrn);
          Mrr[:,:,iTop] = np.matmul(np.linalg.inv(IVpRn),TzeroRn);

        if (transition[iTop]==2):

          uZero = uPlus[iTop+1];
          uP = uPlus[iTop];
          uM = uMinus[iTop];

          factA = (uP**2 + uM**2) / 2 - uZero**2;
          factB = (uP**2 - uM**2) / 2;

          D1p = factA *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * chemRough2[iTop+1]);
          D1m = factA *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * chemRough2[iTop+1]);

          D3p = factA *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * chemRough2[iTop+1]);
          D3m = factA *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * chemRough2[iTop+1]);

          D2p = factB *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * magRough2[iTop+1]);
          D2m = factB *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * magRough2[iTop+1]);

          D4p = factB *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * magRough2[iTop+1]);
          D4m = factB *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * magRough2[iTop+1]);
        
          #CorrectTill Now

          # n -> r
          RzeroUnr= np.zeros((2,2),dtype=complex);
          IVnr= np.zeros((2,2),dtype=complex);
          IVpNr= np.zeros((2,2),dtype=complex);
          TzeroNr = np.zeros((2,2),dtype=complex);

          RzeroUnr[0,0] = (D1p + D2p) / ((uP + uZero)**2) + (D3p - D4p) / ((uM + uZero)**2);
          RzeroUnr[1,1] = RzeroUnr[0,0];
          RzeroUnr[0,1] = -1j * ((D1p + D2p) / (uP + uZero)**2 - (D3p - D4p) / (uM + uZero)**2);
          RzeroUnr[1,0] = -RzeroUnr[0,1];
          RzeroUnr = -0.5 * RzeroUnr;

          IVnr[0,0] = (D1m + D2m) / (uP**2 - uZero**2) + (D3m - D4m) / (uM**2 - uZero**2);
          IVnr[1,1] = IVnr[0,0];
          IVnr[0,1] = -1j * ((D1m + D2m) / (uP**2 - uZero**2) - (D3m - D4m) / (uM**2 - uZero**2));
          IVnr[1,0] = -IVnr[0,1];
          IVnr = 0.5 * IVnr;

          IVpNr[0,0] = (D1m + D2m) / (uP**2 - uZero**2);
          IVpNr[1,1] = (D3m - D4m) / (uM**2 - uZero**2);

          TzeroNr[0,0] = uZero / (uP + uZero);
          TzeroNr[0,1] = -1j * (uZero / (uP + uZero));
          TzeroNr[1,0] = uZero / (uM + uZero);
          TzeroNr[1,1] = 1j * (uZero / (uM + uZero));

          Mtr[:,:,iTop] = np.matmul(np.linalg.inv(IVnr), RzeroUnr);
          Mrr[:,:,iTop] = np.matmul(np.linalg.inv(IVpNr), TzeroNr);
        
          # r -> n

          RzeroUrn= np.zeros((2,2),dtype=complex);
          IVrn = np.zeros((2,2),dtype=complex);
          IVpRn = np.zeros((2,2),dtype=complex);
          TzeroRn= np.zeros((2,2),dtype=complex);

          RzeroUrn[0,0] = (D1p + D2p) / (uP + uZero)**2;
          RzeroUrn[1,1] = (D3p - D4p) / (uM + uZero)**2;

          IVrn = IVpNr;

          IVpRn = IVnr;

          TzeroRn[0,0] = 2 * uP / (uP + uZero);
          TzeroRn[0,1] = 2 * uM / (uM + uZero);
          TzeroRn[1,0] = 1j * 2 * uP / (uP + uZero);
          TzeroRn[1,1] = -1j * 2 * uM / (uM + uZero);

          Mrt[:,:,iTop] = np.matmul(np.linalg.inv(IVrn), RzeroUrn);
          Mtt[:,:,iTop] = np.matmul(np.linalg.inv(IVpRn),TzeroRn);

        if (transition[iTop]==3):

          uPup = uPlus[iTop];
          uPdw = uPlus[iTop+1];
          uMup = uMinus[iTop];
          uMdw = uMinus[iTop+1];

          factA = (uPdw**2 + uMdw**2) / 2 - (uPup**2 + uMup**2) / 2;
          factB = (uPdw**2 - uMdw**2) / 2 - (uPup**2 - uMup**2) / 2;

          D5p = factA * np.exp(-kZero**2 / 2 * ((uPdw + uPup)**2) * chemRough2[iTop+1]);
          D5m = factA * np.exp(-kZero**2 / 2 * ((uPdw - uPup)**2) * chemRough2[iTop+1]);

          D7p = factA * np.exp(-kZero**2 / 2 * ((uMdw + uMup)**2) * chemRough2[iTop+1]);
          D7m = factA * np.exp(-kZero**2 / 2 * ((uMdw - uMup)**2) * chemRough2[iTop+1]);

          D6p = factB * np.exp(-kZero**2 / 2 * ((uPdw + uPup)**2) * magRough2[iTop+1]);
          D6m = factB * np.exp(-kZero**2 / 2 * ((uPdw - uPup)**2) * magRough2[iTop+1]);

          D8p = factB * np.exp(-kZero**2 / 2 * ((uMdw + uMup)**2) * magRough2[iTop+1]);
          D8m = factB * np.exp(-kZero**2 / 2 * ((uMdw - uMup)**2) * magRough2[iTop+1]);
        
          RzeroU= np.zeros((2,2),dtype=complex);
          IV= np.zeros((2,2),dtype=complex);
          IVp= np.zeros((2,2),dtype=complex);
          Tzero1= np.zeros((2,2),dtype=complex);
          Tzero2= np.zeros((2,2),dtype=complex);

          RzeroU[0,0] = -(D5p + D6p) / (uPdw + uPup)**2;
          RzeroU[1,1] = -(D7p - D8p) / (uMdw + uMup)**2;

          IV[0,0] = (D5m + D6m) / (uPdw**2 - uPup**2);
          IV[1,1] = (D7m - D8m) / (uMdw**2 - uMup**2);
        
          IVp = IV; #switching upper and lower layers switches signs in D's and the u-sums

          Tzero1[0,0] = 2 * uPup / (uPdw + uPup);
          Tzero1[1,1] = 2 * uMup / (uMdw + uMup);
          Tzero2[0,0] = 2 * uPdw / (uPdw + uPup);
          Tzero2[1,1] = 2 * uMdw / (uMdw + uMup);

          Mrt[:,:,iTop] = np.matmul(np.linalg.inv(IV), RzeroU);
          Mtt[:,:,iTop] = np.matmul(np.linalg.inv(IVp) , Tzero1);

          Mtr[:,:,iTop] = -Mrt[:,:,iTop];
          Mrr[:,:,iTop] = np.matmul(np.linalg.inv(IV) , Tzero2);

        # multiply all M-matrices by the inverse F's

        Mtt[:,:,iTop] = np.matmul(Mtt[:,:,iTop] , invFmatrix[:,:,iTop]);
        # Mtr = Mtr;
        Mrt[:,:,iTop] = np.matmul(invFmatrix[:,:,iTop], np.matmul(Mrt[:,:,iTop], invFmatrix[:,:,iTop]));
        Mrr[:,:,iTop] = np.matmul(invFmatrix[:,:,iTop], Mrr[:,:,iTop]);

      # W-MATRICES
      Wtt= np.zeros((2,2,nLayers),dtype=complex); 
      Wtr= np.zeros((2,2,nLayers),dtype=complex);
      Wrt= np.zeros((2,2,nLayers),dtype=complex);
      Wrr = np.zeros((2,2,nLayers),dtype=complex);
    
      Wtt[:,:,0] = np.identity(2);
      Wrr[:,:,0] = np.identity(2);

      y1= np.arange(1,nLayers);

      for iBot in y1:
        iTop = iBot - 1;
        
        Wtt[:,:,iBot] = np.matmul(Mtt[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) - np.matmul(Wtr[:,:,iTop] , Mrt[:,:,iTop])), Wtt[:,:,iTop]));
        Wtr[:,:,iBot] = Mtr [:,:,iTop] + np.matmul(Mtt[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) -np.matmul(Wtr[:,:,iTop],Mrt[:,:,iTop])),np.matmul(Wtr[:,:,iTop], Mrr[:,:,iTop])));
        Wrr[:,:,iBot] = np.matmul(Wrr[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) - np.matmul(Mrt[:,:,iTop] , Wtr[:,:,iTop])), Mrr[:,:,iTop]));
        Wrt[:,:,iBot] = Wrt[:,:,iTop] + np.matmul(Wrr[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) -np.matmul(Mrt[:,:,iTop],Wtr[:,:,iTop])),np.matmul(Mrt[:,:,iTop], Wtt[:,:,iTop])));
        
      WrtEnd = Wrt[:,:,-1];

      intPlus[iTheta] = sum(np.square(abs(np.matmul(WrtEnd,tZeroPlus))));
      intMinus[iTheta] = sum(np.square(abs(np.matmul(WrtEnd,tZeroMinus))));
      asymmRatio[iTheta]= (intPlus[iTheta]-intMinus[iTheta])/(intPlus[iTheta]+intMinus[iTheta]);
      intdiff[iTheta]= (intPlus[iTheta]-intMinus[iTheta]);

    import matplotlib.pyplot as plt
    if (bob=='d'):
       plt.plot(thetaRangeDeg, asymmRatio, 'k-',linewidth=2);
       
       plt.savefig('./images/fig.png',transparent=False)
       #plt.show()
       plt.close()
    else:
       plt.semilogy(thetaRangeDeg,(intPlus+intMinus)/2,'k-',LineWidth=2);
     
       plt.savefig('./images/fig.png',transparent=False)
       #plt.show()
       plt.close()
   






    

# bob = 'd';  # 'd' means plots asym ratio, if empty plots plus/minus curve
# farbe = 'k';
# thetaRangeDeg = np.arange(0.0,34.001, 0.2);

# import math
# pi=math.pi
# lambda1 = 12398.42 / beamEnergy; # in Angstrom
# kZero = 2*pi / lambda1; # in Angstrom^-1
# tZeroPlus = 1 /math.sqrt(2) * np.array((1,1j),dtype=complex);
# tZeroMinus = 1 / math.sqrt(2) * np.array((1,-1j),dtype=complex);
# RZERO = 2.81794E-5; # classical electron radius in Angstrom
# NA = 6.022E23; # Avogadro's constant

# # transform thetaRange
# thetaRangeRad = thetaRangeDeg * pi / 180;
# QzRange = 4 * pi / lambda1 * np.sin(thetaRangeRad);

# #read out properties
# nLayers = len(properties);
# thick = properties[:,2];
# chemRough2 = np.square(properties[:,3]);
# magRough2 = np.square(properties[:,4]);
# rho = properties[:,5];
# molW = properties[:,6];
# fc = properties[:,7:9];
# fm = properties[:,9:11];

# #sorted

# # calculate more properties
# rhoAtom = NA * rho/ molW; # in Angstrom^-3
# delta = lambda1**2 / (2 * pi) * RZERO * rhoAtom * fc[:,0];
# delta.astype('complex');
# beta = lambda1**2 / (2 * pi) * RZERO * rhoAtom * fc[:,1];
# beta.astype('complex');
# refrIdx = 1 - delta + 1j * beta;

# A = (-fc[:,0] + 1j*fc[:,1])* magicConstant1 ;
# B = (fm[:,0] + 1j*fm[:,1]) * magicConstant1* magicConstant2 ;

# # find which layers are magnetic
# magIdx = np.where(properties[:,1]==1)[0];
# nonMagIdx = np.where(properties[:,1]==0)[0];

# nonMagIdx.shape

# intPlus= np.zeros((len(thetaRangeRad)));
# intMinus= np.zeros((len(thetaRangeRad))); 
# asymmRatio = np.zeros((len(thetaRangeRad)));
# intdiff = np.zeros((len(thetaRangeRad)));
# iTheta = -1;

# for theta in thetaRangeRad:
#   iTheta = iTheta+1;
#   #calculate uPlus, uMinus...
#   uPlus= np.zeros((nLayers),dtype=complex);
#   uMinus= np.zeros((nLayers),dtype=complex);
#   uZero = np.zeros((nLayers),dtype=complex);
  
#   # ... non-magnetic
#   if nonMagIdx.any():
#     Kemcho=np.sin(theta)**2 - 2 * delta[nonMagIdx] + 2j* beta[nonMagIdx];
#     Kemcho=Kemcho.astype('complex');
#     uPlus[nonMagIdx] = np.sqrt(Kemcho);
#     uMinus[nonMagIdx] = uPlus[nonMagIdx];
#   #... magnetic
#   if magIdx.any():
#     Kemcho1= np.sin(theta)**2 + 4 * pi / (kZero**2)  * RZERO * rhoAtom[magIdx] * (A[magIdx] + B[magIdx]);
#     Kemcho1=Kemcho1.astype('complex');
#     uPlus[magIdx] = np.sqrt(Kemcho1);
#     Kemcho2= np.sin(theta)**2 + 4 * pi / (kZero**2) * RZERO * rhoAtom[magIdx]*(A[magIdx] - B[magIdx]);
#     Kemcho2=Kemcho2.astype('complex');
#     uMinus[magIdx] = np.sqrt(Kemcho2);
    
#   invFmatrix = np.zeros((2,2,nLayers-1),dtype=complex);
#   x1=1j * kZero * np.multiply(uPlus[0:len(uPlus)-1],thick[0:len(thick)-1]);
#   x1=x1.astype('complex');
#   invFmatrix[0,0,:]=np.exp(x1);
#   x2=1j * kZero * np.multiply(uMinus[0:len(uMinus)-1],thick[0:len(thick)-1]);
#   x2=x2.astype('complex');
#   invFmatrix[1,1,:] = np.exp(x2); 

#   Mtt= np.zeros((2,2,nLayers-1),dtype=complex);
#   Mtr= np.zeros((2,2,nLayers-1),dtype=complex);
#   Mrt= np.zeros((2,2,nLayers-1),dtype=complex);
#   Mrr= np.zeros((2,2,nLayers-1),dtype=complex);

#   transition = (properties[1:len(properties),1] == 1) + 2*(properties[0:len(properties)-1,1] == 1);
#   x3=np.arange(0,nLayers-1)

#   for iTop in x3:
#     # iTop 1 is between layer 1 and 2 etc.
#     # iBot = iTop + 1;

#     if (transition[iTop]==0):
#       uZeroTop = uPlus[iTop];
#       uZeroBot = uPlus[iTop+1];

#       Mtt[:,:,iTop] = (2 * uZeroTop / (uZeroTop + uZeroBot) * np.exp( (kZero**2) / 2 * ((uZeroTop - uZeroBot)**2) * chemRough2[iTop+1]) ) * np.identity(2);
#       Mrr[:,:,iTop] = (2 * uZeroBot / (uZeroTop + uZeroBot) * np.exp( kZero**2 / 2 * ((uZeroTop - uZeroBot)**2) * chemRough2[iTop+1]) ) * np.identity(2);
#       Mtr[:,:,iTop] = ((uZeroBot - uZeroTop) / (uZeroTop + uZeroBot) * np.exp( -2 * (kZero**2) * uZeroTop * uZeroBot * chemRough2[iTop+1]) ) *np.identity(2);
#       Mrt[:,:,iTop] = -Mtr[:,:,iTop];


#     if (transition[iTop]==1):
#       uZero = uPlus[iTop];
#       uP = uPlus[iTop+1];
#       uM = uMinus[iTop+1];

#       factA = (uP**2 + uM**2) / 2 - uZero**2;
#       factB = (uP**2 - uM**2) / 2;

#       D1p = factA *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * chemRough2[iTop+1]);
#       D1m = factA *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * chemRough2[iTop+1]);

#       D3p = factA *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * chemRough2[iTop+1]);
#       D3m = factA *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * chemRough2[iTop+1]);

#       D2p = factB *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * magRough2[iTop+1]);
#       D2m = factB *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * magRough2[iTop+1]);

#       D4p = factB *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * magRough2[iTop+1]);
#       D4m = factB *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * magRough2[iTop+1]);
      
#       #CorrectTill Now

#       # n -> r
#       RzeroUnr= np.zeros((2,2),dtype=complex);
#       IVnr= np.zeros((2,2),dtype=complex);
#       IVpNr= np.zeros((2,2),dtype=complex);
#       TzeroNr = np.zeros((2,2),dtype=complex);

#       RzeroUnr[0,0] = (D1p + D2p) / ((uP + uZero)**2) + (D3p - D4p) / ((uM + uZero)**2);
#       RzeroUnr[1,1] = RzeroUnr[0,0];
#       RzeroUnr[0,1] = -1j * ((D1p + D2p) / (uP + uZero)**2 - (D3p - D4p) / (uM + uZero)**2);
#       RzeroUnr[1,0] = -RzeroUnr[0,1];
#       RzeroUnr = -0.5 * RzeroUnr;

#       IVnr[0,0] = (D1m + D2m) / (uP**2 - uZero**2) + (D3m - D4m) / (uM**2 - uZero**2);
#       IVnr[1,1] = IVnr[0,0];
#       IVnr[0,1] = -1j * ((D1m + D2m) / (uP**2 - uZero**2) - (D3m - D4m) / (uM**2 - uZero**2));
#       IVnr[1,0] = -IVnr[0,1];
#       IVnr = 0.5 * IVnr;

#       IVpNr[0,0] = (D1m + D2m) / (uP**2 - uZero**2);
#       IVpNr[1,1] = (D3m - D4m) / (uM**2 - uZero**2);

#       TzeroNr[0,0] = uZero / (uP + uZero);
#       TzeroNr[0,1] = -1j * (uZero / (uP + uZero));
#       TzeroNr[1,0] = uZero / (uM + uZero);
#       TzeroNr[1,1] = 1j * (uZero / (uM + uZero));

#       Mrt[:,:,iTop] = np.matmul(np.linalg.inv(IVnr), RzeroUnr);
#       Mtt[:,:,iTop] = np.matmul(np.linalg.inv(IVpNr), TzeroNr);
      
#       # r -> n

#       RzeroUrn= np.zeros((2,2),dtype=complex);
#       IVrn = np.zeros((2,2),dtype=complex);
#       IVpRn = np.zeros((2,2),dtype=complex);
#       TzeroRn= np.zeros((2,2),dtype=complex);

#       RzeroUrn[0,0] = (D1p + D2p) / (uP + uZero)**2;
#       RzeroUrn[1,1] = (D3p - D4p) / (uM + uZero)**2;

#       IVrn = IVpNr;

#       IVpRn = IVnr;

#       TzeroRn[0,0] = 2 * uP / (uP + uZero);
#       TzeroRn[0,1] = 2 * uM / (uM + uZero);
#       TzeroRn[1,0] = 1j * 2 * uP / (uP + uZero);
#       TzeroRn[1,1] = -1j * 2 * uM / (uM + uZero);

#       Mtr[:,:,iTop] = np.matmul(np.linalg.inv(IVrn), RzeroUrn);
#       Mrr[:,:,iTop] = np.matmul(np.linalg.inv(IVpRn),TzeroRn);

#     if (transition[iTop]==2):

#       uZero = uPlus[iTop+1];
#       uP = uPlus[iTop];
#       uM = uMinus[iTop];

#       factA = (uP**2 + uM**2) / 2 - uZero**2;
#       factB = (uP**2 - uM**2) / 2;

#       D1p = factA *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * chemRough2[iTop+1]);
#       D1m = factA *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * chemRough2[iTop+1]);

#       D3p = factA *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * chemRough2[iTop+1]);
#       D3m = factA *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * chemRough2[iTop+1]);

#       D2p = factB *np.exp(-kZero**2 / 2 * ((uP + uZero)**2) * magRough2[iTop+1]);
#       D2m = factB *np.exp(-kZero**2 / 2 * ((uP - uZero)**2) * magRough2[iTop+1]);

#       D4p = factB *np.exp(-kZero**2 / 2 * ((uM + uZero)**2) * magRough2[iTop+1]);
#       D4m = factB *np.exp(-kZero**2 / 2 * ((uM - uZero)**2) * magRough2[iTop+1]);
      
#       #CorrectTill Now

#       # n -> r
#       RzeroUnr= np.zeros((2,2),dtype=complex);
#       IVnr= np.zeros((2,2),dtype=complex);
#       IVpNr= np.zeros((2,2),dtype=complex);
#       TzeroNr = np.zeros((2,2),dtype=complex);

#       RzeroUnr[0,0] = (D1p + D2p) / ((uP + uZero)**2) + (D3p - D4p) / ((uM + uZero)**2);
#       RzeroUnr[1,1] = RzeroUnr[0,0];
#       RzeroUnr[0,1] = -1j * ((D1p + D2p) / (uP + uZero)**2 - (D3p - D4p) / (uM + uZero)**2);
#       RzeroUnr[1,0] = -RzeroUnr[0,1];
#       RzeroUnr = -0.5 * RzeroUnr;

#       IVnr[0,0] = (D1m + D2m) / (uP**2 - uZero**2) + (D3m - D4m) / (uM**2 - uZero**2);
#       IVnr[1,1] = IVnr[0,0];
#       IVnr[0,1] = -1j * ((D1m + D2m) / (uP**2 - uZero**2) - (D3m - D4m) / (uM**2 - uZero**2));
#       IVnr[1,0] = -IVnr[0,1];
#       IVnr = 0.5 * IVnr;

#       IVpNr[0,0] = (D1m + D2m) / (uP**2 - uZero**2);
#       IVpNr[1,1] = (D3m - D4m) / (uM**2 - uZero**2);

#       TzeroNr[0,0] = uZero / (uP + uZero);
#       TzeroNr[0,1] = -1j * (uZero / (uP + uZero));
#       TzeroNr[1,0] = uZero / (uM + uZero);
#       TzeroNr[1,1] = 1j * (uZero / (uM + uZero));

#       Mtr[:,:,iTop] = np.matmul(np.linalg.inv(IVnr), RzeroUnr);
#       Mrr[:,:,iTop] = np.matmul(np.linalg.inv(IVpNr), TzeroNr);
      
#       # r -> n

#       RzeroUrn= np.zeros((2,2),dtype=complex);
#       IVrn = np.zeros((2,2),dtype=complex);
#       IVpRn = np.zeros((2,2),dtype=complex);
#       TzeroRn= np.zeros((2,2),dtype=complex);

#       RzeroUrn[0,0] = (D1p + D2p) / (uP + uZero)**2;
#       RzeroUrn[1,1] = (D3p - D4p) / (uM + uZero)**2;

#       IVrn = IVpNr;

#       IVpRn = IVnr;

#       TzeroRn[0,0] = 2 * uP / (uP + uZero);
#       TzeroRn[0,1] = 2 * uM / (uM + uZero);
#       TzeroRn[1,0] = 1j * 2 * uP / (uP + uZero);
#       TzeroRn[1,1] = -1j * 2 * uM / (uM + uZero);

#       Mrt[:,:,iTop] = np.matmul(np.linalg.inv(IVrn), RzeroUrn);
#       Mtt[:,:,iTop] = np.matmul(np.linalg.inv(IVpRn),TzeroRn);

#     if (transition[iTop]==3):

#       uPup = uPlus[iTop];
#       uPdw = uPlus[iTop+1];
#       uMup = uMinus[iTop];
#       uMdw = uMinus[iTop+1];

#       factA = (uPdw**2 + uMdw**2) / 2 - (uPup**2 + uMup**2) / 2;
#       factB = (uPdw**2 - uMdw**2) / 2 - (uPup**2 - uMup**2) / 2;

#       D5p = factA * np.exp(-kZero**2 / 2 * ((uPdw + uPup)**2) * chemRough2[iTop+1]);
#       D5m = factA * np.exp(-kZero**2 / 2 * ((uPdw - uPup)**2) * chemRough2[iTop+1]);

#       D7p = factA * np.exp(-kZero**2 / 2 * ((uMdw + uMup)**2) * chemRough2[iTop+1]);
#       D7m = factA * np.exp(-kZero**2 / 2 * ((uMdw - uMup)**2) * chemRough2[iTop+1]);

#       D6p = factB * np.exp(-kZero**2 / 2 * ((uPdw + uPup)**2) * magRough2[iTop+1]);
#       D6m = factB * np.exp(-kZero**2 / 2 * ((uPdw - uPup)**2) * magRough2[iTop+1]);

#       D8p = factB * np.exp(-kZero**2 / 2 * ((uMdw + uMup)**2) * magRough2[iTop+1]);
#       D8m = factB * np.exp(-kZero**2 / 2 * ((uMdw - uMup)**2) * magRough2[iTop+1]);
      
#       RzeroU= np.zeros((2,2),dtype=complex);
#       IV= np.zeros((2,2),dtype=complex);
#       IVp= np.zeros((2,2),dtype=complex);
#       Tzero1= np.zeros((2,2),dtype=complex);
#       Tzero2= np.zeros((2,2),dtype=complex);

#       RzeroU[0,0] = -(D5p + D6p) / (uPdw + uPup)**2;
#       RzeroU[1,1] = -(D7p - D8p) / (uMdw + uMup)**2;

#       IV[0,0] = (D5m + D6m) / (uPdw**2 - uPup**2);
#       IV[1,1] = (D7m - D8m) / (uMdw**2 - uMup**2);
      
#       IVp = IV; #switching upper and lower layers switches signs in D's and the u-sums

#       Tzero1[0,0] = 2 * uPup / (uPdw + uPup);
#       Tzero1[1,1] = 2 * uMup / (uMdw + uMup);
#       Tzero2[0,0] = 2 * uPdw / (uPdw + uPup);
#       Tzero2[1,1] = 2 * uMdw / (uMdw + uMup);

#       Mrt[:,:,iTop] = np.matmul(np.linalg.inv(IV), RzeroU);
#       Mtt[:,:,iTop] = np.matmul(np.linalg.inv(IVp) , Tzero1);

#       Mtr[:,:,iTop] = -Mrt[:,:,iTop];
#       Mrr[:,:,iTop] = np.matmul(np.linalg.inv(IV) , Tzero2);

#     # multiply all M-matrices by the inverse F's

#     Mtt[:,:,iTop] = np.matmul(Mtt[:,:,iTop] , invFmatrix[:,:,iTop]);
#     # Mtr = Mtr;
#     Mrt[:,:,iTop] = np.matmul(invFmatrix[:,:,iTop], np.matmul(Mrt[:,:,iTop], invFmatrix[:,:,iTop]));
#     Mrr[:,:,iTop] = np.matmul(invFmatrix[:,:,iTop], Mrr[:,:,iTop]);

#   # W-MATRICES
#   Wtt= np.zeros((2,2,nLayers),dtype=complex); 
#   Wtr= np.zeros((2,2,nLayers),dtype=complex);
#   Wrt= np.zeros((2,2,nLayers),dtype=complex);
#   Wrr = np.zeros((2,2,nLayers),dtype=complex);
  
#   Wtt[:,:,0] = np.identity(2);
#   Wrr[:,:,0] = np.identity(2);

#   y1= np.arange(1,nLayers);

#   for iBot in y1:
#     iTop = iBot - 1;
    
#     Wtt[:,:,iBot] = np.matmul(Mtt[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) - np.matmul(Wtr[:,:,iTop] , Mrt[:,:,iTop])), Wtt[:,:,iTop]));
#     Wtr[:,:,iBot] = Mtr [:,:,iTop] + np.matmul(Mtt[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) -np.matmul(Wtr[:,:,iTop],Mrt[:,:,iTop])),np.matmul(Wtr[:,:,iTop], Mrr[:,:,iTop])));
#     Wrr[:,:,iBot] = np.matmul(Wrr[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) - np.matmul(Mrt[:,:,iTop] , Wtr[:,:,iTop])), Mrr[:,:,iTop]));
#     Wrt[:,:,iBot] = Wrt[:,:,iTop] + np.matmul(Wrr[:,:,iTop], np.matmul(np.linalg.inv(np.identity(2) -np.matmul(Mrt[:,:,iTop],Wtr[:,:,iTop])),np.matmul(Mrt[:,:,iTop], Wtt[:,:,iTop])));
    
#   WrtEnd = Wrt[:,:,-1];

#   intPlus[iTheta] = sum(np.square(abs(np.matmul(WrtEnd,tZeroPlus))));
#   intMinus[iTheta] = sum(np.square(abs(np.matmul(WrtEnd,tZeroMinus))));
#   asymmRatio[iTheta]= (intPlus[iTheta]-intMinus[iTheta])/(intPlus[iTheta]+intMinus[iTheta]);
#   intdiff[iTheta]= (intPlus[iTheta]-intMinus[iTheta]);

# import matplotlib.pyplot as plt
# if (bob=='d'):
#    plt.plot(thetaRangeDeg, asymmRatio, 'k-',linewidth=2);
#    plt.savefig('./images/fig.png')
#    plt.show()
# else:
#    plt.semilogy(thetaRangeDeg,(intPlus+intMinus)/2,'k-',LineWidth=2);
#    plt.savefig('./images/fig.png')
#    plt.show()




def home(request):


   calci();
   import base64

   with open('images/fig.png', "rb") as image_file:
       image_data = base64.b64encode(image_file.read()).decode('utf-8')
   
   nnn["image"] = image_data
   
   

   return render(request, 'home.html', nnn)

def add(request):
    be=float(request.POST['beamEnergy'])
    v1=float(request.POST['array[1]'])
    v2=float(request.POST['array[2]'])
    v3=float(request.POST['array[3]'])
    v4=float(request.POST['array[4]'])
    v5=float(request.POST['array[5]'])
    v6=float(request.POST['array[6]'])
    v7=float(request.POST['array[7]'])
    v8=float(request.POST['array[8]'])
    v9=float(request.POST['array[9]'])
    v10=float(request.POST['array[10]'])
    v11=float(request.POST['array[11]'])
    v12=float(request.POST['array[12]'])
    v13=float(request.POST['array[13]'])
    v14=float(request.POST['array[14]'])
    b=request.POST['graphtype']

    global magicConstant1
    global beamEnergy
    global magicConstant2
    global bob
    beamEnergy=be
    magicConstant1=[v1,v2,v3,v4,v5,v6,v7]
    magicConstant1=[v8,v9,v10,v11,v12,v13,v14]
    bob=b
    
    calci();
    import base64

    with open('images/fig.png', "rb") as image_file:
       image_data = base64.b64encode(image_file.read()).decode('utf-8')
   
    nnn["image"] = image_data
    
    return render(request, 'result.html', nnn) 



    
