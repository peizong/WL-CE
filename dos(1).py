#!/usr/bin/python
import numpy as np
from scipy.interpolate import splev,splrep,interp1d, spline
import matplotlib.pyplot as plt
import matplotlib

class ppdos:
  """post processing for one dos file """
  def __init__(self, gamma, dos_for_N_atoms, name_of_dos_file,skip_header_file):
    self.gamma=gamma
    self.N=dos_for_N_atoms
    self.filename=name_of_dos_file
    self.skip_header=skip_header_file
    self.data=[]
    self.sro_spline=[]
    self.sro_can=[]
    self.g_spline=[]
    self.U_can=[]
    self.T_mcan=[]
    self.specific_heat_can=[]
    self.specific_heat_mcan=[]
    self.binder_E_can=[]
    self.binder_sro_can=[]
  def read_data(self):
    dat= np.genfromtxt(self.filename,delimiter="	",skip_header=self.skip_header)
    #print dat
    dat[:,0]=dat[:,0]*self.N
    nonzero=(dat==0).sum(1)
    dat_nz=dat[nonzero==0,:]
    dat_nz[:,1]=dat_nz[:,1]-np.min(dat_nz[:,1])
    if dat_nz[0,0]>dat_nz[1,0]:
      dat_nz[:,0]=dat_nz[:,0][::-1]
      dat_nz[:,1]=dat_nz[:,1][::-1]
    self.data=dat_nz
    print self.data
  def cal_sro_spline(self,Ncol):
    self.read_data()
    x_new =np.linspace(np.min(self.data[:,0]), np.max(self.data[:,0]), num=0.5*len(self.data), endpoint=True)
    mid=splrep(self.data[:,0],self.data[:,Ncol])
    y_new=splev(x_new,mid)
    self.sro_spline= np.transpose(np.array([x_new,y_new]))
  def plot_sro(self):
    self.cal_sro_spline(2)
    sro=self.sro_spline
    plt.plot(1000*sro[:,0]/self.N,sro[:,1],linewidth=3.0)
    plt.xlabel('E/meV per atom')
    plt.ylabel('Short Range Order')
    plt.axis([1000*np.min(sro[:,0])/self.N,1000*np.max(sro[:,0])/self.N,np.min(sro[:,1]),np.max(sro[:,1])])
  def cal_sro_can(self,T,Ncol):
    self.cal_sro_spline(Ncol)
    sro=self.sro_spline
    self.cal_g_spline()
    g=self.g_spline
    k_B=1.38e-23/1.602e-19
    beta=1/(k_B*T)
    Z,avg_sro=0,0
    dE=g[1,0]-g[0,0]
#    print sro
#    print g
    for i in range(0,len(sro)):
      Z +=dE*np.exp(g[i,1]-beta*g[i,0]) # change i,0
      avg_sro +=sro[i,1]*dE*np.exp(g[i,1]-beta*g[i,0])
    avg_sro=avg_sro/Z
    self.sro_can= avg_sro
  def plot_sro_can(self,T):
    for i in range(0,3): #ternary
      self.cal_sro_can(T,2+i)
      plt.plot(T,self.sro_can,linewidth=3.0)
    plt.xlabel('T/K')
    plt.ylabel('Short Range Order')
    plt.axis([100,2000, -0.0,0.5]) #np.min(self.sro_can),np.max(self.sro_can)])
    plt.legend([str(self.N)+'_'+str(self.gamma)], loc='best')
  def cal_binder_sro_can(self,T):
    self.cal_sro_spline()
    sro=self.sro_spline
    self.cal_g_spline()
    g=self.g_spline
    #g[:,0]=g[:,0]*self.N
    k_B=1.38e-23/1.602e-19
    beta=1/(k_B*T)
    Z,avg_sro2,avg_sro4=0,0,0
    dE=g[1,0]-g[0,0]
    for i in range(0,len(sro)):
      Z +=dE*np.exp(g[i,1]-beta*g[i,0])
      avg_sro2 +=sro[i,1]**2*dE*np.exp(g[i,1]-beta*g[i,0])
      avg_sro4 +=sro[i,1]**4*dE*np.exp(g[i,1]-beta*g[i,0])
    avg_sro2=avg_sro2/Z
    avg_sro4=avg_sro4/Z
    self.binder_sro_can= 1-avg_sro4/(3*avg_sro2**2)
  def plot_binder_sro_can(self,T):
    self.cal_binder_sro_can(T)
    plt.plot(T,self.binder_sro_can,linewidth=3.0)
    plt.xlabel('T/K')
    plt.ylabel('Short Range Order')
    #plt.axis([600,1500, 0,1]) #np.min(self.sro_can),np.max(self.sro_can)])
    plt.legend([str(self.N)+'_'+str(self.gamma)], loc='best')
  def cal_g_spline(self):
    self.read_data()
    x_new =np.linspace(np.min(self.data[:,0]), np.max(self.data[:,0]), num=0.5*len(self.data), endpoint=True)
    mid=splrep(self.data[:,0],self.data[:,1]) 
    y_new=splev(x_new,mid)
    print mid, y_new
    self.g_spline= np.transpose(np.array([x_new,y_new]))                              
  def plot_g(self):
    self.cal_g_spline()
    g=self.g_spline
#    if self.skip_header==1: g[:,0]=g[:,0]+0.5
    plt.plot(1000*g[:,0]/self.N,g[:,1],linewidth=3.0)
    plt.xlabel('E/meV per atom')
    plt.ylabel('ln[g(E)]')
    #plt.axis([np.min(g[:,0])/self.N,np.max(g[:,0])/self.N,0,55])
    plt.axis([1000*np.min(g[:,0])/self.N,1000*np.max(g[:,0])/self.N,0,np.max(g[:,1])])
  def cal_U_can(self,T):
    self.cal_g_spline()
    g=self.g_spline
    k_B=1.38e-23/1.602e-19
    beta=1/(k_B*T)
    Z,avg_E,avg_E2=0,0,0
    dE=g[1,0]-g[0,0]
    for i in range(0,len(g)):
      Z +=dE*np.exp(g[i,1]-beta*g[i,0])
      avg_E +=g[i,0]*dE*np.exp(g[i,1]-beta*g[i,0])
    avg_E=avg_E/Z
    self.U_can= avg_E
  def plot_UT_can(self,T):
    self.cal_U_can(T)
    plt.plot(1000*self.U_can/self.N,T,linewidth=3.0)
    plt.xlabel('U/meV per atom')
    plt.ylabel('T/K')
    plt.axis([1000*np.min(self.data[:,0])/self.N, 1000*np.max(self.data[:,0])/self.N, 600, 1500])
    plt.legend([str(self.N)+'_'+str(self.gamma)], loc='best')
    #plt.show() 
  def cal_specific_heat_can(self,T):
    self.cal_g_spline()
    g=self.g_spline
    print g
    k_B=1.38e-23/1.602e-19
    beta=1/(k_B*T)
    Z,avg_E,avg_E2=0,0,0
    dE=g[1,0]-g[0,0]
    for i in range(0,len(g)):
      Z +=dE*np.exp(g[i,1]-beta*g[i,0])
      avg_E +=g[i,0]*dE*np.exp(g[i,1]-beta*g[i,0])
      avg_E2 += g[i,0]**2*dE*np.exp(g[i,1]-beta*g[i,0])
    avg_E=avg_E/Z
    avg_E2=avg_E2/Z
    self.specific_heat_can= (avg_E2-avg_E**2)/k_B**2/T**2
  def cal_binder_E_can(self,T):
    self.cal_g_spline()
    g=self.g_spline
    k_B=1.38e-23/1.602e-19
    beta=1/(k_B*T)
    Z,avg_E2,avg_E4=0,0,0
    dE=g[1,0]-g[0,0]
    for i in range(0,len(g)):
      Z +=dE*np.exp(g[i,1]-beta*g[i,0])
      avg_E2 += g[i,0]**2*dE*np.exp(g[i,1]-beta*g[i,0])
      avg_E4 += g[i,0]**4*dE*np.exp(g[i,1]-beta*g[i,0])
    avg_E2=avg_E2/Z
    avg_E4=avg_E4/Z
    self.binder_E_can= 1-avg_E4/(3*avg_E2**2)
  def plot_binder_E(self,T):   
    self.cal_binder_E_can(T)        
    plt.plot(T,self.binder_E_can,linewidth=3.0)   
    plt.xlabel('T/K')          
    plt.ylabel('V$_N$')          
    plt.legend([str(self.N)+'_'+str(self.gamma)], loc='best')
  def plot_specific_heat_can(self,T):   
    self.cal_specific_heat_can(T)   
    for i in range(0,len(T)):
      print T[i]-148,"	",self.specific_heat_can[i]
    #print T-218 #118
    plt.plot(T,self.specific_heat_can,linewidth=3.0)   
    plt.xlabel('T/K')          
    plt.ylabel('Cv/a.u.')          
    plt.legend([str(self.N)+'_'+str(self.gamma)], loc='best')
    #plt.axis([np.min(self.data[:,0]), np.max(self.data[:,0]), -0.0, 1500])
    plt.show() 
  def cal_T_mcan(self):
    self.cal_g_spline()
    g=self.g_spline
    k_B=1.38e-23/1.602e-19
    dE=g[1,0]-g[0,0]
    dS_dE=k_B*np.diff(g[:,1])/dE
    self.T_mcan=1.0/dS_dE
  def plot_UT_mcan(self):
    self.cal_T_mcan()
    #plt.plot(self.g_spline[:,0],self.g_spline[:,1])
    plt.plot(self.g_spline[0:len(self.g_spline)-1,0],self.T_mcan,linewidth=3.0)

#caseA3=ppdos(1E-5,128,'128/dos.dat-41eci',1)
#caseA3.cal_g_spline()
#plt.subplot(211)
#plt.plot(caseA3.g_spline[:-1,0],np.diff(caseA3.g_spline[:,1]))
#plt.subplot(212)
#plt.plot(caseA3.g_spline[:,0],caseA3.g_spline[:,1])
#plt.show()

def plot_test_gamma():
        font={'family' : 'normal',
#              'weight' : 'bold',
              'size'   : 22}
        matplotlib.rc('font',**font)
	#create cases:
	T=np.linspace(600,1500,1000)
        caseC3=ppdos(1E-3,64,'64/dos.dat-1E-3',7)
        caseC4=ppdos(1E-4,64,'64/dos.dat-1E-4',9)
	caseC5=ppdos(1E-5,64,'64/dos.dat-1E-5',16)
        caseC6=ppdos(1E-5,64,'64/dos.dat-1E-6',9)
        caseC7=ppdos(1E-5,64,'64/dos.dat-2E-7',9)
	#figure 1: density of states
	plt.subplot(221)
	caseC3.plot_g()
	caseC4.plot_g()
	caseC5.plot_g()
        caseC6.plot_g()
        caseC7.plot_g()
	plt.legend(['64_1E-3','64_1E-4','64_1E-5','64_1E-6','64_2E-7'],loc='best')
	plt.subplot(222)
	caseC3.plot_UT_can(T)
	caseC4.plot_UT_can(T)
	caseC5.plot_UT_can(T)
        caseC6.plot_UT_can(T)
        caseC7.plot_UT_can(T)
	plt.legend(['64_1E-3','64_1E-4','64_1E-5','64_1E-6','64_2E-7'],loc='best')
	#plt.text(-0.065,615,r'(b)')
	#plt.show()

	#figure 2: transition temperature
	plt.subplot(223)
	caseC3.plot_specific_heat_can(T)
	caseC4.plot_specific_heat_can(T)
	caseC5.plot_specific_heat_can(T)                             
        caseC6.plot_specific_heat_can(T) 
        caseC7.plot_specific_heat_can(T) 
	plt.legend(['64_1E-3','64_1E-4','64_1E-5','64_1E-6','64_2E-7'],loc='best')
	plt.subplot(224)
	TC3=T[np.where(caseC3.specific_heat_can == caseC3.specific_heat_can.max())]
	TC4=T[np.where(caseC4.specific_heat_can == caseC4.specific_heat_can.max())]
	TC5=T[np.where(caseC5.specific_heat_can == caseC5.specific_heat_can.max())]
        TC6=T[np.where(caseC6.specific_heat_can == caseC6.specific_heat_can.max())]
        TC7=T[np.where(caseC7.specific_heat_can == caseC7.specific_heat_can.max())]
        gamma=np.array([1E-3,1E-4,1E-5,1E-6,2E-7])
	TcC=np.array([TC3,TC4,TC5,TC6,TC7])
        print TcC
	plt.plot(gamma,TcC)
	plt.scatter(gamma,TcC,alpha=1,s=100)
	plt.xscale('log')
	plt.legend(['#atom=64'],loc='best')
        plt.xlabel('gamma')
        plt.ylabel('T$_c$/K')
	plt.show()

def plot_test_supercell():
        font={'family' : 'normal',
#              'weight' : 'bold',
              'size'   : 22}
        matplotlib.rc('font',**font)
        #create cases:
        T=np.linspace(600,1500,1000)
        caseA5=ppdos(1E-5,48,'48/dos.dat-1E-5',8)
        caseB5=ppdos(1E-5,54,'54/dos.dat-1E-5',8)
        caseC5=ppdos(1E-5,64,'64/dos.dat-1E-5',16)
        caseD5=ppdos(1E-5,128,'128/dos.dat-1E-5',1)
        #figure 1: density of states
        plt.subplot(221)
        caseA5.plot_g()
        caseB5.plot_g()
        caseC5.plot_g()
        caseD5.plot_g()
        plt.legend(['48','54','64','128'],loc='best')
        plt.subplot(222)
        caseA5.plot_UT_can(T)
        caseB5.plot_UT_can(T)
        caseC5.plot_UT_can(T)
        caseD5.plot_UT_can(T)
        plt.legend(['48','54','64','128'],loc='best')
        
        #figure 2: transition temperature
        plt.subplot(223)
        caseA5.plot_specific_heat_can(T)
        caseB5.plot_specific_heat_can(T)
        caseC5.plot_specific_heat_can(T)
        caseD5.plot_specific_heat_can(T)
        plt.legend(['48','54','64','128'],loc='best')
        plt.subplot(224)
        N_atoms5=np.array([48,54,64,128])
        TA5=T[np.where(caseA5.specific_heat_can == caseA5.specific_heat_can.max())]
        TB5=T[np.where(caseB5.specific_heat_can == caseB5.specific_heat_can.max())]
        TC5=T[np.where(caseC5.specific_heat_can == caseC5.specific_heat_can.max())]
        TD5=T[np.where(caseD5.specific_heat_can == caseD5.specific_heat_can.max())]
        Tc5=np.array([TA5,TB5,TC5,TD5])
        print Tc5
        plt.plot(N_atoms5,Tc5)
        plt.scatter(N_atoms5,Tc5,alpha=1, s=100, marker='x', color='g')
        plt.xlabel('# atoms')
        plt.legend(['1E-5'],loc='best')
        plt.ylabel('T$_c$/K')
        plt.show()
def plot_test_CV():
        font={'family' : 'normal',
#              'weight' : 'bold',
              'size'   : 22}
        matplotlib.rc('font',**font)
        #create cases:
        T=np.linspace(600,1500,1000)
        caseD5=ppdos(1E-5,128,'128/dos.dat-1E-5',1)
        caseD52=ppdos(1E-5,128,'128/dos.dat-15eci',1)
        caseD53=ppdos(1E-5,128,'128/dos.dat-41eci-6',1)
        #figure 1: density of states
        plt.subplot(221)
        #caseD5.read_data()
        #caseD5.data[:,0]=caseD5.data[:,0]-2
        caseD52.plot_g()
        caseD5.plot_g()
        #caseD5.data[:,0]=caseD5.data[:,0]-2
        caseD53.plot_g()
        plt.legend(['15Js,CV=7meV','16Js,CV=15meV','41Js,CV=5meV'],loc='best')
        plt.axis([-95,-60,0,55])
        plt.subplot(222)
        caseD5.plot_UT_can(T)
        caseD52.plot_UT_can(T)
        caseD53.plot_UT_can(T)
        plt.legend(['15Js,CV=7meV','16Js,CV=15meV','41Js,CV=5meV'],loc='best')
 
        #figure 2: transition temperature
        plt.subplot(223)
        caseD5.plot_specific_heat_can(T)
        caseD52.plot_specific_heat_can(T)
        caseD53.plot_specific_heat_can(T)
        plt.legend(['15Js,CV=7meV','16Js,CV=15meV','41Js,CV=5meV'],loc='best')
        plt.subplot(224)
        TD5=T[np.where(caseD5.specific_heat_can == caseD5.specific_heat_can.max())]
        TD52=T[np.where(caseD52.specific_heat_can == caseD52.specific_heat_can.max())]
        TD53=T[np.where(caseD53.specific_heat_can == caseD53.specific_heat_can.max())]
        Tc_CV=np.array([TD52,TD5,TD53])
        print Tc_CV
        #CV=np.array([15,7,5])
        eci_n=np.array([15,16,41])
        #plt.plot(CV,Tc_CV)
        plt.plot(eci_n,Tc_CV)
        plt.scatter(eci_n,Tc_CV,alpha=1, s=100, marker='o', color='g')
        #plt.xlabel('Cross-Validation Score/meV')          
        plt.xlabel('Number of ECIs') 
        plt.ylabel('Transition Temperature /K') 
        plt.legend(['1E-5,128'],loc='best')
        plt.axis([10,45,830,930])
        plt.show()
def plot_sro():
        font={'family' : 'normal',
#              'weight' : 'bold',
              'size'   : 22}
        matplotlib.rc('font',**font) 
        T=np.linspace(600,1500,1000)
        case16=ppdos('16eci',128,'128/sro.dat-16eci',1)
        case15=ppdos('15eci',128,'128/sro.dat-15eci',1)
        case41=ppdos('41eci',128,'128/sro.dat-41eci-4',1)
        plt.subplot(121)
        case16.plot_sro()
        case15.plot_sro()
        case41.plot_sro()
        plt.legend(['16Js,CV=15meV','15Js,CV=7meV','41Js,CV=5meV'],loc='best')
        plt.subplot(122)
        case16.plot_sro_can(T)
        case15.plot_sro_can(T)
        case41.plot_sro_can(T)
        plt.legend(['16Js,CV=15meV','15Js,CV=7meV','41Js,CV=5meV'],loc='best')
        #plt.subplot(133)
        #case16.plot_binder_sro_can(T)
        #case15.plot_binder_sro_can(T)
        #case41.plot_binder_sro_can(T)
        plt.show()
def plot_binder_cumulant():
        font={'family' : 'normal',
#              'weight' : 'bold',
              'size'   : 22}
        matplotlib.rc('font',**font)
        T=np.linspace(600,2000,1000)
        plt.subplot(211)
        caseA5=ppdos(1E-5,48,'48/dos.dat-1E-5',8)
        caseB5=ppdos(1E-5,54,'54/dos.dat-1E-5',8)
        caseC5=ppdos(1E-5,64,'64/dos.dat-1E-5',16)
        caseD5=ppdos(1E-5,128,'128/dos.dat-1E-5',1)
        caseA5.plot_binder_E(T)
        caseB5.plot_binder_E(T)
        caseC5.plot_binder_E(T)
        caseD5.plot_binder_E(T)
        plt.legend(['48','54','64','128'],loc='best')
        plt.subplot(212)
        T=np.linspace(600,2000,1000)
        sroA=ppdos(1E-5,48,'48/sro.dat',1)
        sroB=ppdos(1E-5,54,'54/sro.dat',1)
        sroC=ppdos(1E-5,64,'64/sro.dat',1)
        sroD=ppdos(1E-5,128,'128/sro.dat-16eci',1)
        sroA.plot_binder_sro_can(T)
        sroB.plot_binder_sro_can(T)
        sroC.plot_binder_sro_can(T)
        sroD.plot_binder_sro_can(T)
        plt.legend(['48','54','64','128'],loc='best')
        plt.show()
def plot_sro3():
        font={'family' : 'normal',
#              'weight' : 'bold',
              'size'   : 22}
        matplotlib.rc('font',**font) 
        T=np.linspace(100,2000,1000)
        case=ppdos('41eci',108,'sro.dat',1)
        plt.subplot(121)
        case.plot_sro()
        plt.legend(['NiCoCr'],loc='best')
        plt.subplot(122)
        case.plot_sro_can(T)
        plt.legend(['Ni','Co','Cr'],loc='best')
        #plt.subplot(133)
        #case16.plot_binder_sro_can(T)
        #case15.plot_binder_sro_can(T)
        #case41.plot_binder_sro_can(T)
        plt.show()
#plot_test_gamma()
plt.subplot(221)
plot_test_supercell()
plt.subplot(222)
plot_test_CV()
plt.subplot(223)
plot_sro3()
plt.subplot(224)
plot_binder_cumulant()
#-----------ternary-------------------
##caseNiCoCr=ppdos(1E-3,240,'dos345/3-body/cdos.dat',1)
#caseNiCoCr=ppdos(1E-3,240,'dos345/15eci-25mev/cdos.dat',1)
##caseNiCoCr.read_data()
#T=np.linspace(200,2000,1000)
#caseNiCoCr.plot_specific_heat_can(T)
##caseNiCoCr.plot_UT_mcan()
#------------end----------------------
plt.show()


