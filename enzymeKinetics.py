import matplotlib.pyplot as plt
import math as math
import numpy as np
from scipy.optimize import curve_fit

class Probe(object):
    def __init__ (self, source,absorption = [], seconds = [], MIN = [], deltaE = []):
        self.source=source
        self.seconds=seconds
        self.MIN=MIN
        self.deltaE=deltaE
        self.absorption=absorption

#parse absorbtion file to ABSORPTION[] and SECONDS[]
    def parse(self):
        SECONDS=[]
        ABSORPTION=[]
        with open (self.source,'r') as verd:
        #read the whole file and split lines, or each character will be a line :(
            verd_data = verd.read().splitlines()
            i=0
            for line in verd_data:
                #Iteration 89 :  #DATA
                if (i > 89):
                    slash=line.find('\t')
                    #line looks like '55,000000\t0,03247
                    SECONDS.append(float(line[0:slash].replace(',','.')))
                    ABSORPTION.append(float(line[(slash + 1): (len(line))].replace(',','.')))
                i=i+1
       # print ('Time series :\n ', SECONDS, 'with length ', len(SECONDS))
       # print ('Absorption values:\n ', ABSORPTION, 'with length ', len(ABSORPTION))
        self.absorption= ABSORPTION
        self.seconds= SECONDS
   
    def set_dE(self): 
        self.deltaE=[self.absorption[0]]
        i=1
        for absorp in self.absorption[1:]:
            currentdE=absorp
            self.deltaE.append(currentdE)
         #   print ( 'dE with absorbtion ', absorp , ' at ', self.seconds[i], ' is equall ', self.deltaE[i])
            if (self.seconds[i]==60):
                break
            i=i+1
        
    def plot_dE_dT(self,probe):
        #60 Seconds
        plt.plot(self.MIN[0:13],self.deltaE,label=probe)

    def secToMin(self):
        for sec in self.seconds:
            self.MIN.append(sec/60)
    
#fit line plot and output optimized []        
def fitLine(x,y,name):
    #Polynomial coefficients, highest power first
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),label=name)
    #return a ( y=ax+b )
    return np.polyfit(x, y, 1)

def computeV0(dEdt):
    #405 nm
    ETA=14305
    #10mm=1cm
    d=1
    velocity=(dEdt)/(ETA*d)
    return velocity


def f(S,V_max,Km):
    return (V_max*S)/(Km+S)



print('Michailis-Menten Kinetics. Author Anastasia Grekova')
#wait
input("Press Enter to continue...")
sources = {'unverd':'/home/anastasia/Documents/unverd.asc',
        '2x':'/home/anastasia/Documents/2fach.asc',
        '4x':'/home/anastasia/Documents/4fach.asc',
        '10x':'/home/anastasia/Documents/10fach.asc',
        '20x':'/home/anastasia/Documents/20fach.asc',
        '40x':'/home/anastasia/Documents/40fach.asc',
        '100x':'/home/anastasia/Documents/100fach.asc'}
print('Source files are ', sources)
#wait
input("Press Enter to continue...")
#[S] der Messung
verdReihe=[20,10,5,2,1,0.5,0.02]
SUBSTRATE = [x/100 for x in verdReihe]

print ('[S] in mmol/L ', SUBSTRATE)
data=sources
for source in sources.keys():
    data[source]=Probe(sources[source])
#Absorption plot
for name,probe in data.items():
    probe.parse()
    probe.set_dE()
    probe.secToMin()
    probe.plot_dE_dT(name)

plt.xlabel('time in min.')
plt.ylabel('absorption')
plt.title("dE/dt")
plt.legend()
plt.show()

#Fit Lines to linear function and plot

timeTill1Min=[0.0, 0.08333333333333333, 0.16666666666666666, 0.25, 0.3333333333333333, 0.4166666666666667, 0.5, 0.5833333333333334, 0.6666666666666666, 0.75, 0.8333333333333334, 0.9166666666666666, 1.0]

dEdt={}

for name,probe in data.items():
    probe_lin=fitLine(timeTill1Min, probe.deltaE, name)
    dEdt[name]=probe_lin[0]
input("dE/dt as coeff from normalized plots. Press Enter to continue...")
print (dEdt)
plt.xlabel('time in min.')
plt.ylabel('absorption')
plt.title("normalized dE/dt")
plt.legend()    
plt.show()

#Substrat concentrarion / v0 Plot
v0={}
v0_quot={}
#compute Vo  ans 1/v0
for name,data in dEdt.items():
    data_velocity= computeV0(data)
    #in MUMOL/L*min
    v0[name]=data_velocity*1000000
    velocity_quot=1/(data_velocity*1000000)
    v0_quot[name]=velocity_quot
input("V0 in MUmol/l*min.Press Enter to continue...")
print (v0)
input("1/V0 in L*min/MUmol. Press Enter to continue...")
print (v0_quot)


#compute curve fitting for all [S]
v0_fit=v0.copy()
v0_fit_list=list(v0_fit.values())
#for function optimization it is important, that values are observed from lower to higher, so we have to reverse list to 100x,40x,20x,10x,4x,2x,unv

v0_fit_list.reverse()
v0_fit_list = [0] + v0_fit_list
#mmol/L 
S_fit =SUBSTRATE.copy()
S_fit.reverse()
S_fit = [0] + S_fit
#popt - best-fit values, pcov error matrix
popt, pcov = curve_fit(f,S_fit,v0_fit_list)
input('Curve fit parameters. Press ENTER to continue')
print('Best fit values are Km ', popt[1] ,' mmol/L and Vmax  ', popt[0], ' MUmol/min*L with an error', pcov)
y_fit=[]
for x in S_fit:
    y_fit.append(f(x,*popt))
plt.plot(S_fit,y_fit)
#plot observed V0 for different [S]
v0_names = list(v0.keys())
v0_data= list(v0.values())
for i,txt in enumerate( v0_names):
    plt.scatter(SUBSTRATE[i],v0_data[i])
    plt.annotate(txt, (SUBSTRATE[i],v0_data[i]))
plt.xlabel('concentration in mmol/L')
plt.ylabel('v0 in MUmol/L*min')
plt.title('Michaelis-Menten Diagram')
plt.show()

reversedSubstrate = [1/x for x in SUBSTRATE]

#plot Lineveawer-Burk 1/V0 L*min/MUmol 1/[S] L/mmol
reversedv0=v0_quot
reversedv0_list= list(reversedv0.values())
reversedv0_names= list(reversedv0.keys())

for i,txt in enumerate(reversedv0_names):
    plt.scatter(reversedSubstrate[i],reversedv0_list[i])
    plt.annotate(txt, (reversedSubstrate[i],reversedv0_list[i]))
#y=ax+b [a,b]
fit_Burk_line = fitLine(reversedSubstrate, reversedv0_list, 'regression')
a_Burk=fit_Burk_line[0]
b_Burk=fit_Burk_line[1]
kM_Burk=a_Burk/b_Burk
v0_max_Burk= 1/b_Burk
input('BURK DIAGRAMM - press ENTER')
print ('BURK DIAGRAMM, a= ', a_Burk, ' b= ', b_Burk)
print ('Km in mmol/L', kM_Burk, 'V max in MUmol*L/min', v0_max_Burk)
plt.xlabel('1/concentration in L/mmol')
plt.ylabel('1/v0 in L*min/MUmol')
plt.title('Lineweaver-Burk Diagram')
plt.show()

#Eadie-Hofstee v0 and v0/[s]
#mmol/L
substrateEadie = SUBSTRATE
v0S={}
i=0

for label,data in v0.items():
    v0S[label]=data/substrateEadie[i]
    plt.scatter(v0S[label],data)
    plt.annotate(label,(v0S[label],data))
    i=i+1
#y=ax+b [a,b])
input('Press ENTER to continue')
print('V0/[S]', v0S)
fit_Eadie_line = fitLine(list(v0S.values()), list(v0.values()), 'regression')
a_Eadie=fit_Eadie_line[0]
b_Eadie=fit_Eadie_line[1]
kM_Eadie=-a_Eadie
v0_max_Eadie=b_Eadie

input('EADIE DIAGRAMM - press ENTER')
print ('EADIE DIAGRAMM, a= ', a_Eadie, ' b= ', b_Eadie)
print ('Km in mmol/L', kM_Eadie, 'V max in MUmol*L/min', v0_max_Eadie)

plt.xlabel('v0/[S] in 1*10^-3/min')
plt.ylabel('v0 in MUmol/L*min')
plt.title('Eadie/Hofstee Diagram')
plt.show()

