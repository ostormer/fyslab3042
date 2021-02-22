# TFY41xx Fysikk vaaren 2021.
#
# Programmet tar utgangspunkt i hoeyden til de 8 festepunktene.
# Deretter beregnes baneformen y(x) ved hjelp av 7 tredjegradspolynomer, 
# et for hvert intervall mellom to festepunkter, slik at baade banen y, 
# dens stigningstall y' = dy/dx og dens andrederiverte
# y'' = d2y/dx2 er kontinuerlige i de 6 indre festepunktene.
# I tillegg velges null krumning (andrederivert) 
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# De ulike banene er satt opp med tanke paa at kula skal 
# (1) fullfoere hele banen selv om den taper noe mekanisk energi underveis;
# (2) rulle rent, uten aa gli ("slure").

#Importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

#Horisontal avstand mellom festepunktene er 0.200 m
h = 0.200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
#Skriv inn y-verdiene til banens 8 festepunkter i tabellen yfast.
#Et vilkaarlig eksempel:
yfast = np.asarray([0.350,0.315,0.299,0.306,0.281,0.236,0.168,0.200])
#Erstatt med egne tallverdier avlest i tracker.
#Programmet beregner de 7 tredjegradspolynomene, et
#for hvert intervall mellom to festepunkter,
#med funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')
#Funksjonen cs kan naa brukes til aa regne ut y(x), y'(x) og y''(x)
#for en vilkaarlig horisontal posisjon x, eventuelt for mange horisontale
#posisjoner lagret i en tabell:
#cs(x)   tilsvarer y(x)
#cs(x,1) tilsvarer y'(x)
#cs(x,2) tilsvarer y''(x)
#Her lager vi en tabell med x-verdier mellom 0 og 1.4 m
xmin = 0.000
xmax = 1.401
dx = 0.001
x = np.arange(xmin, xmax, dx)   
#Funksjonen arange returnerer verdier paa det "halvaapne" intervallet
#[xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
#dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400, 
#dvs x blir en tabell med 1401 elementer
Nx = len(x)
y = cs(x)       #y=tabell med 1401 verdier for y(x)
dy = cs(x,1)    #dy=tabell med 1401 verdier for y'(x)
d2y = cs(x,2)   #d2y=tabell med 1401 verdier for y''(x)

#Plotteeksempel: Banens form y(x)
baneform = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,y,xfast,yfast,'*')
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$y(x)$ (m)',fontsize=20)
plt.ylim(0.10,0.40)
plt.grid()
plt.show()
#Figurer kan lagres i det formatet du foretrekker:
#baneform.savefig("baneform.pdf", bbox_inches='tight')
#baneform.savefig("baneform.png", bbox_inches='tight')
#baneform.savefig("baneform.eps", bbox_inches='tight')