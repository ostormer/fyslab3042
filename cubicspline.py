#Importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
import numpy as np
from scipy.interpolate import CubicSpline

# Banens startpunkter (t=0):
t, x, y = np.genfromtxt("data/bane.csv", delimiter=",", skip_header=2, unpack=True)

# Festepunktene til banen - 8 stk:
xFast = x
yFast = y

# Beregner tredjegradspolynomene av banen.
# Gir y(x), y'(x) og y''(x) 
# av cs(x), cs(x,1) og cs(x,2)
cs = CubicSpline(xFast, yFast, bc_type='natural')

# Tabell med x-verdier mellom 0 og 1.4 m 
xMin = 0.000
xMax = 1.401
dx = 0.001
x = np.arange(xMin, xMax, dx)

Nx = len(x)     # Lengden av x
y = cs(x)       # y(x) med 1401 verdier
dy = cs(x,1)    # dy(x) med 1401 verdier 
d2y = cs(x,2)   # d2y(x) med 1401 verdier

# TIME:
time_steps = 10000

# Konstanter

# Gravitasjonskonstanten:
g = 9.81                      
# Konstant knyttet til Treghetsmoment:
c = 2/5
# Kulas masse i kg:
M = 0.031
# X_0, første X-posisjon til kula:
x0 = x[0]
# Y_0, første Y-posisjon til kula:
y0 = cs(x0)
# Startfart: > 0 (Euler):
v0 = 0.001
# Starttid:                   
t0 = 0    
# Tidssteg:                      
dt = 1.15/time_steps                   
# Array med x-posisjoner: 
x_t = [x0] * (time_steps+1)
# Array med hastigheter:               
v_t = [v0] * (time_steps+1)
# Array med medgått tid:
t_list = np.linspace(t0, dt*time_steps, (time_steps+1))

# Beregner akselerasjon med hensyn paa x
def a(x):
    return -(5*g/7)*np.sin(np.arctan(cs(x,1)))


# Eulers metode for å beregne X, Y og v
# av t, samt tidsstegene
last_time_step = 0
for n in range(time_steps):
    x_t[n + 1] = x_t[n] + v_t[n] * dt
    current_a = a(x_t[n])
    v_t[n + 1] = v_t[n] + current_a * dt
    if x_t[n+1] >= 1.4:
        last_time_step = n+1
        break
if last_time_step==0:
    raise("Error")

x_t = x_t[:last_time_step+1]
v_t = v_t[:last_time_step+1]
t_list = t_list[:last_time_step+1]

# Beregner N og F med hensyn på X
v_x = np.sqrt((2*g*(y0-y))/(1+c))
k = d2y/(1 + dy**2)**(3/2)
a_normal = v_x**2*k
N = M*(g*np.cos(np.arctan(dy)) + a_normal)
F = (2*M*g*np.sin(np.arctan(dy)))/7

# Metode for sammenlikning
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx 

# Plotter eksperimentelt resultat (x(t)):
def printXofTCompare():
    t, x, y, v = np.genfromtxt("data/10.csv", delimiter=",", skip_header=2, unpack=True)
    compare_y = []
    t_list2 = np.array(t_list)
    x_t2 = np.array(x_t)
    for t_value in t:
        compare_y.append(abs(x[np.where(t == t_value)]-x_t2[find_nearest(t_list2, t_value)]))
    compare_y = [X[0] for X in compare_y]

    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(t, x, color="green", label="Eksperimentell posisjon")
    plt.plot(t_list, x_t, color="blue", label="Numerisk posisjon")
    plt.plot(t, compare_y, color="red", label="Differanse")
    plt.fill_between(t, compare_y, 0, color="red")
    plt.legend()
    plt.xlabel("Tid t - (s)", fontsize = 28)
    plt.ylabel("Posisjon x - (m)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/XofTComp")
    #plt.show()

# Plotter X av t
def printXofT():
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(t_list, x_t)
    plt.xlabel("Tid t - (s)", fontsize = 28)
    plt.ylabel("Posisjon x - (m)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/XofT")
    # plt.show()

# Plotter V av t
def printVofT():
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(t_list, v_t)
    plt.xlabel("Tid t - (s)", fontsize = 28)
    plt.ylabel("Hastighet v - (m/s)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/VofT")
    # plt.show()

# Plotter V av x
def printVofX():
    t, x, y, v = np.genfromtxt("data/10.csv", delimiter=",", skip_header=2, unpack=True)
    v.append(1.36598)
    x.append(1.4)
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(x_t, v_t, label="Numerisk hastighet")
    plt.scatter(x, v, marker="o", color="red", label="Eksperimentell hastighet")
    plt.legend(fontsize=22)
    plt.xlabel("Posisjon x - (m)", fontsize = 28)
    plt.ylabel("Hastighet v - (m/s)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/VofX")
    plt.show()

# Plotter F av X
def printFofX():
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(x, F)
    plt.xlabel("Posisjon x - (m)", fontsize = 28)
    plt.ylabel("Friksjonskraft f - (N)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/FofX")
    # plt.show()

# Plotter N av X
def printNofX():
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(x, N)
    plt.xlabel("Posisjon x - (m)", fontsize = 28)
    plt.ylabel("Normalkraft N - (N)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/NofX")
    # plt.show()

# Plotter Y av X - Banens form
def printYofX():
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(x, y)
    plt.xlabel("Posisjon x - (m)", fontsize = 28)
    plt.ylabel("Posisjon y - (m)", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/YofX")
    # plt.show()

# Plotter F/N av X
def printFNofX():
    plt.figure(figsize=(12,7), facecolor="w", edgecolor="w")
    plt.plot(x, abs(F/N))
    plt.xlabel("Posisjon x -(m)", fontsize = 28)
    plt.ylabel("[f/N]", fontsize = 28)
    plt.tick_params(axis='x', labelsize=26)
    plt.tick_params(axis='y', labelsize=26)
    plt.subplots_adjust(bottom=0.2)
    plt.grid()
    plt.savefig("figurer/FNofX")
    # plt.show()


if __name__ == "__main__":
    printVofX()
    printXofT()
    printVofT()
    printFofX()
    printNofX()
    printYofX()
    printFNofX()
    printXofTCompare()
    
    print("Sluttfart: ", v_t[-1])
    print("Maksfart: ", max(v_t))
    print("Tid ved høyeste hastighet: ", t_list[v_t.index(max(v_t))])
    print("Medgått tid: ", t_list[-1])
    
    