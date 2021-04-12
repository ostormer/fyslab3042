import numpy as np

# Samler 
end_speeds = []

# Iterer over hver av tracker-analysene, beregner sluttfart i hvert tilfelle.
# Tracker utelater siste frame, vi beregner sluttfart i denne:
for i in range(1,11):
    filename = "data/" + str(i) + ".csv"
    t, x, y, v = np.genfromtxt(filename, delimiter=",", skip_header=2, unpack=True)
    dist = np.sqrt((x[-1] - x[-2])**2 + (y[-1] - y[-2])**2)         # Avstand mellom siste frames  
    velocity = dist/(t[-1]-t[-2])                                   # Avstand delt på tidsendring: Sluttfart
    end_speeds.append(velocity)

# Beregner snitt, standardavvik og standardfeil fra sluttfartene i hvert eksempel
mean = sum(end_speeds)/len(end_speeds)

# Standardavvik fra formler i teori
standard_diviation = np.sqrt((1/(len(end_speeds)-1)) * sum([(x-mean)**2 for x in end_speeds]))

# Standardfeil fra formler i teori
standard_error = standard_diviation/np.sqrt(len(end_speeds))

# Printer resultatet
print("Mean", mean)
print("Standard diviation", standard_diviation)
print("Standard error: ", standard_error)
print(end_speeds)