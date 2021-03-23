import numpy as np

end_speeds = []

for i in range(1,10):
    filename = "data/" + str(i) + ".csv"
    t, x, y, v = np.genfromtxt(filename, delimiter=",", skip_header=2, unpack=True)
    dist = np.sqrt((x[-1] - x[-2])**2 + (y[-1] - y[-2])**2)
    print("dist", dist)
    velocity = dist/(t[-1]-t[-2])
    print("velocity_" + str(i) + "=", velocity)
    print("t", t[-2], t[-1])
    print("x", x[-2], x[-1])
    print("y", y[-2], x[-1])
    end_speeds.append(velocity)

print(end_speeds)
#we can now use the collected end_speeds to calculate mean, standard diviation and standard error.

mean = sum(end_speeds)/len(end_speeds)

standard_diviation = np.sqrt((1/(len(end_speeds)-1)) * sum([(x-mean)**2 for x in end_speeds]))

standard_error = standard_diviation/np.sqrt(len(end_speeds))

print("Mean", mean)
print("Standard diviation", standard_diviation)
print("Standard error: ", standard_error)
