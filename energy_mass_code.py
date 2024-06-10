#ME 104 Final Project Gear Selector
import numpy as np
import matplotlib.pyplot as plt
#%%

pitches = [20,32,48,64]
fos = 1.2
#plastic, brass,steel
frictions = [.2, .15, .1]
p_ses = [12*10**(-6), 14*10**-6, 8*10**(-6)]
densities = [.0513, .307, .284]

min_face_w = 1/8

l = 6

pla_strength = 60*10**6

sizes = [5/8/2,3/4/2,1/2,1.5/2,2/2,2.5/2]

valid_combinations = []

# Iterate through all possible combinations
for r1 in sizes:
    for r2 in sizes:
        for r3 in sizes:
            for r4 in sizes:
                for r5 in sizes:
                    for r6 in sizes:
                        for r7 in sizes:
                            result = (l * r2 * r4 * r6) / (r1 * r3 * r5 * r7)
                            if result < 0.0985 and result > 0.08:
                                valid_combinations.append((l, r1, r2, r3, r4, r5, r6, r7, result))
      
#%%                    
shaft_mu = .001
shaft_r = .25
pressure_angle = 0.349066
      

def calc_effs(valid_combinations, friction,p):
    effs = []
    for i in valid_combinations:
        n1 = (1 - friction*(pressure_angle + np.pi/(2*i[2]*p)
                            )) / (1 - friction*(pressure_angle - np.pi/(2*i[1]*p)))
        n2 = (1 - friction*(pressure_angle + np.pi/(2*i[4]*p)
                            )) / (1 - friction*(pressure_angle - np.pi/(2*i[3]*p)))
        n3 = (1 - friction*(pressure_angle + np.pi/(2*i[6]*p)
                            )) / (1 - friction*(pressure_angle - np.pi/(2*i[5]*p)))
        n4 = (1 - (shaft_r / i[1]) * shaft_mu)
        n5 = (1 - (shaft_r / i[2]) * shaft_mu)
        n6 = (1 - (shaft_r / i[3]) * shaft_mu)
        n7 = (1 - (shaft_r / i[4]) * shaft_mu)
        n8 = (1 - (shaft_r / i[5]) * shaft_mu)
        n9 = (1 - (shaft_r / i[6]) * shaft_mu)
        n10 = (1 - (shaft_r / i[7]) * shaft_mu)
        eff = n1*n2*n3*n4*n5*n6*n7*n8*n9*n10
        effs.append(eff)
    return effs

steel_effs24 = calc_effs(valid_combinations,frictions[2],pitches[0]) 
steel_effs32 = calc_effs(valid_combinations,frictions[2],pitches[1]) 
steel_effs48 = calc_effs(valid_combinations,frictions[2],pitches[2]) 
steel_effs64 = calc_effs(valid_combinations,frictions[2],pitches[3])

#%%
def gearMass(T,N,p_se):
    m = 4*T*np.pi * fos* p_se*N**(7/8)
    if m < 0:
        print(m)
    return m

def gearWidth(T,N,se):
    w = 16*T *fos * pitches[0]**2 / (se * N*(N-11)**1/8)
    if w < 0:
        print(w)
    return w

in_to_m = 1/39.37
Nm_to_lbin = 8.851

def calculateMasses(valid_combinations,effs,mat,p):
    masses = []
    j = 0
    for i in valid_combinations:
        m1 = gearMass((10*9.8* l/ effs[j]) * in_to_m * Nm_to_lbin, 2*i[1]*p,p_ses[mat])
        m2 = gearMass(10*9.8*l*i[2]/ (effs[j]*i[1]) * in_to_m * Nm_to_lbin, 2*i[2]*p,p_ses[mat])
        m3 = gearMass(10*9.8*l*i[2]/ (effs[j]*i[1]) * in_to_m * Nm_to_lbin, 2*i[3]*p,p_ses[mat])
        m4 = gearMass(10*9.8*l*i[2]*i[4]/ (effs[j]*i[1]*i[3]) * in_to_m * Nm_to_lbin, 2*i[4]*p,p_ses[mat])
        m5 = gearMass(10*9.8*l*i[2]*i[4]/ (effs[j]*i[1]*i[3]) * in_to_m * Nm_to_lbin, 2*i[5]*p,p_ses[mat])
        m6 = gearMass(10*9.8*l*i[2]*i[4] * i[6]/ (effs[j]*i[1]*i[3]* i[5]) * in_to_m * Nm_to_lbin, 2*i[5]*p,p_ses[mat])
        m7 = gearMass(10*9.8*l*i[2]*i[4] * i[6]/ (effs[j]*i[1]*i[3]* i[5]) * in_to_m * Nm_to_lbin, 2*i[5]*p,p_ses[mat])
        mm1 = np.pi * i[1]**2 * densities[mat] * min_face_w
        mm2 = np.pi * i[2]**2 * densities[mat] * min_face_w
        mm3 = np.pi * i[3]**2 * densities[mat] * min_face_w
        mm4 = np.pi * i[4]**2 * densities[mat] * min_face_w
        mm5 = np.pi * i[5]**2 * densities[mat] * min_face_w
        mm6 = np.pi * i[6]**2 * densities[mat] * min_face_w
        mm7 = np.pi * i[7]**2 * densities[mat] * min_face_w
        mass = max(m1,mm1) +  max(m2,mm2) +  max(m3,mm3) +  max(m4,mm4) + max(m5,mm5)+ max(m6,mm6) + max(m7,mm7)
        masses.append(mass)
        j+= 1
    return masses
    

steel_masses24 = calculateMasses(valid_combinations,steel_effs24, 2,pitches[0])
steel_masses32 = calculateMasses(valid_combinations,steel_effs32, 2,pitches[1])
steel_masses48 = calculateMasses(valid_combinations,steel_effs48, 2,pitches[2])
steel_masses64 = calculateMasses(valid_combinations,steel_effs64, 2,pitches[3])

#%%


plt.scatter(steel_masses24, steel_effs24, label='Pitch 24', alpha=0.7)
plt.scatter(steel_masses32, steel_effs32, label='Pitch 32', alpha=0.7)
plt.scatter(steel_masses48, steel_effs48, label='Pitch 48', alpha=0.7)
plt.scatter(steel_masses64, steel_effs64, label='Pitch 64', alpha=0.7)

plt.title('Efficiency vs Mass')
plt.xlabel('Mass')
plt.ylabel('Efficiency')
plt.legend()
plt.grid(True)
plt.show()
#%%
# Step 1: Calculate the average mass and efficiency for steel
norm_mass = np.linalg.norm(steel_masses24)
norm_eff = np.linalg.norm(steel_effs24)

# Step 2: Normalize the masses and efficiencies
normalized_masses = [mass / norm_mass for mass in steel_masses24]
normalized_effs = [eff / norm_eff for eff in steel_effs24]

# Step 3: Calculate the combined score
scores = []
for i in range(len(normalized_masses)):
    score = 0.55 * normalized_effs[i] - 0.45 * normalized_masses[i]
    scores.append(score)

# Step 4: Find the best combination
best_index = scores.index(max(scores))
best_combination = valid_combinations[best_index]


# Step 5: Print the best combination
print("Best Steel Combination: ", best_combination)
print("Best combination efficency", steel_effs24[best_index])
print("Best combination mass", steel_masses24[best_index])

Se_steel = 100000
PLA_s = 7250

def calcTorques(i,eff):
    t1 = (10*9.8* l/ eff) * in_to_m * Nm_to_lbin
    t2 = 10*9.8*l*i[2]/ (eff*i[1]) * in_to_m * Nm_to_lbin
    t3 = 10*9.8*l*i[2]/ (eff*i[1]) * in_to_m * Nm_to_lbin
    t4 = 10*9.8*l*i[2]*i[4]/ (eff*i[1]*i[3]) * in_to_m * Nm_to_lbin
    t5 = 10*9.8*l*i[2]*i[4]/ (eff*i[1]*i[3]) * in_to_m * Nm_to_lbin
    t6 = 10*9.8*l*i[2]*i[4] * i[6]/ (eff*i[1]*i[3]* i[5]) * in_to_m * Nm_to_lbin
    t7 = 10*9.8*l*i[2]*i[4] * i[6]/ (eff*i[1]*i[3]* i[5]) * in_to_m * Nm_to_lbin
    return t1,t2,t3,t4,t5,t6,t7

def calculateWidths(i, eff,p):
    w1 = gearWidth((10*9.8* l/ eff) * in_to_m * Nm_to_lbin, 2*i[1]*p,Se_steel)
    w2 = gearWidth(10*9.8*l*i[2]/ (eff*i[1]) * in_to_m * Nm_to_lbin,2*i[2]*p, Se_steel)
    w3 = gearWidth(10*9.8*l*i[2]/ (eff*i[1]) * in_to_m * Nm_to_lbin, 2*i[3]*p,Se_steel)
    w4 = gearWidth(10*9.8*l*i[2]*i[4]/ (eff*i[1]*i[3]) * in_to_m * Nm_to_lbin, 2*i[4]*p,Se_steel)
    w5 = gearWidth(10*9.8*l*i[2]*i[4]/ (eff*i[1]*i[3]) * in_to_m * Nm_to_lbin, 2*i[5]*p,Se_steel)
    w6 = gearWidth(10*9.8*l*i[2]*i[4] * i[6]/ (eff*i[1]*i[3]* i[5]) * in_to_m * Nm_to_lbin, 2*i[6]*p,Se_steel)
    w7 = gearWidth(10*9.8*l*i[2]*i[4] * i[6]/ (eff*i[1]*i[3]* i[5]) * in_to_m * Nm_to_lbin, 2*i[7]*p,Se_steel)
    return w1,w2,w3,w4,w5,w6,w7

def calcGearForces(t,r):
    f1 = t[0]/1.25
    f2 = f1
    f3 = t[2]/r[3]
    f4 = f3
    f5 = t[4]/r[5]
    f6 = f5
    f7= 1*9.8 / 4.448
    return f1,f2,f3,f4,f5,f6,f7

def calcShaftForces(f):
    s2 = f[1] - f[2]
    s3 = f[3] - f[4]
    s4 = f[5] - f[6]
    return s2,s3,s4

widths = calculateWidths(valid_combinations[best_index], steel_effs24[best_index], 24)

print(f"Calculated gear widths for the best combination: {widths}")

L = 2.5

def calcTorsion(torques):
    t1 = torques[2]*2/(np.pi*.125**3)
    t2 = torques[4]*2/(np.pi*.125**3)
    t3 = torques[6]*2/(np.pi*.125**3)
    return t1,t2,t3

def calcBending(loads):
    b1 = loads[0]*L/(np.pi*0.125**3)
    b2 = loads[1]*L/(np.pi*0.125**3)
    b3 = loads[2]*L/(np.pi*0.125**3)
    return b1,b2,b3

#%%
torques = calcTorques(best_combination,steel_effs24[best_index])
forces = calcGearForces(torques, best_combination)
shafts = calcShaftForces(forces)
shaftT = calcTorsion(torques)
shaftB = calcBending(shafts)

print("Torques: ", torques)
print("Gear Forces: ", forces)
print("Shaft Forces: ", shafts)
print("Shaft torsion: ", shaftT)
print("Shaft bending: ", shaftB)



