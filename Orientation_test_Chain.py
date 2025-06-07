### Purpose: This Script will compute Hermans orientation factor (between e2e vector and direction of pull) for systems undergoing deformation ####
### Syntax: python orientation.py (filename).dump ======== orientation.csv ###
### Author: Taofeek Tejuosho###
### Date: July 2024 ###

import math
import numpy as np
import csv

import time
#Creating a function that calculates the average of the stress autocorrelation from the 3 components
start_time = time.time()
def readTimesteps(inputFile, timestepsCount=None):
    counter = 0  # Initialize counter here
    if not timestepsCount:
        timestepsCount = float('inf')  # Set to a very large number if not provided
    
    timesteps = {}
    with open(inputFile) as f:
        line = f.readline()
        while line != '':
            if line.split()[0] == "ITEM:":
                timestep = int(f.readline().split()[0])
                timesteps[timestep] = {}
                filler = f.readline()
                atoms = int(f.readline().split()[0])
                timesteps[timestep]["atoms"] = atoms
                filler = f.readline()
                x = f.readline()
                x_pos = float(x.split()[1])
                x_neg = float(x.split()[0])
                timesteps[timestep]["-x"] = x_neg
                timesteps[timestep]["x"] = x_pos
                y = f.readline()
                y_pos = float(y.split()[1])
                y_neg = float(y.split()[0])
                timesteps[timestep]["-y"] = y_neg
                timesteps[timestep]["y"] = y_pos
                z = f.readline()
                z_pos = float(z.split()[1])
                z_neg = float(z.split()[0])
                timesteps[timestep]["-z"] = z_neg
                timesteps[timestep]["z"] = z_pos
                headers = f.readline()
                timesteps[timestep]["headers"] = headers
                timesteps[timestep]["mols"] = {}
                for i in range(atoms):
                    # line = f.readline()
                    # id = int(line.split()[0])
                    # mol = int(line.split()[1])
                    # type = int(line.split()[2])
                    # q = int(line.split()[3])
                    # xs = float(line.split()[4]) * (timesteps[timestep]['x'] - timesteps[timestep]['-x'])
                    # ys = float(line.split()[5]) * (timesteps[timestep]['y'] - timesteps[timestep]['-y'])
                    # zs = float(line.split()[6]) * (timesteps[timestep]['z'] - timesteps[timestep]['-z'])
                    # ix = int(line.split()[7])
                    # iy = int(line.split()[8])
                    # iz = int(line.split()[9])
                    line = f.readline()
                    id = int(line.split()[0])
                    mol = int(line.split()[1])
                    type = int(line.split()[2])
                    xs = float(line.split()[3]) * (timesteps[timestep]['x'] - timesteps[timestep]['-x'])
                    ys = float(line.split()[4]) * (timesteps[timestep]['y'] - timesteps[timestep]['-y'])
                    zs = float(line.split()[5]) * (timesteps[timestep]['z'] - timesteps[timestep]['-z'])
                    ix = int(line.split()[6])
                    iy = int(line.split()[7])
                    iz = int(line.split()[8])
                    xbox = (timesteps[timestep]['x'] - timesteps[timestep]['-x'])
                    ybox = (timesteps[timestep]['y'] - timesteps[timestep]['-y'])
                    zbox = (timesteps[timestep]['z'] - timesteps[timestep]['-z'])
                    if xs > xbox:
                        xs = xs - xbox
                    elif xs < 0:
                        xs = xs + xbox
                    if ys > ybox:
                        ys = ys - ybox
                    elif ys < 0:
                        ys = ys + ybox
                    if zs > zbox:
                        zs = zs - zbox
                    elif zs < 0:
                        zs = zs + zbox
                    # xs = xbox * float(line.split()[4]) + ix * xbox
                    # ys = ybox * float(line.split()[5]) + iy * ybox
                    # zs = zbox * float(line.split()[6]) + iz * zbox
                    xs = xbox * float(line.split()[3]) + ix * xbox
                    ys = ybox * float(line.split()[4]) + iy * ybox
                    zs = zbox * float(line.split()[5]) + iz * zbox
                    if mol in timesteps[timestep]["mols"]:
                        timesteps[timestep]["mols"][mol]["atoms"].append((id, type, xs, ys, zs, ix, iy, iz))
                    else:
                        timesteps[timestep]["mols"][mol] = {}
                        timesteps[timestep]["mols"][mol]["atoms"] = []
                        timesteps[timestep]["mols"][mol]["atoms"].append((id, type, xs, ys, zs, ix, iy, iz))

                for mol in timesteps[timestep]["mols"]:
                    x = 0
                    y = 0
                    z = 0
                    length = len(timesteps[timestep]["mols"][mol]["atoms"])
                    timesteps[timestep]["mols"][mol]["length"] = length
                    for i in range(length):
                        x = x + timesteps[timestep]["mols"][mol]["atoms"][i][2] / length
                        y = y + timesteps[timestep]["mols"][mol]["atoms"][i][3] / length
                        z = z + timesteps[timestep]["mols"][mol]["atoms"][i][4] / length
                    timesteps[timestep]["mols"][mol]["xcm"] = x
                    timesteps[timestep]["mols"][mol]["ycm"] = y
                    timesteps[timestep]["mols"][mol]["zcm"] = z
            line = f.readline()
            counter = counter + 1
            if counter == timestepsCount:
                return timesteps
    return timesteps

def Rg2Calculator(timesteps):
    for timestep in timesteps:
        for mol in timesteps[timestep]["mols"]:
            length = timesteps[timestep]["mols"][mol]["length"]
            rg2 = 0
            for i in range(length):
                x = timesteps[timestep]["mols"][mol]["atoms"][i][2]
                y = timesteps[timestep]["mols"][mol]["atoms"][i][3]
                z = timesteps[timestep]["mols"][mol]["atoms"][i][4]
                x = x - timesteps[timestep]["mols"][mol]["xcm"]
                y = y - timesteps[timestep]["mols"][mol]["ycm"]
                z = z - timesteps[timestep]["mols"][mol]["zcm"]
                x2 = x ** 2
                y2 = y ** 2
                z2 = z ** 2
                rg2 = ((x2 + y2 + z2) / length) + rg2
            timesteps[timestep]["mols"][mol]["rg2"] = rg2
    return timesteps

def writeRg2ToFile(timesteps, outFile):
    with open(outFile, 'w') as f:
        f.write("rg2, molNumber, timestep\n")
        for timestep in timesteps:
            for mol in timesteps[timestep]["mols"]:
                rg2 = str(timesteps[timestep]["mols"][mol]["rg2"])
                molNumber = str(mol)
                timestepWrite = str(timestep)
                f.write(rg2 + ", " + molNumber + ", " + timestepWrite + "\n")

def writeRg2AverageToFile(timesteps, outFile, frames=1):
    with open(outFile, 'w') as f:
        f.write("rg2\n")
        mols = {}
        frameCount = 0
        for timestep in timesteps:
            if frameCount < frames:
                for mol in timesteps[timestep]["mols"]:
                    if mol not in mols:
                        mols[mol] = timesteps[timestep]["mols"][mol]["rg2"]
                    else:
                        mols[mol] = mols[mol] + timesteps[timestep]["mols"][mol]["rg2"]
                frameCount = frameCount + 1
            else:
                continue
        sum = 0
        for mol in mols:
            sum = sum + (mols[mol] / frames)
        f.write('# ' + str(frames) + '\n')
        f.write('polymer number, rg squared \n')
        for mol in mols:
            f.write(str(mol) + ', ' + str(mols[mol] / frames) + '\n')
        sum = sum / len(mols)
        rg_square_root = math.sqrt(sum)
        f.write('------------------AVERAGE Square radius of gyration----------------\n')
        f.write(str(sum) + '\n')
        f.write('--------------Rg------------\n')
        f.write(str(rg_square_root))


#x - direction
# def calculate_orientation_order_parameter(timesteps, min_chain_length, max_chain_length):
#     coslist = []
#     strain = []
#     initial_xbox = timesteps[list(timesteps.keys())[0]]['x'] - timesteps[list(timesteps.keys())[0]]['-x']
#     for timestep in timesteps:
#         xbox = timesteps[timestep]['x'] - timesteps[timestep]['-x']
#         strain_value = (xbox - initial_xbox) / initial_xbox
#         strain.append(strain_value)
#         for mol in timesteps[timestep]['mols']:
#             chain_length = timesteps[timestep]['mols'][mol]['length']
#             if min_chain_length <= chain_length <= max_chain_length:
#                 atoms = timesteps[timestep]['mols'][mol]['atoms']
#                 x1, y1, z1 = atoms[0][2], atoms[0][3], atoms[0][4]
#                 x2, y2, z2 = atoms[len(atoms) - 1][2], atoms[len(atoms) - 1][3], atoms[len(atoms) - 1][4]
#                 #print(x2,y2,z2)
#                 ABx = x2 - x1
#                 ABy = y2 - y1
#                 ABz = z2 - z1
#                 AB = math.sqrt(ABx**2 + ABy**2 + ABz**2)
#                 cosT = ABx*xbox / (AB * xbox)
#                 cosT2 = cosT * cosT
#                 coslist.append(cosT2)

#     tframes = len(timesteps)
#     num_chains = len(coslist) // tframes  # Number of chains per frame
#     f = [0] * tframes
#     for i in range(tframes):
#         mean = np.mean(coslist[i * num_chains: (i + 1) * num_chains])
#         f[i] = ((3 * mean) - 1) / 2

#     return strain, f
    
# z - direction

def calculate_orientation_order_parameter(timesteps, min_chain_length, max_chain_length):
    coslist = []
    strain = []
    initial_zbox = timesteps[list(timesteps.keys())[0]]['z'] - timesteps[list(timesteps.keys())[0]]['-z']
    for timestep in timesteps:
        zbox = timesteps[timestep]['z'] - timesteps[timestep]['-z']
        strain_value = (zbox - initial_zbox) / initial_zbox
        strain.append(strain_value)
        for mol in timesteps[timestep]['mols']:
            chain_length = timesteps[timestep]['mols'][mol]['length']
            if min_chain_length <= chain_length <= max_chain_length:
                atoms = timesteps[timestep]['mols'][mol]['atoms']
                x1, y1, z1 = atoms[0][2], atoms[0][3], atoms[0][4]
                x2, y2, z2 = atoms[len(atoms) - 1][2], atoms[len(atoms) - 1][3], atoms[len(atoms) - 1][4]
                #print(x2,y2,z2)
                ABx = x2 - x1
                ABy = y2 - y1
                ABz = z2 - z1
                AB = math.sqrt(ABx**2 + ABy**2 + ABz**2)
                cosT = ABz*zbox / (AB * zbox)
                cosT2 = cosT * cosT
                coslist.append(cosT2)

    tframes = len(timesteps)
    num_chains = len(coslist) // tframes  # Number of chains per frame
    f = [0] * tframes
    for i in range(tframes):
        mean = np.mean(coslist[i * num_chains: (i + 1) * num_chains])
        f[i] = ((3 * mean) - 1) / 2

    return strain, f
    
    


def write_orientation_to_csv(strain, timesteps, order_parameter, output_file):
    """
    Writes strain, timesteps, and order parameter to a CSV file.
    """
    with open(output_file, mode='w') as file:
        writer = csv.writer(file)
        # Write the header row
        writer.writerow(["Timestep", "Strain", "Orientation Parameter"])
        # Write data rows
        for ts, s, o in zip(timesteps, strain, order_parameter):
            writer.writerow([ts, s, o])
    print("Data has been saved to:", output_file)

# Example usage
input_file = 'production.dump'
timesteps = readTimesteps(input_file, 10000)  # assuming readTimesteps is defined as in your provided code
min_chain_length = 345
max_chain_length = 378
strain, order_parameter = calculate_orientation_order_parameter(timesteps, min_chain_length, max_chain_length)

output_file = 'orientation_data_1_0_360_W=10.csv'
write_orientation_to_csv(strain, list(timesteps.keys()), order_parameter, output_file)


end_time = time.time()
elapsed_time_seconds = end_time - start_time
elapsed_time_minutes = elapsed_time_seconds/60
print(f"Execution time: {elapsed_time_minutes:.2f} minutes")