# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


#This script includes contributions of multiple authors; an original version was written by Mark Stevens, 
#and it was modified in the research group of Lisa Hall and later by the Taofeek Tejuosho from the Janani Sampath group

import math
import numpy as np
import csv
import time

start_time = time.time()

def readTimesteps(inputFile, timestepsCount=None):
    if timestepsCount:
        counter = 0
    else:
        timestepsCount = 5
    timesteps = {}
    with open(inputFile) as f:
        line = f.readline()
        while line != '':
            if line.split()[0] == "ITEM:":
                timestep = int(f.readline().split()[0])
                timesteps[timestep] = {}
                f.readline()
                atoms = int(f.readline().split()[0])
                timesteps[timestep]["atoms"] = atoms
                f.readline()
                x = f.readline()
                timesteps[timestep]["-x"] = float(x.split()[0])
                timesteps[timestep]["x"] = float(x.split()[1])
                y = f.readline()
                timesteps[timestep]["-y"] = float(y.split()[0])
                timesteps[timestep]["y"] = float(y.split()[1])
                z = f.readline()
                timesteps[timestep]["-z"] = float(z.split()[0])
                timesteps[timestep]["z"] = float(z.split()[1])
                timesteps[timestep]["headers"] = f.readline()
                timesteps[timestep]["mols"] = {}

                for i in range(atoms):
                    line = f.readline()
                    tokens = line.split()
                    id = int(tokens[0])
                    mol = int(tokens[1])
                    type = int(tokens[2])
                    q = int(tokens[3])
                    xs = float(tokens[4])
                    ys = float(tokens[5])
                    zs = float(tokens[6])
                    ix = int(tokens[7])
                    iy = int(tokens[8])
                    iz = int(tokens[9])
                    xbox = timesteps[timestep]['x'] - timesteps[timestep]['-x']
                    ybox = timesteps[timestep]['y'] - timesteps[timestep]['-y']
                    zbox = timesteps[timestep]['z'] - timesteps[timestep]['-z']
                    x = xs * xbox + ix * xbox
                    y = ys * ybox + iy * ybox
                    z = zs * zbox + iz * zbox
                    atom_data = (id, type, q, x, y, z, ix, iy, iz)

                    if mol not in timesteps[timestep]["mols"]:
                        timesteps[timestep]["mols"][mol] = {"atoms": []}
                    timesteps[timestep]["mols"][mol]["atoms"].append(atom_data)

                for mol in timesteps[timestep]["mols"]:
                    atoms = timesteps[timestep]["mols"][mol]["atoms"]
                    length = len(atoms)
                    xcm = sum([atom[3] for atom in atoms]) / length
                    ycm = sum([atom[4] for atom in atoms]) / length
                    zcm = sum([atom[5] for atom in atoms]) / length
                    timesteps[timestep]["mols"][mol]["xcm"] = xcm
                    timesteps[timestep]["mols"][mol]["ycm"] = ycm
                    timesteps[timestep]["mols"][mol]["zcm"] = zcm
                    timesteps[timestep]["mols"][mol]["length"] = length
            line = f.readline()
            counter += 1
            if counter == timestepsCount:
                return timesteps
    return timesteps

def calculate_structure_factor_anisotropic_by_frame_range(timesteps, q_values, start_frame=0, end_frame=None, min_length=None, max_length=None):
    S_q_parallel = np.zeros_like(q_values, dtype=np.float64)  # q || z (parallel)
    S_q_perp = np.zeros_like(q_values, dtype=np.float64)      # q ⟂ z (in x-y plane)
    count_parallel = 0
    count_perp = 0

    timestep_keys = sorted(timesteps.keys())
    total_frames = len(timestep_keys)

    if end_frame is None or end_frame > total_frames:
        end_frame = total_frames

    for frame_idx in range(start_frame, end_frame):
        timestep = timestep_keys[frame_idx]

        for mol in timesteps[timestep]["mols"]:
            N = timesteps[timestep]["mols"][mol]["length"]
            if (min_length and N < min_length) or (max_length and N > max_length):
                continue

            positions = np.array([(atom[3], atom[4], atom[5]) for atom in timesteps[timestep]["mols"][mol]["atoms"]])
            x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]

            for i, q in enumerate(q_values):
                # PARALLEL to deformation (along z)
                phase_z = np.exp(1j * q * z)
                Sz = (1 / N) * np.abs(np.sum(phase_z)) ** 2
                S_q_parallel[i] += Sz
                count_parallel += 1

                # PERPENDICULAR to deformation (x-y plane, radial)
                qx = q / np.sqrt(2)
                qy = q / np.sqrt(2)
                phase_perp = np.exp(1j * (qx * x + qy * y))
                S_perp = (1 / N) * np.abs(np.sum(phase_perp)) ** 2
                S_q_perp[i] += S_perp
                count_perp += 1

    if count_parallel > 0:
        S_q_parallel /= count_parallel
    if count_perp > 0:
        S_q_perp /= count_perp

    return q_values, S_q_perp, S_q_parallel

def write_structure_factor_to_csv(q_values, S_q_perp, S_q_parallel, outFile):
    with open(outFile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["q", "S(q)_perp_xy", "S(q)_parallel_z"])
        for q, sq_perp, sq_parallel in zip(q_values, S_q_perp, S_q_parallel):
            writer.writerow([q, sq_perp, sq_parallel])

# MAIN EXECUTION
input_file = "production.dump"
output_file = "sofk_anisotropic_frames_t=51M.csv"
q_values = np.logspace(-3, 3, 200)  # q from 0.01 to 100

# Read full trajectory
timesteps = readTimesteps(input_file, timestepsCount=2000)

# Frame range (not timestep values!)
start_frame = 1029    # 6th frame
end_frame = 1032 # up to (but not including) the 21st frame
min_length = 345
max_length = 378

q_values, S_q_perp, S_q_parallel = calculate_structure_factor_anisotropic_by_frame_range(
    timesteps,
    q_values,
    start_frame=start_frame,
    end_frame=end_frame,
    min_length=min_length,
    max_length=max_length
)

write_structure_factor_to_csv(q_values, S_q_perp, S_q_parallel, output_file)

end_time = time.time()
print(f"Execution time: {(end_time - start_time)/60:.2f} minutes")
print("DONE calculation")
