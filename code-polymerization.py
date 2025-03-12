import numpy as np
import os
import random
def process_pdb(input_file, output_file, r_offset, x_offset, y_offset, z_offset):
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as f1:
            for l in f.readlines():
                s = int(l[23:26])
                r = int(s + r_offset)
                new_r = str(r)
                while len(new_r) != 3:
                    new_r = ' '+ new_r
                i = float(l[30:38])
                x = float(i + x_offset)
                x = round(x, 3)
                new_x = "{:0.3f}".format(x)
                while len(new_x) != 8:
    	            new_x = ' '+ new_x
                j = float(l[38:46])
                y = float(j + y_offset)
                y = round(y, 3)
                new_y = "{:0.3f}".format(y)
                while len(new_y) != 8:
    	            new_y = ' '+ new_y
                k = float(l[46:54])
                z = float(k + z_offset)
                z = round(z, 3)
                new_z = "{:0.3f}".format(z)
                while len(new_z) != 7:
    	            new_z = ' '+ new_z
                new_line = f"{l[:23] + new_r + '    ' + new_x + new_y + new_z + l[53:]}"
                f1.writelines(new_line)

process_pdb('res1.pdb', 'res1-rotated.pdb', 0, 0, 0, 0)

process_pdb('res2-temp.pdb', 'res2.pdb', 1, -0.073, -8.000, 0.966)

for i in range(2, 35):
    input_file = f'res{i}.pdb'
    output_file = f'res{i+1}.pdb'
    process_pdb(input_file, output_file, 1, -0.073, -8.000, 0.966)

process_pdb('res35-temp.pdb', 'res35.pdb', 34, -2.482, -272.000, 32.844)

############# concatenate files ###############
def concatenate_pdbs(input_files, output_file):
    with open(output_file, 'w') as outfile:
        for file in input_files:
            with open(file, 'r') as infile:
                for line in infile:
                    outfile.write(line)

# List of input files to concatenate
#input_files = ['res1.pdb', 'res2.pdb', 'res3.pdb', ..., 'res10.pdb']

# Generate the list of input files dynamically
input_files = [f'res{i}.pdb' for i in range(1, 36)]

# Output file name
output_file = '35mer.pdb'

# Concatenate the input files into the output file
concatenate_pdbs(input_files, output_file)

#####################################################################
########### rotate each monomer in a particular direction ###########
#####################################################################
# Function to calculate the centroid of a set of points
def calculate_centroid(points):
    return np.mean(points, axis=0)

# Function to rotate a point around the origin
def rotate_point(point, angle, axis):
    angle = np.radians(angle)
    if axis.lower() == 'x':
        rotation_matrix = np.array([[1, 0, 0],
                                    [0, np.cos(angle), -np.sin(angle)],
                                    [0, np.sin(angle), np.cos(angle)]])
    elif axis.lower() == 'y':
        rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                                    [0, 1, 0],
                                    [-np.sin(angle), 0, np.cos(angle)]])
    elif axis.lower() == 'z':
        rotation_matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                                    [np.sin(angle), np.cos(angle), 0],
                                    [0, 0, 1]])
    else:
        raise ValueError("Axis must be 'x', 'y', or 'z'")
    
    return np.dot(rotation_matrix, point)

# Function to read PDB file and extract atom coordinates
def read_pdb(filename):
    atoms = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append([x, y, z])
    return np.array(atoms)

# Function to write PDB file with modified coordinates
def write_pdb(filename, lines, new_coords):
    with open(filename, 'w') as file:
        atom_idx = 0
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x, y, z = new_coords[atom_idx]
                new_line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
                file.write(new_line)
                atom_idx += 1
            else:
                file.write(line)

# Define the initial rotation angle and increment
initial_angle = 20
angle_increment = 15
max_distance = 15.0

previous_centroid = None
previous_coords = None

# Define the range of file names
for i in range(2, 36):
    input_filename = f"res{i}.pdb"
    output_filename = f"res{i}-rotated.pdb"

    # Calculate rotation angle for this iteration
    rotation_angle = (initial_angle + (i - 2) * angle_increment) % 360

    # Read input PDB file
    atoms = read_pdb(input_filename)
    with open(input_filename, 'r') as f:
        lines = f.readlines()

    # Determine if this file pair should be rotated by 180 degrees
    rotate_180 = random.choice([True, False])

    if rotate_180:
        # Rotate coordinates
        atoms = [rotate_point(atom, 180, 'y') for atom in atoms]

    rotated_coords = np.array(atoms)

    # Calculate the centroid of the rotated coordinates
    centroid = np.mean(rotated_coords, axis=0)

    # Adjust position to ensure centroids are not more than 10 Å apart
    if previous_centroid is not None:
        direction_vector = centroid - previous_centroid
        distance = np.linalg.norm(direction_vector)
        if distance > max_distance:
            adjustment_vector = (direction_vector / distance) * (distance - max_distance)
            rotated_coords -= adjustment_vector

        # Ensure no overlap by slightly translating if needed
        while previous_coords is not None and np.any(np.linalg.norm(rotated_coords[:, None] - previous_coords, axis=2) < 1.0):
            rotated_coords += np.random.uniform(-1.0, 1.0, rotated_coords.shape)

    # Write rotated and adjusted coordinates to output file
    write_pdb(output_filename, lines, rotated_coords)

    # Update previous centroid and coordinates
    previous_centroid = np.mean(rotated_coords, axis=0)
    previous_coords = rotated_coords
############# concatenate files ################
def concatenate_pdbs(input_files, output_file):
    with open(output_file, 'w') as outfile:
        for file in input_files:
            with open(file, 'r') as infile:
                for line in infile:
                    outfile.write(line)

# List of input files to concatenate
#input_files = ['res1.pdb', 'res2.pdb', 'res3.pdb', ..., 'res10.pdb']

# Generate the list of input files dynamically
input_files = [f'res{i}-rotated.pdb' for i in range(1, 36)]

# Output file name
output_file = '35mer-rotated.pdb'

# Concatenate the input files into the output file
concatenate_pdbs(input_files, output_file)
############# delete files ###############
# Define the directory containing the files
directory = "path"

# Define the files you want to delete
files_to_delete = [f"res{i}.pdb" for i in range(1, 36) if i != 1]
more_files_to_delete = [f"res{i}-rotated.pdb" for i in range(1, 36) if i != 1]

# Iterate through the files and delete them
for filename in files_to_delete:
    filepath = os.path.join(directory, filename)
    if os.path.exists(filepath):
        os.remove(filepath)
        print(f"Deleted: {filepath}")
    else:
        print(f"File not found: {filepath}")
# Iterate through the files and delete them
for filename in more_files_to_delete:
    filepath = os.path.join(directory, filename)
    if os.path.exists(filepath):
        os.remove(filepath)
        print(f"Deleted: {filepath}")
    else:
        print(f"File not found: {filepath}")
