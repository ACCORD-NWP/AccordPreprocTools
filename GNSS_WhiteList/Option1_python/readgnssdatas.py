# script to read all the gnss observations in a folder
# it reads all the files in a given folder and creates a file called list1
# if it is called again it adds ALL the stations to the list again
# if folder contain other files or formats weird behaviours should be expected
import os

def read_file(filename):
    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        
        data_sets = []
        current_set = []
        
        for line in file:
            line = line.strip()  # Strip any leading/trailing whitespace
            current_set.append(line)
            
            if len(current_set) == 3:
                data_sets.append(current_set)
                current_set = []
        
        return data_sets

def extract_data(data_sets):
    extracted_data = []
    for data_set in data_sets:
        line1 = data_set[0].split()
        line2 = data_set[1].split()
        
        id_value = line1[-1]
        dat1 = line1[3]
        dat2 = line1[4]
        dat3 = line2[2]
        
        extracted_data.append((id_value, dat1, dat2, dat3))
    return extracted_data

def write_to_file(filename, data):
    with open(filename, 'a') as file:  # Use 'a' to append to the file
        for id_value, dat1, dat2, dat3 in data:
            file.write(f"{id_value:<8} {float(dat1):7.2f} {float(dat2):7.2f} {float(dat3):7.2f} 17.0 0.0220 9.01 0.10\n")

def process_files_in_folder(folder_path, output_filename):
    processed_ids = set()
    
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            data_sets = read_file(file_path)
            extracted_data = extract_data(data_sets)
            
            # Filter out entries with duplicate IDs
            unique_data = [entry for entry in extracted_data if entry[0] not in processed_ids]
            
            # Update the set of processed IDs
            processed_ids.update(entry[0] for entry in unique_data)
            
            write_to_file(output_filename, unique_data)

# Example usage
folder_path = '/lustre/utmp/pns/gnssobs1/'  # Replace with the path to your folder
output_filename = 'list20'
process_files_in_folder(folder_path, output_filename)

print("Data has been written to the output file.")

