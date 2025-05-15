# script to read a list of stations and extract the data from daily observations usage from an odb for future statistical analysis

import os

def read_file0(filename):
    with open(filename, 'r') as file:
        
        data_sets = []
        current_set = []
        
        for line in file:
            line = line.strip()  # Strip any leading/trailing whitespace
            current_set.append(line)
            
            if len(current_set) == 1:
                data_sets.append(current_set)
                current_set = []
        
        return data_sets

def read_file(filename):
    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        
        data_sets = []
        current_set = []
        
        for line in file:
            line = line.strip()  # Strip any leading/trailing whitespace
            current_set.append(line)
            
            if len(current_set) == 1:
                data_sets.append(current_set)
                current_set = []
        
        return data_sets

def extract_data(data_sets):
    extracted_data = []
    for data_set in data_sets:

        line1 = data_set[0].split()

        
        id_value = line1[0]
        dat1 = line1[4]
        dat2 = line1[5]
        #dat3 = 9
        
        extracted_data.append((id_value, dat1, dat2))
    return extracted_data

def write_to_file(filename, data,id1):
    with open(filename+'.st', 'a') as file:  # Use 'a' to append to the file
        for id_value, dat1, dat2 in data:

            if id_value == id1 :
                file.write(f"{id_value:<8} {float(dat1):9.4f} {float(dat2):10.6f} \n")
                

def process_files_in_folder(folder_path, output_filename,id1):
    processed_ids = set()
    
    for filename in os.listdir(folder_path):
        if filename.lower().startswith('gnss'):  # Add this condition
            file_path = os.path.join(folder_path, filename)
            #print(file_path)
            if os.path.isfile(file_path):
                data_sets = read_file(file_path)
                extracted_data = extract_data(data_sets)
                
                # Filter out entries with duplicate IDs
                #unique_data = [entry for entry in extracted_data if entry[0] not in processed_ids]
                
                # Update the set of processed IDs
                #processed_ids.update(entry[0] for entry in unique_data)
                
                write_to_file(output_filename, extracted_data,id1)

# Example usage
folder_path = '/lustre/utmp/pns/odbgnss/exp1/'  # Replace with the path to your folder

datas=read_file0('/lustre/utmp/pns/odbgnss/exp1/listuniq')
ids=extract_data(datas)
os.chdir('/lustre/utmp/pns/calc_wl_gnss/t02')
for item in ids:
    id1=item[0]
    output_filename = id1[1:-1]
    process_files_in_folder(folder_path, output_filename,id1)



