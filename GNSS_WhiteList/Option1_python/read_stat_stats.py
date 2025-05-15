# script to read all the gnss observations in a folder
# it reads all the files in a given folder and creates a file called list1
# if it is called again it adds ALL the stations to the list again
#if folder contain other files or formats weird behaviours should be expected
import os
import statistics



def read_file(filename):
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

def extract_data(data_sets):
    extracted_data = []
    for data_set in data_sets:
        #print(data_set)
        line1 = data_set[0].split()

        
        id_value = line1[0]
        dat1 = line1[1]
        dat2 = line1[2]
        #dat3 = 9
        
        extracted_data.append((id_value, dat1, dat2))
    return extracted_data

def write_to_file(output_filename, filename, mean, std_dev,skew):
    with open(output_filename, 'a') as file:  # Use 'a' to append to the file
        file.write(f"{filename:<8} {float(mean):9.4f} {float(std_dev):10.6f}  {float(skew):10.6f}\n")
            

def process_files_in_folder(folder_path, output_filename):

    
    for filename in sorted(os.listdir(folder_path)):
        if filename.lower().endswith('st'):  # Add this condition
            file_path = os.path.join(folder_path, filename)
            #print(file_path)
            if os.path.isfile(file_path):
                data_sets = read_file(file_path)
                extracted_data = extract_data(data_sets)
                
                               
                fg_dep = [float(item[2]) for item in extracted_data]
                if fg_dep:
                    mean = statistics.mean(fg_dep) * 1000
                    median = statistics.median(fg_dep) * 1000                  
                    if len(fg_dep) > 1:
                        std_dev = statistics.stdev(fg_dep) * 1000
                    else:
                        std_dev = 99999  #  'N/A', or any other 
                    skew = 3 * (mean - median) / std_dev    
                        # Write results to output file
                    write_to_file(output_filename, filename, mean, std_dev,skew)
                




# Example usage
#folder_path = '/pred/pns/bell'  # Replace with the path to your folder
#data_sets = read_file(folder_path)
#extracted_data = extract_data(data_sets)
#fg_dep = [float(item[2]) for item in extracted_data]


#mean = statistics.mean(fg_dep)*1000
#std_dev = statistics.stdev(fg_dep)*1000

#print(mean, std_dev)
folder_path='/lustre/utmp/pns/calc_wl_gnss/t02'
output_filename='stats'
os.chdir(folder_path)
process_files_in_folder(folder_path, output_filename)
