#Script to pick one processing center when there is multiple for one station
#The statistics in the file 'stats' must be in alphabetical order!!!
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

def extract_data2(data_sets):
    extracted_data = []
    for data_set in data_sets:
        #print(data_set)
        line1 = data_set[0].split()

        
        id_value = line1[0]
        dat1 = line1[1]
        dat2 = line1[2]
        dat3 = line1[3]        
        #dat3 = 9
        
        extracted_data.append((id_value, dat1, dat2,dat3))
    return extracted_data

def write_to_file(output_filename, lat, lon, alt, filename, mean, std_dev):
    with open(output_filename, 'a') as file:  # Use 'a' to append to the file
        file.write(f"{filename:<8} {float(lat):5.2f} {float(lon):5.2f} {float(alt):7.2f} {15.0:4.1f}{mean:8.4f} {float(std_dev):6.2f} {0.1:4.2f} \n")
            

def process_files_in_folder(folder_path, output_filename):

    
    for filename in sorted(os.listdir(folder_path)):
        if filename.lower().endswith('st'):  # Add this condition
            file_path = os.path.join(folder_path, filename)
            #print(file_path)
            if os.path.isfile(file_path):
                data_sets = read_file(file_path)
                extracted_data = extract_data(data_sets)
                
                               
                fg_dep = [float(item[2]) for item in extracted_data]
                skew = 3 * (mean - median) / std
                if fg_dep:
                    mean = statistics.mean(fg_dep) * 1000
                    median = statistics.median(data) * 1000
                    skew = 3 * (mean - median) / std
                    if len(fg_dep) > 1:
                        std_dev = statistics.stdev(fg_dep) * 1000
                    else:
                        std_dev = 99999  #  'N/A', or any other 
                    skew = 3 * (mean - median) / std    
                        # Write results to output file
                    write_to_file(output_filename, filename, mean, std_dev,skew)
                
def sort():
    data_sets = read_file('stats')
    extracted_data = extract_data(data_sets)
    ids = [item[0] for item in extracted_data]
    print (ids)
    mean=[float(item[1]) for item in extracted_data]
    sd = [item[2] for item in extracted_data]
    #it is the comparing the first four characters but the name has a quotation mark... so 1:4 instead of 0:3
    cut_id = [id_string[0:4] for id_string in ids]
    print(cut_id)
    
    data_sets2 = read_file ('/lustre/utmp/pns/odbgnss/exp1/listuniq') 
    extracted_data_2 = extract_data2(data_sets2)
    lat = [item[1] for item in extracted_data_2]
    lon = [item[2] for item in extracted_data_2]
    alt = [item[3] for item in extracted_data_2]
    i = 0
    while i < len(ids)-1 :
        min_sd = sd[i]
        min_id = ids[i]
        min_mean=mean[i]
        if cut_id[i] == cut_id[i + 1]:
            #Initialize with first element of the sequence

            
            # checking if the following ids are from the same station
            while i < len(ids) - 1 and cut_id[i] == cut_id[i + 1]:
                # Compare current minimum with next element
                if sd[i + 1] < min_sd and float(sd[i + 1]) < 99999:
                    min_sd = sd[i + 1]
                    min_id = ids[i + 1]
                    min_mean = mean[i + 1]
                i += 1
            
            #Print the ID with minimum sd value for this sequence
            print(f"For ID prefix {cut_id[i]}, minimum is: {min_id} with value {min_sd}")
        if float(min_sd) < 99999: 
              
            write_to_file ('wlist',lat[i],lon[i],alt[i],min_id[0:8],min_mean/1000,min_sd)
        i += 1
    if cut_id[i] !=  cut_id[i-1]:    #this is needed to add the last element of the list that it is not taken account in the loop 
        write_to_file ('wlist',lat[i],lon[i],alt[i],ids[i][0:8],mean[i]/1000,sd[i]  )      
    return     
    


os.chdir('/lustre/utmp/pns/calc_wl_gnss/t02')
sort()

