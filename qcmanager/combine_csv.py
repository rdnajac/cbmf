#!/usr/bin/env python3
#
## CSV utils

# function to comibine multuple csv files into one csv file
# the function looks in multiple folders for files of the same name
# files are as follows:
#     */performance_metrics.csv
#     */flagstat.csv
#     */summary.csv
#
# make the new first column, the folder name and add that to the csv file
# make sure we only add the header once

import csv
import os
import sys

def combine_csv_files(meta_folder):
    file_names = ['performance_metrics.csv', 'flagstat.csv', 'summary.csv']
    combined_files = {file_name: f"combined_{file_name}" for file_name in file_names}
    header_written = {file_name: False for file_name in file_names}
    
    # Initialize writers for each combined output file
    writers = {}
    output_files = {}
    for file_name in file_names:
        output_file_path = os.path.join(meta_folder, combined_files[file_name])
        output_files[file_name] = open(output_file_path, 'w', newline='')
        writers[file_name] = csv.writer(output_files[file_name])
    
    # Iterate through each subfolder in the meta-folder
    for subfolder in os.listdir(meta_folder):
        subfolder_path = os.path.join(meta_folder, subfolder)
        if os.path.isdir(subfolder_path):
            for file_name in file_names:
                file_path = os.path.join(subfolder_path, file_name)
                if os.path.isfile(file_path):
                    with open(file_path, 'r') as in_csv:
                        reader = csv.reader(in_csv)
                        header = next(reader)
                        
                        if not header_written[file_name]:
                            writers[file_name].writerow(['Folder'] + header)
                            header_written[file_name] = True
                        
                        for row in reader:
                            writers[file_name].writerow([subfolder] + row)
    
    # Close all output files
    for file_name in file_names:
        output_files[file_name].close()

# # Example usage
# meta_folder = 'metafolder'
# combine_csv_files()
# use argv
if __name__ == '__main__':
    meta_folder = sys.argv[1]
    combine_csv_files(meta_folder)
