import csv
import numpy as np


# def read_csv_to_matrix(file_path):
#     with open(file_path, newline='', encoding='utf-8') as csvfile:
#         reader = csv.reader(csvfile)
#         col_names = next(reader)[1:]
#         index = []
#         data = []
#         for row in reader:
#             try:
#                 index.append(row[0])
#                 data.append(list(map(float, row[1:])))
#             except ValueError:
#                 print(f"non-numeric value found in row {row[0]}, skiping this row")
#     return np.array(data), col_names, index

def read_csv_to_matrix(file_path):
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        # Use DictReader instead of Reader for memory efficiency
        reader = csv.DictReader(csvfile)
        # Get the column names from the first row
        col_names = reader.fieldnames
        # Initialize empty lists for the data and the index
        data = []
        index = []
        for row in reader:
            # Append the first column as index
            index.append(row[col_names[0]])
            # Append the rest of the columns as data
            data.append([float(row[col]) for col in col_names[1:]])
    # Convert the data list to a matrix using numpy
    data = np.array(data)
    return data, col_names[1:], index

def read_ressource_to_matrix(file_path, delim = '\t'):
    with open(file_path, newline='', encoding='utf-8') as csvfile:
        # Use DictReader instead of Reader for memory efficiency
        # and pass the delimiter as '\t'
        reader = csv.DictReader(csvfile, delimiter=delim)
        # Get the column names from the first row
        col_names = reader.fieldnames
        # Initialize empty lists for the data 
        data = []
        for row in reader:
            # Append the columns as data
            data.append([row[col] for col in col_names])
    # Convert the data list to a matrix using numpy
    data = np.array(data)
    return data, col_names