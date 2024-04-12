import numpy as np
import csv

def readcsv8(infile): # 8-colume csv ignoring first 3 lines ("infile"=open('filename'))
    arr = np.zeros((1,8))
    csv_reader = csv.reader(infile)
    next(csv_reader) # Skip 4 title lines in our csv file
    next(csv_reader)
    next(csv_reader)
    for row in csv_reader:
        arr = np.append(arr, [[row[0],row[1],row[2],row[3],
                               row[4],row[5],row[6],row[7]]], axis=0)
    arr = np.delete(arr, (0), axis=0)
    infile.close()
    return arr
