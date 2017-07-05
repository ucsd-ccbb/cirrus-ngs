__author__ = 'Guorong Xu<g1xu@ucsd.edu>'

import os
import sys
import subprocess
from multiprocessing import Pool

def upload_file(input_file):

    print input_file

    s3_location = ""
    subprocess.call(["aws", "s3", "cp", input_file, s3_location])

if __name__ == '__main__':
    work_dir = sys.argv[1]
    file_extension = sys.argv[2]

    file_list = []

    for dirpath, directories, filenames in os.walk(work_dir):
        for filename in filenames:
            if filename.endswith(file_extension):
                input_file = os.path.join(dirpath, filename)
                file_list.append(input_file)

    pool = Pool(4)
    pool.map(upload_file, file_list)
    pool.close()
    pool.join()