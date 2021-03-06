import wget
import os
import sys

_, output_dir, num_files, full_url = sys.argv

if not os.path.isdir(output_dir):
    
    os.mkdir(output_dir)

for i in range(int(num_files)):
    
    # fix problems with numbers of different digit lengths later in a notebook
    new_fName = os.path.join(output_dir, "file_" + str(i+1) + ".tsv")
        
    # need try/except because if a protein has no domains, there is no file generated by HMMER, so wget will throw an error
    if not os.path.isfile(new_fName):
    
        try:
            url = full_url + str(i+1) + "/score?format=tsv"
            filename = wget.download(url, out=new_fName)
        except:
            pass