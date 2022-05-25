import os

ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
prot_path = ROOT + '/DYNAMIC/Moving_protein_in_water/FRAMES/'
entete="CRYST1  150.000  150.000  150.000  90.00  90.00  90.00 P 1           1 " # dimension of the box

files = os.listdir(prot_path)
for file in files:
    if file.endswith(('.pdb')):
        print(file)

        #read
        frame = open(prot_path + file.strip('/'), 'r')
        lines=frame.readlines()
        frame.close()

        #write
        new_file=open(prot_path + file.strip('/'), 'w+')
        for line in lines :
            if line.startswith('CRYST1'):
                change = "" + line.replace(line,entete) +"\n"
                print(change)
                new_file.write(change)
            elif line.strip('\n') !='ENDMDL' :
                new_file.write(line)
                
        new_file.close()
