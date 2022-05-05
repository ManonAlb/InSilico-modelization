import requests
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
import os


def center_water_and_protein(pdbid, pdb_title):

    ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    protein_path = ROOT + '/DYNAMIC/Fixed_protein_in_water/'
    water_path = protein_path + 'FRAMES/'
    files = os.listdir(water_path)

    #CENTER EACH BOX

    for file in files :
         if file.startswith('{}'.format(pdb_title)) and file.endswith(('.pdb')):

                ppdb = PandasPdb()
                ppdb.read_pdb(water_path + '/' + file)

                y_lenght = ppdb.df['HETATM'].loc[:,'y_coord'].max() - ppdb.df['HETATM'].loc[:,'y_coord'].min()
                x_lenght = ppdb.df['HETATM'].loc[:,'x_coord'].max() - ppdb.df['HETATM'].loc[:,'x_coord'].min()
                z_lenght = ppdb.df['HETATM'].loc[:,'z_coord'].max() - ppdb.df['HETATM'].loc[:,'z_coord'].min()

                ppdb.df['HETATM'].loc[:,'x_coord'] = (ppdb.df['HETATM'].loc[:,'x_coord']) -ppdb.df['HETATM'].loc[:,'x_coord'].max() + x_lenght/2

                ppdb.df['HETATM'].loc[:,'y_coord'] = ppdb.df['HETATM'].loc[:,'y_coord'] -ppdb.df['HETATM'].loc[:,'y_coord'].max() + y_lenght/2

                ppdb.df['HETATM'].loc[:,'z_coord'] = ppdb.df['HETATM'].loc[:,'z_coord'] -ppdb.df['HETATM'].loc[:,'z_coord'].max() + z_lenght/2

                ppdb.to_pdb( water_path + '/{}'.format(file),
                records=None,
                gz=False,
                append_newline=True) #save pdb

     #CENTER THE PROTEIN

    pdb_path = ROOT + '/DYNAMIC/Fixed_protein_in_water' + '/' +  f"{pdbid}.pdb"
    pdb_url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    r = requests.get(pdb_url)
    r.raise_for_status()

    with open(pdb_path, "wb") as f:
            f.write(r.content)

    translated_protein = open( protein_path + 'translated_prot.pdb', 'w' )

    with open(pdb_path, 'r') as ubiq :
        for lines in ubiq :
            if lines.startswith('ATOM'):
                translated_protein.write(lines)
        translated_protein.close()


    ppdb = PandasPdb()
    ppdb.read_pdb('{}translated_prot.pdb'.format(protein_path))

    y_lenght = ppdb.df['ATOM'].loc[:,'y_coord'].max() - ppdb.df['ATOM'].loc[:,'y_coord'].min()
    x_lenght = ppdb.df['ATOM'].loc[:,'x_coord'].max() - ppdb.df['ATOM'].loc[:,'x_coord'].min()
    z_lenght = ppdb.df['ATOM'].loc[:,'z_coord'].max() - ppdb.df['ATOM'].loc[:,'z_coord'].min()


    ppdb.df['ATOM'].loc[:,'x_coord'] = (ppdb.df['ATOM'].loc[:,'x_coord']) -ppdb.df['ATOM'].loc[:,'x_coord'].max() + x_lenght/2

    ppdb.df['ATOM'].loc[:,'y_coord'] = ppdb.df['ATOM'].loc[:,'y_coord'] -ppdb.df['ATOM'].loc[:,'y_coord'].max() + y_lenght/2

    ppdb.df['ATOM'].loc[:,'z_coord'] = ppdb.df['ATOM'].loc[:,'z_coord'] -ppdb.df['ATOM'].loc[:,'z_coord'].max() + z_lenght/2

    ppdb.to_pdb('{}/translated_prot.pdb'.format(protein_path),
                    records=None,
                    gz=False,
                    append_newline=True) #save pdb
    return('complete')

center_water_and_protein('1UBQ', 'Water_box')



#ADDITION


def add_protein_still_to_water(protein_name):


    ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    name = 'translated_prot.pdb'
    protein_path = ROOT + '/DYNAMIC/Fixed_protein_in_water/'
    water_path = protein_path + 'FRAMES/'
    head = "CRYST1  150.000  150.000  150.000  90.00  90.00  90.00 P 1           1 " # dimension of the box
    files = os.listdir(water_path)


    with open( protein_path + name, 'r', encoding='utf-8-sig') as prot :
                protein_pdb = prot.read()

    for file in files:

        if file.startswith(('Water_box')) and file.endswith(('.pdb')):

            output = open( water_path  + protein_name+ '_' + file.strip('/'), 'w')

            #Remove 'END' from water pdb
            frame = open(water_path + file.strip('/'), 'r')
            lines=frame.readlines()
            frame.close()

            with open(water_path + file.strip('/'), 'w') as w :
                for line in lines:
                    if line.startswith("END"):
                        print(line)
                    else :
                        w.write(line)
                w.close()

#MERGE BOTH
            with open(water_path + file.strip('/'), 'r') as f:
                water_pdb=f.read()

            output.write( water_pdb + protein_pdb)

    return('complete')


add_protein_still_to_water( 'ubiquitin' )
