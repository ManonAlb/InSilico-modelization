from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
import os



def clean_protein(pdbid):

    ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))

    protein_path = ROOT + '/DYNAMIC/Fixed_protein_in_water' + '/' +  "{}.pdb".format(pdbid)
    clean_protein = open( ROOT + '/DYNAMIC/Fixed_protein_in_water' + '/' + 'clean_prot.pdb', 'w' )

    with open(protein_path, 'r') as prot :
        for lines in prot :
            if lines.startswith('ATOM'):
                clean_protein.write(lines)
    return('cleaned')

def translated_protein(pdbid):

    ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    clean_protein_path = ROOT + '/DYNAMIC/Fixed_protein_in_water' + '/'
    ppdb = PandasPdb()
    ppdb.read_pdb('{}clean_prot.pdb'.format(clean_protein_path))

    y_lenght = ppdb.df['ATOM'].loc[:,'y_coord'].max() - ppdb.df['ATOM'].loc[:,'y_coord'].min()
    x_lenght = ppdb.df['ATOM'].loc[:,'x_coord'].max() - ppdb.df['ATOM'].loc[:,'x_coord'].min()
    z_lenght = ppdb.df['ATOM'].loc[:,'z_coord'].max() - ppdb.df['ATOM'].loc[:,'z_coord'].min()


    ppdb.df['ATOM'].loc[:,'x_coord'] = (ppdb.df['ATOM'].loc[:,'x_coord']) -ppdb.df['ATOM'].loc[:,'x_coord'].max() + x_lenght/2

    ppdb.df['ATOM'].loc[:,'y_coord'] = ppdb.df['ATOM'].loc[:,'y_coord'] -ppdb.df['ATOM'].loc[:,'y_coord'].max() + y_lenght/2

    ppdb.df['ATOM'].loc[:,'z_coord'] = ppdb.df['ATOM'].loc[:,'z_coord'] -ppdb.df['ATOM'].loc[:,'z_coord'].max() + z_lenght/2

    ppdb.to_pdb('{}translated_{}.pdb'.format(clean_protein_path, pdbid),
                records=None,
                gz=False,
                append_newline=True) #save pdb
    return('protein centered')

def center_water( pdb_title):

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

                ppdb.to_pdb( water_path + '/translated_{}'.format(file),
                records=None,
                gz=False,
                append_newline=True) #save pdb


    return('complete')

def add_protein_still_to_water(protein_id):


    ROOT = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    name = 'translated_{}.pdb'.format(protein_id)
    protein_path = ROOT + '/DYNAMIC/Fixed_protein_in_water/'
    water_path = protein_path + 'FRAMES/'
    head = "CRYST1  150.000  150.000  150.000  90.00  90.00  90.00 P 1           1 " # dimension of the box
    files = os.listdir(water_path)


    with open( protein_path + name, 'r', encoding='utf-8-sig') as prot :
                protein_pdb = prot.read()

    for file in files:

        if file.startswith(('Water_box')) and file.endswith(('.pdb')):

            output = open( water_path  + protein_id+ '_' + file.strip('/'), 'w')

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



clean_protein('1ubq')
translated_protein('1ubq')
#center_water('Water_box')
add_protein_still_to_water( '1ubq' )
