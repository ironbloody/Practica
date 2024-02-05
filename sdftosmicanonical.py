import sys
from rdkit import Chem

# Funcion para convertir el archivo sdf en smile
def convertir(file_name):
    # Reemplazar por el nombre de su archivo
    sdf_file = Chem.SDMolSupplier("platinum.sdf")
    output = open('smiles.txt', "w")
    for mol in sdf_file:
        if mol is not None:
            # Solo se busca que sean canonicos y no isomericos
            smi = Chem.MolToSmiles(mol, isomericSmiles=False)
            output.write(f"{smi}\n")
    output.close()
    
    try:
        # Web de donde se obtienen los smiles
        url ="https://cactus.nci.nih.gov/chemical/structure/" + smi+"/iupac_name" 
        ans = urlopen(url).read().decode('utf8')
        return ans
    except:
        return 'Molecula no disponible'
        
if __name__ == "__main__":
    convertir(sys.argv[0])

