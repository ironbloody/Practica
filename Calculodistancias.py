import numpy as np
from rdkit import Chem
import os


# obtencion de distancias
def obtener_distancias(mol, conformer_number, output_file):
    # elimina hidrogenos sin vecinos 
    mol = Chem.RemoveHs(mol, implicitOnly=False)

    # Obtiene las coordenadas de los atomos
    coords = mol.GetConformer().GetPositions()
    distancias_por_tipo = {}

    # Agrupa distancias por el tipo de enlace
    for i in range(mol.GetNumAtoms()):
        # GetSymbol obtiene la letra del atomo
        simbolo_i = mol.GetAtomWithIdx(i).GetSymbol()

        # Condicion para que no tome los enlaces C - C y O - H
        for j in range(i+1, mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetSymbol() == 'H' and mol.GetAtomWithIdx(j).GetSymbol() == 'H' or mol.GetAtomWithIdx(i).GetSymbol() == 'O' and mol.GetAtomWithIdx(j).GetSymbol() == 'H':
                continue 
            # Condicion para que la distancia sea de menos 1.9 Angstrom
            distancia = np.linalg.norm(coords[i] - coords[j])
            if distancia < 1.9:
                simbolo_j = mol.GetAtomWithIdx(j).GetSymbol()
                tipo_enlace = f"{simbolo_i}-{simbolo_j}"
                
                # Condicion para no aparezcan los mismos enlaces una y otra vez y ademas agruparlos por tipo de enlace
                if tipo_enlace not in distancias_por_tipo:
                    distancias_por_tipo[tipo_enlace] = []
                
                distancias_por_tipo[tipo_enlace].append((simbolo_i, simbolo_j, distancia))

    # For que imprime las distancias
    output_file.write(f"Conformer {conformer_number} - carpeta {carpeta}:\n")
    for tipo_enlace, distancias in distancias_por_tipo.items():
        output_file.write(f"Bonds {tipo_enlace} links:\n")
        for simbolo_i, simbolo_j, distancia in distancias:
            output_file.write(f"  {simbolo_i} - {simbolo_j}: {distancia:.2f} Å\n")

# Carga del archivo
sdf_file = "moltiverse_1.sdf"
suppl = Chem.SDMolSupplier(sdf_file)

# Archivo txt, en este se podran ver todas las distancias por molecula y conformeros
archivo_salida = 'output.txt'
if os.path.exists(archivo_salida):
    os.remove(archivo_salida)
# Carpteas donde el programa buscara el archivo SDF
carpetas = [f"moltiverse_{i}" for i in range(1, 11)]

with open('output.txt', 'w') as output_file:
    # Abre cada una de las carpetas
    for carpeta in carpetas:
        carpeta_actual = os.path.join(".", carpeta)
        
        # y busca los archivos sdf de cada uno
        for archivo_sdf in os.listdir(carpeta_actual):
            if archivo_sdf.endswith(".sdf"):
                archivo_sdf_path = os.path.join(carpeta_actual, archivo_sdf)
                
                # Se carga el archivo SDF
                suppl = Chem.SDMolSupplier(archivo_sdf_path)
                
                # For que iterara sobre cada conformero en el archivo SDF
                for idx, mol in enumerate(suppl):
                    if mol is not None:
                        print(f"CONFORMERO {idx} - Carpeta: {carpeta}, Archivo: {archivo_sdf}")
                        obtener_distancias(mol, idx, output_file)
                    else:
                        output_file.write(f"No se pudo cargar la molécula desde {archivo_sdf_path}.\n")
