import os
import numpy as np
from rdkit import Chem
import matplotlib.pyplot as plt

# Funcion la cual es la principal para obtener distancias
def obtener_distancias(mol, idx):
    # Elimina hidrógenos sin vecinos, esto es porque daba un warning pero no se soluciono 
    mol = Chem.RemoveHs(mol, implicitOnly=False)

    # Se obtiene las coordenadas de los atomos
    coords = mol.GetConformer().GetPositions()
    distancias_por_tipo = {}

    # Se agrupan las distancias por tipo de enlace
    for i in range(mol.GetNumAtoms()):
        simbolo_i = mol.GetAtomWithIdx(i).GetSymbol()

        # Condicion para que no tome los enlaces C - C y O - H
        for j in range(i+1, mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetSymbol() == 'H' and mol.GetAtomWithIdx(j).GetSymbol() == 'H' or mol.GetAtomWithIdx(i).GetSymbol() == 'O' and mol.GetAtomWithIdx(j).GetSymbol() == 'H':
                continue 
            # Condicion para que la distancia sea de menos 1.9 Å
            distancia = np.linalg.norm(coords[i] - coords[j])
            if distancia < 1.9:
                simbolo_j = mol.GetAtomWithIdx(j).GetSymbol()
                tipo_enlace = f"{simbolo_i}-{simbolo_j}"
                
                # Condición para no aparecer los mismos enlaces una y otra vez y agruparlos por 
                # tipo de enlace 
                if tipo_enlace not in distancias_por_tipo:
                    distancias_por_tipo[tipo_enlace] = []
                
                distancias_por_tipo[tipo_enlace].append(distancia)

    return distancias_por_tipo

# Si existe output.txt lo elimina, si no se sobreescribe 
archivo_salida = 'output.txt'
if os.path.exists(archivo_salida):
    os.remove(archivo_salida)

# Variable que obtiene el nombre de las carpetas
carpetas = [f"moltiverse_{i}" for i in range(1, 11)]

# Almacenar distancias por tipo de enlace y molecula
distancias_molecula = {carpeta: {} for carpeta in carpetas}

for carpeta in carpetas:
    carpeta_actual = os.path.join(".", carpeta)

    # For que iterar sobre cada archivo SDF en la carpeta actual
    for archivo_sdf in os.listdir(carpeta_actual):
        if archivo_sdf.endswith(".sdf"):
            archivo_sdf_path = os.path.join(carpeta_actual, archivo_sdf)

            # Se carga el archivo SDF
            suppl = Chem.SDMolSupplier(archivo_sdf_path)

            # For que iterar sobre cada molécula en el archivo SDF
            for idx, mol in enumerate(suppl):
                if mol is not None:
                    distancias_tipo = obtener_distancias(mol, idx)

                    # Se almacenan las distancias
                    for tipo_enlace, distancias in distancias_tipo.items():
                        if tipo_enlace not in distancias_molecula[carpeta]:
                            distancias_molecula[carpeta][tipo_enlace] = []

                        distancias_molecula[carpeta][tipo_enlace].extend(distancias)
                        
# Creacion del boxplot
Filas = 5
Columnas = 2

fig, axs = plt.subplots(Filas, Columnas, figsize=(12, 3 * Filas), sharex=True)

# output que almacena las distancias
output_file = 'Distancias_boxplot.txt'

for i in range(Filas):
    for j in range(Columnas):
        index = i * Columnas + j

        if index < len(carpetas):
            carpeta = carpetas[index]
            distancias_por_tipo = distancias_molecula[carpeta]

            tipos_enlace = list(distancias_por_tipo.keys())
            distancias_organizadas = list(distancias_por_tipo.values())

            color = plt.cm.viridis(index / len(carpetas))

            axs[i, j].boxplot(distancias_organizadas, labels=tipos_enlace, vert=False,
                              boxprops=dict(color='black'),
                              medianprops=dict(color='black'))

            axs[i, j].set_title(f'{carpeta}')

            # Se guardan las distancias en el archivo txt output
            with open(output_file, 'a') as f:
                f.write(f'\nResults for Carpeta {carpeta}:\n')
                for tipo_enlace, distancias in distancias_por_tipo.items():
                    f.write(f'{tipo_enlace}: {distancias}\n')

# Ejex X e Y
for ax in axs[-1, :]:
    ax.set_xlabel('Distancia de enlace (Å)')

for ax in axs[:, 0]:
    ax.set_ylabel('Tipo de Enlace')

# Se guarda la figura en un png en caso de volver a utilizarla
plt.tight_layout()
plt.savefig('boxplot_figure.png')
plt.show()
