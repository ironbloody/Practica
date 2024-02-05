from Bio.PDB import PDBParser, Superimposer
import matplotlib.pyplot as plt
import os
import re

# funcion que calcula el RMSD
def calcular_RMSD(structure):
    valores_RMSD = []
    models = list(structure.get_models())
    num_models = len(models)

    for i in range(num_models - 1):
        modelo_ref = models[i]
        atomos_ref = [atom for atom in modelo_ref.get_atoms()]

        for j in range(i + 1, num_models):
            modelo = models[j]
            atomo = [atom for atom in modelo.get_atoms()]
            
            # Superpone los atomos
            superimposer = Superimposer()
            superimposer.set_atoms(atomos_ref, atomo)
            superimposer.apply(atomo)
            
            # Calcula RMSD
            rmsd = superimposer.rms
            valores_RMSD.append(rmsd)

    return valores_RMSD

if __name__ == "__main__":
    # directorio donde estan las carpetetas moltiverse, cambiar por el destino de su carpeta
    directorio_raiz = '/home/kaladin/Descargas/data'
    patron_archivos = re.compile(r'moltiverse_(\d+)\.pdb')

    archivos_abiertos = set()
    pdb_files = []
    
    # For hecho para recorrer los archivos y ademas que no abra archivos ya abiertots
    for root, dirs, files in os.walk(directorio_raiz):
        for archivo in files:
            match = patron_archivos.match(archivo)
            if match:
                num_archivo = int(match.group(1))
                if 1 <= num_archivo <= 10:
                    nuevo_nombre = f"moltiverse_{num_archivo}.pdb"
                    if nuevo_nombre not in archivos_abiertos:
                        pdb_files.append(os.path.join(root, nuevo_nombre))
                        archivos_abiertos.add(nuevo_nombre)
                        print(f"Archivo: {nuevo_nombre}")
    
    # Creacion grilla
    filas = 5
    columnas = 2
    total_subplots = 5 * 2
    num_archivos = len(pdb_files)
    fig, axs = plt.subplots(filas, columnas, figsize=(12, 3 * filas), sharex=True, sharey=True)
   
    # Rellenar grilla
    for i, pdb_file in enumerate(pdb_files):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)

        fila = i // columnas
        columna = i % columnas

        if fila < filas and columna < columnas:
            ax = axs[fila, columna]

            valores_RMSD = calcular_RMSD(structure)
            color = plt.cm.viridis(i / len(pdb_files))
            ax.hist(valores_RMSD, bins=20, color=color, edgecolor='black', alpha=0.7)
            ax.set_xlabel('RMSD (Å)')
            ax.set_ylabel('Frecuencia')
            ax.set_title(f'RMSD - {os.path.basename(pdb_file)}')

            if valores_RMSD:
                average_rmsd = sum(valores_RMSD) / len(valores_RMSD)
                ax.axvline(average_rmsd, color='red', linestyle='dashed', linewidth=2, label=f'Promedio: {average_rmsd:.3f} Å')
                ax.legend()

    plt.tight_layout()
    plt.show()

