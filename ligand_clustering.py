import argparse
import os
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def parse_args():
    parser = argparse.ArgumentParser(description="Cluster ligand PDB files and select the lowest energy conformations.")
    parser.add_argument("file_path", type=str, help="Path to the directory containing PDB files.")
    parser.add_argument("--optimal_k", type=int, default=None, help="Optimal number of clusters (K). If not provided, it will be determined using the elbow method.")
    parser.add_argument("--output_dir", type=str, default="output", help="Directory to save the output SVG files.")
    return parser.parse_args()

def read_ligand_coord_from_pdb(file_path, pdb_names):
    """ Load atomic coordinates and energy values from multiple PDB files. """
    parser = PDBParser(QUIET=True)
    data, energies = [], []
    pdb_files = []

    for pdb_name in pdb_names:
        coords = []
        energy = None
        file_full_path = os.path.join(file_path, pdb_name)
        with open(file_full_path, 'r') as f:
            for line in f:
                if "REMARK VINA RESULT" in line:
                    energy = float(line.split()[3])
                if "HETATM" in line:
                    coords.append([float(line[31:39]), float(line[39:47]), float(line[47:55])])
        
        if coords and energy is not None:
            data.append(np.array(coords))  # Keep coordinates as a 2D array
            energies.append(energy)
            pdb_files.append(pdb_name)

    # Flatten each coordinate set into one row per ligand
    flat_data = [coord.reshape(-1) for coord in data]
    return {'x': flat_data, 'y': energies, 'pdb_names': pdb_files}


def calculate_rmsd(coord1, coord2):
    """ Calculate Root Mean Square Deviation between two sets of coordinates. """
    coord1, coord2 = np.array(coord1), np.array(coord2)
    if coord1.shape != coord2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    return np.sqrt(np.mean(np.sum((coord1 - coord2)**2, axis=1)))

def k_rmsd_clustering(rmsd_matrix, k):
    """ Perform K-RMSD clustering. """
    clustering = KMeans(n_clusters=k)
    labels = clustering.fit_predict(rmsd_matrix)
    return labels

def kmeans_best_k(file_path, pdb_names, output_dir):
    data = read_ligand_coord_from_pdb(file_path, pdb_names)
    x = data['x']
    if len(x) < 1:
        print("Insufficient data to perform clustering.")
        return
    
    distortions = []
    max_k = min(len(x), 10)  # Set max_k to be the minimum of the number of ligands or 10
    K = range(1, max_k + 1)  # Ensure k does not exceed the number of samples
    
    for k in K:
        rmsd_matrix = squareform(pdist(np.array(x), metric=lambda u, v: calculate_rmsd(u.reshape(-1, 3), v.reshape(-1, 3))))
        model = KMeans(n_clusters=k).fit(rmsd_matrix)
        distortions.append(np.sqrt(model.inertia_))

    print(distortions)
    plt.plot(K, distortions, 'bx-')
    plt.xlabel('Number of clusters (K)')
    plt.ylabel('SSE')
    plt.savefig(os.path.join(output_dir, "elbow_plot.svg"))
    plt.show()



def cluster(file_path, pdb_names, k, output_dir):
    data = read_ligand_coord_from_pdb(file_path, pdb_names)
    x = np.array(data['x'])
    y = data['y']

    # 计算每对 PDB 文件之间的 RMSD 矩阵
    rmsd_matrix = squareform(pdist(x, metric=lambda u, v: calculate_rmsd(u.reshape(-1, 3), v.reshape(-1, 3))))

    # 基于 RMSD 矩阵进行 K-RMSD 聚类
    labels = k_rmsd_clustering(rmsd_matrix, k)
    
    pca = PCA(n_components=2)
    x_new = pca.fit_transform(rmsd_matrix)
    
    conformations = {}
    for i in range(len(labels)):
        if labels[i] in conformations:
            if y[i] < conformations[labels[i]]["energy"]:
                conformations[labels[i]] = {"pdb_name": data["pdb_names"][i], "energy": y[i]}
        else:
            conformations[labels[i]] = {"pdb_name": data["pdb_names"][i], "energy": y[i]}
    
    plt.scatter(x_new[:, 0], x_new[:, 1], c=labels)
    plt.xlabel("PCA1")
    plt.ylabel("PCA2")
    plt.savefig(os.path.join(output_dir, "cluster_plot.svg"))
    plt.show()
    
    return conformations

def main():
    args = parse_args()
    file_path = args.file_path
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pdb_names = [f for f in os.listdir(file_path) if f.endswith('.pdb')]
    
    if not pdb_names:
        print("No PDB files found in the specified directory.")
        return
    
    if args.optimal_k is None:
        kmeans_best_k(file_path, pdb_names, output_dir)
        optimal_k = int(input("Please enter the optimal number of clusters (K) based on the elbow plot: "))
    else:
        optimal_k = args.optimal_k
    
    conformations = cluster(file_path, pdb_names, optimal_k, output_dir)
    print("Selected conformations with lowest energy for each cluster:")
    for cluster_id, conformation in conformations.items():
        print(f"Cluster {cluster_id}: PDB Name = {conformation['pdb_name']}, Energy = {conformation['energy']}")

if __name__ == "__main__":
    main()

