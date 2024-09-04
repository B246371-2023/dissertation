import json
import treescore
import pandas as pd
import toytree
import numpy as np
from scipy.stats import describe

# 输入文件路径
sequence_dataset_path = "/home/s2530615/fold_tree/CR/sequence_dataset.csv"
tree_file_paths = ["/home/s2530615/fold_tree/CR/sequences.aln.fst.nwk.rooted.final"]

# 读取Uniprot数据集
uniprot_df = pd.read_csv(sequence_dataset_path)
scores = {}
stats = {}

for t in tree_file_paths:
    tree = toytree.tree(t, format=0)
    lineages = treescore.make_lineages(uniprot_df)
    tree = treescore.label_leaves(tree, lineages)
    overlap = treescore.getTaxOverlap(tree.treenode)
    taxscore = tree.treenode.score

    lengths = np.array([node.dist for node in tree.treenode.traverse()])
    lengths /= np.sum(lengths)

    treescore.getTaxOverlap_root(tree.treenode)
    root_score = treescore.sum_rootscore(tree.treenode)

    distances = np.array([node.get_distance(tree.treenode) for node in tree.treenode.get_leaves()])
    distances_norm = distances / np.mean(distances)

    scores[t] = {
        'score': taxscore,
        'stats': describe(lengths),
        'ultrametricity': describe(distances),
        'ultrametricity_norm': describe(distances_norm),
        'root_score': root_score
    }

# 输出文件路径
output_path = "/home/s2530615/fold_tree/test/output/treescores_seq_tree.json"
with open(output_path, 'w') as snakeout:
    snakeout.write(json.dumps(scores))

print(f"TaxScore计算完成，结果已保存到 {output_path}")

