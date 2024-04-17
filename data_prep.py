import argparse
import random
import numpy as np
import pandas as pd
from tqdm import tqdm


def txt2csv(data_path, dataset_str):
    df = pd.read_csv(
        data_path+dataset_str+'.txt', sep='\t', header=None, 
        skiprows=1, names=['source', 'target', 'sign'])
    df.to_csv(data_path+dataset_str+'_edges.csv', header=False, index=False, sep=',')

    node_ids = np.unique(df['source'].tolist() + df['target'].tolist())
    num_nodes = node_ids.size

    # write the number of nodes and edges to the first line
    f = open(data_path+dataset_str+'_edges.csv', mode='r+', encoding='utf-8')
    original_content = f.read()
    f.seek(0)
    f.write(f"# {str(num_nodes)} # {str(df.shape[0])}\n")
    f.write(original_content)
    f.close()


def parse_youtube_pokec(data_path, dataset_str):
    if dataset_str == "youtube":
        file_name = "com-youtube.ungraph"
        num_edges = 2987624
        skip_line = 4
    if dataset_str == "pokec":
        file_name = "soc-pokec-relationships"
        num_edges = 30622564
        skip_line = 0
    edge_ids = list(range(num_edges))
    num_positive_edges = int(0.7 * num_edges)
    positive_edge_ids = set(random.sample(edge_ids, num_positive_edges))
    source, target, sign = [], [], []
    data = open(data_path+file_name+'.txt', 'r')
    lines = data.readlines()
    data.close()
    # skip the header
    for id, line in enumerate(tqdm(lines[skip_line:], desc="Processing lines", unit="lines")):
        s, t = map(int, line.split('\t'))
        source.append(s)
        target.append(t)
        if id in positive_edge_ids:
            sign.append(1)
        else:
            sign.append(-1)

    node_ids = np.unique(source + target)
    num_nodes = node_ids.size
    # create new node IDs
    node_dict = {node: i for i, node in enumerate(node_ids)}

    f = open(data_path+dataset_str+'_edges.csv', mode='w', encoding='utf-8')
    f.write("# " + str(num_nodes) + " # " + str(num_edges) + "\n")
    for s, t, r in tqdm(zip(source, target, sign), desc="Writing to file", total=len(source), unit="edges"):
        edge = f"{node_dict[s]},{node_dict[t]},{r}"
        f.write(edge + '\n')
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=20159, help='random seed')
    parser.add_argument('--dataset', type=str, default='all', help='dataset: wikielections, epinions, wikipolitics, youtube, pokec')

    opt = parser.parse_args()

    if opt.seed is not None:
        np.random.seed(opt.seed)

    data_path = "./data/"
    opt.dataset = opt.dataset.lower()
    if opt.dataset == "all":
        txt2csv(data_path, "wikielections")
        txt2csv(data_path, "epinions")
        txt2csv(data_path, "wikipolitics")
        parse_youtube_pokec(data_path, "youtube")
        parse_youtube_pokec(data_path, "pokec")
    else:
        if opt.dataset in ['youtube', 'pokec']:
            parse_youtube_pokec(data_path, opt.dataset)
        else:
            txt2csv(data_path, opt.dataset)
