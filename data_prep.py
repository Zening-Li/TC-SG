import os
import random
import argparse
import numpy as np
import pandas as pd


def make_directed_to_undirected(df: pd.DataFrame):
    df["min_node"] = df[["src", "dst"]].min(axis=1)
    df["max_node"] = df[["src", "dst"]].max(axis=1)

    df.drop(["src", "dst"], axis=1, inplace=True)
    df.drop_duplicates(subset=["min_node", "max_node"], keep="last", inplace=True)

    df.rename(columns={"min_node": "src", "max_node": "dst"}, inplace=True)

    df.sort_values(by=["src", "dst"], inplace=True)
    return df


def write_edges_to_file(df: pd.DataFrame, output_path: str):
    nodes = pd.unique(df[["src", "dst"]].values.ravel("K"))
    node_dict = {node: i for i, node in enumerate(nodes)}

    df["src"] = df["src"].map(node_dict)
    df["dst"] = df["dst"].map(node_dict)

    with open(output_path, mode="w", encoding="utf-8") as f:
        f.write(f"# {nodes.shape[0]} # {df.shape[0]}\n")
        df.to_csv(f, columns=["src", "dst", "sign"], header=False, index=False, sep=",")


def parse_scale_data(data_path: str, dataset: str, pos_ratio: float = 0.7):
    if dataset == "youtube":
        file_path = os.path.join(data_path, "com-youtube.ungraph.txt")
        skiprows = 4
        sep_char = "\t"
    elif dataset == "pokec":
        file_path = os.path.join(data_path, "soc-pokec-relationships.txt")
        skiprows = 0
        sep_char = "\t"
    elif dataset == "dbpedia":
        file_path = os.path.join(data_path, "out.dbpedia-link")
        skiprows = 2
        sep_char = "\s+"
    else:
        raise ValueError(f"Dataset {dataset} is not supported")
    output_path = os.path.join(data_path, f"{dataset}_undirect.csv")

    df = pd.read_csv(
        file_path,
        names=["src", "dst"],
        sep=sep_char,
        skiprows=skiprows,
        index_col=False,
    )

    # make the graph undirected
    df = make_directed_to_undirected(df)

    df["sign"] = np.where(np.random.rand(df.shape[0]) < pos_ratio, 1, -1)

    write_edges_to_file(df, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20159, help="random seed")
    parser.add_argument(
        "--dataset",
        type=str,
        default="wiki-vote",
        help="dataset: wiki-vote, epinions, wikisigned, youtube, pokec, dbpedia",
    )
    parser.add_argument("--pos_ratio", type=float, default=0.7, help="positive ratio")

    opt = parser.parse_args()

    if opt.seed is not None:
        random.seed(opt.seed)
        np.random.seed(opt.seed)

    data_path = "./data/"
    opt.dataset = opt.dataset.lower()
    if opt.dataset in ["youtube", "pokec", "dbpedia"]:
        parse_scale_data(data_path, opt.dataset, opt.pos_ratio)
    else:
        input_path = os.path.join(data_path, opt.dataset + ".txt")

        if opt.dataset == "wiki-vote":
            input_path = os.path.join(data_path, "wikielections.txt")
        if opt.dataset == "wikisigned":
            input_path = os.path.join(data_path, "wikipolitics.txt")

        output_path = os.path.join(data_path, opt.dataset + "_undirect.csv")
        df = pd.read_csv(
            input_path,
            sep="\t",
            header=None,
            skiprows=1,
            names=["src", "dst", "sign"],
        )

        df = make_directed_to_undirected(df)

        write_edges_to_file(df, output_path)
