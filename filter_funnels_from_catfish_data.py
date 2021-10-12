import networkx as nx
from os import mkdir
from os import listdir
import argparse
import re


header_regex = re.compile('# graph number = ([0-9]*) name = (.*)')
edge_regex = re.compile('(\d*) (\d*) (\d*\.\d*)')  # noqa


def enumerate_graphs(graph_file):
    def read_next_graph(f):
        header_line = f.readline()

        if header_line == '':
            return None

        m = header_regex.match(header_line)
        if m is None:
            raise Exception('Misformed graph header line.')

        line = f.readline()

        graph = nx.DiGraph()

        while not line == '':
            last_pos = f.tell()
            line = f.readline()

            if line == '':
                break
            elif line[0] == '#':
                f.seek(last_pos)
                break

            list = line.split()

            u = int(list[0])
            v = int(list[1])

            graph.add_edge(u, v)

        return graph

    with open(graph_file) as f:
        while True:
            graph_data = read_next_graph(f)
            if graph_data is None:
                break
            else:
                yield graph_data


def read_instances(graph_file):
    for graphdata in enumerate_graphs(graph_file):
        yield graphdata


def is_forking_vertex(G, v):
    return G.out_degree(v) > 1


def is_reached_by_merging_vertex(G, v, reached_by_merging_vertex):
    in_neighbors = list(G.predecessors(v))

    if len(in_neighbors) == 1:
        return reached_by_merging_vertex[in_neighbors[0]]

    return len(in_neighbors) > 1


def is_funnel(G):  # G is assumed to be a DAG

    reached_by_merging_vertex = {v: False for v in G.nodes}

    for v in nx.topological_sort(G):
        reached_by_merging_vertex[v] =\
            is_reached_by_merging_vertex(G, v, reached_by_merging_vertex)
        if reached_by_merging_vertex[v] and is_forking_vertex(G, v):
            return False

    return True


def process_file(graph_file):
    """Walk over file, writing instances that are funnels to new graph and
    truth files in the output_dir"""
    # Iterate over every graph-instance inside the input file
    for graph in read_instances(graph_file):
        print(is_funnel(graph))
        print(graph.edges())


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str)
    parser.add_argument("output_dir", type=str)
    args = parser.parse_args()

    try:
        mkdir(args.output_dir)
    except FileExistsError:
        pass

    files = [x for x in listdir(args.input_dir) if x[-5:] == "graph"]
    # there are files 1 through 100, both .graph and .truth
    for graph_file in files[0:1]:
        truth_file = graph_file.split(".")[0] + ".truth"
        print(graph_file)
        print(truth_file)
        process_file(args.input_dir + "/" + graph_file)
