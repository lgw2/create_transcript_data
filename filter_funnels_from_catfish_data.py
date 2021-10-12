import networkx as nx
from os import mkdir
from os import listdir
import argparse
import re


header_regex = re.compile('# graph number = ([0-9]*) name = (.*)')
edge_regex = re.compile('(\d*) (\d*) (\d*\.\d*)')  # noqa


def enumerate_graphs(graph_file):
    def read_next_graph(f):
        file_lines = []
        header_line = f.readline()
        file_lines.append(header_line)

        if header_line == '':
            return None

        m = header_regex.match(header_line)
        if m is None:
            raise Exception('Misformed graph header line.')
        (graph_number, graph_name) = (m.group(1), m.group(2))

        line = f.readline()
        file_lines.append(line)

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
            file_lines.append(line)

        return graph, graph_name, graph_number, file_lines

    with open(graph_file) as f:
        while True:
            graph_data = read_next_graph(f)
            if graph_data is None:
                break
            else:
                yield graph_data


def enumerate_decompositions(decomposition_file):
    def read_next_decomposition(f):
        header_line = f.readline()
        file_lines = []
        file_lines.append(header_line)

        if header_line == '':
            return None

        m = header_regex.match(header_line)
        if m is None:
            raise Exception('Misformed graph header line.')
        (graph_number, graph_name) = (m.group(1), m.group(2))

        path_decomposition = []
        line = header_line
        while not line == '':
            last_pos = f.tell()
            line = f.readline()

            if line == '':
                break
            elif line[0] == '#':
                f.seek(last_pos)
                break

            l = line.split()  # noqa
            l = list(map(lambda x: int(x), l))  # noqa

            path_decomposition.append((l[0], l[1:]))
            file_lines.append(line)

        return (graph_name, graph_number, path_decomposition, file_lines)

    with open(decomposition_file) as f:
        while True:
            decomposition = read_next_decomposition(f)
            if decomposition is None:
                break
            else:
                yield decomposition


def read_instances(graph_file, truth_file):
    for graphdata, truthdata in zip(enumerate_graphs(graph_file),
                                    enumerate_decompositions(
                                    truth_file)):
        yield (graphdata, truthdata)


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


def process_file(in_graph_file, in_truth_file, out_graph_file, out_truth_file):
    """Walk over file, writing instances that are funnels to new graph and
    truth files in the output_dir"""
    # Iterate over every graph-instance inside the input file
    graph_out = open(out_graph_file, "w")
    truth_out = open(out_truth_file, "w")
    for graph, truth in read_instances(in_graph_file, in_truth_file):
        if not is_funnel(graph[0]):
            for line in graph[3]:
                graph_out.write(line)
            for line in truth[3]:
                truth_out.write(line)
    graph_out.close()
    truth_out.close()


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
    for graph_file in files:
        truth_file = graph_file.split(".")[0] + ".truth"
        print(f"Processing {graph_file}")
        process_file(args.input_dir + "/" + graph_file,
                     args.input_dir + "/" + truth_file,
                     args.output_dir + "/" + graph_file,
                     args.output_dir + "/" + truth_file)
