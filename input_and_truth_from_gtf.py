import networkx as nx
from os import mkdir
import argparse
import numpy as np
from collections import defaultdict


SOURCE = (0, 0)
SINK = (-1, -1)


parser = argparse.ArgumentParser()
parser.add_argument("input_filename", type=str)
parser.add_argument('--simulate_cov', dest='simulate_cov',
                    default=False, action='store_true')
parser.add_argument('--catfish_format', dest='catfish_format',
                    default=False, action='store_true')
parser.add_argument('--seed', type=int, required=False)
parser.add_argument('--filter-funnels', dest='filter_funnels',
                    default=False, action='store_true')
args = parser.parse_args()
filename = args.input_filename
# assume input filename is a .gtf file
dir_to_write = filename.split('.gtf')[0]
input_gtf = open(filename, 'r').read().split('\n')
input_gtf = list(filter(lambda line: not line.startswith('#') and
                        line, input_gtf))


def process_attribute(att_item):
    att_item = att_item.strip()
    return (att_item.split(' ')[0], att_item.split('"')[-2])


def process_attributes(att_line):
    return {key: value for (key, value) in
            map(lambda att_item: process_attribute(att_item),
                att_line.split(';')[:-1])}


def process_line(line):
    line = line.split('\t')
    attributes = process_attributes(line[-1])

    return {
        'seqname': line[0],
        'feature': line[2],
        'start': int(line[3]),
        'end': int(line[4]),
        'strand': line[6],
        'transcript_id':
        attributes.get('transcript_id', None) if
        'transcript_id' in attributes else None,
        # One represents the unweighted case
        'cov': float(attributes.get('cov', None)) if 'cov' in attributes else 1
    }


data = list(map(lambda line: process_line(line), input_gtf))

# This notebook only works with the data on the positive strand and
# only exons and transcripts
data = list(filter(lambda datum: datum['strand'] == '+' and
                   (datum['feature'] == 'transcript' or
                    datum['feature'] == 'exon'), data))


# Group by seqname, and then by transcript_id
def group_data(data):
    sequences = dict()

    for datum in data:

        if datum['seqname'] not in sequences:
            sequences[datum['seqname']] = dict()
        sequence = sequences[datum['seqname']]

        if datum['transcript_id'] not in sequence:
            sequence[datum['transcript_id']] = {
                'exons': list()
            }

        transcript = sequence[datum['transcript_id']]

        if datum['feature'] == 'transcript':
            transcript['cov'] = datum['cov']
        elif datum['feature'] == 'exon':
            transcript['exons'].append((datum['start'], datum['end']))

    # Sort exons
    for sequence in sequences.values():
        for transcript in sequence.values():
            transcript['exons'].sort(key=lambda exon: exon[0])

    return sequences


data = group_data(data)

# at this point, let's add cov values.
if args.seed:
    np.random.seed(args.seed)
if args.simulate_cov:
    for chromosome in data:
        for transcript in data[chromosome]:
            data[chromosome][transcript]['cov'] =\
                1000 * np.random.lognormal(-4, 4)

# if coverages are non-integer, round
for chromosome in data:
    for transcript in data[chromosome]:
        data[chromosome][transcript]['cov'] =\
            round(data[chromosome][transcript]['cov'])


def add_pseudo_exons(sequence):
    # Extract exon_endpoints and sorts them according to its position and in
    # case of ties it puts first the starting positions
    exons_endpoints = list()
    for transcript_id in sequence:
        transcript = sequence[transcript_id]
        transcript_exons_endpoints = list()
        for i, exon in enumerate(transcript['exons']):
            transcript_exons_endpoints.append({'position': exon[0],
                                               'transcript_id': transcript_id,
                                               'exon_index': i,
                                               'start_point': True})
            transcript_exons_endpoints.append({'position': exon[1],
                                               'transcript_id': transcript_id,
                                               'exon_index': i,
                                               'start_point': False})
        exons_endpoints.append(transcript_exons_endpoints)

    exon_endpoint_pos_list = [item for sublist in exons_endpoints for item in
                              sublist]
    exon_endpoint_pos_list.sort(key=lambda x: [x['position'],
                                               not(x['start_point'])])

    # Mark the exons when overlapping with another
    active_exons = dict()
    for exon_endpoint in exon_endpoint_pos_list:
        if exon_endpoint['start_point']:
            exon_endpoint['starting_points'] = list()
            active_exons[(exon_endpoint['transcript_id'],
                          exon_endpoint['exon_index'])] = exon_endpoint

            for key in active_exons:
                active_exons[key]['starting_points'].append(
                    exon_endpoint['position'])
        else:
            for key in active_exons:
                active_exons[key]['starting_points'].append(
                    exon_endpoint['position']+1)
            del active_exons[(exon_endpoint['transcript_id'],
                              exon_endpoint['exon_index'])]

    # Compute pseudo-exons
    for exon_endpoint in exon_endpoint_pos_list:
        if exon_endpoint['start_point']:
            exon_endpoint['pseudo_exons'] = list()
            previous_value = exon_endpoint['starting_points'][0]
            for i in range(1, len(exon_endpoint['starting_points'])):
                if previous_value != exon_endpoint['starting_points'][i]:
                    exon_endpoint['pseudo_exons'].append(
                        (previous_value,
                         exon_endpoint['starting_points'][i]-1))
                    previous_value = exon_endpoint['starting_points'][i]

    # Add them to the sequence
    for transcript in sequence.values():
        transcript['pseudo_exons'] = list()

    for exon_endpoint in exon_endpoint_pos_list:
        if exon_endpoint['start_point']:
            transcript = sequence[exon_endpoint['transcript_id']]
            for pseudo_exon in exon_endpoint['pseudo_exons']:
                transcript['pseudo_exons'].append(pseudo_exon)

    for transcript in sequence.values():
        transcript['pseudo_exons'].sort(key=lambda p_exon: p_exon[0])

    return sequence


for sequence in data:
    add_pseudo_exons(data[sequence])


def build_splicing_graphs(sequence):

    # Build vertex and edge sets
    vertices = set([item for sublist in
                    map(lambda transcript: transcript['pseudo_exons'],
                        sequence.values()) for item in sublist])
    edges = dict()
    transcripts_starting_at = dict()
    for vertex in vertices:
        transcripts_starting_at[vertex] = list()

    for (transcript_id, transcript) in sequence.items():
        # Consecutive pseudo exons in pseudo_exons are linked by an edge
        cov = transcript['cov']
        pseudo_exons = transcript['pseudo_exons']
        transcripts_starting_at[pseudo_exons[0]].append({
            'transcript_id': transcript_id,
            'pseudo_exons': pseudo_exons,
            'cov': cov
        })

        for i in range(len(pseudo_exons)-1):
            current_p_exon = pseudo_exons[i]
            next_p_exon = pseudo_exons[i+1]
            edge = (current_p_exon, next_p_exon)
            if edge not in edges:
                edges[edge] = (current_p_exon, next_p_exon, {'cov': cov})
            else:
                edges[edge] = (edges[edge][0], edges[edge][1],
                               {'cov': edges[edge][2]['cov'] + cov})

    # Build graph to find weakly connected components
    G = nx.DiGraph(edges.values())
    G.update(nodes=vertices)

    components = list()
    for component_v in nx.weakly_connected_components(G):
        component_dict = {'graph': G.subgraph(component_v)}
        transcripts_component = list()

        for vertex in component_v:
            transcripts_component += transcripts_starting_at[vertex]

        component_dict['transcripts'] = transcripts_component

        # add an edge from source to first pseudo exon and from last pseudo
        # exon to sink for this component
        # component_dict['graph'].add_edge(SOURCE,
        components.append(component_dict)

    return components


def check_flow(graph, transcripts):
    sources_and_sinks = set()
    for transcript in transcripts:
        source = transcript["pseudo_exons"][0]
        sink = transcript["pseudo_exons"][-1]
        sources_and_sinks.add(source)
        sources_and_sinks.add(sink)
    weights = nx.get_edge_attributes(graph, "cov")
    for node in graph.nodes():
        if node not in sources_and_sinks:
            out_edges = graph.out_edges(node)
            out_weight = 0
            for edge in out_edges:
                out_weight += weights[edge]
            in_edges = graph.in_edges(node)
            in_weight = 0
            for edge in in_edges:
                in_weight += weights[edge]
            assert out_weight == in_weight


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


# Write in a single file all components in sg format
def store_components_to_sg(filename, components, filter_funnels):
    output_sg = open(filename, 'w')

    graph_number = 0
    for component in components:
        graph = component['graph']
        transcripts = component['transcripts']
        check_flow(graph, transcripts)
        # Only store graph with at least one transcript of length > 1 whose
        # coverage is also greater than 0
        num_transcripts = len(list(filter(lambda transcript:
                                          len(transcript['pseudo_exons']) > 1,
                                          transcripts)))
        total_cov = sum([x["cov"] for x in transcripts])
        if num_transcripts > 0 and total_cov > 0:
            if filter_funnels:
                # check whether this graph is a funnel
                # TODO: do this based on coverage values, not just all of the
                # transcripts
                # build temporary graph using only cov>0 edges to test funnel
                new_graph = nx.DiGraph()
                for edge in graph.edges():
                    cov = graph.edges[edge]['cov']
                    if cov > 0:
                        new_graph.add_edge(edge[0], edge[1])
                if is_funnel(new_graph):
                    continue  # skip this instance
            name = ','.join(list(map(lambda transcript:
                                     transcript['transcript_id'],
                                     transcripts)))
            if args.catfish_format:
                output_sg.write(
                    f'# graph number = {graph_number} name = {name}\n')
            else:
                output_sg.write(
                    f'H # graph number = {graph_number} name = {name}\n')

            if args.catfish_format:
                topo_sort = list(nx.topological_sort(graph))
                sink = len(list(topo_sort)) + 1
                output_sg.write(f'{sink + 1}\n')
                mapping = dict(zip(list(topo_sort),
                                   range(1, sink)))
                graph = nx.relabel_nodes(graph, mapping, copy=True)

            for edge in graph.edges:
                cov = graph.edges[edge]['cov']
                if cov > 0:
                    if args.catfish_format:
                        output_sg.write(
                            f'{edge[0]} {edge[1]} {cov}\n') # noqa
                    else:
                        output_sg.write(
                            f'L\t({edge[0][0]},{edge[0][1]})\t+\t({edge[1][0]},{edge[1][1]})\t+\t{cov}\n') # noqa
            # add edges from source to start pseudo exons and end pseudo
            # exons to sink
            new_edges = defaultdict(int)
            for transcript in component['transcripts']:
                # print(f"{transcript['pseudo_exons'][0]}...{transcript['pseudo_exons'][-1]}, # {transcript['cov']}") # noqa
                if args.catfish_format:
                    new_edges[(0, mapping[transcript['pseudo_exons'][0]])] +=\
                        transcript['cov']
                    new_edges[(mapping[transcript['pseudo_exons'][-1]],
                               sink)] +=\
                        transcript['cov']
                else:
                    new_edges[(SOURCE, transcript['pseudo_exons'][0])] +=\
                        transcript['cov']
                    new_edges[(transcript['pseudo_exons'][-1], SINK)] +=\
                        transcript['cov']
            for edge in new_edges:
                u = edge[0]
                v = edge[1]
                cov = new_edges[edge]
                if cov > 0:
                    if args.catfish_format:
                        output_sg.write(f'{u} {v} {cov}\n') # noqa
                    else:
                        output_sg.write(f'L\t({u[0]},{u[1]})\t+\t({v[0]},{v[1]})\t+\t{cov}\n') # noqa

            graph_number += 1
    output_sg.close()


# Write in a single file all truth graphs in catfish's format
def store_transcripts_to_truth_file(filename, components, filter_funnels):
    output_truth_file = open(filename, 'w')

    graph_number = 0
    for component in components:
        transcripts = component['transcripts']
        graph = component['graph']
        # Only store graph with at least one transcript of length > 1 whose
        # coverage is also greater than 0
        num_transcripts = len(list(filter(lambda transcript:
                                          len(transcript['pseudo_exons']) > 1,
                                          transcripts)))
        total_cov = sum([x["cov"] for x in transcripts])
        if num_transcripts > 0 and total_cov > 0:
            if filter_funnels:
                new_graph = nx.DiGraph()
                for edge in graph.edges():
                    cov = graph.edges[edge]['cov']
                    if cov > 0:
                        new_graph.add_edge(edge[0], edge[1])
                if is_funnel(new_graph):
                    continue  # skip this instance
                # check whether this graph is a funnel
                if is_funnel(new_graph):
                    continue  # skip this instance
            name = ','.join(list(map(lambda transcript:
                                     transcript['transcript_id'],
                                     transcripts)))
            output_truth_file.write(
                f'# graph number = {graph_number} name = {name}\n')

            if args.catfish_format:
                topo_sort = list(nx.topological_sort(graph))
                sink = len(list(topo_sort)) + 1
                mapping = dict(zip(list(topo_sort),
                                   range(1, sink)))
                graph = nx.relabel_nodes(graph, mapping, copy=True)

            for transcript in transcripts:
                # make a temp version of pseudo_exons with source and sink
                temp_pseudo_exons = transcript["pseudo_exons"]
                temp_pseudo_exons.insert(0, SOURCE)
                temp_pseudo_exons.append(SINK)
                cov = transcript['cov']
                # don't write weight 0 paths
                if cov > 0:
                    if args.catfish_format:
                        mid_edges = " ".join(
                            list(map(lambda p_exon:
                                     f"{mapping[p_exon]}",
                                     temp_pseudo_exons[1:-1])))
                        output_truth_file.write(f'{cov} {0} {mid_edges} {sink}\n') # noqa
                    else:
                        output_truth_file.write(f'{cov} {" ".join(list(map(lambda p_exon: f"({p_exon[0]},{p_exon[1]})", temp_pseudo_exons)))}\n') # noqa

            graph_number += 1
    output_truth_file.close()


try:
    mkdir(dir_to_write)
except FileExistsError:
    pass

# "sequence" actually refers to the chromosome I believe
for sequence in data:
    components = build_splicing_graphs(data[sequence])
    if args.catfish_format:
        file_extension = 'graph'
    else:
        file_extension = 'sg'
    store_components_to_sg(f'./{dir_to_write}/{sequence}.{file_extension}',
                           components, args.filter_funnels)
    store_transcripts_to_truth_file(f'./{dir_to_write}/{sequence}.truth',
                                    components, args.filter_funnels)
