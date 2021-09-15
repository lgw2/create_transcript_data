import networkx as nx
from os import mkdir


input_gtf = open('human.gtf', 'r').read().split('\n')
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
        components.append(component_dict)

    return components


# Write in a single file all components in gfa format
def store_components_to_gfa(filename, components):
    output_gfa = open(filename, 'w')

    graph_number = 0
    for component in components:
        graph = component['graph']
        transcripts = component['transcripts']
        # Only store graph with at least one transcript of length > 1
        if len(list(filter(lambda transcript:
                           len(transcript['pseudo_exons']) > 1,
                           transcripts))) > 0:
            name = ','.join(list(map(lambda transcript:
                                     transcript['transcript_id'],
                                     transcripts)))
            output_gfa.write(
                f'H # graph number = {graph_number} name = {name}\n')

            for edge in graph.edges:
                cov = graph.edges[edge]['cov']
                output_gfa.write(
                    f'L\t({edge[0][0]},{edge[0][1]})\t+\t({edge[1][0]},{edge[1][1]})\t+\t{cov}\n') # noqa

            graph_number += 1
    output_gfa.close()


# Write in a single file all truth graphs in catfish's format
def store_transcripts_to_truth_file(filename, components):
    output_truth_file = open(filename, 'w')

    graph_number = 0
    for component in components:
        transcripts = component['transcripts']
        # Only store graph with at least one transcript of length > 1
        if len(list(filter(lambda transcript: len(transcript['pseudo_exons'])
                           > 1, transcripts))) > 0:
            name = ','.join(list(map(lambda transcript:
                                     transcript['transcript_id'],
                                     transcripts)))
            output_truth_file.write(
                f'# graph number = {graph_number} name = {name}\n')

            for transcript in transcripts:
                cov = transcript['cov']
                output_truth_file.write(f'{cov} {" ".join(list(map(lambda p_exon: f"({p_exon[0]},{p_exon[1]})", transcript["pseudo_exons"])))}\n') # noqa

            graph_number += 1
    output_truth_file.close()


try:
    mkdir('./human')
except FileExistsError:
    pass

for sequence in data:
    components = build_splicing_graphs(data[sequence])
    store_components_to_gfa(f'./human/{sequence}.gfa', components)
    store_transcripts_to_truth_file(f'./human/{sequence}.truth', components)
