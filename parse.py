import csv
import cPickle
import random
from graph_tool.all import *


def parse_tsv(infile, outfile):
    adj = {}
    dis_dict = {}
    gene_dict = {}
    with open(infile) as f:
        reader = csv.reader(f, delimiter='\t')
        # skip header
        next(reader, None)
        # header = 'geneId geneName description diseaseId diseaseName score NofPmids NofSnps sources'
        for row in reader:
            gene_id = row[0]
            gene_name = row[1]
            disease_id = row[3]
            disease_name = row[4]
            score = row[5]
            if disease_id not in adj:
                adj[disease_id] = []

            adj[disease_id].append((gene_id, score))

            if disease_id not in dis_dict:
                dis_dict[disease_id] = disease_name

            if gene_id not in gene_dict:
                gene_dict[gene_id] = gene_name

    with open(outfile, 'wb') as f:
        cPickle.dump((adj, dis_dict, gene_dict), f)

    return adj, dis_dict, gene_dict


def build_graph(input_adj):
    # init network
    g = Graph(directed=False)
    disease_vert_ind = {}
    gene_vert_ind = {}

    # edge and vertex properties
    v_name = g.new_vertex_property("string")
    # 0 for disease, 1 for gene
    v_type = g.new_vertex_property("int")
    g.vertex_properties["name"] = v_name
    g.vertex_properties["type"] = v_type
    # scores from disgenet
    e_score = g.new_edge_property("float")
    g.edge_properties["score"] = e_score

    for disease_id, genes in input_adj.iteritems():
        if disease_id not in disease_vert_ind:
            dv = g.add_vertex()
            g.vp.name[dv] = disease_id
            g.vp.type[dv] = 0
            disease_vert_ind[disease_id] = g.vertex_index[dv]
        else:
            dv = g.vertex(disease_vert_ind[disease_id])

        # each gene is (id, score)
        for gene_id, score in genes:
            if gene_id not in gene_vert_ind:
                gv = g.add_vertex()
                g.vp.name[gv] = gene_id
                g.vp.type[gv] = 1
                gene_vert_ind[gene_id] = g.vertex_index[gv]
            else:
                gv = g.vertex(gene_vert_ind[gene_id])

            new_edge = g.add_edge(dv, gv)
            g.ep.score[new_edge] = float(score)

    return g


def check_graph(in_graph, input_adj, disease_dict, gene_dict):
    # checks on the total number of edges and vertices
    assert(in_graph.num_edges() == sum([len(v) for k, v in input_adj.iteritems()]))
    assert(in_graph.num_vertices() == len(gene_dict) + len(disease_dict))

    for i in range(100):
        index = random.randint(1, in_graph.num_vertices())
        vert = in_graph.vertex(index)
        # make sure to select a disease vertex
        while in_graph.vp.type[vert]:
            index = random.randint(1, in_graph.num_vertices())
            vert = in_graph.vertex(index)

        disease_id = in_graph.vp.name[vert]
        # same number of genes connected
        assert(vert.out_degree() == len(input_adj[disease_id]))

        for out_v in vert.out_neighbours():
            assert(in_graph.vp.name[in_graph.vertex(out_v)] in zip(*input_adj[disease_id])[0])


def get_connected(in_graph):
    lcc = label_largest_component(in_graph)
    u = GraphView(in_graph, vfilt=lcc)  # extract the largest component as a graph
    return u


if __name__ == "__main__":
    A, D, G = parse_tsv('data/all_gene_disease_associations.tsv', 'data/disgenet.pickle')
    A, D, G = cPickle.load(open('data/disgenet.pickle'))