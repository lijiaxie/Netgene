import csv
import cPickle
import random
from graph_tool.all import *


def parse_tsv(infile, outfile):
    adjacency = {}
    diseases = {}
    genes = {}
    with open(infile) as f:
        reader = csv.reader(f, delimiter='\t')
        # skip header
        next(reader, None)
        # header = 'geneId geneName description diseaseId diseaseName score NofPmids NofSnps sources'
        for row in reader:
            geneID = row[0]
            geneName = row[1]
            diseaseID = row[3]
            diseaseName = row[4]
            score = row[5]
            if diseaseID not in adjacency:
                adjacency[diseaseID] = []

            adjacency[diseaseID].append((geneID, score))

            if diseaseID not in diseases:
                diseases[diseaseID] = diseaseName

            if geneID not in genes:
                genes[geneID] = geneName

    with open(outfile, 'wb') as f:
        cPickle.dump((adjacency, diseases, genes), f)

    return adjacency, diseases, genes

def build_graph(adjacency):
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

    for diseaseID, genes in adjacency.iteritems():
        if diseaseID not in disease_vert_ind:
            dv = g.add_vertex()
            g.vp.name[dv] = diseaseID
            g.vp.type[dv] = 0
            disease_vert_ind[diseaseID] = g.vertex_index[dv]
        else:
            dv = g.vertex(disease_vert_ind[diseaseID])

        # each gene is (id, score)
        for geneID, score in genes:
            if geneID not in gene_vert_ind:
                gv = g.add_vertex()
                g.vp.name[gv] = geneID
                g.vp.type[gv] = 1
                gene_vert_ind[geneID] = g.vertex_index[gv]
            else:
                gv = g.vertex(gene_vert_ind[geneID])

            e = g.add_edge(dv, gv)
            g.ep.score[e] = float(score)

    return g

def check_graph(in_graph, adjacency, diseases, genes):
    # checks on the total number of edges and vertices
    assert(in_graph.num_edges() == sum([len(v) for k, v in adjacency.iteritems()]))
    assert(in_graph.num_vertices() == len(genes) + len(diseases))

    for i in range(100):
        index = random.randint(1, in_graph.num_vertices())
        vert = in_graph.vertex(index)
        # make sure to select a disease vertex
        while in_graph.vp.type[vert]:
            index = random.randint(1, in_graph.num_vertices())
            vert = in_graph.vertex(index)

        diseaseID = in_graph.vp.name[vert]
        # same number of genes connected
        assert(vert.out_degree() == len(adjacency[diseaseID]))

        for out_v in vert.out_neighbours():
            assert(in_graph.vp.name[in_graph.vertex(out_v)] in zip(*adjacency[diseaseID])[0])



if __name__ == "__main__":
    A, D, G = parse_tsv('data/all_gene_disease_associations.tsv', 'data/disgenet.pickle')
    A, D, G = cPickle.load(open('data/disgenet.pickle'))