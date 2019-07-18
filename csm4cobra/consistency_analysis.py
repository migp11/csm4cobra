import networkx as nx
import pandas as pd
from cobra.flux_analysis import find_blocked_reactions


class UmFinder:
    def __init__(self, cobra_model, cc_method='fva', report=True):
        self._model = cobra_model

        if report:
            print("===========================================================")
            print("Initializing UmFinder Builder using")
            print("Model: %s" % cobra_model.id)

            print("- Nº of reactions: %i" % len(self._model.reactions))
            print("- Nº of metabolites: %i" % len(self._model.metabolites))

            print("\nChecking network consistency (may take some minutes)")
            print("Finding blocked reaction using method: %s\n" % cc_method)

        self._blocked_reactions = find_blocked_reactions(self.model)
        self._gap_metabolites = UmFinder.find_gap_metabolites(self.model, self.blocked_reactions)
        self._gap_graph = UmFinder.create_gap_graph(self.model, self._gap_metabolites, self._blocked_reactions)

        unconnected_modules = nx.connected_components(self._gap_graph.to_undirected())
        self._unconnected_modules = sorted(unconnected_modules, key=lambda x: len(x), reverse=True)

        if report:
            print("- Nº of blocked reactions: %i" % len(self._blocked_reactions))
            print("- Nº of gap metabolites: %i" % len(self._gap_metabolites))
            print("- Nº of unconnected modules: %i" % len(self.unconnected_modules))

            if len(self.unconnected_modules):
                df_ums = self.unconnected_modules_frame()
                df_biggest_um = df_ums.node_type[df_ums.um_id == 1]
                rxns = df_biggest_um.index[df_biggest_um =='rxn']
                mets = df_biggest_um.index[df_biggest_um == 'met']
                print("- N of reactions in the biggest unconnected module: %i" % len(rxns))
                print("- N of metabolites in the biggest unconnected module: %i" % len(mets))


    @property
    def model(self):
        return self._model

    @property
    def gap_metabolites(self):
        return frozenset(self._gap_metabolites)

    @property
    def gap_graph(self):
        return self._gap_graph

    @property
    def blocked_reactions(self):
        return frozenset(self._blocked_reactions)

    @property
    def unconnected_modules(self):
        return self._unconnected_modules

    def update(self):
        self._blocked_reactions = find_blocked_reactions(self.model)
        self._gap_metabolites = UmFinder.find_gap_metabolites(self.model, self.blocked_reactions)
        self._gap_graph = UmFinder.create_gap_graph(self.model, self._gap_metabolites, self._blocked_reactions)
        self._unconnected_modules = nx.connected_component_subgraphs(self._gap_graph.to_undirected())

        self._unconnected_modules = sorted(self._unconnected_modules, key=lambda x: len(x), reverse=True)

    def unconnected_module_subgraphs(self):
        for um in self.unconnected_modules:
            yield self.gap_graph.subgraph(um)

    def unconnected_modules_frame(self):
        columns = ['node_id', 'node_type', 'um_id']
        data = {}
        counter = 0
        for i, um in enumerate(self.unconnected_modules):
            for e in um:
                if e in self.gap_metabolites:
                    e_type = 'met'
                elif e in self.blocked_reactions:
                    e_type = 'rxn'
                else:
                    e_type = None
                data[counter] = (e, e_type, i+1)
                counter += 1

        return pd.DataFrame.from_dict(data, orient='index', columns=columns)

    @staticmethod
    def find_gap_metabolites(model, blocked_reactions):
        gap_metabolites = []
        for m in model.metabolites:
            reactions = set([r.id for r in m.reactions])
            if reactions.issubset(blocked_reactions):
                gap_metabolites.append(m.id)

        return gap_metabolites

    @staticmethod
    def create_metabolic_graph(cobra_model, directed=True, reactions=None, rev_rxn_label='reversible'):

        graph = nx.DiGraph()
        if not directed:
            graph = nx.Graph()

        if not reactions:
            reactions = cobra_model.reactions

        if not hasattr(reactions[0], 'id'):
            reactions = [cobra_model.reactions.get_by_id(r) for r in reactions]

        for r in reactions:
            graph.add_node(r.id, label=r.id, text=r.id, node_class="rxn", node_id=r.id)
            for m in r.metabolites:
                if m.id not in graph.nodes():
                    graph.add_node(m.id, label=m.id, text=m.id, node_class="met", node_id=m.id)

                (tail, head) = (r.id, m.id)
                if r.get_coefficient(m) < 0:
                    (tail, head) = (m.id, r.id)

                graph.add_edge(tail, head)
                graph[tail][head][rev_rxn_label] = r.lower_bound < 0

        return graph

    @staticmethod
    def create_gap_graph(model, gap_metabolites, blocked_reactions):

        if hasattr(gap_metabolites[0], 'id'):
            gap_metabolites = [m.id for m in gap_metabolites]

        if hasattr(blocked_reactions[0], 'id'):
            blocked_reactions = [r.id for r in blocked_reactions]

        graph = UmFinder.create_metabolic_graph(model)
        gap_graph = graph.subgraph(gap_metabolites + blocked_reactions)

        return gap_graph




def main():
    model = cobra.test.create_test_model('ecoli')
    um_finder = UmFinder(model)

