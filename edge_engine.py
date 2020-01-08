
ACTION_MODES = ("activation", "binding", "catalysis", "expression", "inhibition", "ptmod", "reaction")

EVIDENCE_COLUMNS = ("neighborhood", "fusion", "cooccurence", "coexpression", "experimental", "database", "textmining")

NON_DIRECTIONAL = ("binding", "expression", "reaction")


class SdblEdgeProperty:
    __slots__ = ["score", "arrowtype", "direction", "penwidth"]

    def __init__(self, score=0, arrowtype="none", direction="forward",
            penwidth=1.0):
        self.score = score
        self.arrowtype = arrowtype
        self.direction = direction
        self.penwidth = penwidth

    def __iadd__(self, other):
        self.score = max(self.score, other.score)

        if self.direction != other.direction:
            self.direction = "both"

        return self


class SdblEdgeEngine:
    NON_DIRECTIONAL = ("binding", "expression", "reaction")

    def __init__(self, db_results):
        self.edges = dict()
        self.db_results = db_results
        weights = [r[-1] for r in db_results]
        self.min_weight = min(weights)
        self.distance = max(weights) - self.min_weight

        if self.distance == 0:
            self.distance = 1

    def edge_property_factory(self, db_result):
        mode = db_result[2]
        score = db_result[-1]
        penwidth = ((float(score) - self.min_weight) / self.distance) + 1

        if mode in ACTION_MODES:
            if mode in NON_DIRECTIONAL:
                direction = "none"
            else:
                k = tuple(sorted(db_result[:2]))
                o = tuple(db_result[:2])

                if db_result[5] == 1:
                    if k != o:
                        direction = "back"
                    else:
                        direction = "forward"
                else:
                    if k != o:
                        direction = "forward"
                    else:
                        direction = "back"

            if mode == "inhibition":
                arrowtype = "tee"
            else:
                arrowtype = "normal"

        else:
            direction = "none"
            arrowtype = "normal"

        return SdblEdgeProperty(score, arrowtype, direction, penwidth)

    def generate_edges(self):
        for res in self.db_results:
            key = tuple(sorted([res[0], res[1]]) + [res[2]])

            if key not in self.edges:
                self.edges[key] = self.edge_property_factory(res)
            else:
                self.edges[key] += self.edge_property_factory(res)

    def __iter__(self):
        if len(self.edges) == 0:
            self.generate_edges()
        return ((k, v) for k, v in self.edges.items())

