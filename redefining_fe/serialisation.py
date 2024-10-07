import json
import redefining_fe
import redefining_fe.cells
import redefining_fe.groups
import sympy as sp


class FETripleEncoder(json.JSONEncoder):

    # seen_objects = {}

    def __init__(self, *args, **kwargs):
        print("Creating new")
        self.seen_objects = {}
        super(FETripleEncoder, self).__init__(*args, **kwargs)

    def default(self, o):
        o_dict = {}

        if isinstance(o, sp.core.containers.Tuple) or isinstance(o, sp.Expr):
            return sp.srepr(o)

        print(o)
        print(o.dict_id())
        if o.dict_id() in self.seen_objects.keys():
            print("Cache hit")
            return o.dict_id()

        if hasattr(o, "_to_dict"):
            o_dict = o._to_dict()
            self.seen_objects[list(o_dict.keys())[0]] = o_dict
            return o_dict

        return super().default(o)

    def obj_id(self, o):
        return str(o.__class__) + str(id(o))
