import json
# import redefining_fe.cells
# import redefining_fe.groups
import sympy as sp
import re


class FETripleEncoder(json.JSONEncoder):

    def __init__(self, *args, **kwargs):
        self.seen_objects = {}
        super(FETripleEncoder, self).__init__(*args, **kwargs)

    def default(self, o):
        o_dict = {}

        if isinstance(o, sp.core.containers.Tuple) or isinstance(o, sp.Expr):
            return sp.srepr(o)

        # print(o)
        # print(o.dict_id())
        if o.dict_id() in self.seen_objects.keys():
            # print("Cache hit")
            return o.dict_id()

        if hasattr(o, "_to_dict"):
            o_dict = o._to_dict()
            if isinstance(o_dict, dict):
                self.seen_objects[list(o_dict.keys())[0]] = o_dict
            return o_dict

        return super().default(o)

    def encode(self, o):
        dict_str = super(FETripleEncoder, self).encode(o)
        print("encoded")
        for seen_obj in self.seen_objects.keys():
            matches = [(m.start(), m.end()) for m in re.finditer(seen_obj, dict_str)]
            if len(matches) > 1:
                print(seen_obj)
                # need to replace whole dict found
                for (s, e) in matches:
                    if dict_str[e+1] == ":":
                        s2, e2 = self.bracket_matching(dict_str[e:])
                        found_dict = dict_str[e+s2:e+e2+1]
                        dict_str = dict_str[:e+s2] + seen_obj + dict_str[e+e2:] + ", {" + seen_obj + ":" + found_dict + "}"
        dict_str = "{" + dict_str + "}"
        return dict_str

    def bracket_matching(self, dict_str):
        bracket_counter = 0
        print(dict_str)
        print("BRACK")
        for m in re.finditer("{|}", dict_str):
            if m.group() == "{":
                if bracket_counter == 0:
                    start_val = m.start()
                bracket_counter += 1
            elif m.group() == "}":
                bracket_counter -= 1
                if bracket_counter == 0:
                    return start_val, m.end()
        raise ValueError("End of parsing reached without matching brackets")

    def obj_id(self, o):
        return str(o.__class__) + str(id(o))
