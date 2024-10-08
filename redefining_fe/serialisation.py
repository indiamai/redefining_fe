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

        if o.dict_id() in self.seen_objects.keys():
            return o.dict_id()

        if hasattr(o, "_to_dict"):
            o_dict = o._to_dict()
            if isinstance(o_dict, dict):
                self.seen_objects[list(o_dict.keys())[0]] = o_dict
            return o_dict

        return super().default(o)

    def encode(self, o):
        dict_str = super(FETripleEncoder, self).encode(o) + ", {"
        for seen_obj in self.seen_objects.keys():
            matches = [(m.start(), m.end()) for m in re.finditer(seen_obj, dict_str)]
            if len(matches) > 1:
                remove_s = None
                for (s, e) in matches:
                    if dict_str[e+1] == ":":
                        # this will only occur once
                        remove_s = s
                if remove_s:
                    start_index = dict_str[:remove_s].rindex("{")
                    s2, e2 = self.bracket_matching(dict_str[start_index:],)
                    found_dict = dict_str[start_index+s2:start_index+e2]
                    dict_str = dict_str[:start_index+s2] + seen_obj + dict_str[start_index+e2:] + found_dict + ","
        new_dict_str = "{" + dict_str + "}}"
        return new_dict_str

    def bracket_matching(self, dict_str):
        bracket_counter = 0
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
