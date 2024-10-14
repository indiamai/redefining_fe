import json
# from redefining_fe import *
from redefining_fe.cells import Point
from redefining_fe.groups import GroupRepresentation
from redefining_fe.triples import ElementTriple
import sympy as sp
from sympy.combinatorics import Permutation
import re


class FETripleEncoder(json.JSONEncoder):
    # need to think more about how to handle sympy/perm groups

    def __init__(self, *args, **kwargs):
        self.seen_objects = {}
        super(FETripleEncoder, self).__init__(*args, **kwargs)

    def default(self, o):
        o_dict = {}
        print(o)
        print(self.seen_objects.keys())
        if isinstance(o, sp.core.containers.Tuple) or isinstance(o, sp.Expr):
            return sp.srepr(o)
        # if isinstance(o, PermutationGroup):
        #     return [g for g in o.elements]
        if isinstance(o, Permutation):
            return o.array_form

        if o.dict_id() in self.seen_objects.keys():
            return o.dict_id()

        if hasattr(o, "_to_dict"):
            o_dict = o._to_dict()
            if isinstance(o_dict, dict):
                self.seen_objects[list(o_dict.keys())[0]] = o_dict
            return o_dict

        return super().default(o)

    def encode(self, o):
        # return super(FETripleEncoder, self).encode(o)
        dict_str = "{ \"encoded_obj\":" + super(FETripleEncoder, self).encode(o) + ", \"encoded_obj_dict\": {"
        # return dict_str
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
                    s2, e2 = bracket_matching(dict_str[start_index:],)
                    found_dict = dict_str[start_index+s2:start_index+e2]
                    dict_str = dict_str[:start_index+s2] + seen_obj + dict_str[start_index+e2:] + found_dict + ","
        if dict_str[-1] == ",":
            dict_str = dict_str[:-1]
        new_dict_str = dict_str + "}}"
        return new_dict_str

    def obj_id(self, o):
        return str(o.__class__) + str(id(o))


class FETripleDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, dct):
        print(type(dct))
        print(dct.keys())
        for key in dct.keys():
            if key.startswith("Triple"):
                print("Triple")
                print(key)
                print(dct[key].keys())
                return ElementTriple._from_dict(dct[key])
            elif key.startswith("PointKernel"):
                pass
            elif key.startswith("Point"):
                return Point._from_dict(dct[key])
            elif key.startswith("GroupRep"):
                return GroupRepresentation._from_dict(dct[key])
        return dct

    def decode(self, o_str):
        # print(o_str)
        print(o_str.find("encoded_obj"))
        obj_base = o_str.find("encoded_obj")
        obj_s, obj_e = bracket_matching(o_str[obj_base:])
        dict_base = o_str.find("encoded_obj")
        dict_s, dict_e = bracket_matching(o_str[dict_base:])
        obj = o_str[obj_base + obj_s: obj_base + obj_e]
        obj_dict = o_str[dict_base + dict_s: dict_base + dict_e]
        print(obj)
        print(obj_dict)


def bracket_matching(dict_str):
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
