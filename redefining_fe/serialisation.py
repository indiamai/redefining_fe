import json
# from redefining_fe import *
from redefining_fe.cells import Point, Edge
from redefining_fe.groups import GroupRepresentation, GroupMemberRep
from redefining_fe.triples import ElementTriple
import sympy as sp
from sympy.combinatorics import Permutation
import re


class FETripleEncoder(json.JSONEncoder):

    def __init__(self, *args, **kwargs):
        self.seen_objects = {}
        super(FETripleEncoder, self).__init__(*args, **kwargs)

    def default(self, o):
        o_dict = {}
        if isinstance(o, sp.core.containers.Tuple) or isinstance(o, sp.Expr):
            return {"sympy": sp.srepr(o)}
        if isinstance(o, Permutation):
            return o.array_form

        if o.dict_id() in self.seen_objects.keys():
            return {"stored": o.dict_id()}

        if hasattr(o, "_to_dict"):
            o_dict = o._to_dict()
            if isinstance(o_dict, dict):
                self.seen_objects[list(o_dict.keys())[0]] = o_dict
            return o_dict

        return super().default(o)

    def encode(self, o):
        # return super(FETripleEncoder, self).encode(o)
        dict_str = "{ \"encoded_obj\":" + super(FETripleEncoder, self).encode(o) + ", \"encoded_obj_dict\": {"
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
                    s2, e2 = bracket_matching(dict_str[start_index:], with_brackets=False)
                    found_dict = dict_str[start_index+s2:start_index+e2]
                    dict_str = dict_str[:start_index+s2] + "\"stored\": \"" + seen_obj + "\"" + dict_str[start_index+e2:] + found_dict + ","
        if dict_str[-1] == ",":
            dict_str = dict_str[:-1]
        new_dict_str = dict_str + "}}"
        return new_dict_str

    def obj_id(self, o):
        return str(o.__class__) + str(id(o))


class FETripleDecoder(json.JSONDecoder):
    def __init__(self, *args, **kwargs):
        self.decoded_obj_dict = {}
        self.store_res = False
        super().__init__(object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, dct):
        res = dct
        for key in dct.keys():
            if key.startswith("Triple"):
                print("Triple")
                res = ElementTriple._from_dict(dct[key])
            elif key.startswith("PointKernel"):
                raise NotImplementedError
            elif key.startswith("Point"):
                res = Point._from_dict(dct[key])
            elif key.startswith("Edge"):
                res = Edge._from_dict(dct[key])
            elif key.startswith("GroupRep"):
                res = GroupRepresentation._from_dict(dct[key])
            elif key.startswith("GroupMemberRep"):
                res = GroupMemberRep._from_dict(dct[key])
            elif key == "sympy":
                res = sp.parse_expr(dct[key])
            elif key == "stored":
                res = self.decoded_obj_dict[dct[key]]
            if self.store_res:
                # need to think about order of serialization/how to handle this!
                self.decoded_obj_dict[key] = res
        return res

    def decode(self, o_str):
        # extract encoded object
        obj_base = o_str.find("encoded_obj")
        obj_s, obj_e = bracket_matching(o_str[obj_base:])

        # extract reference dictionary
        dict_base = o_str.find("encoded_obj_dict")
        dict_s, dict_e = bracket_matching(o_str[dict_base:])
        obj = o_str[obj_base + obj_s: obj_base + obj_e]
        obj_dict = o_str[dict_base + dict_s: dict_base + dict_e]
        # print("obj", obj)
        print("obj_dict", obj_dict)
        print(re.match("'", obj_dict))
        if obj_dict != "":
            self.store_res = True
            print("obj_dict_decoded", super(FETripleDecoder, self).decode(obj_dict))
            self.store_res = False
        return super(FETripleDecoder, self).decode(obj)


def bracket_matching(dict_str, with_brackets=True):
    # extracts the whole entry of a dictionary - either a matched set of brackets or
    #  the text before a comma if a comma is seen before brackets.
    # with bracket argument allows you to control if the text is returned with the surrounding brackets
    bracket_counter = 0
    modifier = 0
    if not with_brackets:
        modifier = 1

    for m in re.finditer("{|}|,", dict_str):
        if m.group() == "{":
            if bracket_counter == 0:
                start_val = m.start() + modifier
            bracket_counter += 1
        elif m.group() == "}":
            bracket_counter -= 1
            if bracket_counter == 0:
                return start_val, m.end() - modifier
        elif m.group() == ",":
            # no brackets found
            if bracket_counter == 0:
                return dict_str.find(":") + 1, m.start()
    raise ValueError("End of parsing reached without matching brackets")
