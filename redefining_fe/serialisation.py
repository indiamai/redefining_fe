import json
from redefining_fe import *
from redefining_fe.spaces.polynomial_spaces import RestrictedPolynomialSpace, ConstructedPolynomialSpace
from redefining_fe.spaces.element_sobolev_spaces import ElementSobolevSpace
from redefining_fe.spaces.interpolation_spaces import InterpolationSpace
from redefining_fe.traces import Trace
from redefining_fe.triples import ImmersedDOFs
import sympy as sp
from sympy.combinatorics import Permutation
import re


class ElementSerialiser():

    def __init__(self):
        self.obj_id_counter = {}
        self.seen_objs = {}
        self.obj_storage = {}

        self.obj_types = {"Cell": Point,
                          "Edge": Edge,
                          "Triple": ElementTriple,
                          "Group": GroupRepresentation,
                          "SobolevSpace": ElementSobolevSpace,
                          "InterpolationSpace": InterpolationSpace,
                          "PolynomialSpace": PolynomialSpace,
                          "RestrictedPolynomialSpace": RestrictedPolynomialSpace,
                          "ConstructedPolynomialSpace": ConstructedPolynomialSpace,
                          "DOF": DOF,
                          "ImmersedDOF": ImmersedDOFs,
                          "DOFGen": DOFGenerator,
                          "Delta": DeltaPairing,
                          "L2Inner": L2InnerProd,
                          "PolynomialKernel": PolynomialKernel,
                          "PointKernel": PointKernel,
                          "Trace": Trace
                          }

    def encode(self, obj):
        base_obj = self.encode_traverse(obj)
        self.obj_storage["encoded_obj"] = base_obj
        return json.dumps(self.obj_storage, indent=2)

    def decode(self, obj_str):
        obj_dict = json.loads(obj_str)
        obj = self.decode_traverse(obj_dict["encoded_obj"], obj_dict)
        return obj

    def encode_traverse(self, obj, path=[]):
        obj_dict = {}
        if isinstance(obj, list) or isinstance(obj, tuple):
            res_array = [{} for i in range(len(obj))]
            for i in range(len(obj)):
                dfs_res = self.encode_traverse(obj[i], path + [i])
                res_array[i] = dfs_res
            if isinstance(obj, tuple):
                return tuple(res_array)
            return res_array

        if obj in self.seen_objs.keys():
            print("repeat hit")
            return self.seen_objs[obj]["id"]

        if hasattr(obj, "_to_dict"):
            for (key, val) in obj._to_dict().items():
                obj_dict[key] = self.encode_traverse(val, path + [key])
            obj_id = self.get_id(obj)
            self.store_obj(obj, obj.dict_id(), obj_id, obj_dict, path)
            return obj.dict_id() + " " + str(obj_id)

        if isinstance(obj, sp.core.containers.Tuple) or isinstance(obj, sp.Expr):
            return "Sympy " + sp.srepr(obj)
        # print(type(obj))
        return obj

    def get_id(self, obj):
        obj_name = obj.dict_id()
        if obj_name in self.obj_id_counter.keys():
            obj_id = self.obj_id_counter[obj_name]
            self.obj_id_counter[obj_name] += 1
        else:
            obj_id = 0
            self.obj_id_counter[obj_name] = 1
        return obj_id

    def store_obj(self, obj, name, obj_id, obj_dict, path):
        self.seen_objs[obj] = {"id": name + " " + str(obj_id), "path": path, "dict": obj_dict}
        if name in self.obj_storage.keys():
            self.obj_storage[name][obj_id] = obj_dict
        else:
            self.obj_storage[name] = {obj_id: obj_dict}

    def decode_traverse(self, obj, obj_dict):
        # better way to identify if soemthing is an obj string
        if isinstance(obj, str):
            split_str = obj.split(" ")
            if split_str[0] in self.obj_types.keys():
                name, obj_id = split_str[0], split_str[1]
                sub_dict = obj_dict[name][obj_id]
                for (key, value) in sub_dict.items():
                    sub_dict[key] = self.decode_traverse(value, obj_dict)
                return self.obj_types[name]._from_dict(sub_dict)
            elif split_str[0] == "Sympy":
                return sp.parse_expr(" ".join(split_str[1:]))
            else:
                return obj
        elif isinstance(obj, list) or isinstance(obj, tuple):
            res_array = [0 for i in range(len(obj))]
            for i in range(len(obj)):
                dfs_res = self.decode_traverse(obj[i], obj_dict)
                res_array[i] = dfs_res
            if isinstance(obj, tuple):
                return tuple(res_array)
            return res_array

        return obj





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
