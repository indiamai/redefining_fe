import json
import redefining_fe
import redefining_fe.cells
import redefining_fe.groups
import sympy


class FETripleEncoder(json.JSONEncoder):

    seen_objects = {}

    # def __init__(self, *args):
    #     self.seen_objects = {}
    #     super(self, FETripleEncoder).__init__(*args)

    def default(self, o):
        o_dict = {}
        print("seen objs", self.seen_objects)
        if str(o) in self.seen_objects.keys():
            return self.seen_objects[str(o)]
        if isinstance(o, redefining_fe.triples.ElementTriple):
            o_dict = {self.obj_id(o): json.dumps(o.cell, cls=FETripleEncoder)}
        elif isinstance(o, redefining_fe.cells.Point):
            o_dict = {self.obj_id(o): {"dim": o.dimension,
                                       "group": json.dumps(o.group, cls=FETripleEncoder),
                                       "edges": [json.dumps(c, cls=FETripleEncoder) for c in o.connections]}}
        elif isinstance(o, redefining_fe.groups.GroupRepresentation) or isinstance(o, redefining_fe.groups.GroupMemberRep):
            o_dict = {self.obj_id(o): str(o)}
        elif isinstance(o, redefining_fe.cells.Edge):
            print(type(o.attachment))
            print(o.attachment)
            o_dict = {self.obj_id(o): {"attachment": json.dumps(o.attachment, cls=FETripleEncoder),
                                       "point": json.dumps(o.point, cls=FETripleEncoder),
                                       "orientation": json.dumps(o.o, cls=FETripleEncoder)}}
        elif isinstance(o, sympy.core.containers.Tuple):
            print("SYMPY")
            return {self.obj_id(o): str(o)}
        self.seen_objects[self.obj_id(o)] = "ref_" + self.obj_id(o)
        print(type(o).__bases__)
        if len(o_dict.keys()) == 0:
            return super().default(o)
        return o_dict
    
    def obj_id(self, o):
        return str(o.__class__) + str(id(o))
    # def serialise(self, obj):
    #     obj_dict = obj.__dict__
    #     for val in obj_dict.values():
    #         if val not in self.seen_objects.keys():
