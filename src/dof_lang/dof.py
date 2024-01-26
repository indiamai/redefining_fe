

class DegreeOfFreedom():

    def __init__(self):
        self.functional = None

        # immersion details
        self.output_level = 1
        self.initial_level = 0

        self.entity_association = (0, 0)

    def __call__(self, g):
        # this function applies g to all relevant points

        # changes self.entity association
        # then also apply to the functional?
