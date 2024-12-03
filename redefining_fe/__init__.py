
from redefining_fe.cells import Point, Edge, n_sided_polygon, make_tetrahedron, constructCellComplex
from redefining_fe.groups import r, S1, S2, S3, D4, Z3, Z4, C3, C4, S4, A4, tri_C3, tet_edges, tet_faces, sq_edges, GroupRepresentation, PermutationSetRepresentation, get_cyc_group, get_sym_group
from redefining_fe.dof import DeltaPairing, DOF, L2InnerProd, MyTestFunction, PointKernel, PolynomialKernel
from redefining_fe.triples import ElementTriple, DOFGenerator, immerse
from redefining_fe.traces import TrH1, TrGrad, TrHess, TrHCurl, TrHDiv

from redefining_fe.spaces.element_sobolev_spaces import CellH1, CellL2, CellHDiv, CellHCurl, CellH2
from redefining_fe.spaces.polynomial_spaces import P0, P1, P2, P3, Q2, PolynomialSpace
from redefining_fe.spaces.interpolation_spaces import C0, L2, H1, HDiv
