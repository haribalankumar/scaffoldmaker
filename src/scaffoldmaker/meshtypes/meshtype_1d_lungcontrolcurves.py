"""
Generates a 1-D  mesh using control curves
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.zinc.node import Node

class MeshType_1d_lungcontrolcurves(Scaffold_base):
    '''
    Generates a 3-D lung mesh with variable numbers
    of elements around
    '''

    posteriorPathDefaultScaffoldPackages = {
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 4
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [ [  -12.7,  -11.4, 1.31 ], [ -0.3,  -1.5, 1.6 ], [ -24.0,  -6.0, -12.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [ -13.0, -12.6, 7.0 ], [ -0.2,  -0.45, 2.5 ], [ -22.0,  -4.0,  -8.0 ], [  0.0,  0.0, 0.0 ] ],
                [ [ -13.0, -12.3, 11.4 ], [ 0.017, -0.15, 2.20 ], [ -10.0,  20.0,   8.0 ], [ 0.0,  0.0, 0.0 ] ],
                [ [ -13.5, -11.4, 15.2 ], [ -0.15, 0.8, 2.15 ], [  -5.0,   4.0,  29.0 ], [   0.0,  0.0, 0.0 ] ],
                [ [ -13.4,  -9.21, 17.3 ], [ -0.067, 2.06, 0.37 ], [  -2.0,  10.0,  22.0 ], [ 0.0,  0.0, 0.0 ] ] ])
            } ),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-96.0, -45.0, 282.21], [3.0, 15.0, -16.0], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-99.5, -45.6, 180.1], [-7, -2.3, -44.45], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-99.5, -63.6, 107.9], [-2.017, -10.15, -32.20], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-106.0, -87.4, 52.2], [0.05, -18.0, -25.15], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-100.4, -133.21, 27.0], [7.32, -36.06, -3.85], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 1
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,   0.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [  -2.0,  11.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 0.0, 0.0 ] ] ] )
            } )
        }

    anteriorPathDefaultScaffoldPackages = {
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 4
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
                [[-15.1, -2.7, 12.1], [0.45, -0.6, 1.12], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]],
                [[-14.8, -3.44, 13.9], [0.252, -0.61, 0.62], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                [[-13.9, -4.7, 15.9], [0.4, -1.04, 0.6], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                [[-13.6, -6.2, 17.2], [0.03, -1.37, 0.85], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                [[-13.4,  -9.21, 17.3], [0.394, -1.71, -0.64 ], [ -24.0, -6.0, -12.0 ], [ 0.0,  0.0, 0.0 ]] ])
            } ),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-110.1, -175.7, 138.1], [0.45, 6, 11.12], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]],
                    [[-98.8, -169.44, 114.9], [7.452, 13.61, -13.12], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-83.3, -153.1, 101.9], [-2.4, 11.14, -18.2], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-87.6, -154.2, 57.2], [1.13, 4.37, -21.0], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-99.4, -137.21, 27.3], [-6.1, 28.71, -8.64], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]]])
        }),
        'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 1
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
                [ [   0.0,   10.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 0.0, 0.0, 0.0 ] ],
                [ [  -2.0,  21.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 0.0, 0.0 ] ] ] )
            } )
        }

    accesslobePathDefaultScaffoldPackages = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[ -13.0, -12.3, 11.4 ], [ 1.6, 0.7, -0.34], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-11.4, -9.9, 12.4], [-0.13, 1.16, 0.55], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-12.5, -7.75, 13.6], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-13.6, -6.0, 14.5], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-13.9, -4.7, 15.9], [-0.3, 0.82, 1.1], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-101.0, -63.2, 105.4], [-13.6, 1.2, -1.14], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-114.4, -85.1, 109.4], [3.81, -13.16, 3.15], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-109.5, -103.25, 112.6], [4.46, -20.86, 3.5], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-96.6, -129.0, 111.5], [6.1, -11.15, -6.13], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-83.1, -152.7, 100.2], [4.1, -6.22, -9.1], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 1
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
                    [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
        })
    }

    lateralsidePathDefaultScaffoldPackages = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-99.0, -61.8, 103.4 ], [14, 7.7, -5], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-78.9, -66.7, 99.2], [21.2, -8.1, 0.3], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-61.0, -90.0, 100.7], [11, -30, -3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-62.0, -127.5, 99.6], [ -5.92, -27.75, 0.66], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-81.1, -152.7, 101.9], [-8.3, -8.84, -1.0], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-98.0, -46.3, 284.4], [-17.9, -16.2, 1.8], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-130.9, -89., 234.2], [-2.3, 1.5, 1.07], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-111.0, -115.0, 181.7], [14.6, -6.6, 19.3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-16.0, -6.5, 14.6], [0.73, 0.85, 0.61], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-13.9, -4.7, 15.9], [0.94, 1.04, -0.16], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 1
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
                    [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
        })
    }


    basemedialPathDefaultScaffoldPackages = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-12.7, -11.4, 1.31], [1.8, 1.62, -0.18], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-11.0, -9.6, 5.2], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-13.2, -7.8, 9], [-0.6, 1.1, 1.63], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-14, -5.81, 11.9], [-0.34, 1.26, 1.3], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-15.1, -2.7, 12.1], [-0.76, 1.64, 0.13], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-98.7, -46.4, 284.31], [-18.0, -16.2, 1.8], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-130.0, -89.7, 234.2], [7.71, -30.22, -38.05], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-111.2, -115.8, 181], [14.67, -6.5, -19.63], [0.0, -10.0, 0.0], [0.0, 0.0, 0.0]],
                    [[-99, -142.81, 142.5], [-5.94, -17.26, -11.3], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
                    [[-113.1, -172.7, 138.1], [7.6, -16.4, 1.3], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 1
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
                    [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
        })
    }

    baselateralPathDefaultScaffoldPackages = {
        'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-12.7, -11.4, 1.31], [-2.6, 0.24, -1.133], [-24.0, -6.0, -12.0], [0.0, -1.0, 0.0]],
                    [[-17.4, -9.2, 2.26], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-18.2, -6.4, 5.5], [0.2, 1.05, 1.3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-17.4, -3.9, 8.6], [1.1, 0.93, 1.35], [-5.0, 4.0, 29.0], [0.0, 1.0, 0.0]],
                    [[-15.1, -2.7, 12.1], [1.61, 0.74, 1.0], [-2.0, 10.0, 22.0], [5.0, 0.0, 0.0]]])
        }),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 4
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[-89.7, -45.4, 283.1], [39.0, -5.24, -0.533], [-24.0, -6.0, -12.0], [0.0, -1.0, 0.0]],
                    [[-47.4, -55.2, 263.0], [32.73, -25.51, -16.38], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-35.2, -109.4, 211.5], [-4.5, -2.55, -16.5], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
                    [[-58.7, -149.0, 156.3], [-32.0, -32.3, -32.5], [-5.0, 4.0, 29.0], [0.0, 1.0, 0.0]],
                    [[-111.1, -175.7, 139.1], [-16.1, -7.4, -10.0], [-2.0, 10.0, 22.0], [5.0, 0.0, 0.0]]])
        }),
        'Human 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 1
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
                    [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
                    [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
        })
    }


    @staticmethod
    def getName():
        return '1D Lung control curve'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Pig 1',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Mouse 1' in parameterSetName:
            posteriorPathOption = cls.posteriorPathDefaultScaffoldPackages['Mouse 1']
            anteriorPathOption = cls.anteriorPathDefaultScaffoldPackages['Mouse 1']
            accesslobePathOption = cls.accesslobePathDefaultScaffoldPackages['Mouse 1']
            lateralsidePathOption = cls.lateralsidePathDefaultScaffoldPackages['Mouse 1']
            basemedialPathOption = cls.basemedialPathDefaultScaffoldPackages['Mouse 1']
            baselateralPathOption = cls.baselateralPathDefaultScaffoldPackages['Mouse 1']

        if 'Human 1' in parameterSetName:
            posteriorPathOption = cls.posteriorPathDefaultScaffoldPackages['Human 1']
            anteriorPathOption = cls.anteriorPathDefaultScaffoldPackages['Human 1']
            lateralsidePathOption = cls.lateralsidePathDefaultScaffoldPackages['Human 1']
            accesslobePathOption = cls.accesslobePathDefaultScaffoldPackages['Human 1']
            basemedialPathOption = cls.basemedialPathDefaultScaffoldPackages['Human 1']
            baselateralPathOption = cls.baselateralPathDefaultScaffoldPackages['Human 1']
        else:
            posteriorPathOption = cls.posteriorPathDefaultScaffoldPackages['Mouse 1']
            anteriorPathOption = cls.anteriorPathDefaultScaffoldPackages['Mouse 1']
            accesslobePathOption = cls.accesslobePathDefaultScaffoldPackages['Mouse 1']
            lateralsidePathOption = cls.lateralsidePathDefaultScaffoldPackages['Mouse 1']
            basemedialPathOption = cls.basemedialPathDefaultScaffoldPackages['Mouse 1']
            baselateralPathOption = cls.baselateralPathDefaultScaffoldPackages['Mouse 1']


        options = {
            'Posterior path' : copy.deepcopy(posteriorPathOption),
            'Anterior path': copy.deepcopy(anteriorPathOption),
            'Accesslobe path': copy.deepcopy(accesslobePathOption),
            'Lateralside path': copy.deepcopy(lateralsidePathOption),
            'Basemedial path': copy.deepcopy(basemedialPathOption),
            'Baselateral path': copy.deepcopy(baselateralPathOption),
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Posterior path',
            'Anterior path',
            'Accesslobe path',
            'Lateralside path',
            'Basemedial path',
            'Baselateral path',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Posterior path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Anterior path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Accesslobe path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Lateralside path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Basemedial path':
            return [ MeshType_1d_path1 ]
        if optionName == 'Baselateral path':
            return [ MeshType_1d_path1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Posterior path':
            return list(cls.posteriorPathDefaultScaffoldPackages.keys())
        if optionName == 'Anterior path':
            return list(cls.anteriorPathDefaultScaffoldPackages.keys())
        if optionName == 'Accesslobe path':
            return list(cls.accesslobePathDefaultScaffoldPackages.keys())
        if optionName == 'Lateralside path':
            return list(cls.lateralsidePathDefaultScaffoldPackages.keys())
        if optionName == 'Basemedial path':
            return list(cls.basemedialPathDefaultScaffoldPackages.keys())
        if optionName == 'Baselateral path':
            return list(cls.baselateralPathDefaultScaffoldPackages.keys())
        assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
            cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
            'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
        return scaffoldType.getParameterSetNames()

    @classmethod
    def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
        '''
        :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
        :return: ScaffoldPackage.
        '''
        if parameterSetName:
            assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
                'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
                ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
        if optionName == 'Posterior path':
            if not parameterSetName:
                parameterSetName = list(cls.posteriorPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.posteriorPathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Anterior path':
            if not parameterSetName:
                parameterSetName = list(cls.anteriorPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.anteriorPathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Accesslobe path':
            if not parameterSetName:
                parameterSetName = list(cls.accesslobePathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.accesslobePathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Lateralside path':
            if not parameterSetName:
                parameterSetName = list(cls.lateralsidePathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.lateralsidePathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Basemedial path':
            if not parameterSetName:
                parameterSetName = list(cls.basemedialPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.basemedialPathDefaultScaffoldPackages[parameterSetName])
        if optionName == 'Baselateral path':
            if not parameterSetName:
                parameterSetName = list(cls.baselateralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.baselateralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Posterior path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Posterior path'):
            options['Posterior path'] = cls.getOptionScaffoldPackage('Posterior path', MeshType_1d_path1)
        if not options['Anterior path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Anterior path'):
            options['Anterior path'] = cls.getOptionScaffoldPackage('Anterior path', MeshType_1d_path1)
        if not options['Accesslobe path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Accesslobe path'):
            options['Accesslobe path'] = cls.getOptionScaffoldPackage('Accesslobe path', MeshType_1d_path1)
        if not options['Lateralside path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Lateral path'):
            options['Lateralside path'] = cls.getOptionScaffoldPackage('Lateralside path', MeshType_1d_path1)
        if not options['Basemedial path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Basemedial path'):
            options['Basemedial path'] = cls.getOptionScaffoldPackage('Basemedial path', MeshType_1d_path1)
        if not options['Baselateral path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Baselateral path'):
            options['Baselateral path'] = cls.getOptionScaffoldPackage('Baselateral path', MeshType_1d_path1)
        for key in [
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1

    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: annotationGroups
        """
        posteriorPath = options['Posterior path']
        anteriorPath = options['Anterior path']
        accesslobePath = options['Accesslobe path']
        lateralsidePath = options['Lateralside path']
        basemedialPath = options['Basemedial path']
        baselateralPath = options['Baselateral path']

        elementsCountAlong = 4

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        zero = [0,0,0]

        ## Posterior path
        tmpRegion = region.createRegion()
        posteriorPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        tmpRegion = region.createRegion()
        anteriorPath.generate(tmpRegion)
        cxAnterior, cd1Anterior, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        tmpRegion = region.createRegion()
        accesslobePath.generate(tmpRegion)
        cxAccesslobe, cd1Accesslobe, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        tmpRegion = region.createRegion()
        lateralsidePath.generate(tmpRegion)
        cxLateralside, cd1Lateralside, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        tmpRegion = region.createRegion()
        basemedialPath.generate(tmpRegion)
        cxBasemedial, cd1Basemedial, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        tmpRegion = region.createRegion()
        baselateralPath.generate(tmpRegion)
        cxBaselateral, cd1Baselateral, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        # del tmpRegion

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm, components_count=3)
        cache = fm.createFieldcache()

        #################
        # Create nodes
        #################
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        nodeIdentifier = 1
        # x = [ 0.0, 0.0, 0.0 ]
        # dx_ds1 = [ length/elementsCount, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 1.0, 0.0 ]
        d2x_ds1ds2 = [ 0.0, 0.0, 0.0 ]
        for n in range(elementsCountAlong+1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(elementsCountAlong + 1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxAnterior[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Anterior[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

        #accessory lobe edge
        for n in range(elementsCountAlong +1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxAccesslobe[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Accesslobe[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

        #lateralside edge (not a physiological edge- for convenience only)
        for n in range(elementsCountAlong+1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxLateralside[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Lateralside[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

        # base medial edge
        for n in range(elementsCountAlong + 1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxBasemedial[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Basemedial[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1

        # base lateral  edge
        for n in range(elementsCountAlong + 1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxBaselateral[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Baselateral[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, d2x_ds1ds2)
            nodeIdentifier = nodeIdentifier + 1


        #################
        # Create elements
        #################
        mesh = fm.findMeshByDimension(1)
        cubicHermiteBasis = fm.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(cubicHermiteBasis)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1
        for e in range(elementsCountAlong):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1, e + 2 ])
            elementIdentifier = elementIdentifier + 1

        e0 = e+2
        for e in range(elementsCountAlong):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1 + e0, e + 2 + e0])
            elementIdentifier = elementIdentifier + 1

        e0 = e + 2 + e0
        for e in range(elementsCountAlong):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1 + e0, e + 2 + e0])
            elementIdentifier = elementIdentifier + 1

        e0 = e + 2 + e0
        for e in range(elementsCountAlong):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1 + e0, e + 2 + e0])
            elementIdentifier = elementIdentifier + 1

        e0 = e + 2 + e0
        for e in range(elementsCountAlong):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1 + e0, e + 2 + e0])
            elementIdentifier = elementIdentifier + 1

        e0 = e + 2 + e0
        for e in range(elementsCountAlong):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1 + e0, e + 2 + e0])
            elementIdentifier = elementIdentifier + 1

        fm.endChange()

        return []

def extractPathParametersFromRegion(region):
    '''
    Returns parameters of all nodes in region in identifier order.
    Assumes nodes in region have field coordinates (1 to 3 components).
    Currently limited to nodes with exactly value, d_ds1, d_ds2, d2_ds12,
    same as path 1 scaffold.
    :return: cx, cd1, cd2, cd12 (all padded with zeroes to 3 components)
    '''
    fm = region.getFieldmodule()
    coordinates = fm.findFieldByName('coordinates').castFiniteElement()
    componentsCount = coordinates.getNumberOfComponents()
    assert componentsCount in [1, 2,
                              3], 'extractPathParametersFromRegion.  Invalid coordinates number of components'
    cache = fm.createFieldcache()
    cx = []
    cd1 = []
    cd2 = []
    cd12 = []
    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodeIter = nodes.createNodeiterator()
    node = nodeIter.next()
    while node.isValid():
        cache.setNode(node)
        result, x = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, componentsCount)
        result, d1 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, componentsCount)
        result, d2 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, componentsCount)
        result, d12 = coordinates.getNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, componentsCount)
        for c in range(componentsCount, 3):
            x.append(0.0)
            d1.append(0.0)
            d2.append(0.0)
            d12.append(0.0)

        cx.append(x)
        cd1.append(d1)
        cd2.append(d2)
        cd12.append(d12)
        node = nodeIter.next()
    return cx, cd1, cd2, cd12
