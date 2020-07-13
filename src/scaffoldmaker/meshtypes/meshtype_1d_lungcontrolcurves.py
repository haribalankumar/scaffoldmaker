"""
Generates a 1-D  mesh using control curves
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractxyzPathParametersFromRegion
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
# from opencmiss.utils.zinc.general import ChangeManager

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

    curvePathDefaultScaffoldPackages = {
        'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 20
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3], [
                    [[-2.4, -11, 1.72], [-0.2, -0.45, 2.5], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-3 , -12.7, 7], [0.017, -0.15, 2.20], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[ -3.9, -12.7, 11.4], [-0.15, 0.8, 2.15], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[ -3.8 , -11.4, 14.9], [-0.067, 2.06, 0.37], [-2.0, 10.0, 22.0], [0, 0, 0]],
                    [[ -4,  -8.7,  17.5], [0.03, -1.37, 0.85], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[ -3.3, -5.1, 16], [0.4, -1.04, 0.6], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[ -4.4 , -3.3, 13.1], [0.252, -0.61, 0.62], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[ -5.7, -2.7,  11.2], [0.45, -0.6, 1.12], [-2.0, 10.0, 22.0], [0, 0, 0]],
                    [[ -6.2, -3.2,  9.82], [-0.13, 1.16, 0.55], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[ -1.6,  -10.2, 11.7], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[ -1.8, -8.4, 13], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[ -4,  -5.54,  13],  [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[ -6.6, -10.3, 13.4], [21.2, -8.1, 0.3], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[  -6.5, -6.6, 13.3], [11, -30, -3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-5.4, -4.6,  13.2], [-5.92, -27.75, 0.66], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-1, -9.2, 4.6], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-3.8, -8.1, 6.9], [-0.6, 1.1, 1.63], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[ -5.6, -5.6, 9], [-0.34, 1.26, 1.3], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-5.5, -10.3, 1.5], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-7.65, -7.8, 3.8], [0.2, 1.05, 1.3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-8, -5.3, 6.8], [1.1, 0.93, 1.35], [-5.0, 4.0, 29.0],[0,0,0] ] ] )
            } ),
        'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
            'scaffoldSettings': {
                'Coordinate dimensions': 3,
                'Length': 1.0,
                'Number of elements': 20
            },
            'meshEdits': exnodeStringFromNodeValues(
                [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3], [
                    [[-12.7, -11.4, 1.31], [-0.3, -1.5, 1.6], [-24.0, -6.0, -12.0], [0, 0, 0]],
                    [[-13.0, -12.6, 7.0], [-0.2, -0.45, 2.5], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-13.0, -12.3, 11.4], [0.017, -0.15, 2.20], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-13.5, -11.4, 15.2], [-0.15, 0.8, 2.15], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-13.4, -9.21, 17.3], [-0.067, 2.06, 0.37], [-2.0, 10.0, 22.0], [0, 0, 0]],
                    [[-13.6, -6.2, 17.2], [0.03, -1.37, 0.85], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-13.9, -4.7, 15.9], [0.4, -1.04, 0.6], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-14.8, -3.44, 13.9], [0.252, -0.61, 0.62], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-15.1, -2.7, 12.1], [0.45, -0.6, 1.12], [-2.0, 10.0, 22.0], [0, 0, 0]],
                    [[-11.4, -9.9, 12.4], [-0.13, 1.16, 0.55], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-12.5, -7.75, 13.6], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-13.6, -6.0, 14.5], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-78.9, -66.7, 99.2], [21.2, -8.1, 0.3], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-61.0, -90.0, 100.7], [11, -30, -3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-62.0, -127.5, 99.6], [-5.92, -27.75, 0.66], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-11.0, -9.6, 5.2], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-13.2, -7.8, 9], [-0.6, 1.1, 1.63], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-14, -5.81, 11.9], [-0.34, 1.26, 1.3], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-17.4, -9.2, 2.26], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-18.2, -6.4, 5.5], [0.2, 1.05, 1.3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-17.4, -3.9, 8.6], [1.1, 0.93, 1.35], [-5.0, 4.0, 29.0], [0, 0, 0] ] ] )
            } ),
        }

    @staticmethod
    def getName():
        return '1D Lung control curve'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Pig 1',
            'Mouse 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Pig 1' in parameterSetName:
            curvePathOption = cls.curvePathDefaultScaffoldPackages['Pig 1']
        else:
            curvePathOption = cls.curvePathDefaultScaffoldPackages['Mouse 1']

        options = {
            'Curve path' : copy.deepcopy(curvePathOption),
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Curve path',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Curve path':
            return [ MeshType_1d_path1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Curve path':
            return list(cls.curvePathDefaultScaffoldPackages.keys())
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
        if optionName == 'Curve path':
            if not parameterSetName:
                parameterSetName = list(cls.curvePathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.curvePathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Curve path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Curve path'):
            options['Curve path'] = cls.getOptionScaffoldPackage('Curve path', MeshType_1d_path1)
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
        curvePath = options['Curve path']

        elementsCountAlong = 20

        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        zero = [0,0,0]

        ## Posterior path
        tmpRegion = region.createRegion()
        curvePath.generate(tmpRegion)
        cx, cd1, cd2, cd3 = extractxyzPathParametersFromRegion(tmpRegion)

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
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)

        nodeIdentifier = 1
        dx_ds2 = [ 0.0, 1.0, 0.0 ]
        d2x_ds1ds2 = [ 0.0, 0.0, 0.0 ]
        for n in range(elementsCountAlong):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            print('cx=',cx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d2x_ds1ds2)
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
        # for e in range(elementsCountAlong):
        for e in range(4):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            element.setNodesByIdentifier(eft, [ e + 1, e + 2 ])
            elementIdentifier = elementIdentifier + 1

        fm.endChange()

        return []











    #
    #
    # posteriorPathDefaultScaffoldPackages = {
    #     'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings' : {
    #             'Coordinate dimensions' : 3,
    #             'Length' : 1.0,
    #             'Number of elements' : 4
    #             },
    #         'meshEdits' : exnodeStringFromNodeValues(
    #             [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
    #             [ [  -12.7,  -11.4, 1.31 ], [ -0.3,  -1.5, 1.6 ], [ -24.0,  -6.0, -12.0 ], [ 0.0, 0.0, 0.0 ] ],
    #             [ [ -13.0, -12.6, 7.0 ], [ -0.2,  -0.45, 2.5 ], [ -22.0,  -4.0,  -8.0 ], [  0.0,  0.0, 0.0 ] ],
    #             [ [ -13.0, -12.3, 11.4 ], [ 0.017, -0.15, 2.20 ], [ -10.0,  20.0,   8.0 ], [ 0.0,  0.0, 0.0 ] ],
    #             [ [ -13.5, -11.4, 15.2 ], [ -0.15, 0.8, 2.15 ], [  -5.0,   4.0,  29.0 ], [   0.0,  0.0, 0.0 ] ],
    #             [ [ -13.4,  -9.21, 17.3 ], [ -0.067, 2.06, 0.37 ], [  -2.0,  10.0,  22.0 ], [ 0.0,  0.0, 0.0 ] ] ])
    #         } ),
    #     'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-96.0, -45.0, 282.21], [3.0, 15.0, -16.0], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-99.5, -45.6, 180.1], [-7, -2.3, -44.45], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-99.5, -63.6, 107.9], [-2.017, -10.15, -32.20], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-106.0, -87.4, 52.2], [0.05, -18.0, -25.15], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-100.4, -133.21, 27.0], [7.32, -36.06, -3.85], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings' : {
    #             'Coordinate dimensions' : 3,
    #             'Length' : 1.0,
    #             'Number of elements' : 1
    #             },
    #         'meshEdits' : exnodeStringFromNodeValues(
    #             [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
    #             [ [   0.0,   0.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 0.0, 0.0, 0.0 ] ],
    #             [ [  -2.0,  11.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 0.0, 0.0 ] ] ] )
    #         } )
    #     }
    #
    # anteriorPathDefaultScaffoldPackages = {
    #     'Mouse 1' : ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings' : {
    #             'Coordinate dimensions' : 3,
    #             'Length' : 1.0,
    #             'Number of elements' : 4
    #             },
    #         'meshEdits' : exnodeStringFromNodeValues(
    #             [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2  ], [
    #             [[-15.1, -2.7, 12.1], [0.45, -0.6, 1.12], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]],
    #             [[-14.8, -3.44, 13.9], [0.252, -0.61, 0.62], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #             [[-13.9, -4.7, 15.9], [0.4, -1.04, 0.6], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #             [[-13.6, -6.2, 17.2], [0.03, -1.37, 0.85], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #             [[-13.4,  -9.21, 17.3], [0.394, -1.71, -0.64 ], [ -24.0, -6.0, -12.0 ], [ 0.0,  0.0, 0.0 ]] ])
    #         } ),
    #     'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-110.1, -175.7, 138.1], [0.45, 6, 11.12], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]],
    #                 [[-98.8, -169.44, 114.9], [7.452, 13.61, -13.12], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-83.3, -153.1, 101.9], [-2.4, 11.14, -18.2], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-87.6, -154.2, 57.2], [1.13, 4.37, -21.0], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-99.4, -137.21, 27.3], [-6.1, 28.71, -8.64], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Human 1' : ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings' : {
    #             'Coordinate dimensions' : 3,
    #             'Length' : 1.0,
    #             'Number of elements' : 1
    #             },
    #         'meshEdits' : exnodeStringFromNodeValues(
    #             [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2 ], [
    #             [ [   0.0,   10.0,  0.0 ], [  6.0, 12.0,  -2.0 ], [ 2.0,  1.0,  2.0 ], [ 0.0, 0.0, 0.0 ] ],
    #             [ [  -2.0,  21.0, -3.0 ], [ -8.0,  4.0,   9.0 ], [ 2.0,  2.0,  1.0 ], [ 0.0, 0.0, 0.0 ] ] ] )
    #         } )
    #     }
    #
    # accesslobePathDefaultScaffoldPackages = {
    #     'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[ -13.0, -12.3, 11.4 ], [ 1.6, 0.7, -0.34], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-11.4, -9.9, 12.4], [-0.13, 1.16, 0.55], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-12.5, -7.75, 13.6], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-13.6, -6.0, 14.5], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-13.9, -4.7, 15.9], [-0.3, 0.82, 1.1], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-101.0, -63.2, 105.4], [-13.6, 1.2, -1.14], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-114.4, -85.1, 109.4], [3.81, -13.16, 3.15], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-109.5, -103.25, 112.6], [4.46, -20.86, 3.5], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-96.6, -129.0, 111.5], [6.1, -11.15, -6.13], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-83.1, -152.7, 100.2], [4.1, -6.22, -9.1], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Human 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 1
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
    #                 [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
    #     })
    # }
    #
    # lateralsidePathDefaultScaffoldPackages = {
    #     'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-99.0, -61.8, 103.4 ], [14, 7.7, -5], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-78.9, -66.7, 99.2], [21.2, -8.1, 0.3], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-61.0, -90.0, 100.7], [11, -30, -3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-62.0, -127.5, 99.6], [ -5.92, -27.75, 0.66], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-81.1, -152.7, 101.9], [-8.3, -8.84, -1.0], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-98.0, -46.3, 284.4], [-17.9, -16.2, 1.8], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-130.9, -89., 234.2], [-2.3, 1.5, 1.07], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-111.0, -115.0, 181.7], [14.6, -6.6, 19.3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-16.0, -6.5, 14.6], [0.73, 0.85, 0.61], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-13.9, -4.7, 15.9], [0.94, 1.04, -0.16], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Human 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 1
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
    #                 [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
    #     })
    # }
    #
    #
    # basemedialPathDefaultScaffoldPackages = {
    #     'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-12.7, -11.4, 1.31], [1.8, 1.62, -0.18], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-11.0, -9.6, 5.2], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-13.2, -7.8, 9], [-0.6, 1.1, 1.63], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-14, -5.81, 11.9], [-0.34, 1.26, 1.3], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-15.1, -2.7, 12.1], [-0.76, 1.64, 0.13], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-98.7, -46.4, 284.31], [-18.0, -16.2, 1.8], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
    #                 [[-130.0, -89.7, 234.2], [7.71, -30.22, -38.05], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-111.2, -115.8, 181], [14.67, -6.5, -19.63], [0.0, -10.0, 0.0], [0.0, 0.0, 0.0]],
    #                 [[-99, -142.81, 142.5], [-5.94, -17.26, -11.3], [-5.0, 4.0, 29.0], [0.0, 0.0, 0.0]],
    #                 [[-113.1, -172.7, 138.1], [7.6, -16.4, 1.3], [-2.0, 10.0, 22.0], [0.0, 0.0, 0.0]]])
    #     }),
    #     'Human 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 1
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
    #                 [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
    #     })
    # }
    #
    # baselateralPathDefaultScaffoldPackages = {
    #     'Mouse 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-12.7, -11.4, 1.31], [-2.6, 0.24, -1.133], [-24.0, -6.0, -12.0], [0.0, -1.0, 0.0]],
    #                 [[-17.4, -9.2, 2.26], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-18.2, -6.4, 5.5], [0.2, 1.05, 1.3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-17.4, -3.9, 8.6], [1.1, 0.93, 1.35], [-5.0, 4.0, 29.0], [0.0, 1.0, 0.0]],
    #                 [[-15.1, -2.7, 12.1], [1.61, 0.74, 1.0], [-2.0, 10.0, 22.0], [5.0, 0.0, 0.0]]])
    #     }),
    #     'Pig 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 4
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[-89.7, -45.4, 283.1], [39.0, -5.24, -0.533], [-24.0, -6.0, -12.0], [0.0, -1.0, 0.0]],
    #                 [[-47.4, -55.2, 263.0], [32.73, -25.51, -16.38], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
    #                 [[-35.2, -109.4, 211.5], [-4.5, -2.55, -16.5], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
    #                 [[-58.7, -149.0, 156.3], [-32.0, -32.3, -32.5], [-5.0, 4.0, 29.0], [0.0, 1.0, 0.0]],
    #                 [[-111.1, -175.7, 139.1], [-16.1, -7.4, -10.0], [-2.0, 10.0, 22.0], [5.0, 0.0, 0.0]]])
    #     }),
    #     'Human 1': ScaffoldPackage(MeshType_1d_path1, {
    #         'scaffoldSettings': {
    #             'Coordinate dimensions': 3,
    #             'Length': 1.0,
    #             'Number of elements': 1
    #         },
    #         'meshEdits': exnodeStringFromNodeValues(
    #             [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], [
    #                 [[0.0, 0.0, 0.0], [6.0, 12.0, -2.0], [2.0, 1.0, 2.0], [6.0, 0.0, 3.0]],
    #                 [[-2.0, 11.0, -3.0], [-8.0, 4.0, 9.0], [2.0, 2.0, 1.0], [0.0, 1.0, 2.0]]])
    #     })
    # }
    #
    #

    # @staticmethod
    # def getName():
    #     return '1D Lung control curve'
    #
    # @staticmethod
    # def getParameterSetNames():
    #     return [
    #         'Default',
    #         'Human 1',
    #         'Pig 1',
    #         'Mouse 1']
    #
    # @classmethod
    # def getDefaultOptions(cls, parameterSetName='Default'):
    #     if 'Mouse 1' in parameterSetName:
    #         posteriorPathOption = cls.posteriorPathDefaultScaffoldPackages['Mouse 1']
    #         anteriorPathOption = cls.anteriorPathDefaultScaffoldPackages['Mouse 1']
    #         accesslobePathOption = cls.accesslobePathDefaultScaffoldPackages['Mouse 1']
    #         lateralsidePathOption = cls.lateralsidePathDefaultScaffoldPackages['Mouse 1']
    #         basemedialPathOption = cls.basemedialPathDefaultScaffoldPackages['Mouse 1']
    #         baselateralPathOption = cls.baselateralPathDefaultScaffoldPackages['Mouse 1']
    #
    #     if 'Human 1' in parameterSetName:
    #         posteriorPathOption = cls.posteriorPathDefaultScaffoldPackages['Human 1']
    #         anteriorPathOption = cls.anteriorPathDefaultScaffoldPackages['Human 1']
    #         lateralsidePathOption = cls.lateralsidePathDefaultScaffoldPackages['Human 1']
    #         accesslobePathOption = cls.accesslobePathDefaultScaffoldPackages['Human 1']
    #         basemedialPathOption = cls.basemedialPathDefaultScaffoldPackages['Human 1']
    #         baselateralPathOption = cls.baselateralPathDefaultScaffoldPackages['Human 1']
    #     else:
    #         posteriorPathOption = cls.posteriorPathDefaultScaffoldPackages['Mouse 1']
    #         anteriorPathOption = cls.anteriorPathDefaultScaffoldPackages['Mouse 1']
    #         accesslobePathOption = cls.accesslobePathDefaultScaffoldPackages['Mouse 1']
    #         lateralsidePathOption = cls.lateralsidePathDefaultScaffoldPackages['Mouse 1']
    #         basemedialPathOption = cls.basemedialPathDefaultScaffoldPackages['Mouse 1']
    #         baselateralPathOption = cls.baselateralPathDefaultScaffoldPackages['Mouse 1']
    #
    #     options = {
    #         'Posterior path' : copy.deepcopy(posteriorPathOption),
    #         'Anterior path': copy.deepcopy(anteriorPathOption),
    #         'Accesslobe path': copy.deepcopy(accesslobePathOption),
    #         'Lateralside path': copy.deepcopy(lateralsidePathOption),
    #         'Basemedial path': copy.deepcopy(basemedialPathOption),
    #         'Baselateral path': copy.deepcopy(baselateralPathOption),
    #         'Refine' : False,
    #         'Refine number of elements around' : 1,
    #         'Refine number of elements along' : 1,
    #         'Refine number of elements through wall' : 1
    #         }
    #     return options
    #
    # @staticmethod
    # def getOrderedOptionNames():
    #     return [
    #         'Posterior path',
    #         'Anterior path',
    #         'Accesslobe path',
    #         'Lateralside path',
    #         'Basemedial path',
    #         'Baselateral path',
    #         'Refine',
    #         'Refine number of elements around',
    #         'Refine number of elements along',
    #         'Refine number of elements through wall' ]
    #
    # @classmethod
    # def getOptionValidScaffoldTypes(cls, optionName):
    #     if optionName == 'Posterior path':
    #         return [ MeshType_1d_path1 ]
    #     if optionName == 'Anterior path':
    #         return [ MeshType_1d_path1 ]
    #     if optionName == 'Accesslobe path':
    #         return [ MeshType_1d_path1 ]
    #     if optionName == 'Lateralside path':
    #         return [ MeshType_1d_path1 ]
    #     if optionName == 'Basemedial path':
    #         return [ MeshType_1d_path1 ]
    #     if optionName == 'Baselateral path':
    #         return [ MeshType_1d_path1 ]
    #     return []
    #
    # @classmethod
    # def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
    #     if optionName == 'Posterior path':
    #         return list(cls.posteriorPathDefaultScaffoldPackages.keys())
    #     if optionName == 'Anterior path':
    #         return list(cls.anteriorPathDefaultScaffoldPackages.keys())
    #     if optionName == 'Accesslobe path':
    #         return list(cls.accesslobePathDefaultScaffoldPackages.keys())
    #     if optionName == 'Lateralside path':
    #         return list(cls.lateralsidePathDefaultScaffoldPackages.keys())
    #     if optionName == 'Basemedial path':
    #         return list(cls.basemedialPathDefaultScaffoldPackages.keys())
    #     if optionName == 'Baselateral path':
    #         return list(cls.baselateralPathDefaultScaffoldPackages.keys())
    #     assert scaffoldType in cls.getOptionValidScaffoldTypes(optionName), \
    #         cls.__name__ + '.getOptionScaffoldTypeParameterSetNames.  ' + \
    #         'Invalid option \'' + optionName + '\' scaffold type ' + scaffoldType.getName()
    #     return scaffoldType.getParameterSetNames()
    #
    # @classmethod
    # def getOptionScaffoldPackage(cls, optionName, scaffoldType, parameterSetName=None):
    #     '''
    #     :param parameterSetName:  Name of valid parameter set for option Scaffold, or None for default.
    #     :return: ScaffoldPackage.
    #     '''
    #     if parameterSetName:
    #         assert parameterSetName in cls.getOptionScaffoldTypeParameterSetNames(optionName, scaffoldType), \
    #             'Invalid parameter set ' + str(parameterSetName) + ' for scaffold ' + str(scaffoldType.getName()) + \
    #             ' in option ' + str(optionName) + ' of scaffold ' + cls.getName()
    #     if optionName == 'Posterior path':
    #         if not parameterSetName:
    #             parameterSetName = list(cls.posteriorPathDefaultScaffoldPackages.keys())[0]
    #         return copy.deepcopy(cls.posteriorPathDefaultScaffoldPackages[parameterSetName])
    #     if optionName == 'Anterior path':
    #         if not parameterSetName:
    #             parameterSetName = list(cls.anteriorPathDefaultScaffoldPackages.keys())[0]
    #         return copy.deepcopy(cls.anteriorPathDefaultScaffoldPackages[parameterSetName])
    #     if optionName == 'Accesslobe path':
    #         if not parameterSetName:
    #             parameterSetName = list(cls.accesslobePathDefaultScaffoldPackages.keys())[0]
    #         return copy.deepcopy(cls.accesslobePathDefaultScaffoldPackages[parameterSetName])
    #     if optionName == 'Lateralside path':
    #         if not parameterSetName:
    #             parameterSetName = list(cls.lateralsidePathDefaultScaffoldPackages.keys())[0]
    #         return copy.deepcopy(cls.lateralsidePathDefaultScaffoldPackages[parameterSetName])
    #     if optionName == 'Basemedial path':
    #         if not parameterSetName:
    #             parameterSetName = list(cls.basemedialPathDefaultScaffoldPackages.keys())[0]
    #         return copy.deepcopy(cls.basemedialPathDefaultScaffoldPackages[parameterSetName])
    #     if optionName == 'Baselateral path':
    #         if not parameterSetName:
    #             parameterSetName = list(cls.baselateralPathDefaultScaffoldPackages.keys())[0]
    #         return copy.deepcopy(cls.baselateralPathDefaultScaffoldPackages[parameterSetName])
    #     assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'
    #
    # @classmethod
    # def checkOptions(cls, options):
    #     if not options['Posterior path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Posterior path'):
    #         options['Posterior path'] = cls.getOptionScaffoldPackage('Posterior path', MeshType_1d_path1)
    #     if not options['Anterior path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Anterior path'):
    #         options['Anterior path'] = cls.getOptionScaffoldPackage('Anterior path', MeshType_1d_path1)
    #     if not options['Accesslobe path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Accesslobe path'):
    #         options['Accesslobe path'] = cls.getOptionScaffoldPackage('Accesslobe path', MeshType_1d_path1)
    #     if not options['Lateralside path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Lateral path'):
    #         options['Lateralside path'] = cls.getOptionScaffoldPackage('Lateralside path', MeshType_1d_path1)
    #     if not options['Basemedial path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Basemedial path'):
    #         options['Basemedial path'] = cls.getOptionScaffoldPackage('Basemedial path', MeshType_1d_path1)
    #     if not options['Baselateral path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Baselateral path'):
    #         options['Baselateral path'] = cls.getOptionScaffoldPackage('Baselateral path', MeshType_1d_path1)
    #     for key in [
    #         'Refine number of elements around',
    #         'Refine number of elements along',
    #         'Refine number of elements through wall']:
    #         if options[key] < 1:
    #             options[key] = 1
    #
    # @classmethod
    # def generateBaseMesh(cls, region, options):
    #     """
    #     Generate the base tricubic Hermite mesh. See also generateMesh().
    #     :param region: Zinc region to define model in. Must be empty.
    #     :param options: Dict containing options. See getDefaultOptions().
    #     :return: annotationGroups
    #     """
    #     posteriorPath = options['Posterior path']
    #     anteriorPath = options['Anterior path']
    #     accesslobePath = options['Accesslobe path']
    #     lateralsidePath = options['Lateralside path']
    #     basemedialPath = options['Basemedial path']
    #     baselateralPath = options['Baselateral path']
    #
    #     elementsCountAlong = 4
    #
    #     firstNodeIdentifier = 1
    #     firstElementIdentifier = 1
    #




# class GetNode:
#     '''
#     Stores a tree of 1-D bifurcating curves through nested children.
#     '''
#
#     def __init__(self, x, d1, d2):
#         '''
#         :param x: coordinates
#         :param d1: Parent/primary coordinate derivative
#         :param r: Parent/primary radius
#         '''
#         self._x = x
#         self._d1 = [ d1 ]  # list to support versions
#         self._d2 = [ d2 ]  # list to support versions
#
# class ControlCurves:
#     '''
#     Class for generating 1D curves and converting to Zinc model.
#     '''
#
#     def __init__(self, number_of_curves, options):
#         '''
#         '''
#         self._curvecount = number_of_curves
#         self._options = options
#         # rootDirection = [ 0.0, 0.0, rootLength ]
#         # self._rootNode = TreeNode([ 0.0, 0.0, 0.0 ], rootDirection, rootRadius)
#         self._rootNode.addCurve(self._createNodeTree(1, [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ]))
#
#     def _createNodeTree(self, curve, x1, d1, d2):
#         '''
#         Create node with specified x1, d1, d2 and recursively add nodes until generationCount.
#         :param forkNormal: Unit direction normal to d1 and child branches.
#         :return: Top node of tree.
#         '''
#         node = GetNode(x1, d1, d2)
#         if curve < self.number_of_curves:
#             node.addCurve(self._createNodeTree(curve))
#         return node
#
#     def getRootNode(self):
#         return self._rootNode
#
#     def generateZincModel(self, region, nextNodeIdentifier=1, nextElementIdentifier=1):
#         '''
#         Generate Zinc nodes and elements in region to represent tree.
#         :return: Final nextNodeIdentifier, nextElementIdentifier.
#         '''
#         self._fieldmodule = region.getFieldmodule()
#         self._nodes = self._fieldmodule.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
#         self._mesh1d = self._fieldmodule.findMeshByDimension(1)
#         self._cubicHermiteBasis = self._fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
#         self._linearBasis = self._fieldmodule.createElementbasis(1, Elementbasis.FUNCTION_TYPE_LINEAR_LAGRANGE)
#         self._nodetemplates = {}  # indexed by (d1VersionsCount, rVersionsCount)
#         self._elementtemplates = {}  # indexed by start (d1Version, rVersion)
#         with ChangeManager(self._fieldmodule):
#             self._coordinates = findOrCreateFieldCoordinates(self._fieldmodule)
#             # self._radius = findOrCreateFieldFiniteElement(self._fieldmodule, "radius", components_count=1, managed=True)
#             self._fieldcache = self._fieldmodule.createFieldcache()
#             parentNode = None
#             self._generateZincModelTree(self._rootNode, parentNode, nextNodeIdentifier, nextElementIdentifier)
#         return nextNodeIdentifier, nextElementIdentifier
#
#     def _getZincNodetemplate(self, d1VersionsCount, d2VersionsCount):
#         '''
#         Get node template for specified numbers of d1 and r versions.
#         :return: Zinc Nodetemplate
#         '''
#         templateId = (d1VersionsCount, d2VersionsCount)
#
#         nodetemplate = self._nodetemplates.get(templateId)
#         if not nodetemplate:
#             assert (d1VersionsCount > 0) and (d2VersionsCount > 0)
#             nodetemplate = self._nodes.createNodetemplate()
#             nodetemplate.defineField(self._coordinates)
#             nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS1, d1VersionsCount)
#
#             # nodetemplate.defineField(self._radius)
#             nodetemplate.defineField(self._coordinates)
#             nodetemplate.setValueNumberOfVersions(self._coordinates, -1, Node.VALUE_LABEL_D_DS2, d2VersionsCount)
#             # nodetemplate.setValueNumberOfVersions(self._radius, -1, Node.VALUE_LABEL_VALUE, rVersionsCount)
#             self._nodetemplates[templateId] = nodetemplate
#         return nodetemplate
#
#     def _getZincElementtemplate(self, d1Version, rVersion):
#         '''
#         Get node template for specified numbers of d1 and r versions.
#         :return: Zinc Elementtemplate
#         '''
#         templateId = (d1Version, rVersion)
#         elementtemplate = self._elementtemplates.get(templateId)
#         if not elementtemplate:
#             assert (d1Version > 0) and (rVersion > 0)
#             elementtemplate = self._mesh1d.createElementtemplate()
#             elementtemplate.setElementShapeType(Element.SHAPE_TYPE_LINE)
#             eftCoordinates = self._mesh1d.createElementfieldtemplate(self._cubicHermiteBasis)
#             eftCoordinates.setTermNodeParameter(2, 1, 1, Node.VALUE_LABEL_D_DS1, d1Version)
#             elementtemplate.defineField(self._coordinates, -1, eftCoordinates)
#             eftRadius = self._mesh1d.createElementfieldtemplate(self._linearBasis)
#             eftRadius.setTermNodeParameter(1, 1, 1, Node.VALUE_LABEL_VALUE, rVersion)
#             elementtemplate.defineField(self._radius, -1, eftRadius)
#             self._elementtemplates[templateId] = elementtemplate
#         return elementtemplate
#
#     def _generateZincModelTree(self, treeNode, parentNode, nextNodeIdentifier, nextElementIdentifier):
#         '''
#         :return: Final nextNodeIdentifier, nextElementIdentifier.
#         '''
#         d1VersionsCount = len(treeNode._d1)
#         d2VersionsCount = len(treeNode._d2)
#         nodetemplate = self._getZincNodetemplate(d1VersionsCount, d2VersionsCount)
#         node = self._nodes.createNode(nextNodeIdentifier, nodetemplate)
#         nextNodeIdentifier += 1
#         self._fieldcache.setNode(node)
#         self._coordinates.setNodeParameters(self._fieldcache, -1, Node.VALUE_LABEL_VALUE, 1, treeNode._x)
#         for i in range(d1VersionsCount):
#             self._coordinates.setNodeParameters(self._fieldcache, -1, Node.VALUE_LABEL_D_DS1, i + 1, treeNode._d1[i])
#         for i in range(d2VersionsCount):
#             self._coordinates.setNodeParameters(self._fieldcache, -1, Node.VALUE_LABEL_D_DS2, i + 1, treeNode._d2[i])
#
#         if parentNode:
#             elementtemplate = self._getZincElementtemplate(treeNode._parent_d1_index + 1, treeNode._parent_r_index + 1)
#             element = self._mesh1d.createElement(nextElementIdentifier, elementtemplate)
#             nextElementIdentifier += 1
#             eftCoordinates = element.getElementfieldtemplate(self._coordinates, -1)
#             eftRadius = element.getElementfieldtemplate(self._radius, -1)
#             # must set node for both efts
#             for eft in ( eftCoordinates, eftRadius ):
#                 element.setNode(eft, 1, parentNode)
#                 element.setNode(eft, 2, node)
#
#         for childTreeNode in treeNode._children:
#             nextNodeIdentifier, nextElementIdentifier = self._generateZincModelTree(childTreeNode, node, nextNodeIdentifier, nextElementIdentifier)
#
#         return nextNodeIdentifier, nextElementIdentifier
#
