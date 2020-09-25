"""
Generates a 3-D lung mesh using control curves
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_path1 import MeshType_1d_path1, extractPathParametersFromRegion
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite

from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.zinc.node import Node

class MeshType_3d_lungmesh(Scaffold_base):
    '''
    Generates a 3-D lung mesh using curves
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
                    [[-13.0, -12.3, 11.4 ], [ -1.28, -3.75, 0.1], [-24.0, -6.0, -12.0], [0.0, 0.0, 0.0]],
                    [[-15.9, -12., 12.2], [-2.3, 1.5, 1.07], [-22.0, -4.0, -8.0], [0.0, 0.0, 0.0]],
                    [[-17.0, -9.0, 13.7], [0.6, 2.6, 1.3], [-10.0, 20.0, 8.0], [0.0, 0.0, 0.0]],
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
        return '3D Lung mesh'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
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
            'Use cross derivatives' : False,
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
            'Use cross derivatives',
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
        useCrossDerivatives = options['Use cross derivatives']

        elementsCountAlong = 4
        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        fm = region.getFieldmodule()
        fm.beginChange()
        cache = fm.createFieldcache()

        # Coordinates field
        coordinates = findOrCreateFieldCoordinates(fm)
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

        mesh = fm.findMeshByDimension(3)
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        zero = [0,0,0]

        ## Posterior path
        tmpRegion = region.createRegion()
        posteriorPath.generate(tmpRegion)
        cx, cd1, cd2, cd12 = \
            (tmpRegion)
        del tmpRegion

        tmpRegion = region.createRegion()
        anteriorPath.generate(tmpRegion)
        cxAnterior, cd1Anterior, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)

        # cx_bladder, cd1_bladder, cd2_bladder, cd12_bladder = extractPathParametersFromRegion(tmpRegion,
        #     [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2], groupName='urinary bladder')
        # # for i in range(len(cx_bladder)):
        # cx, cd1, cd2, cd12 = extractPathParametersFromRegion(tmpRegion,
        #                                                      [Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1,
        #                                                       Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D2_DS1DS2])
        del tmpRegion

        tmpRegion = region.createRegion()
        accesslobePath.generate(tmpRegion)
        cxAccesslobe, cd1Accesslobe, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        del tmpRegion

        tmpRegion = region.createRegion()
        lateralsidePath.generate(tmpRegion)
        cxLateralside, cd1Lateralside, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        del tmpRegion

        tmpRegion = region.createRegion()
        basemedialPath.generate(tmpRegion)
        cxBasemedial, cd1Basemedial, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        del tmpRegion

        tmpRegion = region.createRegion()
        baselateralPath.generate(tmpRegion)
        cxBaselateral, cd1Baselateral, cd2, cd12 = extractPathParametersFromRegion(tmpRegion)
        del tmpRegion

        for n in range(len(cxLateralside)):
            print('cxBaselat=', cxBaselateral[n][0],cxBaselateral[n][1],cxBaselateral[n][2])
            print('cxlatside=', cxLateralside[n][0],cxLateralside[n][1],cxLateralside[n][2])
            print('cxAccesslobe=', cxAccesslobe[n][0],cxAccesslobe[n][1],cxAccesslobe[n][2])

        cxmidmedial = []
        cxmidlateral = []
        for n in range(len(cxLateralside)):
            for j in range(3):
                print('first n,j loop=',n,j)
                cxmidlateral[n][j] = 0.5 * (cxLateralside[n][j] + cxBaselateral[n][j])
                cxmidmedial[n][j] = 0.5 * (cxAccesslobe[n][j] + cxBasemedial[n][j])

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm, components_count=3)
        cache = fm.createFieldcache()

        # _, sPosteriorOffsetderiv2, _, _, _ = \
        #     interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountAlong)

        ##############
        # Create nodes
        ##############
        nodeIdentifier = firstNodeIdentifier

        for n in range(len(cxBaselateral)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxBaselateral[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Baselateral[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(1, len(cxBasemedial)-2):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxBasemedial[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Basemedial[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(0, 0):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1
        for n in range(1, len(cxLateralside)-2):
            print('n=',n)
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxmidlateral[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Lateralside[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1
        for n in range(len(cxLateralside)-1, len(cxLateralside)-1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxAnterior[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd1Anterior[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(1, len(cxAccesslobe)-1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxmidmedial[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Accesslobe[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1


        for n in range(len(cxLateralside)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxLateralside[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Lateralside[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(1, len(cxAccesslobe)-1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxAccesslobe[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Accesslobe[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1


        for n in range(elementsCountAlong-1,elementsCountAlong):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1


        for n in range(elementsCountAlong, elementsCountAlong):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cxAnterior[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1Anterior[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1


        ####################
        # Create elements
        ####################
        elementIdentifier = firstElementIdentifier
        #for edge elements
        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [1, 2, 9, 1, 6, 9, 14]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1

        #form normal elements (mid lung)
        for e2 in range(1, elementsCountAlong-2):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bn1 = 2
            bn2 = bn1 + elementsCountAlong
            bn3 = bn1 + 2*elementsCountAlong
            bn4 = bn2 + 2*elementsCountAlong
            nodeIdentifiers = [bn1, bn1+1, bn3, bn3+1, bn2, bn2+1, bn4, bn4+1]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

        #for edge elements
        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [4, 5, 12, 13, 8, 5, 16, 13]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1


        fm.endChange()

        return []


    @classmethod
    def refineMesh(cls, meshrefinement, options):
        """
        Refine source mesh into separate region, with change of basis.
        :param meshrefinement: MeshRefinement, which knows source and target region.
        :param options: Dict containing options. See getDefaultOptions().
        """
        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountAlong = options['Refine number of elements along']
        refineElementsCountThroughWall = options['Refine number of elements through wall']

        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong,
                                                       refineElementsCountThroughWall)
        return

