"""
Generates a 3-D lung mesh along the control lines
"""

import copy
from scaffoldmaker.meshtypes.meshtype_1d_lungpath1 import MeshType_1d_lungpath1, extractxyzPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import tubemesh
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.zinc.node import Node
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field


class MeshType_3d_lungcontrolcurves(Scaffold_base):
    '''
    Generates a 3-D lung mesh
    '''

    centralPathDefaultScaffoldPackages = {
        'Human 1' : ScaffoldPackage(MeshType_1d_lungpath1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 20,
                'Species' : 'Human'
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3], [
                    [[-2.4, -11, 1.72], [-2.0, -4.0, -8.0], [-0.2, -0.45, 2.5], [10, 1, 0]],
                    [[-3, -12.7, 7], [-1.0, 2.0, 8.0], [0.017, -0.15, 2.20], [0, 0, 0]],
                    [[-3.9, -12.7, 11.4], [-5.0, 4.0, 29.0], [-0.15, 0.8, 2.15], [10, 1, 0]],
                    [[-3.8, -11.4, 14.9], [-2.0, 10.0, 22.0], [-0.067, 2.06, 0.37], [0, 0, 0]],
                    [[-4, -8.7, 17.5], [-22.0, -4.0, -8.0], [0.03, -1.37, 0.85], [0, 0, 0]],
                    [[-3.3, -5.1, 16], [-10.0, 20.0, 8.0],[0.4, -1.04, 0.6],  [0, 0, 0]],
                    [[-4.4, -3.3, 13.1],  [-5.0, 4.0, 29.0], [0.252, -0.61, 0.62],[0, 0, 0]],
                    [[-5.7, -2.7, 11.2], [-2.0, 10.0, 22.0], [0.45, -0.6, 1.12], [0, 0, 0]],
                    [[-6.2, -3.2, 9.82],  [-22.0, -4.0, -8.0], [-0.13, 1.16, 0.55], [0, 0, 0]],
                    [[-1.6, -10.2, 11.7], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-1.8, -8.4, 13], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-4, -5.54, 13], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-6.6, -10.3, 13.4], [21.2, -8.1, 0.3], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-6.5, -6.6, 13.3], [11, -30, -3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-5.4, -4.6, 13.2], [-5.9, -27.8, 0.66], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-1, -9.2, 4.6], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-3.8, -8.1, 6.9], [-0.6, 1.1, 1.63], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-5.6, -5.6, 9], [-0.34, 1.26, 1.3], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-5.5, -10.3, 1.5], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-7.65, -7.8, 3.8], [0.2, 1.05, 1.3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-8, -5.3, 6.8], [1.1, 0.93, 1.35], [-5.0, 4.0, 29.0], [0, 0, 0]]])
        } ),
        'Mouse 1' : ScaffoldPackage(MeshType_1d_lungpath1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 20,
                'Species' : 'Mouse'
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ], [
                    [[-2.4, -11, 1.72], [2.1, 2.7, 2.0], [-0.1, -0.16, 2.1], [10, 1, 0]],
                    [[-3, -12.7, 7], [1.0, 1.7, 1.3], [-0.2, -0.2, 1.70], [0, 0, 0]],
                    [[-3.9, -12.7, 11.4], [1.5, 1.2, 0.4], [-0.15, 0.8, 2.15], [-1, 0, 0]],
                    [[-3.8, -11.4, 14.9], [1.1, 1.8, 0.1], [-0.2, 1.06, 2.5], [0, 0, 0]],
                    [[-4, -8.7, 17.5], [-3.0, 2.0, 1.2], [-1.1, 3.0, 0.7], [0, 0, 0]],
                    [[-3.3, -5.1, 16], [-1.0, 2.0, 1.0], [-0.6, 1.2, -1.6], [0, 0, 0]],
                    [[-4.4, -3.3, 13.1], [-0.08, 1.1, 0.2], [-0.12, 0.61, -1.2], [0, 0, 0]],
                    [[-5.7, -2.7, 11.2], [0.1, 1.6, 0.6], [-0.45, -0.01, -1.12], [0, 0, 0]],
                    [[-6.2, -3.2, 9.82], [-0.7, 2, 0.2], [-0.3, 0.16, -1.3], [0, 0, 0]],
                    [[-1.6, -10.2, 11.7], [0.43, 1.1, 1.0], [-0.4, 0.2, 3.0], [0, 0, 0]],
                    [[-1.8, -8.4, 13], [-0.68, 1.15, 0.613], [-0.3, -0.5, 1.8], [0, 0, 0]],
                    [[-4, -5.54, 13], [-1.0, 2.0, 0.5], [-0.2, -0.5, 2.0], [0, 0, 0]],
                    [[-6.6, -10.3, 13.4], [-0.06, 2.1, 0.4], [1.1, 0.1, 1.4], [0, 0, 0]],
                    [[-6.5, -6.6, 13.3], [0.7, 2.2, 0.3], [1.3, -0.3, 1.2], [0, 0, 0]],
                    [[-5.4, -4.6, 13.2], [0.04, 1.6, 0.1], [0.9, -0.9, 1.8], [0, 0, 0]],
                    [[-1, -9.2, 4.6], [-0.51, 1.22, 2.75], [-0.8, -0.8, 3.5], [0, 0, 0]],
                    [[-3.8, -8.1, 6.9], [-0.6, 1.1, 1.63], [0.4, -0.7, 2.0], [0, 1, 0]],
                    [[-5.6, -5.6, 9], [-0.34, 1.26, 1.3], [0.4, 0.2, 2.0], [0, 1, 0]],
                    [[-5.5, -10.3, 1.5], [-0.73, 1.1, 1.38], [-2.0, -4.0, -8.0], [0, 1, 0]],
                    [[-7.65, -7.8, 3.8], [0.2, 1.05, 1.3], [0.05, -1.2, 2.0], [0, 1, 0]],
                    [[-8, -5.3, 6.8], [1.1, 0.93, 1.35], [0.5, -0.8, 2.7], [0, 1, 0]]])
        } ),
        'Pig 1' : ScaffoldPackage(MeshType_1d_lungpath1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements' : 20,
                'Species' : 'Pig'
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ], [
                    [[-2.4, -11, 1.72], [-0.2, -0.45, 2.5], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-3, -12.7, 7], [0.017, -0.15, 2.20], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-3.9, -12.7, 11.4], [-0.15, 0.8, 2.15], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-3.8, -11.4, 14.9], [-0.067, 2.06, 0.37], [-2.0, 10.0, 22.0], [0, 0, 0]],
                    [[-4, -8.7, 17.5], [0.03, -1.37, 0.85], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-3.3, -5.1, 16], [0.4, -1.04, 0.6], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-4.4, -3.3, 13.1], [0.252, -0.61, 0.62], [-5.0, 4.0, 29.0], [1, 0, 0]],
                    [[-5.7, -2.7, 11.2], [0.45, -0.6, 1.12], [-2.0, 10.0, 22.0], [0, 0, 0]],
                    [[-6.2, -3.2, 9.82], [-0.13, 1.16, 0.55], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-1.6, -10.2, 11.7], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-1.8, -8.4, 13], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-4, -5.54, 13], [-0.68, 1.15, 0.613], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-6.6, -10.3, 13.4], [-0.06, 2.1, 0.3], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-6.5, -6.6, 13.3], [11, -30, -3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-5.4, -4.6, 13.2], [-5.9, -27.8, 0.66], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-1, -9.2, 4.6], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-3.8, -8.1, 6.9], [-0.6, 1.1, 1.63], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-5.6, -5.6, 9], [-0.34, 1.26, 1.3], [-5.0, 4.0, 29.0], [0, 0, 0]],
                    [[-5.5, -10.3, 1.5], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [0, 0, 0]],
                    [[-7.65, -7.8, 3.8], [0.2, 1.05, 1.3], [-10.0, 20.0, 8.0], [0, 0, 0]],
                    [[-8, -5.3, 6.8], [1.1, 0.93, 1.35], [-5.0, 4.0, 29.0], [0, 0, 0]]])
        } ),
        }

    @staticmethod
    def getName():
        return '3DLung from controlcurve'

    @staticmethod
    def getParameterSetNames():
        return [
            'Default',
            'Human 1',
            'Mouse 1',
            'Pig 1']

    @classmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        if 'Human 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        elif 'Mouse 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        elif 'Pig 1' in parameterSetName:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Pig 1']
        else:
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Human 1']
        options = {
            'Central path' : copy.deepcopy(centralPathOption),
            'Number of segments': 20,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        if 'Mouse' in parameterSetName:
            options['Number of segments'] = 10
        elif 'Pig 1' in parameterSetName:
            options['Number of segments'] = 20
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Central path',
            'Number of segments',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Central path':
            return [ MeshType_1d_lungpath1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Central path':
            return list(cls.centralPathDefaultScaffoldPackages.keys())
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
        if optionName == 'Central path':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Central path'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Central path'):
            options['Central path'] = cls.getOptionScaffoldPackage('Central path', MeshType_1d_lungpath1)
        for key in [
            'Number of segments',
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
        centralPath = options['Central path']

        useCrossDerivatives = True
        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        zero = [0,0,0]

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd3 = extractxyzPathParametersFromRegion(tmpRegion)
        # print('extracted central path for lungs')
        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd3[i], '],')
        del tmpRegion

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
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = eftfactory.createEftBasic()

        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        # Create nodes
        # Coordinates field
        nodeIdentifier = 1
        for n in range(len(cx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        ##########################
        # Create regular elements
        ##########################
        elementIdentifier = 1
        nodeIdentifiers = [16, 17, 10, 11, 19, 20, 13, 14]
        element = mesh.createElement(elementIdentifier, elementtemplate)
        element.setNodesByIdentifier(eft, nodeIdentifiers)

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
