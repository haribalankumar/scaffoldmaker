"""
Generates a 3-D lung mesh along the control lines
"""

import copy
import math
from scaffoldmaker.meshtypes.meshtype_1d_lungpath1 import MeshType_1d_lungpath1, extractxyzPathParametersFromRegion
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.zinc_utils import exnodeStringFromNodeValues
from opencmiss.utils.zinc.general import ChangeManager
from opencmiss.zinc.context import Context
from opencmiss.zinc.fieldmodule import Fieldmodule

from opencmiss.zinc.node import Node
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds


class MeshType_3d_lungcontrolcurves(Scaffold_base):
    '''
    Generates a 3-D lung mesh
    '''

    centralPathDefaultScaffoldPackages = {
        'Human 1' : ScaffoldPackage(MeshType_1d_lungpath1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements per curve' : 2,
                'Species' : 'Human'
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3], [
                    [[-2.4, -11, 1.72], [-2.0, -4.0, -8.0], [-0.2, -0.45, 2.5], [1, 1, 0]],
                    [[-3.9, -12.7, 11.4], [-5.0, 4.0, 29.0], [-0.15, 0.8, 2.15], [1, 1, 0]],
                    [[-4, -8.7, 17.5], [-22.0, -4.0, -8.0], [0.03, -1.37, 0.85], [1, 0, 0]],
                    [[-4.4, -3.3, 13.1],  [-5.0, 4.0, 29.0], [0.252, -0.61, 0.62],[1, 0, 0]],
                    [[-6.2, -3.2, 9.82],  [-22.0, -4.0, -8.0], [-0.13, 1.16, 0.55], [1, 0, 0]],
                    [[-1, -9.2, 4.6], [-0.51, 1.22, 2.75], [-22.0, -4.0, -8.0], [1, 0, 0]],
                    [[-5.5, -10.3, 1.5], [-0.73, 1.1, 1.38], [-22.0, -4.0, -8.0], [1, 0, 0]],
                    [[-1.6, -10.2, 11.7], [-0.46, 0.86, 0.5], [-10.0, 20.0, 8.0], [1, 0, 0]],
                    [[-6.6, -10.3, 13.4], [21.2, -8.1, 0.3], [-22.0, -4.0, -8.0], [1, 0, 0]]
                    ])
        } ),
        'Mouse 1' : ScaffoldPackage(MeshType_1d_lungpath1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements per curve' : 2,
                'Species' : 'Mouse'
                },
            'meshEdits' : exnodeStringFromNodeValues(
                [ Node.VALUE_LABEL_VALUE, Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3 ], [
                    [[-2.4, -11, 1.72], [2.8, 1.17, 2.1], [-0.1, -0.16, 2.1], [3, 0.7, 2]],
                    [[-3.9, -12.7, 11.4], [1.8, 1.0, 0.2], [-0.15, 0.8, 2.15], [1, 0, 0]],
                    [[-3.5, -8.8, 17], [-3.0, 2.0, 1.2], [0.18, 2.9, 1.07], [1, 0, 0]],
                    [[-4.4, -3.3, 13.1], [-0.08, 1.1, 0.2], [0.12, -0.61, 1.2], [1, 0, 0]],
                    [[-6.2, -3.2, 9.82], [-0.7, 2, 0.2], [0.3, -0.16, 1.3], [1, 0.8, 0.54]],
                    [[-1, -9.2, 4.6], [-1.01, 0.83, 2.45], [-0.8, -0.8, 3.5], [1.1, 0.06, 1]],
                    [[-6, -10.3, 1.6], [-2.4, 1.6, 1.0], [-0.07, -2.0, 2.7], [1.4, 0.41, 1.2]],
                    [[-1.6, -10.2, 11.7], [0.43, 1.1, 1.0], [-1.1, 0.42, 2.5], [1, 0, 0]],
                    [[-6.6, -11, 12.5], [-0.7, 2.5, 0.74], [0.7, -0.1, 1.8], [1, 0, 0]]])
        } ),
        'Pig 1' : ScaffoldPackage(MeshType_1d_lungpath1, {
            'scaffoldSettings' : {
                'Coordinate dimensions' : 3,
                'Length' : 1.0,
                'Number of elements per curve' : 4,
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
            centralPathOption = cls.centralPathDefaultScaffoldPackages['Mouse 1']
        options = {
            'Control curves' : copy.deepcopy(centralPathOption),
            'Species' : 'Mouse',
            'Number of segments': 20,
            'Number of elements along': 4,
            'Refine' : False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
            }
        if 'Mouse' in parameterSetName:
            options['Number of segments'] = 20
            options['Species'] = 'Mouse'
        elif 'Pig 1' in parameterSetName:
            options['Number of segments'] = 20
            options['Species'] = 'Pig'
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Control curves',
            'Species',
            'Number of segments',
            'Number of elements along',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall' ]

    @classmethod
    def getOptionValidScaffoldTypes(cls, optionName):
        if optionName == 'Control curves':
            return [ MeshType_1d_lungpath1 ]
        return []

    @classmethod
    def getOptionScaffoldTypeParameterSetNames(cls, optionName, scaffoldType):
        if optionName == 'Control curves':
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
        if optionName == 'Control curves':
            if not parameterSetName:
                parameterSetName = list(cls.centralPathDefaultScaffoldPackages.keys())[0]
            return copy.deepcopy(cls.centralPathDefaultScaffoldPackages[parameterSetName])
        assert False, cls.__name__ + '.getOptionScaffoldPackage:  Option ' + optionName + ' is not a scaffold'

    @classmethod
    def checkOptions(cls, options):
        if not options['Control curves'].getScaffoldType() in cls.getOptionValidScaffoldTypes('Control curves'):
            options['Control curves'] = cls.getOptionScaffoldPackage('Control curves', MeshType_1d_lungpath1)
        for key in [
            'Number of segments',
            'Number of elements along',
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
        centralPath = options['Control curves']
        speciestype = options['Species']
        # elementsCountcurve = options['Number of elements per curve']

        elementsCountAlong = options['Number of elements along']
        useCrossDerivatives = True
        firstNodeIdentifier = 1
        firstElementIdentifier = 1

        zero = [0,0,0]

        # Central path
        tmpRegion = region.createRegion()
        centralPath.generate(tmpRegion)
        cx, cd1, cd2, cd3 = extractxyzPathParametersFromRegion(tmpRegion)
        del tmpRegion

        # for i in range(len(cx)):
        #     print(i, '[', cx[i], ',', cd1[i], ',', cd2[i], ',', cd3[i], '],')

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)
        cache = fm.createFieldcache()

        # Coordinates field
        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)

        if speciestype == 'Mouse':  #6 curves for left lobe
        ###-----------------------
            ##Assign nodes to groups
            accessoryedgecx = []
            accessoryedgecd1 = []
            accessoryedgecd2 = []
            accessoryedgecd3 = []

            posteriorlateralcx = []
            posteriorlateralcd1 = []
            posteriorlateralcd2 = []
            posteriorlateralcd3 = []

            basemedialcx = []
            basemedialcd1 = []
            basemedialcd2 = []
            basemedialcd3 = []

            baselateralcx = []
            baselateralcd1 = []
            baselateralcd2 = []
            baselateralcd3 = []

            midmedialcx = []
            midmedialcd1 = []
            midmedialcd2 = []
            midmedialcd3 = []

            midlateralcx = []
            midlateralcd1 = []
            midlateralcd2 = []
            midlateralcd3 = []

            apicaledgecx = []
            apicaledgecd1 = []
            apicaledgecd2 = []

            midmedialapicalcx = []
            midmedialapicalcd1 = []
            midmedialapicalcd2 = []
            midmedialapicalcd3 = []

            midlateralapicalcx = []
            midlateralapicalcd1 = []
            midlateralapicalcd2 = []
            midlateralapicalcd3 = []

            # accessory edge
            accessoryedgecx.append(cx[1])
            accessoryedgecd1.append(cd1[1])
            accessoryedgecx.append(cx[7])
            accessoryedgecd1.append(cd1[7])
            accessoryedgecx.append(cx[3])
            accessoryedgecd1.append(cd1[3])
            if(elementsCountAlong>2):
                accessoryedgecx, accessoryedgecd1, _, _, _ = \
                    interp.sampleCubicHermiteCurves(accessoryedgecx, accessoryedgecd1, elementsCountOut=elementsCountAlong)
                accessoryedgecd1 = interp.smoothCubicHermiteDerivativesLine(accessoryedgecx, accessoryedgecd1,
                                                fixAllDirections=False,
                                                fixStartDerivative=True, fixEndDerivative=True,
                                                magnitudeScalingMode = interp.DerivativeScalingMode.HARMONIC_MEAN)
            for n in range(elementsCountAlong):
                if(n==0):
                    accessoryedgecd2.append(cd2[1])
                    accessoryedgecd3.append(cd3[1])
                if (n == elementsCountAlong - 1):
                    accessoryedgecd2.append(cd2[3])
                    accessoryedgecd3.append(cd3[3])
                else:
                    accessoryedgecd2.append(cd2[7])
                    accessoryedgecd3.append(cd3[7])

            # Posterior Lateral
            #------------------
            tempcx=[]
            tempcd1=[]
            tempcx.append(cx[1])
            tempcd1.append(cd1[1])
            tempcx.append(cx[8])
            tempcd1.append(cd1[8])
            tempcx.append(cx[3])
            tempcd1.append(cd1[3])
            if(elementsCountAlong>2):
                posteriorlateralcx, posteriorlateralcd1, _, _, _ = \
                   interp.sampleCubicHermiteCurves(tempcx, tempcd1, elementsCountOut=elementsCountAlong)
                posteriorlateralcd1 = interp.smoothCubicHermiteDerivativesLine(posteriorlateralcx, posteriorlateralcd1,
                                                                    fixAllDirections=False,
                                                                    fixStartDerivative=True, fixEndDerivative=True,
                                                                    magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
            for n in range(elementsCountAlong):
                if(n==0):
                    posteriorlateralcd2.append(cd2[1])
                    posteriorlateralcd3.append(cd3[1])
                if (n == elementsCountAlong - 1):
                    posteriorlateralcd2.append(cd2[3])
                    posteriorlateralcd3.append(cd3[3])
                else:
                    posteriorlateralcd2.append(cd2[8])
                    posteriorlateralcd3.append(cd3[8])


            ## base medial
            ## -----------
            tempcx=[]
            tempcd1=[]
            tempcx.append(cx[0])
            tempcd1.append(cd1[0])
            tempcx.append(cx[5])
            tempcd1.append(cd1[5])
            tempcx.append(cx[4])
            tempcd1.append(cd1[4])
            if(elementsCountAlong>2):
                basemedialcx, basemedialcd1, _, _, _ = \
                interp.sampleCubicHermiteCurves(tempcx, tempcd1,
                                                elementsCountOut=elementsCountAlong)
                basemedialcd1 = interp.smoothCubicHermiteDerivativesLine(basemedialcx, basemedialcd1,
                                                                       fixAllDirections=False,
                                                                       fixStartDerivative=True, fixEndDerivative=True,
                                                                       magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
            for n in range(elementsCountAlong):
                if(n==0):
                    basemedialcd2.append(cd2[0])
                    basemedialcd3.append(cd3[0])
                if (n == elementsCountAlong - 1):
                    basemedialcd2.append(cd2[4])
                    basemedialcd3.append(cd3[4])
                else:
                    basemedialcd2.append(cd2[5])
                    basemedialcd3.append(cd3[5])


            ### base lateral
            ### ------------
            tempcx=[]
            tempcd1=[]
            tempcx.append(cx[0])
            tempcd1.append(cd1[0])
            tempcx.append(cx[6])
            tempcd1.append(cd1[6])
            tempcx.append(cx[4])
            tempcd1.append(cd1[4])
            if(elementsCountAlong>2):
                baselateralcx, baselateralcd1, _, _, _ = interp.sampleCubicHermiteCurves(tempcx, tempcd1,
                                                            elementsCountOut=elementsCountAlong)
                baselateralcd1 = interp.smoothCubicHermiteDerivativesLine(baselateralcx, baselateralcd1,
                                                                 fixAllDirections=False,
                                                                 fixStartDerivative=True, fixEndDerivative=True,
                                                                 magnitudeScalingMode=interp.DerivativeScalingMode.HARMONIC_MEAN)
                for n in range(elementsCountAlong):
                    if (n == 0):
                        baselateralcd2.append(cd2[0])
                        baselateralcd3.append(cd3[0])
                    if (n == elementsCountAlong - 1):
                        baselateralcd2.append(cd2[4])
                        baselateralcd3.append(cd3[4])
                    else:
                        baselateralcd2.append(cd2[6])
                        baselateralcd3.append(cd3[6])

            ### Apical
            ### ----------
            if(elementsCountAlong<4):
                apicaledgecx.append(cx[2])
                apicaledgecd1.append(cd1[2])
                apicaledgecd2.append(cd2[2])
            else:
                tempcx=[]
                tempcd1=[]
                tempcx.append(cx[0])
                tempcd1.append(cd2[0])
                tempcx.append(cx[1])
                tempcd1.append(cd2[1])
                tempcx.append(cx[2])
                tempcd1.append(cd2[2])
                tempapicalcx, tempapicalcd2, _, _, _ = interp.sampleCubicHermiteCurves(tempcx, tempcd1,
                                                               elementsCountOut=elementsCountAlong)

                apicaledgecx.append(tempapicalcx[3])
                apicaledgecd2.append(tempapicalcd2[3])
                apicaledgecd1.append(cd1[1])
                apicaledgecx.append(cx[2])
                apicaledgecd2.append(cd2[2])
                apicaledgecd1.append(cd1[2])

                tempcx=[]
                tempcd1=[]
                tempcx.append(cx[2])
                tempcd1.append(cd2[2])
                tempcx.append(cx[3])
                tempcd1.append(cd2[3])
                tempcx.append(cx[4])
                tempcd1.append(cd2[4])
                tempapicalcx, tempapicalcd2, _, _, _ = interp.sampleCubicHermiteCurves(tempcx, tempcd1,
                                                               elementsCountOut=elementsCountAlong)
                apicaledgecx.append(tempapicalcx[1])
                apicaledgecd2.append(tempapicalcd2[1])
                apicaledgecd1.append(cd1[3])

            # Create additional nodes
            ###########################
            if(elementsCountAlong>2):
                temp = []
                for i in range(int(elementsCountAlong/2)//2):
                    for n in range(elementsCountAlong+1):
                        xfrac = 1.0/float(i+2)
                        temp = [(accessoryedgecx[n][c]+basemedialcx[n][c])*xfrac for c in range(3)]
                        midmedialcx.append(temp)
                        midmedialcd1.append(basemedialcd1[n])
                        midmedialcd2.append(basemedialcd2[n])
                        midmedialcd3.append(basemedialcd3[n])

                        temp = [(posteriorlateralcx[n][c]+baselateralcx[n][c])*xfrac for c in range(3)]
                        midlateralcx.append(temp)
                        midlateralcd1.append(baselateralcd1[n])
                        midlateralcd2.append(baselateralcd2[n])
                        midlateralcd3.append(baselateralcd3[n])

            # if(elementsCountAlong>3):
            # temp = []
            # for n in range(elementsCountAlong-1):
            #     temp = [(cx[n+9][c]+cx[n+3][c])*0.5 for c in range(3)]
            #     if(n==0):
            #         temp = [(cx[n+9][c]+cx[n+3][c]+cx[n+4][c]+cx[n+10][c])*0.25 for c in range(3)]
            #     if (n == elementsCountAlong-2):
            #         temp = [(cx[n+9][c] + cx[n+3][c] + cx[n+2][c] + cx[n+8][c])*0.25 for c in range(3)]
            #     midmedialapicalcx.append(temp)
            #     midmedialapicalcd1.append(accessoryedgecd1[n+1])
            #     midmedialapicalcd2.append(accessoryedgecd2[n+1])
            #     midmedialapicalcd3.append(accessoryedgecd3[n+1])
            #
            #     temp = [(cx[n+12][c]+cx[n+3][c])/2 for c in range(3)]
            #     midlateralapicalcx.append(temp)
            #     midlateralapicalcd1.append(posteriorlateralcd1[n+1])
            #     midlateralapicalcd2.append(posteriorlateralcd2[n+1])
            #     midlateralapicalcd3.append(posteriorlateralcd3[n+1])


        # # Create nodes
        # ##################
        # # Coordinates field
        # nodeIdentifier = 1
        # for n in range(len(cx)):
        #     node = nodes.createNode(nodeIdentifier, nodetemplate)
        #     cache.setNode(node)
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, cx[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, cd1[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, cd2[n])
        #     coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, cd3[n])
        #     if useCrossDerivatives:
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
        #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        #     nodeIdentifier = nodeIdentifier + 1


        # Create nodes
        ##################
        # Coordinates field
        nodeIdentifier = 1
        for n in range(len(basemedialcx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, basemedialcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, basemedialcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, basemedialcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, basemedialcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(1, len(basemedialcx)-1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, baselateralcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, baselateralcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, baselateralcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, baselateralcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(len(midmedialcx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, midmedialcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, midmedialcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, midmedialcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, midmedialcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(1, len(midlateralcx)-1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, midlateralcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, midlateralcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, midlateralcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, midlateralcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(len(accessoryedgecx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, accessoryedgecx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, accessoryedgecd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, accessoryedgecd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, accessoryedgecd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(1, len(accessoryedgecx)-1):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, posteriorlateralcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, posteriorlateralcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, posteriorlateralcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, posteriorlateralcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(len(midmedialapicalcx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, midmedialapicalcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, midmedialapicalcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, midmedialapicalcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, midmedialapicalcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(len(midlateralapicalcx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, midlateralapicalcx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, midlateralapicalcd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, midlateralapicalcd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, midlateralapicalcd3[n])
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1

        for n in range(len(apicaledgecx)):
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, apicaledgecx[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, apicaledgecd1[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, apicaledgecd2[n])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            nodeIdentifier = nodeIdentifier + 1
        #
        # ##########################
        # # Create elements
        # ##########################
        # elementIdentifier = 1
        #
        # eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        # eft = eftfactory.createEftBasic()
        #
        # elementtemplate = mesh.createElementtemplate()
        # elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        # elementtemplate.defineField(coordinates, -1, eft)
        # elementtemplateX = mesh.createElementtemplate()
        # elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        #
        # #elements between accessory edge and base(medial/lateral)
        # for n in range(elementsCountAlong*int(elementsCountAlong/2)):
        #     if((n)%4==0):  #wedge elements xi3zero
        #         va = n%4
        #         vb = (n%4 + 1)%elementsCountAlong
        #         eft1 = eftfactory.createEftWedgeXi3Zero(va*100, vb*100)
        #         # setEftScaleFactorIds(eft1, [1], [])
        #         elementtemplateX.defineField(coordinates, -1, eft1)
        #         # nodeIdentifiers = [1, 6, 9, 14, 1, 2, 9, 10]
        #         bni1 = 1 + (2*elementsCountAlong)*(n//4)
        #         bni2 = bni1 + elementsCountAlong + 1
        #         bni3 = bni1 + 2*elementsCountAlong
        #         bni4 = bni3 + (elementsCountAlong + 1)
        #         # nodeIdentifiers = [1, 6, 9, 14, 2, 10]
        #         nodeIdentifiers = [bni1, bni2, bni3, bni4, bni1+1, bni3+1]
        #         print('nodes xi3=0:', bni1, bni2, bni3, bni4, bni1+1, bni3+1)
        #         element = mesh.createElement(elementIdentifier, elementtemplateX)
        #         result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
        #         # result2 = element.setScaleFactors(eft1, [-1])
        #         elementIdentifier = elementIdentifier + 1
        #     elif (n>0 and (n+1)%4==0): #wedge elements xi3One
        #         va = n%4
        #         vb = (n%4 + 1)%elementsCountAlong
        #         eft2 = eftfactory.createEftWedgeXi3One(va*100, vb*100)
        #         # setEftScaleFactorIds(eft2, [1], [])
        #         elementtemplateX.defineField(coordinates, -1, eft2)
        #
        #         bni1 = (2*elementsCountAlong)*((n+1)//4)
        #         bni2 = (elementsCountAlong+1)+(2*elementsCountAlong)*((n)//4)
        #         bni3 = bni1 + 2*elementsCountAlong
        #         bni4 = bni1 + (elementsCountAlong + 1)
        #         # nodeIdentifiers = [1, 6, 9, 14, 2, 10]
        #         nodeIdentifiers = [bni1, bni2, bni3, bni4, bni2-1, bni4-1]
        #         print('xi3=1 elem',bni1, bni2, bni3, bni4, bni2-1, bni4-1)
        #         element = mesh.createElement(elementIdentifier, elementtemplateX)
        #         result = element.setNodesByIdentifier(eft2, nodeIdentifiers)
        #         # result2 = element.setScaleFactors(eft2, [-1])
        #         elementIdentifier = elementIdentifier + 1
        #     else:
        #         eft = eftfactory.createEftBasic()
        #         bni1 = n + (elementsCountAlong)*(n//elementsCountAlong+1)+1
        #         bni2 = bni1 + 2*(elementsCountAlong)
        #         bni3 = bni1 - elementsCountAlong
        #         bni4 = bni1 + (elementsCountAlong)
        #         # nodeIdentifiers = [6, 7, 14, 15, 2, 3, 10, 11]
        #         nodeIdentifiers = [bni1,bni1+1,bni2,bni2+1,bni3,bni3+1,bni4,bni4+1]
        #         print('normal elems=',bni1,bni1+1,bni2,bni2+1,bni3,bni3+1,bni4,bni4+1)
        #         element = mesh.createElement(elementIdentifier, elementtemplate)
        #         result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        #         # result2 = element.setScaleFactors(eft, [-1])
        #         elementIdentifier = elementIdentifier + 1
        #
        #
        #     # nodeIdentifiers = [7, 8, 15, 16, 3, 4, 11, 12]
        #     # element = mesh.createElement(elementIdentifier, elementtemplate)
        #     # result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        #     # elementIdentifier = elementIdentifier + 1
        #
        #     # if(n>0 and (n+1)%4==0):  #wedge elements
        #     #     eft1 = eftfactory.createEftNoCrossDerivatives()
        #     #     setEftScaleFactorIds(eft1, [1], [])
        #     #     # nodeIdentifiers = [1, 6, 9, 14, 1, 2, 9, 10]
        #     #     nodeIdentifiers = [1, 6, 9, 14, 2, 10]
        #     #     elementtemplateX.defineField(coordinates, -1, eft1)
        #     #     element = mesh.createElement(elementIdentifier, elementtemplateX)
        #     #     result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
        #     #     result2 = element.setScaleFactors(eft1, [-1])
        #     #     elementIdentifier = elementIdentifier + 1
        #
        # #elements between accessory edge and Apex
        # eft1 = eftfactory.createEftWedgeXi3Zero(1 * 100, 2 * 100)
        # # setEftScaleFactorIds(eft1, [1], [])
        # elementtemplateX.defineField(coordinates, -1, eft1)
        # nodeIdentifiers = [17, 22, 31, 28, 18, 25]
        # element = mesh.createElement(elementIdentifier, elementtemplateX)
        # result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
        # # result2 = element.setScaleFactors(eft1, [-1])
        # elementIdentifier = elementIdentifier + 1
        #
        # eft = eftfactory.createEftBasic()
        # nodeIdentifiers = [22, 23, 28, 29, 18, 19, 25, 26]
        # element = mesh.createElement(elementIdentifier, elementtemplate)
        # result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        # elementIdentifier = elementIdentifier + 1
        #
        # eft = eftfactory.createEftBasic()
        # nodeIdentifiers = [23, 24, 29, 30, 19, 20, 26, 27]
        # element = mesh.createElement(elementIdentifier, elementtemplate)
        # result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        # elementIdentifier = elementIdentifier + 1
        #
        # eft2 = eftfactory.createEftWedgeXi3One(0 * 100, 1 * 100)
        # # setEftScaleFactorIds(eft2, [1], [])
        # elementtemplateX.defineField(coordinates, -1, eft2)
        # nodeIdentifiers = [24, 21, 30, 33, 20, 27]
        # element = mesh.createElement(elementIdentifier, elementtemplateX)
        # result = element.setNodesByIdentifier(eft2, nodeIdentifiers)
        # # result2 = element.setScaleFactors(eft2, [-1])
        # elementIdentifier = elementIdentifier + 1
        #
        # eft1 = eftfactory.createEftWedgeXi3One(1 * 100, 2 * 100)
        # # setEftScaleFactorIds(eft1, [1], [])
        # for n in range(2):
        #     d2Map = (0,-1,0)
        #     remapEftNodeValueLabel(eft1, [(n+1)*2], Node.VALUE_LABEL_D_DS1,  derivativeSignsToExpressionTerms(
        #                                (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3),
        #                                d2Map))
        #     remapEftNodeValueLabel(eft1, [(n+1)*4], Node.VALUE_LABEL_D_DS1, [(Node.VALUE_LABEL_D_DS1, [])])
        # elementtemplateX.defineField(coordinates, -1, eft1)
        # # nodeIdentifiers = [31,28,32,29,25,26]
        # nodeIdentifiers = [28, 29, 31, 32, 25, 26]
        # element = mesh.createElement(elementIdentifier, elementtemplateX)
        # result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
        # # result2 = element.setScaleFactors(eft1, [-1])
        # elementIdentifier = elementIdentifier + 1
        #
        # eft1 = eftfactory.createEftWedgeXi3One(1 * 100, 2 * 100)
        # # setEftScaleFactorIds(eft1, [1], [])
        # elementtemplateX.defineField(coordinates, -1, eft1)
        # # nodeIdentifiers = [33,30,32,29,27,26]
        # nodeIdentifiers = [29, 30, 32, 33, 26, 27]
        # element = mesh.createElement(elementIdentifier, elementtemplateX)
        # result = element.setNodesByIdentifier(eft1, nodeIdentifiers)
        # # result2 = element.setScaleFactors(eft1, [-1])
        # elementIdentifier = elementIdentifier + 1

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


def derivativeSignsToExpressionTerms(valueLabels, signs):
    '''
    Return remap expression terms for summing derivative[i]*sign[i]
    :param valueLabels: List of node value labels to possibly include.
    :param signs: List of 1 (no scaling), -1 (scale by scale factor 1) or 0 (no term).
    '''
    expressionTerms = []
    for i in range(len(valueLabels)):
        if signs[i] is 1:
            expressionTerms.append((valueLabels[i], []))
        elif signs[i] is -1:
            expressionTerms.append((valueLabels[i], [1]))
    return expressionTerms


# [[-2.4, -11, 1.72], [2.8, 1.17, 2.1], [-0.1, -0.16, 2.1], [3, 0.7, 2]],
# [[-3, -12.7, 7], [1.0, 1.7, 1.3], [-0.2, -0.2, 1.70], [1, 0, 0]],
# [[-3.9, -12.7, 11.4], [1.8, 1.0, 0.2], [-0.15, 0.8, 2.15], [1, 0, 0]],
# [[-3.8, -11.4, 15], [1.1, 1.8, 0.1], [-0.2, 1.1, 2.5], [1, 0, 0]],
# [[-3.5, -8.8, 17], [-3.0, 2.0, 1.2], [0.18, 2.9, 1.07], [1, 0, 0]],
# [[-3.5, -5.1, 16.1], [-1.0, 2.0, 1.0], [0.6, -2.2, 1.6], [1, 0, 0]],
# [[-4.4, -3.3, 13.1], [-0.08, 1.1, 0.2], [0.12, -0.61, 1.2], [1, 0, 0]],
# [[-5.7, -2.7, 11.2], [0.1, 1.6, 0.6], [0.45, 0.01, 1.12], [1, 0, 0]],
# [[-6.2, -3.2, 9.82], [-0.7, 2, 0.2], [0.3, -0.16, 1.3], [1, 0.8, 0.54]],
# [[-1.6, -10.2, 11.7], [0.43, 1.1, 1.0], [-1.1, 0.42, 2.5], [1, 0, 0]],
# [[-1.8, -8.3, 12.5], [-0.68, 1.15, 0.613], [-1.0, -0.3, 1.4], [1, 0, 0]],
# [[-4, -5.5, 13], [-1.0, 2.0, 0.5], [-0.2, -0.5, 2.0], [1, 0, 0]],
# [[-6.6, -11, 12.5], [-0.7, 2.5, 0.74], [0.7, -0.1, 1.8], [1, 0, 0]],
# [[-6.6, -7.7, 13.3], [0.7, 2.2, 0.3], [0.7, -0.4, 1.3], [1, 0, 0]],
# [[-5.4, -4.6, 13.2], [1.1, 1.6, -0.15], [0.9, -0.9, 1.8], [1, 0, 0]],
# [[-1, -9.2, 4.6], [-1.01, 0.83, 2.45], [-0.8, -0.8, 3.5], [1.1, 0.06, 1]],
# [[-3.8, -8.1, 6.9], [-1.16, 1.0, 1.7], [0.4, -0.7, 2.0], [1.4, 0.41, 1]],
# [[-5.6, -5.6, 9], [-0.38, 1.8, 0.7], [0.4, 0.2, 2.0], [1.4, 0.53, 1]],
# [[-6, -10.3, 1.6], [-2.4, 1.6, 1.0], [-0.07, -2.0, 2.7], [1.4, 0.41, 1.2]],
# [[-8, -7.1, 4.4], [-0.74, 1.35, 1.9], [0.05, -1.2, 2.0], [1.8, -0.2, 1]],
# [[-7.6, -4.53, 7.7], [0.86, 1.3, 2.7], [0.5, -0.8, 2.7], [1.32, -0.8, 0.83]]])
