"""
Generates a 3-D lung block mesh with variable numbers of elements
up and along
"""

from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils import vector


class MeshType_3d_lungs1(Scaffold_base):
    '''
    classdocs
    '''

    @staticmethod
    def getName():
        return '3D Lungs 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        return {
            'Number of lung elements up': 4,
            'Number of lung elements laterally': 4,
            'Height': 34,
            'Width': 12,
            'Upcurve coefficient1': 1.0,
            'Upcurve coefficient2': 1.0,
            'Medialsurface coefficient': 0.8,
            'Lateralsurface coefficient': 0.8,
            'Tracheal Rotation angle': 0,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements up': 1
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of lung elements up',
            'Number of lung elements laterally',
            'Height',
            'Width',
            'Upcurve coefficient1',
            'Upcurve coefficient2',
            'Medialsurface coefficient',
            'Lateralsurface coefficient',
            'Tracheal Rotation angle',
            'Use cross derivatives',
            'Refine',
            'Refine number of lung elements around',
            'Refine number of lung elements up'
        ]

    @staticmethod
    def checkOptions(options):

        if options['Number of lung elements up'] < 3:
            options['Number of lung elements up'] = 3
        if options['Number of lung elements laterally'] < 3:
            options['Number of lung elements laterally'] = 3


        if options['Upcurve coefficient1'] > 4:
            options['Upcurve coefficient1'] = 4
        if options['Upcurve coefficient2'] > 4:
                options['Upcurve coefficient2'] = 4

        if options['Medialsurface coefficient'] > 0.95:
            options['Medialsurface coefficient'] = 0.95

        if options['Lateralsurface coefficient'] > 0.95:
            options['Lateralsurface coefficient'] = 0.95

        if options['Tracheal Rotation angle'] > 10:
            options['Tracheal Rotation angle'] = 10

        if options['Medialsurface coefficient'] < 0.05:
            options['Medialsurface coefficient'] = 0.05

        if options['Lateralsurface coefficient'] < 0.05:
            options['Lateralsurface coefficient'] = 0.05


    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        lobeid = 1 #Left lobe
        nodeIdentifier = 1
        elementIdentifier = 1

        generateLobeMesh(region,options,lobeid, nodeIdentifier, elementIdentifier)
        return []  # no annotation groups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            return cls.generateBaseMesh(region, options)
        refineElementsCountAround = options['Refine number of lung elements around']
        refineElementsCountAlong = options['Refine number of lung elements up']

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountAlong)
        return meshrefinement.getAnnotationGroups()


def generateLobeMesh(region, options, lobeid, startNodeIdentifier, startElementIdentifier, meshGroups = []):
    '''
    :param vesselMeshGroups: List (over number of vessels) of list of mesh groups to add vessel elements to.
    :return: nextNodeIdentifier, nextElementIdentifier
    :param meshGroups:  Optional list of Zinc MeshGroup for adding new elements to.
    '''

    elementsCountlateral = options['Number of lung elements laterally']
    elementsCountUp = options['Number of lung elements up']
    lungheight = options['Height']
    lungwidth = options['Width']

    posteriorcurveradiuscoeff = options['Upcurve coefficient1']
    anteriorcurveradiuscoeff = options['Upcurve coefficient2']

    medialsurfcoeff = options['Medialsurface coefficient']

    lateralsurfcoeff = options['Lateralsurface coefficient']
    TrachealRotationAngle = options['Tracheal Rotation angle']
    useCrossDerivatives = options['Use cross derivatives']

    fm = region.getFieldmodule()
    fm.beginChange()
    coordinates = findOrCreateFieldCoordinates(fm)

    # -----------------------------------------------
    # hardcoded values (for starters)
    TrachealRotationAngleRadians = -math.pi / 180 * TrachealRotationAngle
    thetaup = math.pi / 3.0
    thetalateral = math.pi / 6.0

    thetaupPosterior = [-math.pi/6, 0, math.pi/5]
    thetaupAnterior = [-math.pi/3, 0, math.pi/8]

    posteriorcurveradius = lungheight / (2 * math.sin(thetaup))
    anteriorcurveradius = 0.6*lungheight / (2 * math.sin(thetaup))

    LMBcentre = 0.7
    RMBcentre = LMBcentre + lungwidth * 0.1

    zlateralcentre = 12.0
    halflungwidth = lungwidth * 0.5
    lungdepth = lungwidth*1.35


    ycentreleftlateral =  halflungwidth

    lateralsurfradius = lateralsurfcoeff * lungwidth

    # -------------end of hardcoded values-----------------

    nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
    nodetemplateApex = nodes.createNodetemplate()
    nodetemplateApex.defineField(coordinates)
    nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
    nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
    nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
    nodetemplateApex.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
    if useCrossDerivatives:
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)
    else:
        nodetemplate = nodetemplateApex

    cache = fm.createFieldcache()

    #################
    # Create nodes
    #################
    nodeIdentifier = startNodeIdentifier
    lRadiansUp = []

    # # Create node location and derivative 2 on left posterior edge upwards - CURVE
    lRadiansUp = []
    lRadiansLateral = []
    x = []
    nx = []
    d1 = []
    nd1 = []
    d2 = []
    nd2 = []

    zero = [0, 0, 0]
    unitvectorz = [0,0,1]
    unitvectorx = [1,0,0]
    unitvectorminusx = [-1,0,0]
    unitvectory = [0,1,0]

    for n2 in range(3):
        radiansUp = -thetaup + (n2) / (2.0) * (2 * thetaup)
        lRadiansUp.append(radiansUp)

        radiansLateral = -thetalateral + (n2) / (2.0) * (2 * thetalateral)
        lRadiansLateral.append(radiansLateral)

    # -----------------------------------------
    ## POSTERIOR edges

    for n2 in range(3):
        radiansUp = thetaupPosterior[n2]
        cosRadiansUp = math.cos(radiansUp)
        sinRadiansUp = math.sin(radiansUp)
        x1 = [LMBcentre, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
              1.2*zlateralcentre + posteriorcurveradius * sinRadiansUp]
        d1 = unitvectorz
        nx.append(x1)
        nd1.append(d1)

    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixStartDerivative=False, fixEndDerivative=False)
    d21 = smoothedd2[0]
    nd2.append(d21)
    d22 = smoothedd2[1]
    nd2.append(d22)
    d23 = smoothedd2[2]
    nd2.append(d23)
    sPosteriorx, sPosteriorderiv2, _, _, _ = \
        interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountUp)

    nx = []
    nd1 = []
    nd2 = []
    for n2 in range(3):
        radiansUp = thetaupPosterior[n2]
        cosRadiansUp = math.cos(radiansUp)
        sinRadiansUp = math.sin(radiansUp)
        x1 = [LMBcentre + 1.5, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
              1.2*zlateralcentre + posteriorcurveradius * sinRadiansUp]
        d1 = unitvectorz
        nx.append(x1)
        nd1.append(d1)

    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixStartDerivative=False, fixEndDerivative=False)
    d21 = smoothedd2[0]
    nd2.append(d21)
    d22 = smoothedd2[1]
    nd2.append(d22)
    d23 = smoothedd2[2]
    nd2.append(d23)
    sPosteriorOffsetx, sPosteriorOffsetderiv2, _, _, _ = \
        interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountUp)

    # -----------------
    ## Anterior edge
    nx = []
    nd1 = []
    nd2 = []
    for n2 in range(3):
        radiansUp = thetaupAnterior[n2]
        cosRadiansUp = math.cos(radiansUp)
        sinRadiansUp = math.sin(radiansUp)
        x1 = [LMBcentre, -ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
              zlateralcentre + anteriorcurveradius * sinRadiansUp]
        d1 = unitvectorz
        nx.append(x1)
        nd1.append(d1)

    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixStartDerivative=False, fixEndDerivative=False)
    d21 = smoothedd2[0]
    nd2.append(d21)
    d22 = smoothedd2[1]
    nd2.append(d22)
    d23 = smoothedd2[2]
    nd2.append(d23)
    sAnteriorx, sAnteriorderiv2, _, _, _ = \
        interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountUp)

    # -------------------------
    ## Anterior edge offset
    ## ------------------------
    nx = []
    nd1 = []
    nd2 = []
    for n2 in range(3):
        radiansUp = thetaupAnterior[n2]
        cosRadiansUp = math.cos(radiansUp)
        sinRadiansUp = math.sin(radiansUp)
        x1 = [LMBcentre + 0.5, -ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
              zlateralcentre + anteriorcurveradius * sinRadiansUp]
        d1 = unitvectorz
        nx.append(x1)
        nd1.append(d1)

    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixStartDerivative=False, fixEndDerivative=False)
    d21 = smoothedd2[0]
    nd2.append(d21)
    d22 = smoothedd2[1]
    nd2.append(d22)
    d23 = smoothedd2[2]
    nd2.append(d23)
    sAnteriorOffsetx, sAnteriorderivOffset2, _, _, _ = \
        interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountUp)

    # # ADD nodes on posterior edge upwards
    # # ------------------------------------
    latUpNodeId = []
    for n2 in range(0, elementsCountUp):
        layerNodeId = []
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sPosteriorx[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sPosteriorderiv2[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        layerNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1
        latUpNodeId.append(layerNodeId)

    # # ADD nodes on anterior edge upwards
    # latUpNodeId = []
    for n2 in range(0, elementsCountUp):
        layerNodeId = []
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sAnteriorx[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sAnteriorderiv2[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        layerNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1
        latUpNodeId.append(layerNodeId)

    # # ADD nodes on posterior edge OFFSET upwards
    # # -------------------------------------------
    latUpNodeId = []
    for n2 in range(0, elementsCountUp):
        layerNodeId = []
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sPosteriorOffsetx[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sPosteriorOffsetderiv2[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        layerNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1
        latUpNodeId.append(layerNodeId)

    # # ADD nodes on anterior edge upwards
    # latUpNodeId = []
    for n2 in range(0, elementsCountUp):
        layerNodeId = []
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sAnteriorOffsetx[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, sAnteriorderivOffset2[n2])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        layerNodeId.append(nodeIdentifier)
        nodeIdentifier = nodeIdentifier + 1
        latUpNodeId.append(layerNodeId)

    # --------------------------------------------
    # # Create nodes on MEDIAL and LATERAL side
    # --------------------------------------------
    xtop = sAnteriorx[0]
    sMedialcoord = []
    sLateralcoord = []

    for n2 in range(0, elementsCountUp):

        #MEDIAL SIDE
        #------------
        nx = []
        nd1 = []
        nd2 = []
        x1 = sPosteriorx[n2]
        x3 = sAnteriorx[n2]
        mediolateral_50 = 0.5*(x1[1] + x3[1])
        mediolateral_40 = 0.6*x1[1] + 0.4*x3[1]
        zlung = 0.5 * (x1[2] + x3[2])
        zratio = -0.1 + 0.8*(abs(zlung-xtop[2])/lungheight)**0.5

        radiansLateral = lRadiansLateral[0]
        cosRadiansLateral = math.cos(radiansLateral)
        sinRadiansLateral = math.sin(radiansLateral)

        d21 = unitvectory
        nx.append(x1)
        nd1.append(d21)

        radiansLateral = lRadiansLateral[1]
        cosRadiansLateral = math.cos(radiansLateral)
        sinRadiansLateral = math.sin(radiansLateral)
        x2 = [x1[0] + 0.05 * lungdepth,
              mediolateral_50, 0.9*zlung]

        d22 = unitvectory
        nx.append(x2)
        nd1.append(d22)

        radiansLateral = lRadiansLateral[2]
        cosRadiansLateral = math.cos(radiansLateral)
        sinRadiansLateral = math.sin(radiansLateral)

        d23 = unitvectory
        nx.append(x3)
        nd1.append(d23)

        smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixStartDerivative=False, fixEndDerivative=False)

        d21 = smoothedd2[0]
        nd2.append(d21)
        d22 = smoothedd2[1]
        nd2.append(d22)
        d23 = smoothedd2[2]
        nd2.append(d23)
        sMedialx, sMedialderiv2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)
        sMedialcoord.append(sMedialx)

        #LATERAL SIDE
        #---------------
        nx = []
        nd1 = []
        nd2 = []

        x1 = sPosteriorOffsetx[n2]
        x3 = sAnteriorOffsetx[n2]
        mediolateral_50 = 0.5 * (x1[1] + x3[1])
        mediolateral_40 = 0.6*x1[1] + 0.4*x3[1]
        mediolateral_30 = 0.7*x1[1] + 0.3*x3[1]
        radiansLateral = lRadiansLateral[0]
        cosRadiansLateral = math.cos(radiansLateral)
        sinRadiansLateral = math.sin(radiansLateral)

        d21 = unitvectory
        nx.append(x1)
        nd1.append(d21)

        radiansLateral = lRadiansLateral[1]
        cosRadiansLateral = math.cos(radiansLateral)
        sinRadiansLateral = math.sin(radiansLateral)
        x2 = [x1[0] + (zratio+0.05)*lungdepth,
              mediolateral_50,
              zlung]
        d22 = unitvectory
        nx.append(x2)
        nd1.append(d22)

        radiansLateral = lRadiansLateral[2]
        cosRadiansLateral = math.cos(radiansLateral)
        sinRadiansLateral = math.sin(radiansLateral)

        d23 = unitvectory
        nx.append(x3)
        nd1.append(d23)

        smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd1, fixStartDerivative=False, fixEndDerivative=False)
        d21 = smoothedd2[0]
        nd2.append(d21)
        d22 = smoothedd2[1]
        nd2.append(d22)
        d23 = smoothedd2[2]
        nd2.append(d23)
        sLateralx, sLateralderiv2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)
        sLateralcoord.append(sLateralx)

    # for n2 in range(1, elementsCountUp):
        for n3 in range(1, elementsCountlateral):
            layerNodeId = []
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sMedialx[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sMedialderiv2[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            layerNodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
            latUpNodeId.append(layerNodeId)

        for n3 in range(1, elementsCountlateral):
            layerNodeId = []
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, sLateralx[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, sLateralderiv2[n3])
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, zero)
            if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
            layerNodeId.append(nodeIdentifier)
            nodeIdentifier = nodeIdentifier + 1
            latUpNodeId.append(layerNodeId)

    ###################
    # Create elements
    ###################
    e0 = startNodeIdentifier
    elementIdentifier = startElementIdentifier

    mesh = fm.findMeshByDimension(3)
    tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    eft = tricubichermite.createEftBasic()

    tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

    # Regular elements
    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
    elementtemplate.defineField(coordinates, -1, eft)

    delQ = 4 * elementsCountUp
    delX = 2 * elementsCountlateral - 2
    delA = 4 * elementsCountUp + (elementsCountlateral - 1)
    delB = 4 * elementsCountUp + (3 * elementsCountlateral - 3)
    delM = 2 * elementsCountUp
    delP = 4 * elementsCountUp + (elementsCountlateral - 2)

    for e2 in range(1, elementsCountUp):
        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [e0 + e2 -1,
                           e0 + delQ + (delX) * (e2 - 1),
                           e0 + e2,
                           e0 + delQ + (delX) * (e2 - 1) + delX,
                           e0 + e2-1 + delM,
                           e0 + delA + (delX) * (e2 - 1),
                           e0 + e2 + delM,
                           e0 + delB + (delX) * (e2 - 1)
                           ]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        for meshGroup in meshGroups:
            meshGroup.addElement(element)

    for e3 in range(1, elementsCountlateral - 1):
        for e2 in range(1, elementsCountUp):
            e22 = e0 + delQ + (delX) * (e2 - 1) + (e3 - 1)
            element = mesh.createElement(elementIdentifier, elementtemplate)
            nodeIdentifiers = [e22,
                               e22 + 1,
                               e22 + delX,
                               e22 + delX + 1,
                               e22 + (elementsCountlateral - 1),
                               e22 + (elementsCountlateral),
                               e22 + (3 * elementsCountlateral - 3),
                               e22 + (3 * elementsCountlateral - 3) + 1
                               ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1
            for meshGroup in meshGroups:
                meshGroup.addElement(element)

    for e2 in range(1, elementsCountUp):
        element = mesh.createElement(elementIdentifier, elementtemplate)
        nodeIdentifiers = [e0 + delP + (delX) * (e2 - 1),
                           e0 + (e2-1) + elementsCountUp,
                           e0 + delP + (delX) * (e2),
                           e0 + (e2) + elementsCountUp,
                           e0 + delP + (delX) * (e2 - 1) + (elementsCountlateral - 1),
                           e0 + e2-1 + (3 * elementsCountUp),
                           e0 + delP + (delX) * (e2) + (elementsCountlateral - 1),
                           e0 + e2 + (3 * elementsCountUp)
                           ]
        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
        elementIdentifier = elementIdentifier + 1
        for meshGroup in meshGroups:
            meshGroup.addElement(element)

    #coordinates.smooth(fm.createFieldsmoothing())

    fm.endChange()
    return nodeIdentifier, elementIdentifier

