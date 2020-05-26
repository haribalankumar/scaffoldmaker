"""
Generates a 3-D unit solid sphere mesh with variable numbers of elements
around, up the central axis, and radially.
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
            'Number of elements up': 4,
            'Number of elements laterally': 4,
            'Height': 10,
            'Width': 10,
            'Upcurve coefficient1': 1.0,
            'Upcurve coefficient2': 1.0,
            'Medialsurface coefficient': 0.2,
            'Lateralsurface coefficient': 0.2,
            'Tracheal Rotation angle': 0,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around': 1,
            'Refine number of elements up': 1,
        }

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements up',
            'Number of elements laterally',
            'Height',
            'Width',
            'Upcurve coefficient1',
            'Upcurve coefficient2',
            'Medialsurface coefficient',
            'Lateralsurface coefficient',
            'Tracheal Rotation angle',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements up',
        ]

    @staticmethod
    def checkOptions(options):

        if options['Number of elements up'] < 3:
            options['Number of elements up'] = 3
        if options['Number of elements laterally'] < 3:
            options['Number of elements laterally'] = 3

        if options['Width'] < 2:
            options['Width'] = 2
        if options['Height'] < 2:
            options['Height'] = 2

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


    @classmethod
    def generateBaseMesh(cls, region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """
        elementsCountlateral = options['Number of elements laterally']
        elementsCountUp = options['Number of elements up']
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

        #-----------------------------------------------
        # hardcoded values (for starters)
        TrachealRotationAngleRadians = -math.pi/180*TrachealRotationAngle
        thetaup =  math.pi/3.0
        thetalateral = math.pi/6.0
        posteriorcurveradius = lungheight/(2*math.sin(thetaup))
        anteriorcurveradius = posteriorcurveradius

        LMBcentre = 0.0
        RMBcentre = LMBcentre + lungwidth*0.1

        zlateralcentre = 5.0
        ywidthcentre = 35.0
        halflungwidth = lungwidth*0.5
        lungdepth=halflungwidth

        ycentreleftmedial = LMBcentre - halflungwidth
        ycentreleftlateral = LMBcentre + halflungwidth

        ymidLateralMedial = (ycentreleftmedial+ycentreleftlateral)*0.5

        lateralsurfradius =  lateralsurfcoeff*lungwidth
        internalsurfradius =  lateralsurfradius*2.25

        #-------------end of hardcoded values-----------------

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
        nodeIdentifier = 1
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

        for n2 in range(3):
            radiansUp = -thetaup + (n2) / (2.0) * (2 * thetaup)
            lRadiansUp.append(radiansUp)

            radiansLateral = -thetalateral + (n2) / (2.0) * (2 * thetalateral)
            lRadiansLateral.append(radiansLateral)

        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [LMBcentre, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + posteriorcurveradius * sinRadiansUp]
            d1 = [0, -posteriorcurveradius * sinRadiansUp, posteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sPosteriorx, sPosteriorderiv2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [LMBcentre - 2, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + posteriorcurveradius * sinRadiansUp]
            d1 = [0, -posteriorcurveradius * sinRadiansUp, posteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sPosteriorOffsetx, sPosteriorOffsetderiv2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

        zero = [0, 0, 0]

        # -----------------------------------------
        ## Anterior edge
        ## -------------------------------------
        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [LMBcentre, -lungwidth + ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + anteriorcurveradius * sinRadiansUp]
            d1 = [0, anteriorcurveradius * sinRadiansUp, anteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sAnteriorx, sAnteriorderiv2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

        # -----------------------------------------
        ## Anterior edge offset
        ## ---------------------------------------
        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [LMBcentre - 0.5, -lungwidth + ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + anteriorcurveradius * sinRadiansUp]
            d1 = [0, anteriorcurveradius * sinRadiansUp, anteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sAnteriorOffsetx, sAnteriorderivOffset2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

        # Apply tracheal rotation angle
        # ---------------------------------
        xcen = sPosteriorx[0][0]
        ycen = sPosteriorx[0][1]
        zcen = sPosteriorx[0][2]
        for n2 in range(0, elementsCountUp):
            newx = (sPosteriorx[n2][0] - xcen) * math.cos(TrachealRotationAngleRadians) - \
                   (sPosteriorx[n2][1] - ycen) * math.sin(TrachealRotationAngleRadians)

            newy = (sPosteriorx[n2][0] - xcen) * math.sin(TrachealRotationAngleRadians) + \
                   (sPosteriorx[n2][1] - ycen) * math.cos(TrachealRotationAngleRadians)

            sPosteriorx[n2][0] = newx + xcen
            sPosteriorx[n2][1] = newy + ycen

            newx = (sAnteriorx[n2][0] - xcen) * math.cos(TrachealRotationAngleRadians) - \
                   (sAnteriorx[n2][1] - ycen) * math.sin(TrachealRotationAngleRadians)

            newy = (sAnteriorx[n2][0] - xcen) * math.sin(TrachealRotationAngleRadians) + \
                   (sAnteriorx[n2][1] - ycen) * math.cos(TrachealRotationAngleRadians)

            sAnteriorx[n2][0] = newx + xcen
            sAnteriorx[n2][1] = newy + ycen

            newx = (sPosteriorOffsetx[n2][0] - xcen) * math.cos(TrachealRotationAngleRadians) - \
                   (sPosteriorOffsetx[n2][1] - ycen) * math.sin(TrachealRotationAngleRadians)

            newy = (sPosteriorOffsetx[n2][0] - xcen) * math.sin(TrachealRotationAngleRadians) + \
                   (sPosteriorOffsetx[n2][1] - ycen) * math.cos(TrachealRotationAngleRadians)

            sPosteriorOffsetx[n2][0] = newx + xcen
            sPosteriorOffsetx[n2][1] = newy + ycen

            newx = (sAnteriorOffsetx[n2][0] - xcen) * math.cos(TrachealRotationAngleRadians) - \
                   (sAnteriorOffsetx[n2][1] - ycen) * math.sin(TrachealRotationAngleRadians)

            newy = (sAnteriorOffsetx[n2][0] - xcen) * math.sin(TrachealRotationAngleRadians) + \
                   (sAnteriorOffsetx[n2][1] - ycen) * math.cos(TrachealRotationAngleRadians)

            sAnteriorOffsetx[n2][0] = newx + xcen
            sAnteriorOffsetx[n2][1] = newy + ycen

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
        # # ------------------------------------
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

        # ----------------------------------
        # # Create nodes on MEDIAL side
        # ----------------------------------
        for n2 in range(0, elementsCountUp):

            x1 = sPosteriorx[n2]
            x3 = sAnteriorx[n2]
            halflung = 0.5 * (x1[1] + x3[1])
            zlung = 0.5 * (x1[2] + x3[2])
            zratio = (0.3 - 0.55 * abs(zlung / lungheight - 0.5)) ** 0.5

            radiansLateral = lRadiansLateral[0]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)

            d21 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                   0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral, 0]

            radiansLateral = lRadiansLateral[1]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            x2 = [x1[0] - 0.2 * lungdepth - 0.5 * lungdepth * medialsurfcoeff * (cosRadiansLateral - 1.0),
                  halflung,
                  x1[2]]

            d22 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                   -0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral,
                   0]

            radiansLateral = lRadiansLateral[2]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            d23 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                   0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral, 0]

            nx = [x1, x2, x3]
            nd2 = [d21, d22, d23]
            sMedialx, sMedialderiv2, _, _, _  = \
                interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)

            x1 = sPosteriorOffsetx[n2]
            x3 = sAnteriorOffsetx[n2]
            halflung = 0.5 * (x1[1] + x3[1])

            radiansLateral = lRadiansLateral[0]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            d21 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                   0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                   0]

            radiansLateral = lRadiansLateral[1]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            x2 = [x1[0] - zratio * lungdepth - 0.5 * lungdepth * lateralsurfcoeff * (cosRadiansLateral - 1.0),
                  halflung,
                  x1[2]]
            d22 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                   -0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                   0]

            radiansLateral = lRadiansLateral[2]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            d23 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                   0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                   0]

            nx = [x1, x2, x3]
            nd2 = [d21, d22, d23]
            sLateralx, sLateralderiv2, _, _, _  = \
                interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)

            # Apply tracheal rotation angle
            # ---------------------------------
            for n3 in range(1, elementsCountlateral):
                newx = (sMedialx[n3][0] - xcen) * math.cos(TrachealRotationAngleRadians) - \
                       (sMedialx[n3][1] - ycen) * math.sin(TrachealRotationAngleRadians)

                newy = (sMedialx[n3][0] - xcen) * math.sin(TrachealRotationAngleRadians) + \
                       (sMedialx[n3][1] - ycen) * math.cos(TrachealRotationAngleRadians)

                sMedialx[n3][0] = newx + xcen
                sMedialx[n3][1] = newy + ycen

                newx = (sLateralx[n3][0] - xcen) * math.cos(TrachealRotationAngleRadians) - \
                       (sLateralx[n3][1] - ycen) * math.sin(TrachealRotationAngleRadians)

                newy = (sLateralx[n3][0] - xcen) * math.sin(TrachealRotationAngleRadians) + \
                       (sLateralx[n3][1] - ycen) * math.cos(TrachealRotationAngleRadians)
                sLateralx[n3][0] = newx + xcen
                sLateralx[n3][1] = newy + ycen

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

        e0 = nodeIdentifier - 1
        print('left nodes', e0)

        ###############################################################################################
        # RIGHT LUNG
        lungwidth = lungwidth * 0.8

        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [RMBcentre, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + posteriorcurveradius * sinRadiansUp]
            d1 = [0, -posteriorcurveradius * sinRadiansUp, posteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sPosteriorx, sPosteriorderiv2, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1,
                                                                                 elementsCountOut=elementsCountUp)

        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [RMBcentre + 2, ycentreleftlateral + posteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + posteriorcurveradius * sinRadiansUp]
            d1 = [0, -posteriorcurveradius * sinRadiansUp, posteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sPosteriorOffsetx, sPosteriorOffsetderiv2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

        zero = [0, 0, 0]

        # -----------------------------------------
        ## Anterior edge
        ## -------------------------------------
        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [RMBcentre, -lungwidth + ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
            zlateralcentre + anteriorcurveradius * sinRadiansUp]
            d1 = [0, anteriorcurveradius * sinRadiansUp, anteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sAnteriorx, sAnteriorderiv2, _, _, _ = interp.sampleCubicHermiteCurves(nx, nd1,
                                                                               elementsCountOut=elementsCountUp)

        # -----------------------------------------
        ## Anterior edge offset
        ## ---------------------------------------
        nx = []
        nd1 = []
        for n2 in range(3):
            radiansUp = lRadiansUp[n2]
            cosRadiansUp = math.cos(radiansUp)
            sinRadiansUp = math.sin(radiansUp)
            x1 = [RMBcentre + 0.5, -lungwidth + ycentreleftlateral - anteriorcurveradius * (cosRadiansUp - 1.0),
                zlateralcentre + anteriorcurveradius * sinRadiansUp]
            d1 = [0, anteriorcurveradius * sinRadiansUp, anteriorcurveradius * cosRadiansUp]
            nx.append(x1)
            nd1.append(d1)

        sAnteriorOffsetx, sAnteriorderivOffset2, _, _, _ = \
            interp.sampleCubicHermiteCurves(nx, nd1, elementsCountOut=elementsCountUp)

        # Apply tracheal rotation angle
        # ---------------------------------
        xcen = sPosteriorx[0][0]
        ycen = sPosteriorx[0][1]
        zcen = sPosteriorx[0][2]
        for n2 in range(0, elementsCountUp):
            newx = (sPosteriorx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                   (sPosteriorx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

            newy = (sPosteriorx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                   (sPosteriorx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

            sPosteriorx[n2][0] = newx + xcen
            sPosteriorx[n2][1] = newy + ycen

            newx = (sAnteriorx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                   (sAnteriorx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

            newy = (sAnteriorx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                   (sAnteriorx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

            sAnteriorx[n2][0] = newx + xcen
            sAnteriorx[n2][1] = newy + ycen

            newx = (sPosteriorOffsetx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                   (sPosteriorOffsetx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

            newy = (sPosteriorOffsetx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                   (sPosteriorOffsetx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

            sPosteriorOffsetx[n2][0] = newx + xcen
            sPosteriorOffsetx[n2][1] = newy + ycen

            newx = (sAnteriorOffsetx[n2][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                   (sAnteriorOffsetx[n2][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

            newy = (sAnteriorOffsetx[n2][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                   (sAnteriorOffsetx[n2][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

            sAnteriorOffsetx[n2][0] = newx + xcen
            sAnteriorOffsetx[n2][1] = newy + ycen

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
        # # ------------------------------------
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

        # ----------------------------------
        # # Create nodes on MEDIAL side
        # ----------------------------------
        for n2 in range(0, elementsCountUp):

            x1 = sPosteriorx[n2]
            x3 = sAnteriorx[n2]
            halflung = 0.5 * (x1[1] + x3[1])
            zlung = 0.5 * (x1[2] + x3[2])
            zratio = (0.25 - 0.4 * abs(zlung / lungheight - 0.5)) ** 0.5

            radiansLateral = lRadiansLateral[0]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)

            d21 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                   0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral, 0]

            radiansLateral = lRadiansLateral[1]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            x2 = [x1[0] + 0.2 * lungdepth + 0.5 * lungdepth * medialsurfcoeff * (cosRadiansLateral - 1.0),
                  halflung,
                  x1[2]]

            d22 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                   -0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral,
                   0]

            radiansLateral = lRadiansLateral[2]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            d23 = [0.5 * lungdepth * medialsurfcoeff * sinRadiansLateral,
                   0.5 * lungdepth * medialsurfcoeff * cosRadiansLateral, 0]

            nx = [x1, x2, x3]
            nd2 = [d21, d22, d23]
            sMedialx, sMedialderiv2, _, _, _  = \
                interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)

            x1 = sPosteriorOffsetx[n2]
            x3 = sAnteriorOffsetx[n2]
            halflung = 0.5 * (x1[1] + x3[1])

            radiansLateral = lRadiansLateral[0]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            d21 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                   0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                   0]

            radiansLateral = lRadiansLateral[1]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            x2 = [x1[0] + zratio * lungdepth + 0.5 * lungdepth * lateralsurfcoeff * (cosRadiansLateral - 1.0),
                  halflung,
                  x1[2]]
            d22 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                   -0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                   0]

            radiansLateral = lRadiansLateral[2]
            cosRadiansLateral = math.cos(radiansLateral)
            sinRadiansLateral = math.sin(radiansLateral)
            d23 = [0.5 * lungdepth * lateralsurfcoeff * sinRadiansLateral,
                   0.5 * lungwidth * lateralsurfcoeff * cosRadiansLateral,
                   0]

            nx = [x1, x2, x3]
            nd2 = [d21, d22, d23]
            sLateralx, sLateralderiv2,  _, _, _  = \
                interp.sampleCubicHermiteCurves(nx, nd2, elementsCountOut=elementsCountlateral)

            # Apply tracheal rotation angle
            # ---------------------------------
            for n3 in range(1, elementsCountlateral):
                newx = (sMedialx[n3][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                       (sMedialx[n3][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                newy = (sMedialx[n3][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                       (sMedialx[n3][1] - ycen) * math.cos(-TrachealRotationAngleRadians)

                sMedialx[n3][0] = newx + xcen
                sMedialx[n3][1] = newy + ycen

                newx = (sLateralx[n3][0] - xcen) * math.cos(-TrachealRotationAngleRadians) - \
                       (sLateralx[n3][1] - ycen) * math.sin(-TrachealRotationAngleRadians)

                newy = (sLateralx[n3][0] - xcen) * math.sin(-TrachealRotationAngleRadians) + \
                       (sLateralx[n3][1] - ycen) * math.cos(-TrachealRotationAngleRadians)
                sLateralx[n3][0] = newx + xcen
                sLateralx[n3][1] = newy + ycen

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

        #################
        # Create elements
        #################
        mesh = fm.findMeshByDimension(3)
        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()

        tricubicHermiteBasis = fm.createElementbasis(3, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)

        # Regular elements
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        elementtemplate.defineField(coordinates, -1, eft)

        elementIdentifier = 1

        delQ = 4 * elementsCountUp
        delX = 2 * elementsCountlateral - 2
        delA = 4 * elementsCountUp + (elementsCountlateral - 1)
        delB = 4 * elementsCountUp + (3 * elementsCountlateral - 3)
        delM = 2 * elementsCountUp
        delP = 4 * elementsCountUp + (elementsCountlateral - 2)

        for e2 in range(1, elementsCountUp):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            nodeIdentifiers = [e2,
                               1 + delQ + (delX) * (e2 - 1),
                               e2 + 1,
                               1 + delQ + (delX) * (e2 - 1) + delX,
                               e2 + delM,
                               1 + delA + (delX) * (e2 - 1),
                               e2 + delM + 1,
                               1 + delB + (delX) * (e2 - 1)
                               ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

        for e3 in range(1, elementsCountlateral - 1):
            for e2 in range(1, elementsCountUp):
                e22 = 1 + delQ + (delX) * (e2 - 1) + (e3 - 1)
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

        for e2 in range(1, elementsCountUp):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            e22 = 1 + delQ + (delX) * (e2 - 1)
            nodeIdentifiers = [1 + delP + (delX) * (e2 - 1),
                               (e2) + elementsCountUp,
                               1 + delP + (delX) * (e2),
                               (e2) + elementsCountUp + 1,
                               1 + delP + (delX) * (e2 - 1) + (elementsCountlateral - 1),
                               e2 + (3 * elementsCountUp),
                               1 + delP + (delX) * (e2) + (elementsCountlateral - 1),
                               e2 + (3 * elementsCountUp) + 1
                               ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

        ###########################################################################
        # RIGHT LUNG
        ###########################################################################
        e0 = nodeIdentifier - 1

        for e2 in range(1, elementsCountUp):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            nodeIdentifiers = [e0 + e2,
                               e0 + 1 + delQ + (delX) * (e2 - 1),
                               e0 + e2 + 1,
                               e0 + 1 + delQ + (delX) * (e2 - 1) + delX,
                               e0 + e2 + delM,
                               e0 + 1 + delA + (delX) * (e2 - 1),
                               e0 + e2 + delM + 1,
                               e0 + 1 + delB + (delX) * (e2 - 1)
                               ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

        for e3 in range(1, elementsCountlateral - 1):
            for e2 in range(1, elementsCountUp):
                e22 = e0 + 1 + delQ + (delX) * (e2 - 1) + (e3 - 1)
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

        for e2 in range(1, elementsCountUp):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            e22 = 1 + delQ + (delX) * (e2 - 1)
            nodeIdentifiers = [e0 + 1 + delP + (delX) * (e2 - 1),
                               e0 + (e2) + elementsCountUp,
                               e0 + 1 + delP + (delX) * (e2),
                               e0 + (e2) + elementsCountUp + 1,
                               e0 + 1 + delP + (delX) * (e2 - 1) + (elementsCountlateral - 1),
                               e0 + e2 + (3 * elementsCountUp),
                               e0 + 1 + delP + (delX) * (e2) + (elementsCountlateral - 1),
                               e0 + e2 + (3 * elementsCountUp) + 1
                               ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

        fm.endChange()
        # return annotationGroups

    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """
        if not options['Refine']:
            cls.generateBaseMesh(region, options)
            return

        refineElementsCountAround = options['Refine number of elements around']
        refineElementsCountUp = options['Refine number of elements up']
        refineElementsCountRadial = options['Refine number of elements radial']

        baseRegion = region.createRegion()
        cls.generateBaseMesh(baseRegion, options)

        meshrefinement = MeshRefinement(baseRegion, region)
        meshrefinement.refineAllElementsCubeStandard3d(refineElementsCountAround, refineElementsCountUp,
                                                       refineElementsCountRadial)