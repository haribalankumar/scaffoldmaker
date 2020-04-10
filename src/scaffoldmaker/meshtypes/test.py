"""
Generates a 3-D unit tube mesh with variable numbers of elements around, along and
through wall, plus variable wall thickness for unit diameter.
"""

from __future__ import division
import math
import copy
import numpy as np
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.meshrefinement import MeshRefinement
from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes
from scaffoldmaker.utils.geometry import createCirclePoints, getCircleProjectionAxes

from scaffoldmaker.utils import vector

class MeshType_3d_airwaybifurcation1(Scaffold_base):
    '''
    classdocs
    '''
    @staticmethod
    def getName():
        return '3D Airway 1'

    @staticmethod
    def getDefaultOptions(parameterSetName='Default'):
        options = {
            'Number of elements around' : 8,
            'Number of elements along' : 8,
            'Number of elements through wall' : 1,
            'Unit scale': 1.0,
            'Outlet': False,
            'Wall thickness': 0.1,
            'Ostium diameter': 1.0,
            'Daughter length': 0.5,
            'Branch angle 1 degrees': 0.0,
            'Branch angle 1 spread degrees': 0.0,
            'Branch angle 2 degrees': 0.0,
            'Branch end length factor': 1.0,
            'Number of elements along branch',6,
            'Use cross derivatives': False,
            'Refine': False,
            'Refine number of elements around' : 1,
            'Refine number of elements along' : 1,
            'Refine number of elements through wall' : 1
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements around',
            'Number of elements along',
            'Number of elements through wall',
            'Wall thickness',
            'Unit scale',
            'Outlet',
            'Ostium diameter',
            'Daughter length',
            'Branch angle 1 degrees',
            'Branch angle 1 spread degrees',
            'Branch angle 2 degrees',
            'Branch end length factor',
            'Number of elements along branch',
            'Use cross derivatives',
            'Refine',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall'
        ]

    @staticmethod
    def checkOptions(options):
        for key in [
            'Number of elements along',
            'Number of elements through wall',
            'Refine number of elements around',
            'Refine number of elements along',
            'Refine number of elements through wall']:
            if options[key] < 1:
                options[key] = 1
        if (options['Wall thickness'] < 0.0) :
            options['Wall thickness'] = 0.0
        elif (options['Wall thickness'] > 0.5) :
            options['Wall thickness'] = 0.5


    @staticmethod
    def generateBaseMesh(region, options):
        """
        Generate the base tricubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: None
        """

        #elementsCountAround = options['Number of elements around']
        #elementsCountAlong = options['Number of elements along']
        #elementsCountThroughWall = options['Number of elements through wall']

        wallThickness = options['Wall thickness']
        useCrossDerivatives = options['Use cross derivatives']

        vesselAngle1Radians = math.radians(options['Branch angle 1 degrees'])
        vesselAngle1SpreadRadians = math.radians(options['Branch angle 1 spread degrees'])
        vesselAngle2Radians = math.radians(options['Branch angle 2 degrees'])
        Branchlength = math.radians(options['Daughter length'])

        elementsalongBranch = options['Number of elements along branch']
        useCubicHermiteThroughVesselWall = True
        useCubicHermiteThroughOstiumWall = True

        elementsCountAround = 8
        elementsCountAlong = 8
        elementsCountThroughWall = 1
        daughterlength = 0.5

        fm = region.getFieldmodule()
        fm.beginChange()
        coordinates = findOrCreateFieldCoordinates(fm)

        nodes = fm.findNodesetByFieldDomainType(Field.DOMAIN_TYPE_NODES)
        nodetemplate = nodes.createNodetemplate()
        nodetemplate.defineField(coordinates)

        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

        mesh = fm.findMeshByDimension(3)
        tricubichermite = eftfactory_tricubichermite(mesh, useCrossDerivatives)
        eft = tricubichermite.createEftBasic()
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

        # create nodes
        nodeIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        wallThicknessPerElement = wallThickness/elementsCountThroughWall
        x = [ 0.0, 0.0, 0.0 ]
        x1 = [ 0.0, 0.0, 0.0 ]
        dx_ds1 = [ 0.0, 0.0, 0.0 ]
        dx_ds2 = [ 0.0, 0.0, 1.0 / elementsCountAlong ]
        dx_ds3 = [ 0.0, 0.0, 0.0 ]
        zero = [ 0.0, 0.0, 0.0 ]
        axis1 = [ 0.0, -0.5, -0.866 ]
        startPointsx = []
        startPointsd1 = []
        startPointsd2 = []
        startPointsd3 = []
        startNodeId = []
        endPointsx = []
        #endPointsd1,endPointsd2,endPointsd3, endNodeId, endDerivativesMap =  None, None, None, None, None
        #startDerivativesMap =  None

        for n3 in range(elementsCountThroughWall + 1):
            radius = 0.5 + wallThickness*(n3/elementsCountThroughWall - 1.0)
            for n2 in range(elementsCountAlong + 1):
                x[2] = n2 / elementsCountAlong
                for n1 in range(elementsCountAround):
                    radiansAround = n1*radiansPerElementAround
                    cosRadiansAround = math.cos(radiansAround)
                    sinRadiansAround = math.sin(radiansAround)
                    x[0] = radius*cosRadiansAround
                    x[1] = radius*sinRadiansAround
                    dx_ds1[0] = radiansPerElementAround*radius*-sinRadiansAround
                    dx_ds1[1] = radiansPerElementAround*radius*cosRadiansAround
                    dx_ds3[0] = wallThicknessPerElement*cosRadiansAround
                    dx_ds3[1] = wallThicknessPerElement*sinRadiansAround
                    node = nodes.createNode(nodeIdentifier, nodetemplate)

                    # edit nodes
                    if nodeIdentifier==109 or nodeIdentifier==110 or nodeIdentifier==111 or nodeIdentifier==112 or nodeIdentifier==117 or nodeIdentifier==120 or nodeIdentifier==125 or nodeIdentifier==126 or nodeIdentifier==127 or nodeIdentifier==128:
                        startPointsx.append(x)
                        startPointsd1.append(dx_ds1)
                        startPointsd2.append(dx_ds2)
                        startPointsd3.append(dx_ds3)
                        startNodeId.append(nodeIdentifier)
                        x1[0] = x[0] + axis1[0]*daughterlength
                        x1[1] = x[1] + axis1[1]*daughterlength
                        x1[2] = x[2] + axis1[2]*daughterlength

                        endPointsx.append(x1)

                    cache.setNode(node)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
                    coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, dx_ds3)
                    if useCrossDerivatives:
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
                    nodeIdentifier = nodeIdentifier + 1

        # create elements
        elementIdentifier = 1
        now = (elementsCountAlong + 1)*elementsCountAround
        for e3 in range(elementsCountThroughWall):
            for e2 in range(elementsCountAlong):
                for e1 in range(elementsCountAround):
                    if elementIdentifier == 37 or elementIdentifier == 38 or elementIdentifier == 39 or elementIdentifier == 45 or elementIdentifier == 46 or elementIdentifier == 47:
                        elementIdentifier = elementIdentifier + 1
                    else:
                        element = mesh.createElement(elementIdentifier, elementtemplate)
                        bni11 = e3*now + e2*elementsCountAround + e1 + 1
                        bni12 = e3*now + e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                        bni21 = e3*now + (e2 + 1)*elementsCountAround + e1 + 1
                        bni22 = e3*now + (e2 + 1)*elementsCountAround + (e1 + 1)%elementsCountAround + 1
                        nodeIdentifiers = [ bni11, bni12, bni21, bni22, bni11 + now, bni12 + now, bni21 + now, bni22 + now ]
                        result = element.setNodesByIdentifier(eft, nodeIdentifiers)
                        elementIdentifier = elementIdentifier + 1

        unitScale = options['Unit scale']
        ostiumRadius = 0.5 * unitScale * options['Ostium diameter']
        scale = 1.1*(ostiumRadius*2.0)

        #nodeIdentifier, elementIdentifier = createAnnulusMesh3d(
        #    nodes, mesh, nodeIdentifier, elementIdentifier,
        #    startPointsx, startPointsd1, startPointsd2, startPointsd3, startNodeId, startDerivativesMap,
        #    endPointsx, endPointsd1, endPointsd2, endPointsd3, endNodeId, endDerivativesMap,
        #    forceMidLinearXi3 = not useCubicHermiteThroughVesselWall,
        #    elementsCountRadial = elementsalongBranch, meshGroups = [])

        fm.endChange()

    #USAGE
#       nodeIdentifier, elementIdentifier, (lpvox, lpvod1, lpvod2, lpvod3, lpvoNodeId, lpvoPositions) = \
#            generateOstiumMesh(region, lpvOstiumSettings, laTrackSurface, lpvOstiumPosition, lpvOstiumDirection,
#                              nodeIdentifier, elementIdentifier, vesselMeshGroups)


    @classmethod
    def generateMesh(cls, region, options):
        """
        Generate base or refined mesh.
        :param region: Zinc region to create mesh in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        """

        if not options['Refine']:
            return cls.generateBaseMesh(region, options)

        baseRegion = region.createRegion()
        baseAnnotationGroups = cls.generateBaseMesh(baseRegion, options)
        meshrefinement = MeshRefinement(baseRegion, region, baseAnnotationGroups)
        cls.refineMesh(meshrefinement, options)
        return meshrefinement.getAnnotationGroups()

=========================================================================================
nodesOnTrackSurfacex[0][0] = -0.35355339059327384
nodesOnTrackSurfacex[0][1] = -0.35355339059327373
nodesOnTrackSurfacex[0][2] = 0.5

nodesOnTrackSurfacex[1][0] = -0
nodesOnTrackSurfacex[1][1] = -0.5
nodesOnTrackSurfacex[1][2] = 0.5

nodesOnTrackSurfacex[2][0] = 0.3535533905932737
nodesOnTrackSurfacex[2][1] = -0.35355
nodesOnTrackSurfacex[2][2] = 0.5

nodesOnTrackSurfacex[3][0] = -0.35355339059327384
nodesOnTrackSurfacex[3][1] = -0.35355339059327373
nodesOnTrackSurfacex[3][2] = 0.625

nodesOnTrackSurfacex[4][0] = -0
nodesOnTrackSurfacex[4][1] = -0.5
nodesOnTrackSurfacex[4][2] = 0.625

nodesOnTrackSurfacex[5][0] = 0.3535533905932737
nodesOnTrackSurfacex[5][1] = -0.35355339059327384
nodesOnTrackSurfacex[5][2] = 0.625

nodesOnTrackSurfacex[6][0] = -0.3535533905932737
nodesOnTrackSurfacex[6][1] = -0.35355339059327384
nodesOnTrackSurfacex[6][2] = 0.75

nodesOnTrackSurfacex[7][0] = 0
nodesOnTrackSurfacex[7][1] = -0.5
nodesOnTrackSurfacex[7][2] = 0.75

nodesOnTrackSurfacex[8][0] = 0.3535533905932737
nodesOnTrackSurfacex[8][1] = -0.35355339059327384
nodesOnTrackSurfacex[8][2] = 0.75