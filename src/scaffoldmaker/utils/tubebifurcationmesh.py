'''
Utility function for generating tubular mesh from a central line
using a segment profile.
'''
from __future__ import division
import math
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates, findOrCreateFieldTextureCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector

def TubeBifurcationMesh(region,
    xList, d1List, radius1list,
    elementsCountAround, elementsCountAlong,
    annotationGroups, annotationArray,
    firstNodeIdentifier, firstElementIdentifier,
    useCrossDerivatives):
    """
    :param xList: coordinates of centerline points.
    :param d1List: derivatives along axis of segment.
    :param radius1List: derivatives along axis of segment.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :return coordinates and derivatives of warped points.
    """
    xWarpedList = []
    d1WarpedList = []
    d2WarpedList = []
    smoothd2WarpedList = []
    d3WarpedUnitList = []

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [ 0.0, 0.0, 0.0 ]

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

    mesh = fm.findMeshByDimension(3)

    eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
    eft = eftfactory.createEftBasic()

    radiansPerElementAround = 2.0 * math.pi / elementsCountAround
    x = [0.0, 0.0, 0.0]
    dx_ds1 = [0.0, 0.0, 0.0]
    dx_ds2 = [0.0, 0.0, 1.0 / elementsCountAlong]
    zero = [0.0, 0.0, 0.0]
    radius = 0.5
    for n2 in range(elementsCountAlong + 1):
        x[2] = n2 / elementsCountAlong
        for n1 in range(elementsCountAround):
            radiansAround = n1 * radiansPerElementAround
            cosRadiansAround = math.cos(radiansAround)
            sinRadiansAround = math.sin(radiansAround)
            x[0] = radius * cosRadiansAround
            x[1] = radius * sinRadiansAround
            dx_ds1[0] = radiansPerElementAround * radius * -sinRadiansAround
            dx_ds1[1] = radiansPerElementAround * radius * cosRadiansAround
            node = nodes.createNode(nodeIdentifier, nodetemplate)
            cache.setNode(node)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)

    fm.endChange()


    return xWarpedList, d1WarpedList, smoothd2WarpedList, d3WarpedUnitList
