'''
Utility function for generating tubular mesh from a central line
using a segment profile.
'''
from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis, Elementfieldtemplate
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node
from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base

from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear
from scaffoldmaker.utils.eftfactory_tricubichermite import eftfactory_tricubichermite
from scaffoldmaker.utils.eft_utils import remapEftNodeValueLabel, setEftScaleFactorIds

from scaffoldmaker.utils.geometry import createCirclePoints
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils import matrix
from scaffoldmaker.utils import vector


def derivativeSignsToExpressionTerms(valueLabels, signs):
    '''
    Return remap expression terms for summing derivative[i]*sign[i]
    :param valueLabels: List of node value labels to possibly include.
    :param signs: List of 1 (no scaling), -1 (scale by scale factor 1) or 0 (no term).
    '''
    expressionTerms = []
    for i in range(len(valueLabels)):
        if signs[i] is 1:
            expressionTerms.append( ( valueLabels[i], [] ) )
        elif signs[i] is -1:
            expressionTerms.append( ( valueLabels[i], [1] ) )
    return expressionTerms


def createjunctionAirwaySurfaceSegmentPoints(
        xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
        d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
        d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
        segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
        parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
        sxparent, sxDaugh1, sxDaugh2,
        sd1parent, sd1Daugh1, sd1Daugh2,
        sd2parent, sd2Daugh1, sd2Daugh2,
        elementsCountAround, elementsCountAlongSegment, nSegment):

    """
    Warps points in Airway segment to account for bending and twisting
    along central path defined by nodes sx and derivatives sd1 and sd2.
    :param x1List: coordinates of segment points.
    :param d1List: derivatives around axis of segment.
    :param d2List: derivatives along axis of segment.
    :param segmentAxis: axis perpendicular to segment plane.
    :param sx: coordinates of points on central path.
    :param sd1: derivatives of points along central path.
    :param sd2: derivatives representing cross axes.
    groups along the segment.
    :return coordinates and derivatives of warped points.
    """
    xjunctionOuterList = []
    xjunctionInnerList = []

    d1junctionOuterList = []
    d2junctionOuterList = []
    zero = [0,0,0]

    d1junctionInnerList = []
    d2junctionInnerList = []

    n2 = elementsCountAlongSegment * nSegment + elementsCountAlongSegment
    n2t = elementsCountAround*(elementsCountAlongSegment+1)
    lastring = (elementsCountAlongSegment) * elementsCountAround

    #JUNCTION OUTER
    #################
    nx = []
    nd2 = []
    xParentAlongSegment = xParentWarpedList[lastring+1]
    xDaugh1AlongSegment = xDaugh1WarpedList[2]
    xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    xjunctionOuter1 = [0.35*xmid[j]+0.65*xDaugh1AlongSegment[j] for j in range(3)]
    xjunctionOuterList.append(xjunctionOuter1)

    #deriv calc next
    nx.append(xParentAlongSegment)
    nx.append(xjunctionOuter1)
    nx.append(xDaugh1AlongSegment)
    nd2.append(d2ParentWarpedList[lastring+1])
    nd2.append(zero)
    nd2.append(d2Daugh1WarpedList[2])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2junctionOuterList.append(smoothedd2[1])

    nx = []
    nd2 = []
    xParentAlongSegment = xParentWarpedList[lastring+3]
    xDaugh2AlongSegment = xDaugh2WarpedList[2]
    xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    xjunctionOuter2 = [0.35*xmid[j]+0.65*xDaugh2AlongSegment[j] for j in range(3)]
    xjunctionOuterList.append(xjunctionOuter2)
    #deriv calc next
    nx.append(xParentAlongSegment)
    nx.append(xjunctionOuter2)
    nx.append(xDaugh2AlongSegment)
    nd2.append(d2ParentWarpedList[lastring+3])
    nd2.append(zero)
    nd2.append(d2Daugh2WarpedList[2])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2junctionOuterList.append(smoothedd2[1])

    nx = []
    nd2 = []
    xDaugh1AlongSegment = xDaugh1WarpedList[0]
    xDaugh2AlongSegment = xDaugh2WarpedList[0]
    xjunctionOuter3 = [(7*(xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])+2*sxparent[n2][j])/16.0 for j in range(3)]
    xjunctionOuterList.append(xjunctionOuter3)
    nx.append(xDaugh1AlongSegment)
    nx.append(xjunctionOuter3)
    nx.append(xDaugh2AlongSegment)
    nd2.append(d1ParentWarpedList[0])
    nd2.append(zero)
    d2temp = [-1.0 * d2Daugh2WarpedList[0][j] for j in range(3)]
    nd2.append(d2temp)
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2junctionOuterList.append(smoothedd2[1])

    #JUNCTION INNER
    ##################
    xParentAlongSegment = xParentWarpedList[lastring]
    xDaugh1AlongSegment = xDaugh1WarpedList[1]
    xDaugh2AlongSegment = xDaugh2WarpedList[3]
    xjunctionInner1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    xjunctionInnerList.append(xjunctionInner1)

    xParentAlongSegment = xParentWarpedList[lastring+2]
    xDaugh1AlongSegment = xDaugh1WarpedList[3]
    xDaugh2AlongSegment = xDaugh2WarpedList[1]
    xjunctionInner2 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    xjunctionInnerList.append(xjunctionInner2)

    ##DERIVS
    #-----------
    nx = []
    nd2 = []
    nx.append(xDaugh2WarpedList[elementsCountAround-3])
    nx.append(xjunctionInner2)
    nx.append(xjunctionOuter1)
    nx.append(xjunctionInner1)
    nx.append(xDaugh2WarpedList[elementsCountAround-1])
    d2temp = [-1.0 * d2Daugh2WarpedList[elementsCountAround-3][j] for j in range(3)]
    nd2.append(d2temp)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(d2Daugh2WarpedList[elementsCountAround-1])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    d1junctionOuterList.append(d2temp)
    d2temp = [-1.0 * smoothedd2[3][j] for j in range(3)]
    d1junctionInnerList.append(d2temp)
    d1junctionInnerList.append(smoothedd2[1])

    ##DERIVS
    nx = []
    nd2 = []
    nx.append(xDaugh1WarpedList[elementsCountAround-1])
    nx.append(xjunctionInner2)
    nx.append(xjunctionOuter2)
    nx.append(xjunctionInner1)
    nx.append(xDaugh1WarpedList[elementsCountAround-3])
    d2temp = [-1.0 * d2Daugh1WarpedList[elementsCountAround-1][j] for j in range(3)]
    nd2.append(d2temp)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(d2Daugh1WarpedList[elementsCountAround-3])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    d1junctionOuterList.append(smoothedd2[2])

    ##DERIVS
    nx = []
    nd2 = []
    nx.append(xParentWarpedList[n2t-2])
    nx.append(xjunctionInner2)
    nx.append(xjunctionOuter3)
    nx.append(xjunctionInner1)
    nx.append(xParentWarpedList[n2t-4])

    #nd2.append(d2ParentWarpedList[n2t-2])
    d2temp = [(-1.0 * d2ParentWarpedList[n2t-2][j]) for j in range(3)]
    nd2.append(d2temp)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(d2ParentWarpedList[n2t-4])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    d1junctionOuterList.append(d2temp)
    d2temp = [-1.0 * smoothedd2[3][j] for j in range(3)]
    d2junctionInnerList.append(d2temp)
    d2temp = [1.0 * smoothedd2[1][j] for j in range(3)]
    d2junctionInnerList.append(d2temp)

    return xjunctionOuterList, xjunctionInnerList, \
           d1junctionOuterList, d1junctionInnerList, \
           d2junctionOuterList, d2junctionInnerList



def createjunctionAirwaySegmentPoints(
        xParentWarpedList, xDaugh1WarpedList, xDaugh2WarpedList,
        d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList,
        d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList,
        segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
        sxparent, sxDaugh1, sxDaugh2,
        sd1parent, sd1Daugh1, sd1Daugh2,
        sd2parent, sd2Daugh1, sd2Daugh2,
        elementsCountAround, elementsCountAlongSegment, nSegment):

    """
    Warps points in Airway segment to account for bending and twisting
    along central path defined by nodes sx and derivatives sd1 and sd2.
    :param x1List: coordinates of segment points.
    :param d1List: derivatives around axis of segment.
    :param d2List: derivatives along axis of segment.
    :param segmentAxis: axis perpendicular to segment plane.
    :param sx: coordinates of points on central path.
    :param sd1: derivatives of points along central path.
    :param sd2: derivatives representing cross axes.
    groups along the segment.
    :return coordinates and derivatives of warped points.
    """
    xjunctionOuterList = []
    xjunctionInnerList = []

    d1junctionOuterList = []
    d2junctionOuterList = []
    zero = [0,0,0]

    d3junctionOuterList = []
    d1junctionInnerList = []
    d2junctionInnerList = []
    d3junctionInnerList = []

    n2 = elementsCountAlongSegment * nSegment + elementsCountAlongSegment
    n2t = elementsCountAround*(elementsCountAlongSegment+1)
    lastring = (elementsCountAlongSegment) * elementsCountAround

    # #JUNCTION OUTER
    # #################
    # nx = []
    # nd2 = []
    # xParentAlongSegment = xParentWarpedList[9]
    # xDaugh1AlongSegment = xDaugh1WarpedList[2]
    # xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    # xjunctionOuter1 = [(xmid[j]+xDaugh1AlongSegment[j])/2.0 for j in range(3)]
    # xjunctionOuterList.append(xjunctionOuter1)
    # ###deriv calc next
    # nx.append(xParentAlongSegment)
    # nx.append(xjunctionOuter1)
    # nx.append(xDaugh1AlongSegment)
    # nd2.append(d2ParentWarpedList[9])
    # nd2.append(zero)
    # nd2.append(d2Daugh1WarpedList[2])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2junctionOuterList.append(smoothedd2[1])
    #
    # nx = []
    # nd2 = []
    # xParentAlongSegment = xParentWarpedList[11]
    # xDaugh2AlongSegment = xDaugh2WarpedList[2]
    # xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    # xjunctionOuter2 = [(xmid[j]+xDaugh2AlongSegment[j])/2.0 for j in range(3)]
    # xjunctionOuterList.append(xjunctionOuter2)
    # ###deriv calc next
    # nx.append(xParentAlongSegment)
    # nx.append(xjunctionOuter2)
    # nx.append(xDaugh2AlongSegment)
    # nd2.append(d2ParentWarpedList[11])
    # nd2.append(zero)
    # nd2.append(d2Daugh2WarpedList[2])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2junctionOuterList.append(smoothedd2[1])
    #
    # nx = []
    # nd2 = []
    # xDaugh1AlongSegment = xDaugh1WarpedList[0]
    # xDaugh2AlongSegment = xDaugh2WarpedList[0]
    # xjunctionOuter3 = [(2*(xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])+sxparent[n2][j])/5.0 for j in range(3)]
    # xjunctionOuterList.append(xjunctionOuter3)
    # nx.append(xDaugh1AlongSegment)
    # nx.append(xjunctionOuter3)
    # nx.append(xDaugh2AlongSegment)
    # nd2.append(d1ParentWarpedList[0])
    # nd2.append(zero)
    # d2temp = [-1.0 * d2Daugh2WarpedList[0][j] for j in range(3)]
    # nd2.append(d2temp)
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2junctionOuterList.append(smoothedd2[1])
    #
    # ###JUNCTION INNER
    # ##################
    # xParentAlongSegment = xParentWarpedList[8]
    # xDaugh1AlongSegment = xDaugh1WarpedList[1]
    # xDaugh2AlongSegment = xDaugh2WarpedList[1]
    # xjunctionInner1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    # xjunctionInnerList.append(xjunctionInner1)
    #
    # xParentAlongSegment = xParentWarpedList[10]
    # xDaugh1AlongSegment = xDaugh1WarpedList[3]
    # xDaugh2AlongSegment = xDaugh2WarpedList[3]
    # xjunctionInner2 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    # xjunctionInnerList.append(xjunctionInner2)
    #
    # ###DERIVS
    # nx = []
    # nd2 = []
    # nx.append(xDaugh2WarpedList[elementsCountAround-1])
    # nx.append(xjunctionInner2)
    # nx.append(xjunctionOuter1)
    # nx.append(xjunctionInner1)
    # nx.append(xDaugh2WarpedList[elementsCountAround-2])
    # d2temp = [-1.0 * d2Daugh2WarpedList[elementsCountAround-1][j] for j in range(3)]
    # nd2.append(d2temp)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(d2Daugh2WarpedList[elementsCountAround-2])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d1junctionOuterList.append(smoothedd2[2])
    # d2junctionInnerList.append(smoothedd2[3])
    # d1junctionInnerList.append(smoothedd2[1])
    # d1junctionInnerList.append(smoothedd2[3])
    #
    # ###DERIVS
    # nx = []
    # nd2 = []
    # nx.append(xDaugh1WarpedList[elementsCountAround-1])
    # nx.append(xjunctionInner2)
    # nx.append(xjunctionOuter2)
    # nx.append(xjunctionInner1)
    # nx.append(xDaugh1WarpedList[elementsCountAround-2])
    # d2temp = [-1.0 * d2Daugh1WarpedList[elementsCountAround-1][j] for j in range(3)]
    # nd2.append(d2temp)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(d2Daugh1WarpedList[elementsCountAround-2])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d1junctionOuterList.append(smoothedd2[2])
    # d2temp = [-1.0 * smoothedd2[1][j] for j in range(3)]
    # d2junctionInnerList.append(d2temp)
    #
    # ###DERIVS
    # nx = []
    # nd2 = []
    # nx.append(xParentWarpedList[n2t-2])
    # nx.append(xjunctionInner2)
    # nx.append(xjunctionOuter3)
    # nx.append(xjunctionInner1)
    # nx.append(xParentWarpedList[n2t-4])
    # nd2.append(d2ParentWarpedList[n2t-2])
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(zero)
    # d2temp = [(-1.0 * d2ParentWarpedList[n2t-4][j]) for j in range(3)]
    # nd2.append(d2temp)
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    # d1junctionOuterList.append(d2temp)
    #
    # #JUNCTION OUTER
    # #################
    # nx = []
    # nd2 = []
    # xParentAlongSegment = xParentWarpedList[lastring+1]
    # xDaugh1AlongSegment = xDaugh1WarpedList[2]
    # xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    # xjunctionOuter1 = [0.35*xmid[j]+0.65*xDaugh1AlongSegment[j] for j in range(3)]
    # xjunctionOuterList.append(xjunctionOuter1)
    #
    # #deriv calc next
    # nx.append(xParentAlongSegment)
    # nx.append(xjunctionOuter1)
    # nx.append(xDaugh1AlongSegment)
    # nd2.append(d2ParentWarpedList[lastring+1])
    # nd2.append(zero)
    # nd2.append(d2Daugh1WarpedList[2])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2junctionOuterList.append(smoothedd2[1])
    #
    # nx = []
    # nd2 = []
    # xParentAlongSegment = xParentWarpedList[lastring+3]
    # xDaugh2AlongSegment = xDaugh2WarpedList[2]
    # xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    # xjunctionOuter2 = [0.35*xmid[j]+0.65*xDaugh2AlongSegment[j] for j in range(3)]
    # xjunctionOuterList.append(xjunctionOuter2)
    # #deriv calc next
    # nx.append(xParentAlongSegment)
    # nx.append(xjunctionOuter2)
    # nx.append(xDaugh2AlongSegment)
    # nd2.append(d2ParentWarpedList[lastring+3])
    # nd2.append(zero)
    # nd2.append(d2Daugh2WarpedList[2])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2junctionOuterList.append(smoothedd2[1])
    #
    # nx = []
    # nd2 = []
    # xDaugh1AlongSegment = xDaugh1WarpedList[0]
    # xDaugh2AlongSegment = xDaugh2WarpedList[0]
    # xjunctionOuter3 = [(7*(xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])+2*sxparent[n2][j])/16.0 for j in range(3)]
    # xjunctionOuterList.append(xjunctionOuter3)
    # nx.append(xDaugh1AlongSegment)
    # nx.append(xjunctionOuter3)
    # nx.append(xDaugh2AlongSegment)
    # nd2.append(d1ParentWarpedList[0])
    # nd2.append(zero)
    # d2temp = [-1.0 * d2Daugh2WarpedList[0][j] for j in range(3)]
    # nd2.append(d2temp)
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2junctionOuterList.append(smoothedd2[1])
    #
    # #JUNCTION INNER
    # ##################
    # xParentAlongSegment = xParentWarpedList[lastring]
    # xDaugh1AlongSegment = xDaugh1WarpedList[1]
    # xDaugh2AlongSegment = xDaugh2WarpedList[3]
    # xjunctionInner1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    # xjunctionInnerList.append(xjunctionInner1)
    #
    # xParentAlongSegment = xParentWarpedList[lastring+2]
    # xDaugh1AlongSegment = xDaugh1WarpedList[3]
    # xDaugh2AlongSegment = xDaugh2WarpedList[1]
    # xjunctionInner2 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    # xjunctionInnerList.append(xjunctionInner2)
    #
    # ##DERIVS
    # #-----------
    # nx = []
    # nd2 = []
    # nx.append(xDaugh2WarpedList[elementsCountAround-3])
    # nx.append(xjunctionInner2)
    # nx.append(xjunctionOuter1)
    # nx.append(xjunctionInner1)
    # nx.append(xDaugh2WarpedList[elementsCountAround-1])
    # d2temp = [-1.0 * d2Daugh2WarpedList[elementsCountAround-3][j] for j in range(3)]
    # nd2.append(d2temp)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(d2Daugh2WarpedList[elementsCountAround-1])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    # d1junctionOuterList.append(d2temp)
    # d2junctionInnerList.append(smoothedd2[3])
    # d1junctionInnerList.append(smoothedd2[1])
    # d1junctionInnerList.append(smoothedd2[3])
    #
    # ##DERIVS
    # nx = []
    # nd2 = []
    # nx.append(xDaugh1WarpedList[elementsCountAround - 1])
    # nx.append(xjunctionInner2)
    # nx.append(xjunctionOuter2)
    # nx.append(xjunctionInner1)
    # nx.append(xDaugh1WarpedList[elementsCountAround - 3])
    # d2temp = [-1.0 * d2Daugh1WarpedList[elementsCountAround - 1][j] for j in range(3)]
    # nd2.append(d2temp)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(d2Daugh1WarpedList[elementsCountAround - 3])
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative=True, fixEndDerivative=True)
    # d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    # d1junctionOuterList.append(smoothedd2[2])
    # d2temp = [-1.0 * smoothedd2[1][j] for j in range(3)]
    # d2junctionInnerList.append(d2temp)
    #
    # ##DERIVS
    # nx = []
    # nd2 = []
    # nx.append(xParentWarpedList[n2t-2])
    # nx.append(xjunctionInner2)
    # nx.append(xjunctionOuter3)
    # nx.append(xjunctionInner1)
    # nx.append(xParentWarpedList[n2t-4])
    # nd2.append(d2ParentWarpedList[n2t-2])
    # nd2.append(zero)
    # nd2.append(zero)
    # nd2.append(zero)
    # d2temp = [(-1.0 * d2ParentWarpedList[n2t-4][j]) for j in range(3)]
    # nd2.append(d2temp)
    # smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    # d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    # d1junctionOuterList.append(d2temp)

    #JUNCTION OUTER
    #################
    nx = []
    nd2 = []
    xParentAlongSegment = xParentWarpedList[lastring+1]
    xDaugh1AlongSegment = xDaugh1WarpedList[2]
    xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    xjunctionOuter1 = [0.35*xmid[j]+0.65*xDaugh1AlongSegment[j] for j in range(3)]
    xjunctionOuterList.append(xjunctionOuter1)

    #deriv calc next
    nx.append(xParentAlongSegment)
    nx.append(xjunctionOuter1)
    nx.append(xDaugh1AlongSegment)
    nd2.append(d2ParentWarpedList[lastring+1])
    nd2.append(zero)
    nd2.append(d2Daugh1WarpedList[2])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2junctionOuterList.append(smoothedd2[1])

    nx = []
    nd2 = []
    xParentAlongSegment = xParentWarpedList[lastring+3]
    xDaugh2AlongSegment = xDaugh2WarpedList[2]
    xmid = [(xParentAlongSegment[j]+sxparent[n2][j])/2.0 for j in range(3)]
    xjunctionOuter2 = [0.35*xmid[j]+0.65*xDaugh2AlongSegment[j] for j in range(3)]
    xjunctionOuterList.append(xjunctionOuter2)
    #deriv calc next
    nx.append(xParentAlongSegment)
    nx.append(xjunctionOuter2)
    nx.append(xDaugh2AlongSegment)
    nd2.append(d2ParentWarpedList[lastring+3])
    nd2.append(zero)
    nd2.append(d2Daugh2WarpedList[2])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2junctionOuterList.append(smoothedd2[1])

    nx = []
    nd2 = []
    xDaugh1AlongSegment = xDaugh1WarpedList[0]
    xDaugh2AlongSegment = xDaugh2WarpedList[0]
    xjunctionOuter3 = [(7*(xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])+2*sxparent[n2][j])/16.0 for j in range(3)]
    xjunctionOuterList.append(xjunctionOuter3)
    nx.append(xDaugh1AlongSegment)
    nx.append(xjunctionOuter3)
    nx.append(xDaugh2AlongSegment)
    nd2.append(d1ParentWarpedList[0])
    nd2.append(zero)
    d2temp = [-1.0 * d2Daugh2WarpedList[0][j] for j in range(3)]
    nd2.append(d2temp)
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2junctionOuterList.append(smoothedd2[1])

    #JUNCTION INNER
    ##################
    xParentAlongSegment = xParentWarpedList[lastring]
    xDaugh1AlongSegment = xDaugh1WarpedList[1]
    xDaugh2AlongSegment = xDaugh2WarpedList[3]
    xjunctionInner1 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    xjunctionInnerList.append(xjunctionInner1)

    xParentAlongSegment = xParentWarpedList[lastring+2]
    xDaugh1AlongSegment = xDaugh1WarpedList[3]
    xDaugh2AlongSegment = xDaugh2WarpedList[1]
    xjunctionInner2 = [(xParentAlongSegment[j]+xDaugh1AlongSegment[j]+xDaugh2AlongSegment[j])/3.0 for j in range(3)]
    xjunctionInnerList.append(xjunctionInner2)

    ##DERIVS
    #-----------
    nx = []
    nd2 = []
    nx.append(xDaugh2WarpedList[elementsCountAround-3])
    nx.append(xjunctionInner2)
    nx.append(xjunctionOuter1)
    nx.append(xjunctionInner1)
    nx.append(xDaugh2WarpedList[elementsCountAround-1])
    d2temp = [-1.0 * d2Daugh2WarpedList[elementsCountAround-3][j] for j in range(3)]
    nd2.append(d2temp)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(d2Daugh2WarpedList[elementsCountAround-1])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    d1junctionOuterList.append(d2temp)
    d2temp = [-1.0 * smoothedd2[3][j] for j in range(3)]
    d1junctionInnerList.append(d2temp)
    d2temp = [-1.0 * smoothedd2[1][j] for j in range(3)]
    d1junctionInnerList.append(d2temp)

    ##DERIVS
    nx = []
    nd2 = []
    nx.append(xDaugh1WarpedList[elementsCountAround-1])
    nx.append(xjunctionInner2)
    nx.append(xjunctionOuter2)
    nx.append(xjunctionInner1)
    nx.append(xDaugh1WarpedList[elementsCountAround-3])
    d2temp = [-1.0 * d2Daugh1WarpedList[elementsCountAround-1][j] for j in range(3)]
    nd2.append(d2temp)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(d2Daugh1WarpedList[elementsCountAround-3])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    d1junctionOuterList.append(smoothedd2[2])

    ##DERIVS
    nx = []
    nd2 = []
    nx.append(xParentWarpedList[n2t-2])
    nx.append(xjunctionInner2)
    nx.append(xjunctionOuter3)
    nx.append(xjunctionInner1)
    nx.append(xParentWarpedList[n2t-4])

    #nd2.append(d2ParentWarpedList[n2t-2])
    d2temp = [(-1.0 * d2ParentWarpedList[n2t-2][j]) for j in range(3)]
    nd2.append(d2temp)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(zero)
    nd2.append(d2ParentWarpedList[n2t-4])
    smoothedd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    d2temp = [-1.0 * smoothedd2[2][j] for j in range(3)]
    d1junctionOuterList.append(d2temp)
    d2temp = [-1.0 * smoothedd2[3][j] for j in range(3)]
    d2junctionInnerList.append(d2temp)
    d2temp = [1.0 * smoothedd2[1][j] for j in range(3)]
    d2junctionInnerList.append(d2temp)

    ## Calculate unit d3
    #########################
    for n in range(len(xjunctionOuterList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1junctionOuterList[n]),
                                                       vector.normalise(d2junctionOuterList[n])))
        d3junctionOuterList.append(d3Unit)

    # Calculate unit d3 Daughter2
    #################################
    for n in range(len(xjunctionInnerList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1junctionInnerList[n]),
                                                       vector.normalise(d2junctionInnerList[n])))
        d3junctionInnerList.append(d3Unit)

    return xjunctionOuterList, xjunctionInnerList, \
           d1junctionOuterList, d2junctionOuterList, d3junctionOuterList, \
           d1junctionInnerList, d2junctionInnerList, d3junctionInnerList


def warpAirwaySegmentPoints(x1ListParent, x1ListDaugh1, x1ListDaugh2,
                            d1ListParent, d1ListDaugh1, d1ListDaugh2,
                            d2ListParent, d2ListDaugh1, d2ListDaugh2,
                            segmentAxisParent, segmentAxisDaughter1, segmentAxisDaughter2,
                            parentsegmentLength, daughter1segmentLength, daughter2segmentLength,
                            sxParent, sxDaugh1, sxDaugh2,
                            sd1Parent, sd1Daugh1, sd1Daugh2,
                            sd2Parent, sd2Daugh1, sd2Daugh2,
                            elementsCountAround, elementsCountAlongSegment,
                            nSegment,
                            ParentfaceMidPointZ, Daughter1faceMidPointZ, Daughter2faceMidPointZ):
    """
    Warps points in Airway segment to account for bending and twisting
    along central path defined by nodes sx and derivatives sd1 and sd2.
    :param x1List: coordinates of segment points.
    :param d1List: derivatives around axis of segment.
    :param d2List: derivatives along axis of segment.
    :param segmentAxis: axis perpendicular to segment plane.
    :param sx: coordinates of points on central path.
    :param sd1: derivatives of points along central path.
    :param sd2: derivatives representing cross axes.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :param faceMidPointZ: z-coordinate of midpoint for each element
    groups along the segment.
    :return coordinates and derivatives of warped points.
    """
    x1ParentWarpedList = []
    x1Daugh1WarpedList = []
    x1Daugh2WarpedList = []

    d1ParentWarpedList = []
    d1Daugh1WarpedList = []
    d1Daugh2WarpedList = []

    d2ParentWarpedList = []
    d2Daugh1WarpedList = []
    d2Daugh2WarpedList = []

    #smoothd2WarpedList = []
    d3ParentWarpedUnitList = []
    d3Daugh1WarpedUnitList = []
    d3Daugh2WarpedUnitList = []

    #parent
    #######
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = x1ListParent[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1ListParent[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2ListParent[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        xMid = [0.0, 0.0, ParentfaceMidPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1Parent[n2])
        cp = vector.crossproduct3(segmentAxisParent, unitTangent)
        dp = vector.dotproduct(segmentAxisParent, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxisParent, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
        else: # path tangent parallel to segment axis (z-axis)
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                axisRot = [1.0, 0, 0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                midRot = xMid
        translateMatrix = [sxParent[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2Parent[n2]
                    pt = [midRot[j] + sd2Parent[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2Parent[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            x1ParentWarpedList.append(xTranslate)
            d1ParentWarpedList.append(d1Rot2)
            d2ParentWarpedList.append(d2Rot2)

    #Daughter1
    ##########
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = x1ListDaugh1[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1ListDaugh1[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2ListDaugh1[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]

        xMid = [Daughter1faceMidPointZ[nAlongSegment], 0.0, 0.0]
        #xMid = [Daughter1faceMidPointZ[nAlongSegment], 0.0, ParentfaceMidPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1Daugh1[n2])
        cp = vector.crossproduct3(segmentAxisDaughter1, unitTangent)
        dp = vector.dotproduct(segmentAxisDaughter1, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxisDaughter1, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
        else: # path tangent parallel to segment axis (z-axis)
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                #axisRot = [1.0, 0, 0]
                axisRot = [0.0, 0, 1.0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                midRot = xMid
        translateMatrix = [sxDaugh1[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2Daugh1[n2]
                    pt = [midRot[j] + sd2Daugh1[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2Daugh1[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            x1Daugh1WarpedList.append(xTranslate)
            d1Daugh1WarpedList.append(d1Rot2)
            d2Daugh1WarpedList.append(d2Rot2)

    #Daughter2
    ############
    for nAlongSegment in range(elementsCountAlongSegment + 1):
        n2 = elementsCountAlongSegment * nSegment + nAlongSegment
        xElementAlongSegment = x1ListDaugh2[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d1ElementAlongSegment = d1ListDaugh2[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]
        d2ElementAlongSegment = d2ListDaugh2[elementsCountAround*nAlongSegment: elementsCountAround*(nAlongSegment+1)]

        xMid = [-Daughter2faceMidPointZ[nAlongSegment], 0.0, 0.0]
        #xMid = [-Daughter2faceMidPointZ[nAlongSegment], 0.0, ParentfaceMidPointZ[nAlongSegment]]

        # Rotate to align segment axis with tangent of central line
        unitTangent = vector.normalise(sd1Daugh2[n2])
        cp = vector.crossproduct3(segmentAxisDaughter2, unitTangent)
        dp = vector.dotproduct(segmentAxisDaughter2, unitTangent)
        if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
            print('path tnght not parallel to segment axis')
            axisRot = vector.normalise(cp)
            thetaRot = math.acos(vector.dotproduct(segmentAxisDaughter2, unitTangent))
            rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
            midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
        else: # path tangent parallel to segment axis (z-axis)
            #print('checking rotation angle for daughter 2')
            if dp == -1.0: # path tangent opposite direction to segment axis
                thetaRot = math.pi
                ###axisRot = [1.0, 0, 0]
                axisRot = [0.0, 0, 1.0]
                rotFrame = matrix.getRotationMatrixFromAxisAngle(axisRot, thetaRot)
                midRot = [rotFrame[j][0]*xMid[0] + rotFrame[j][1]*xMid[1] + rotFrame[j][2]*xMid[2] for j in range(3)]
            else: # segment axis in same direction as unit tangent
                #print('segment axis same direction as tngt for daugh2')
                midRot = xMid
        translateMatrix = [sxDaugh2[n2][j] - midRot[j] for j in range(3)]

        for n1 in range(elementsCountAround):
            x = xElementAlongSegment[n1]
            d1 = d1ElementAlongSegment[n1]
            d2 = d2ElementAlongSegment[n1]

            if vector.magnitude(cp)> 0.0: # path tangent not parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)]
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)]
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)]

                if n1 == 0:
                    # Project sd2 onto plane normal to sd1
                    v = sd2Daugh2[n2]
                    pt = [midRot[j] + sd2Daugh2[n2][j] for j in range(3)]
                    dist = vector.dotproduct(v, unitTangent)
                    ptOnPlane = [pt[j] - dist*unitTangent[j] for j in range(3)]
                    newVector = [ptOnPlane[j] - midRot[j] for j in range(3)]
                    # Rotate first point to align with planar projection of sd2
                    firstVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    thetaRot2 = math.acos(vector.dotproduct(vector.normalise(newVector), firstVector))
                    cp2 = vector.crossproduct3(vector.normalise(newVector), firstVector)
                    if vector.magnitude(cp2) > 0.0:
                        cp2 = vector.normalise(cp2)
                        signThetaRot2 = vector.dotproduct(unitTangent, cp2)
                        axisRot2 = unitTangent
                        rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, -signThetaRot2*thetaRot2)
                    else:
                        rotFrame2 = [ [1, 0, 0], [0, 1, 0], [0, 0, 1]]

            else: # path tangent parallel to segment axis
                xRot1 = [rotFrame[j][0]*x[0] + rotFrame[j][1]*x[1] + rotFrame[j][2]*x[2] for j in range(3)] if dp == -1.0 else x
                d1Rot1 = [rotFrame[j][0]*d1[0] + rotFrame[j][1]*d1[1] + rotFrame[j][2]*d1[2] for j in range(3)] if dp == -1.0 else d1
                d2Rot1 = [rotFrame[j][0]*d2[0] + rotFrame[j][1]*d2[1] + rotFrame[j][2]*d2[2] for j in range(3)] if dp == -1.0 else d2

                # Rotate to align start of elementsAround with sd2
                if n1 == 0:
                    v = vector.normalise(sd2Daugh2[n2])
                    startVector = vector.normalise([xRot1[j] - midRot[j] for j in range(3)])
                    axisRot2 = unitTangent
                    thetaRot2 = dp*-math.acos(vector.dotproduct(v, startVector))
                    rotFrame2 = matrix.getRotationMatrixFromAxisAngle(axisRot2, thetaRot2)

            xRot2 = [rotFrame2[j][0]*xRot1[0] + rotFrame2[j][1]*xRot1[1] + rotFrame2[j][2]*xRot1[2] for j in range(3)]
            d1Rot2 = [rotFrame2[j][0]*d1Rot1[0] + rotFrame2[j][1]*d1Rot1[1] + rotFrame2[j][2]*d1Rot1[2] for j in range(3)]
            d2Rot2 = [rotFrame2[j][0]*d2Rot1[0] + rotFrame2[j][1]*d2Rot1[1] + rotFrame2[j][2]*d2Rot1[2] for j in range(3)]
            xTranslate = [xRot2[j] + translateMatrix[j] for j in range(3)]

            x1Daugh2WarpedList.append(xTranslate)
            d1Daugh2WarpedList.append(d1Rot2)
            d2Daugh2WarpedList.append(d2Rot2)

    # # Smooth d2 for segment
    # ########################
    # smoothd2Raw = []
    # for n1 in range(elementsCountAround):
    #     nx = []
    #     nd2 = []
    #     for n2 in range(elementsCountAlongSegment + 1):
    #         n = n2*elementsCountAround + n1
    #         nx.append(xWarpedList[n])
    #         nd2.append(d2WarpedList[n])
    #     smoothd2 = interp.smoothCubicHermiteDerivativesLine(nx, nd2, fixStartDerivative = True, fixEndDerivative = True)
    #     smoothd2Raw.append(smoothd2)
    #
    # # Re-arrange smoothd2
    # for n2 in range(elementsCountAlongSegment + 1):
    #     for n1 in range(elementsCountAround):
    #         smoothd2WarpedList.append(smoothd2Raw[n1][n2])

    # Calculate unit d3 Parent
    ###########################
    for n in range(len(x1ParentWarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1ParentWarpedList[n]),
                                                       vector.normalise(d2ParentWarpedList[n])))
        d3ParentWarpedUnitList.append(d3Unit)

    for n in range(len(d1Daugh1WarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1Daugh1WarpedList[n]),
                                                       vector.normalise(d2Daugh1WarpedList[n])))
        d3Daugh1WarpedUnitList.append(d3Unit)

    for n in range(len(d1Daugh2WarpedList)):
        d3Unit = vector.normalise(vector.crossproduct3(vector.normalise(d1Daugh2WarpedList[n]),
                                                       vector.normalise(d2Daugh2WarpedList[n])))
        d3Daugh2WarpedUnitList.append(d3Unit)


    return x1ParentWarpedList, x1Daugh1WarpedList, x1Daugh2WarpedList, \
           d1ParentWarpedList, d1Daugh1WarpedList, d1Daugh2WarpedList, \
           d2ParentWarpedList, d2Daugh1WarpedList, d2Daugh2WarpedList, \
           d3ParentWarpedUnitList, d3Daugh1WarpedUnitList, d3Daugh2WarpedUnitList


def getAirwayJunctionCoordinatesFromInner(
        xjunctionOuter, d1junctionOuter, d2junctionOuter, d3junctionOuter,
        xjunctionInner, d1junctionInner, d2junctionInner, d3junctionInner,
        wallThicknessList, elementsCountAround,
        elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Generates coordinates from inner to outer surface using coordinates
    and derivatives of inner surface.
    :param xInner: Coordinates on inner surface
    :param d1Inner: Derivatives on inner surface around tube
    :param d2Inner: Derivatives on inner surface along tube
    :param d3Inner: Derivatives on inner surface through wall
    :param wallThicknessList: Wall thickness for each element along tube
    :param elementsCountAround: Number of elements around tube
    :param elementsCountAlong: Number of elements along tube
    :param elementsCountThroughWall: Number of elements through tube wall
    :param transitElementList: stores true if element around is a transition
    element that is between a big and a small element.
    return nodes and derivatives for mesh, and curvature along inner surface.
    """
    curvatureAroundInner = []
    curvatureAlong = []
    curvatureList = []

    xOuterList = []
    d1OuterList = []
    d2OuterList = []
    d3OuterList = []

    xInnerList = []
    d1InnerList = []
    d2InnerList = []
    d3InnerList = []

    for n3 in range(elementsCountThroughWall + 1):
        xi3 = 1 / elementsCountThroughWall * n3

        #junction outer (3 nodes)
        xOuter = []
        for n1 in range(len(xjunctionOuter)):
            wallThickness = wallThicknessList[n1]
            norm = d3junctionOuter[n1]

            # Calculate outer coordinates
            x = [xjunctionOuter[n1][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)

            innerx = xjunctionOuter[n1]
            outerx = xOuter[n1]
            dWall = [wallThickness * c for c in norm]

            # x
            x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
            xOuterList.append(x)

            # dx_ds1
            d1 = d1junctionOuter[n1]
            d1OuterList.append(d1)

            # dx_ds2
            d2 = d2junctionOuter[n1]
            d2OuterList.append(d2)

            # dx_ds3
            d3 = [c * wallThickness / elementsCountThroughWall for c in norm]
            d3OuterList.append(d3)

        #junction inner (2 middles, one front, one back)
        xOuter = []
        for n1 in range(len(xjunctionInner)):
            wallThickness = wallThicknessList[n1]
            norm = d3junctionInner[n1]

            # Calculate outer coordinates
            x = [xjunctionInner[n1][i] + norm[i] * wallThickness for i in range(3)]
            xOuter.append(x)

            innerx = xjunctionInner[n1]
            outerx = xOuter[n1]
            dWall = [wallThickness * c for c in norm]

            # x
            x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
            xInnerList.append(x)

            # dx_ds1
            d1 = d1junctionInner[n1]
            d1InnerList.append(d1)

            # dx_ds2
            d2 = d2junctionOuter[n1]
            d2InnerList.append(d2)

            # dx_ds3
            d3 = [c * wallThickness / elementsCountThroughWall for c in norm]
            d3InnerList.append(d3)

    return xOuterList, d1OuterList, d2OuterList, d3OuterList, \
           xInnerList, d1InnerList, d2InnerList, d3InnerList,\
           curvatureList


def getAirwaySegmentCoordinatesFromInner(
        xParentInner, xDaugh1Inner, xDaugh2Inner,
        d1ParentInner, d1Daugh1Inner, d1Daugh2Inner,
        d2ParentInner, d2Daugh1Inner, d2Daugh2Inner,
        d3ParentInner, d3Daugh1Inner, d3Daugh2Inner,
        wallThicknessList, elementsCountAround,
        elementsCountAlong, elementsCountThroughWall, transitElementList):
    """
    Generates coordinates from inner to outer surface using coordinates
    and derivatives of inner surface.
    :param xInner: Coordinates on inner surface
    :param d1Inner: Derivatives on inner surface around tube
    :param d2Inner: Derivatives on inner surface along tube
    :param d3Inner: Derivatives on inner surface through wall
    :param wallThicknessList: Wall thickness for each element along tube
    :param elementsCountAround: Number of elements around tube
    :param elementsCountAlong: Number of elements along tube
    :param elementsCountThroughWall: Number of elements through tube wall
    :param transitElementList: stores true if element around is a transition
    element that is between a big and a small element.
    return nodes and derivatives for mesh, and curvature along inner surface.
    """

    xOuter = []
    curvatureAroundInner = []
    curvatureAlong = []
    curvatureList = []

    xParentList = []
    d1ParentList = []
    d2ParentList = []
    d3ParentList = []

    xDaughter1List = []
    d1Daughter1List = []
    d2Daughter1List = []
    d3Daughter1List = []

    xDaughter2List = []
    d1Daughter2List = []
    d2Daughter2List = []
    d3Daughter2List = []

    #PARENT
    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3ParentInner[n]

            # Calculate outer coordinates
            x = [xParentInner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)

            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xParentInner[prevIdx], d1ParentInner[prevIdx], xParentInner[n], d1ParentInner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xParentInner[n], d1ParentInner[n], xParentInner[nextIdx], d1ParentInner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xParentInner[n], d2ParentInner[n], xParentInner[n + elementsCountAround], d2ParentInner[n + elementsCountAround], vector.normalise(d3ParentInner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xParentInner[n - elementsCountAround], d2ParentInner[n - elementsCountAround], xParentInner[n], d2ParentInner[n], vector.normalise(d3ParentInner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xParentInner[n - elementsCountAround], d2ParentInner[n - elementsCountAround], xParentInner[n], d2ParentInner[n], vector.normalise(d3ParentInner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xParentInner[n], d2ParentInner[n], xParentInner[n + elementsCountAround], d2ParentInner[n + elementsCountAround], vector.normalise(d3ParentInner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3ParentInner[n]
                innerx = xParentInner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xParentList.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1ParentInner[n]]
                d1ParentList.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xParentInner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2ParentInner[n]]
                d2ParentList.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3ParentList.append(d3)

    xOuter = []
    curvatureAroundInner = []
    curvatureAlong = []
    #DAUGHTER1
    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3Daugh1Inner[n]
            # Calculate outer coordinates
            x = [xDaugh1Inner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)

            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xDaugh1Inner[prevIdx], d1Daugh1Inner[prevIdx], xDaugh1Inner[n], d1Daugh1Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xDaugh1Inner[n], d1Daugh1Inner[n], xDaugh1Inner[nextIdx], d1Daugh1Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n], d2Daugh1Inner[n], xDaugh1Inner[n + elementsCountAround], d2Daugh1Inner[n + elementsCountAround], vector.normalise(d3Daugh1Inner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n - elementsCountAround], d2Daugh1Inner[n - elementsCountAround], xDaugh1Inner[n], d2Daugh1Inner[n], vector.normalise(d3Daugh1Inner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n - elementsCountAround], d2Daugh1Inner[n - elementsCountAround], xDaugh1Inner[n], d2Daugh1Inner[n], vector.normalise(d3Daugh1Inner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xDaugh1Inner[n], d2Daugh1Inner[n], xDaugh1Inner[n + elementsCountAround], d2Daugh1Inner[n + elementsCountAround], vector.normalise(d3Daugh1Inner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3Daugh1Inner[n]
                innerx = xDaugh1Inner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xDaughter1List.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1Daugh1Inner[n]]
                d1Daughter1List.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xDaugh1Inner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2Daugh1Inner[n]]
                d2Daughter1List.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3Daughter1List.append(d3)

    xOuter = []
    curvatureAroundInner = []
    curvatureAlong = []
    #DAUGHTER2
    ###########
    for n2 in range(elementsCountAlong + 1):
        wallThickness = wallThicknessList[n2]
        for n1 in range(elementsCountAround):
            n = n2*elementsCountAround + n1
            norm = d3Daugh2Inner[n]
            # Calculate outer coordinates
            x = [xDaugh2Inner[n][i] + norm[i]*wallThickness for i in range(3)]
            xOuter.append(x)
            # Calculate curvature along elements around
            prevIdx = n - 1 if (n1 != 0) else (n2 + 1)*elementsCountAround - 1
            nextIdx = n + 1 if (n1 < elementsCountAround - 1) else n2*elementsCountAround
            kappam = interp.getCubicHermiteCurvatureSimple(xDaugh2Inner[prevIdx], d1Daugh2Inner[prevIdx], xDaugh2Inner[n], d1Daugh2Inner[n], 1.0)
            kappap = interp.getCubicHermiteCurvatureSimple(xDaugh2Inner[n], d1Daugh2Inner[n], xDaugh2Inner[nextIdx], d1Daugh2Inner[nextIdx], 0.0)
            if not transitElementList[n1] and not transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = 0.5*(kappam + kappap)
            elif transitElementList[n1]:
                curvatureAround = kappam
            elif transitElementList[(n1-1)%elementsCountAround]:
                curvatureAround = kappap
            curvatureAroundInner.append(curvatureAround)

            # Calculate curvature along
            if n2 == 0:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n], d2Daugh2Inner[n], xDaugh2Inner[n + elementsCountAround], d2Daugh2Inner[n + elementsCountAround], vector.normalise(d3Daugh2Inner[n]), 0.0))
            elif n2 == elementsCountAlong:
                curvature = abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n - elementsCountAround], d2Daugh2Inner[n - elementsCountAround], xDaugh2Inner[n], d2Daugh2Inner[n], vector.normalise(d3Daugh2Inner[n]), 1.0))
            else:
                curvature = 0.5*(
                    abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n - elementsCountAround], d2Daugh2Inner[n - elementsCountAround], xDaugh2Inner[n], d2Daugh2Inner[n], vector.normalise(d3Daugh2Inner[n]), 1.0)) +
                    abs(interp.getCubicHermiteCurvature(xDaugh2Inner[n], d2Daugh2Inner[n], xDaugh2Inner[n + elementsCountAround], d2Daugh2Inner[n + elementsCountAround], vector.normalise(d3Daugh2Inner[n]), 0.0)))
            curvatureAlong.append(curvature)

        for n3 in range(elementsCountThroughWall + 1):
            xi3 = 1/elementsCountThroughWall * n3
            for n1 in range(elementsCountAround):
                n = n2*elementsCountAround + n1
                norm = d3Daugh2Inner[n]
                innerx = xDaugh2Inner[n]
                outerx = xOuter[n]
                dWall = [wallThickness*c for c in norm]
                # x
                x = interp.interpolateCubicHermite(innerx, dWall, outerx, dWall, xi3)
                xDaughter2List.append(x)

                # dx_ds1
                factor = 1.0 + wallThickness*xi3 * curvatureAroundInner[n]
                d1 = [ factor*c for c in d1Daugh2Inner[n]]
                d1Daughter2List.append(d1)

                # dx_ds2
                curvature = curvatureAlong[n]
                distance = vector.magnitude([x[i] - xDaugh2Inner[n][i] for i in range(3)])
                factor = 1.0 - curvature*distance
                d2 = [ factor*c for c in d2Daugh2Inner[n]]
                d2Daughter2List.append(d2)
                curvatureList.append(curvature)

                #dx_ds3
                d3 = [c * wallThickness/elementsCountThroughWall for c in norm]
                d3Daughter2List.append(d3)

    return xParentList, d1ParentList, d2ParentList, d3ParentList, \
           xDaughter1List, d1Daughter1List, d2Daughter1List, d3Daughter1List,\
           xDaughter2List, d1Daughter2List, d2Daughter2List, d3Daughter2List,\
           curvatureList


def createAirwaySegmentNodesAndElements\
                (region,
                 xParent, d1Parent, d2Parent, d3Parent,
                 xDaughter1, d1Daughter1, d2Daughter1, d3Daughter1,
                 xDaughter2, d1Daughter2, d2Daughter2, d3Daughter2,
                 xjunctionOuter, xjunctionInner, d1junctionOuter, d1junctionInner,
                 d2junctionOuter, d2junctionInner,
                 d3junctionOuter, d3junctionInner,
                 elementsCountAround, elementsCountAlong, elementsCountThroughWall,
                 firstNodeIdentifier, firstElementIdentifier,
                 useCubicHermiteThroughWall,
                 useCrossDerivatives):

    #    annotationGroups, annotationArray,
    """
    :param xList: coordinates of centerline points.
    :param d1List: derivatives along axis of segment.
    :param radius1List: derivatives along axis of segment.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :return coordinates and derivatives of warped points.
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [0.0, 0.0, 0.0]

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
    if useCubicHermiteThroughWall:
        nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS3, 1)
        if useCrossDerivatives:
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS2DS3, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1)

    mesh = fm.findMeshByDimension(3)

    eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)


    if useCubicHermiteThroughWall:
        eftfactory = eftfactory_tricubichermite(mesh, useCrossDerivatives)
    else:
        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)

    eft1 = eftfactory.createEftBasic()



    elementtemplateStandard = mesh.createElementtemplate()
    elementtemplateStandard.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    elementtemplateX = mesh.createElementtemplate()
    elementtemplateX.setElementShapeType(Element.SHAPE_TYPE_CUBE)

    result = elementtemplate.defineField(coordinates, -1, eft1)


    mapDerivatives = False
    mapEndDerivatives = True
    mapStartDerivatives = False


    # if mapDerivatives:
    #     eft1 = eftfactory.createEftNoCrossDerivatives()
    #     setEftScaleFactorIds(eft1, [1], [])
    #
    #     elementtemplateX.defineField(coordinates, -1, eft1)
    #     elementtemplate = elementtemplateX
    # else:
    #     eft1 = eft
        # elementtemplate = elementtemplate

    # element = mesh.createElement(elementIdentifier, elementtemplate)


    #
    #
    # if mapDerivatives:
    #     result3 = element.setScaleFactors(eft1, [-1.0])


    ###################
    # Create nodes
    ##################

    # Create nodes for Parent
    # Coordinates field
    for n in range(len(xParent)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xParent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Parent[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes for Daughter1
    # Coordinates field
    print('length of parent with thickness = ',len(xDaughter1))

    for n in range(len(xDaughter1)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaughter1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daughter1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daughter1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Daughter1[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes for Daughter2
    # Coordinates field
    for n in range(len(xDaughter2)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaughter2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daughter2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daughter2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3Daughter2[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Coordinates field
    for n in range(len(xjunctionOuter)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xjunctionOuter[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1junctionOuter[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2junctionOuter[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3junctionOuter[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Coordinates field
    for n in range(len(xjunctionInner)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xjunctionInner[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1junctionInner[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2junctionInner[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS3, 1, d3junctionInner[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS2DS3, 1, zero)
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D3_DS1DS2DS3, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    #######################
    # create elements
    ########################

    now = elementsCountAround*(elementsCountThroughWall+1)
    for e2 in range(elementsCountAlong):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                bni11 = e2*now + e3*elementsCountAround + e1 + 1
                bni12 = e2*now + e3*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                bni21 = e2*now + (e3 + 1)*elementsCountAround + e1 + 1
                bni22 = e2*now + (e3 + 1)*elementsCountAround + (e1 + 1) % elementsCountAround + 1
                nodeIdentifiers = [ bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now ]

#               onOpening = e1 > elementsCountAround - 2
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

    now = elementsCountAround * (elementsCountThroughWall + 1)
    offset = (elementsCountAlong+1) * elementsCountAround * (elementsCountThroughWall + 1)
    for e2 in range(elementsCountAlong):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                bni11 = e2 * now + e3 * elementsCountAround + e1 + (offset+1)
                bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + (offset+1)
                bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + (offset+1)
                bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + (offset+1)
                nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

    now = elementsCountAround * (elementsCountThroughWall + 1)
    offset = 2 * (elementsCountAlong+1) * elementsCountAround * (elementsCountThroughWall + 1)
    for e2 in range(elementsCountAlong):
        for e3 in range(elementsCountThroughWall):
            for e1 in range(elementsCountAround):
                bni11 = e2 * now + e3 * elementsCountAround + e1 + (offset+1)
                bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + (offset+1)
                bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + (offset+1)
                bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + (offset+1)
                nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
                element = mesh.createElement(elementIdentifier, elementtemplate)
                element.setNodesByIdentifier(eft1, nodeIdentifiers)
                elementIdentifier = elementIdentifier + 1

    ### FORM JUNCTION ELEMS
    ################################
    # now = (elementsCountThroughWall+1)*elementsCountAround * (elementsCountThroughWall + 1)
    # offset = 3 * (elementsCountAlong+1) * elementsCountAround * (elementsCountThroughWall + 1)
    # for e2 in range(len(xjunctionOuter)+len(xjunctionInner)):
    #     for e3 in range(elementsCountThroughWall):
    #         bni11 = e2 * now + e3 * elementsCountAround + e1 + (offset+1)
    #         bni12 = e2 * now + e3 * elementsCountAround + (e1 + 1) % elementsCountAround + (offset+1)
    #         bni21 = e2 * now + (e3 + 1) * elementsCountAround + e1 + (offset+1)
    #         bni22 = e2 * now + (e3 + 1) * elementsCountAround + (e1 + 1) % elementsCountAround + (offset+1)
    #         nodeIdentifiers = [bni11, bni12, bni11 + now, bni12 + now, bni21, bni22, bni21 + now, bni22 + now]
    #         element = mesh.createElement(elementIdentifier, elementtemplate)
    #         element.setNodesByIdentifier(eft, nodeIdentifiers)
    #         elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [17, 18, 79, 73, 21, 22, 81, 76]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [79, 73, 26, 27, 81, 76, 30, 31]

    # startDerivativesMap[0][1]=(None, (1, 1, 0), None)
    # startDerivativesMap[1][1]=(None, (1, 1, 0), None)

    # for i in range(2):
        # lns = [1, 5] if (i == 0) else [2, 6]
        # for n3 in range(2):
            # derivativesMap = startDerivativesMap[n3][e1]
            # handle different d1 on each side of node
    # if mapDerivatives:
    #     d2Map = (1,1,0)
    #     lns = [1,5]
    #     for n3 in range(2):
    #         ln = [lns[n3]]
    #         remapEftNodeValueLabel(eft1, ln, Node.VALUE_LABEL_D_DS2, \
    #                            derivativeSignsToExpressionTerms(
    #                                (Node.VALUE_LABEL_D_DS1, Node.VALUE_LABEL_D_DS2, Node.VALUE_LABEL_D_DS3), d2Map))
    #         elementtemplateX.defineField(coordinates, -1, eft1)
    #         elementtemplate = elementtemplateX


    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1



    nodeIdentifiers = [18, 19, 73, 80, 22, 23, 76, 82]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [73, 80, 27, 28, 76, 82, 31, 32]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [20, 17, 74, 79, 24, 21, 77, 81]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [74, 79, 51, 52, 77, 81, 55, 56]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [19, 20, 80, 74, 23, 24, 82, 77]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [80, 74, 50, 51, 82, 77, 54, 55]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [75, 79, 25, 26, 78, 81, 29, 30]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [79, 75, 52, 49, 81, 78, 56, 53]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [80, 75, 28, 29, 82, 78, 32, 29]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    nodeIdentifiers = [75, 80, 49, 50, 78, 82, 53, 54]
    element = mesh.createElement(elementIdentifier, elementtemplate)
    element.setNodesByIdentifier(eft1, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return nodeIdentifier, elementIdentifier



def createAirwaySegmentSurfaceNodesAndElements(region,
                                  xParent, d1Parent, d2Parent,
                                  xDaugh1, d1Daugh1, d2Daugh1,
                                  xDaugh2, d1Daugh2, d2Daugh2,
                                  xjunctionOuter, xjunctionInner,
                                               d1junctionOuter, d1junctionInner,
                                               d2junctionOuter, d2junctionInner,
                                  elementsCountAround, elementsCountAlong,
                                  firstNodeIdentifier, firstElementIdentifier,
                                  useCrossDerivatives):
    #    annotationGroups, annotationArray,
    """
    :param xList: coordinates of centerline points.
    :param d1List: derivatives along axis of segment.
    :param radius1List: derivatives along axis of segment.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :return coordinates and derivatives of warped points.
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [0.0, 0.0, 0.0]

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

    #    mesh = fm.findMeshByDimension(3)
    #    eftfactory = eftfactory_bicubichermite(mesh, useCrossDerivatives)
    #    eft = eftfactory.createEftBasic()

    mesh = fm.findMeshByDimension(2)
    bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
    if not useCrossDerivatives:
        for n in range(4):
            eft.setFunctionNumberOfTerms(n * 4 + 4, 0)

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    cache = fm.createFieldcache()
    print('coming to write nodes in create surface node')

    #####################
    # CREATE NODES
    #####################
    # Create nodes - PARENT
    # Coordinates field
    for n in range(len(xParent)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xParent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Parent[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1


    # Create nodes - DauGH1
    # Coordinates field
    for n in range(len(xDaugh1)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaugh1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daugh1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daugh1[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1


    # Create nodes - DAUGH2
    # Coordinates field
    for n in range(len(xDaugh2)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaugh2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daugh2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daugh2[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes - JUNCTION
    nodeIdentifierOuter = []
    nodeIdentifierInner = []

    for n in range(len(xjunctionOuter)):
        nodeIdentifierOuter.append(nodeIdentifier)

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xjunctionOuter[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1junctionOuter[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2junctionOuter[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1


    for n in range(len(xjunctionInner)):
        nodeIdentifierInner.append(nodeIdentifier)

        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xjunctionInner[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1junctionInner[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2junctionInner[n])
        if useCrossDerivatives:
            coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    #####################
    # CREATE ELEMENTS
    #####################
    # # create elements - Parent
    # ##################
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2 * elementsCountAround + e1 + 1
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + 1
            nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    # # create elements - Daughter1
    # ##############################
    node0 = (elementsCountAround) * (elementsCountAlong + 1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2 * elementsCountAround + e1 + node0
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + node0
            nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    # # create elements - Daughter2
    # ##########################3##
    node0 = 2 * (elementsCountAround) * (elementsCountAlong + 1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2 * elementsCountAround + e1 + node0
            bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + node0
            nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1


    # # create elements - junction
    # ##########################3##
    # for e2 in range(len(xjunctionInner)):
    #     nodeinner = nodeIdentifierInner[e2]
    #     bnp1 = (elementsCountAround)*(elementsCountAlong)+1
    #     bnd1 = (elementsCountAround)*(elementsCountAlong+1)+1
    #     bnd2 = 2*(elementsCountAround)*(elementsCountAlong+1)+1
    #
    #     bndparentouter = bnp1 if e2==0 else bnp1+2
    #     bnddaugh1outer = bnd1+1 if e2==0 else bnd1+3
    #     bnddaugh2outer = bnd2+1 if e2==0 else bnd2+3
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [bndparentouter, 10, nodeinner, 37]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [39, nodeinner, 13, bnddaugh1outer]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [12, bndparentouter, 38, nodeinner]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [nodeinner, 38, bnddaugh2outer, 27]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [39, nodeinner, 25, bnddaugh2outer]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #
    #     element = mesh.createElement(elementIdentifier, elementtemplate)
    #     nodeIdentifiers = [nodeinner, 37, bnddaugh1outer, 15]
    #     result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #     elementIdentifier = elementIdentifier + 1
    #

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [9, 10, 40, 37]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1


    # element.setNodesByIdentifier(eftApex1, nodeIdentifiers)
    # # set general linear map coefficients
    # radiansAround = e1 * radiansPerElementAround
    # radiansAroundNext = ((e1 + 1) % elementsCountAround) * radiansPerElementAround
    # scalefactors = [
    #     -1.0,
    #     math.sin(radiansAround), math.cos(radiansAround), radiansPerElementAround,
    #     math.sin(radiansAroundNext), math.cos(radiansAroundNext), radiansPerElementAround
    # ]
    # result = element.setScaleFactors(eftApex1, scalefactors)
    # elementIdentifier = elementIdentifier + 1



    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [40, 37, 14, 15]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1
    #------------------
    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [10, 11, 37, 41]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [37, 41, 15, 16]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1
    #--------------
    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [12, 9, 38, 40]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [38, 40, 27, 28]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1
    #--------------
    # JUNCTION BACK
    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [11, 12, 41, 38]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [41, 38, 26, 27]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1
    #--------------

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [39, 40, 13, 14]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [41, 39, 16, 13]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1
    #------------

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [40, 39, 28, 25]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1

    element = mesh.createElement(elementIdentifier, elementtemplate)
    nodeIdentifiers = [39, 41, 25, 26]
    result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    elementIdentifier = elementIdentifier + 1


    # node0 = 3 * (elementsCountAround) * (elementsCountAlong + 1) + 1
    # for e2 in range(elementsCountAlong):
    #     for e1 in range(elementsCountAround):
    #         element = mesh.createElement(elementIdentifier, elementtemplate)
    #         bni1 = e2 * elementsCountAround + e1 + node0
    #         bni2 = e2 * elementsCountAround + (e1 + 1) % elementsCountAround + node0
    #         nodeIdentifiers = [bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround]
    #         result = element.setNodesByIdentifier(eft, nodeIdentifiers)
    #         elementIdentifier = elementIdentifier + 1

    fm.endChange()

    return nodeIdentifier, elementIdentifier




def createSurfaceNodesAndElements(region,
                                  xParent, d1Parent, d2Parent,
                                  xDaugh1, d1Daugh1, d2Daugh1,
                                  xDaugh2, d1Daugh2, d2Daugh2,
                                  elementsCountAround, elementsCountAlong,
                                  firstNodeIdentifier, firstElementIdentifier,
                                  useCrossDerivatives):
    #    annotationGroups, annotationArray,
    """
    :param xList: coordinates of centerline points.
    :param d1List: derivatives along axis of segment.
    :param radius1List: derivatives along axis of segment.
    :param elementsCountAround: Number of elements around segment.
    :param elementsCountAlongSegment: Number of elements along segment.
    :param nSegment: Segment index along central path.
    :return coordinates and derivatives of warped points.
    """

    nodeIdentifier = firstNodeIdentifier
    elementIdentifier = firstElementIdentifier
    zero = [ 0.0, 0.0, 0.0 ]

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


    #    mesh = fm.findMeshByDimension(3)
#    eftfactory = eftfactory_bicubichermite(mesh, useCrossDerivatives)
#    eft = eftfactory.createEftBasic()

    mesh = fm.findMeshByDimension(2)
    bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
    eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
    if not useCrossDerivatives:
        for n in range(4):
            eft.setFunctionNumberOfTerms(n * 4 + 4, 0)

    elementtemplate = mesh.createElementtemplate()
    elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
    result = elementtemplate.defineField(coordinates, -1, eft)

    cache = fm.createFieldcache()

    # radiansPerElementAround = 2.0 * math.pi / elementsCountAround
    # x = [0.0, 0.0, 0.0]
    # dx_ds1 = [0.0, 0.0, 0.0]
    # dx_ds2 = [0.0, 0.0, 1.0 / elementsCountAlong]
    # zero = [0.0, 0.0, 0.0]

    # #parent
    # radius = radius1list[0]
    # for n2 in range(elementsCountAlong + 1):
    #     x[2] = n2 / elementsCountAlong
    #     for n1 in range(elementsCountAround):
    #         radiansAround = n1 * radiansPerElementAround
    #         cosRadiansAround = math.cos(radiansAround)
    #         sinRadiansAround = math.sin(radiansAround)
    #         x[0] = radius * cosRadiansAround
    #         x[1] = radius * sinRadiansAround
    #         dx_ds1[0] = radiansPerElementAround * radius * -sinRadiansAround
    #         dx_ds1[1] = radiansPerElementAround * radius * cosRadiansAround
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, dx_ds1)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, dx_ds2)
    #         if useCrossDerivatives:
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
    #         nodeIdentifier = nodeIdentifier + 1

    # # Central path
    # cx = [ [ 0.0, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
    # cd1 = [ [ segmentLength, 0.0, 0.0 ], [ segmentLength, 0.0, 0.0 ] ]
    # cd2 = [ [ 0.0, 1.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]
    # cd12 = [ [0.0, 0.0, 0.0 ], [ 0.0, 0.0, 0.0 ] ]
    #
    # # Sample central path
    # sx, sd1, se, sxi, ssf = interp.sampleCubicHermiteCurves(cx, cd1, elementsCountAlong)
    # sd2 = interp.interpolateSampleCubicHermite(cd2, cd12, se, sxi, ssf)[0]

    # ###daughter 1
    # radius = radius1list[1]
    # for n2 in range(elementsCountAlong + 1):
    #     x[2] = n2 / elementsCountAlong
    #     for n1 in range(elementsCountAround):
    #         radiansAround = n1 * radiansPerElementAround
    #         cosRadiansAround = math.cos(radiansAround)
    #         sinRadiansAround = math.sin(radiansAround)
    #         x[0] = radius * cosRadiansAround
    #         x[1] = radius * sinRadiansAround
    #         node = nodes.createNode(nodeIdentifier, nodetemplate)
    #         cache.setNode(node)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, x)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, zero)
    #         coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, zero)
    #         if useCrossDerivatives:
    #             coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
    #         nodeIdentifier = nodeIdentifier + 1

    print('coming to write nodes in create surface node')

    # Create nodes - PARENT
    # Coordinates field
    for n in range(len(xParent)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xParent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Parent[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Parent[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes - DauGH1
    # Coordinates field
    for n in range(len(xDaugh1)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaugh1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daugh1[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daugh1[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1

    # Create nodes - DAUGH2
    # Coordinates field
    for n in range(len(xDaugh2)):
        node = nodes.createNode(nodeIdentifier, nodetemplate)
        cache.setNode(node)
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_VALUE, 1, xDaugh2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS1, 1, d1Daugh2[n])
        coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D_DS2, 1, d2Daugh2[n])
        if useCrossDerivatives:
                coordinates.setNodeParameters(cache, -1, Node.VALUE_LABEL_D2_DS1DS2, 1, zero)
        nodeIdentifier = nodeIdentifier + 1


    # # create elements - Parent
    # ##################
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2*elementsCountAround + e1 + 1
            bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + 1
            nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1


    # # create elements - Daughter1
    # ##############################
    node0 = (elementsCountAround)*(elementsCountAlong+1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2*elementsCountAround + e1 + node0
            bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + node0
            nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1

    # # create elements - Daughter2
    # ##########################3##
    node0 = 2*(elementsCountAround)*(elementsCountAlong+1) + 1
    for e2 in range(elementsCountAlong):
        for e1 in range(elementsCountAround):
            element = mesh.createElement(elementIdentifier, elementtemplate)
            bni1 = e2*elementsCountAround + e1 + node0
            bni2 = e2*elementsCountAround + (e1 + 1)%elementsCountAround + node0
            nodeIdentifiers = [ bni1, bni2, bni1 + elementsCountAround, bni2 + elementsCountAround ]
            result = element.setNodesByIdentifier(eft, nodeIdentifiers)
            elementIdentifier = elementIdentifier + 1


    fm.endChange()

    return nodeIdentifier,elementIdentifier

