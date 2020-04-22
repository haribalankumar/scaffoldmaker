"""
Generates a 3-D bifurcation mesh with variable numbers of elements around and up.
"""

from __future__ import division
import math
import copy
from opencmiss.utils.zinc.field import findOrCreateFieldCoordinates
from opencmiss.zinc.element import Element, Elementbasis
from opencmiss.zinc.field import Field
from opencmiss.zinc.node import Node

from scaffoldmaker.annotation.annotationgroup import AnnotationGroup
#from scaffoldmaker.meshtypes.meshtype_3d_ostium1 import MeshType_3d_ostium1

from scaffoldmaker.meshtypes.scaffold_base import Scaffold_base
from scaffoldmaker.scaffoldpackage import ScaffoldPackage
##from scaffoldmaker.utils.annulusmesh import createAnnulusMesh3d
##from scaffoldmaker.utils.eftfactory_bicubichermitelinear import eftfactory_bicubichermitelinear

from scaffoldmaker.utils.tubebifurcationmesh import CreateTubeBifurcationMesh
from scaffoldmaker.utils import interpolation as interp
from scaffoldmaker.utils.interpolation import getCubicHermiteBasis, smoothCubicHermiteDerivativesLine
from scaffoldmaker.utils import vector

from scaffoldmaker.utils.meshrefinement import MeshRefinement
#from scaffoldmaker.utils.tracksurface import TrackSurface, TrackSurfacePosition, calculate_surface_axes


class MeshType_3d_airwaybifurcation1(Scaffold_base):
    '''
    3-D Airway Bifurcation scaffold.
    '''

    @staticmethod
    def getName():
        return '3D Bifurcation template 1'

    @staticmethod
    def getDefaultOptions(cls, parameterSetName='Default'):
        options = {
            'Number of elements along': 2,
            'Number of elements around': 4,  # should be even
            'Use cross derivatives': False
        }
        return options

    @staticmethod
    def getOrderedOptionNames():
        return [
            'Number of elements along',
            'Number of elements around',
            'Use cross derivatives'
        ]

    @staticmethod
    def checkOptions(cls, options):
        if (options['Number of elements along'] < 4):
            options['Number of elements along'] = 4
        if (options['Number of elements around'] > 8):
            options['Number of elements around'] = 8
            if (options['Number of elements around'] < 4):
                options['Number of elements around'] = 4
        elif (options['Number of elements around'] % 2) == 1:
            options['Number of elements around'] += 1

    @staticmethod
    def generateMesh(region, options):
        '''
        Generate the base bicubic Hermite mesh. See also generateMesh().
        :param region: Zinc region to define model in. Must be empty.
        :param options: Dict containing options. See getDefaultOptions().
        :return: list of AnnotationGroup
        '''
        elementsCountAlong = options['Number of elements along']
        elementsCountAround = options['Number of elements around']

        useCrossDerivatives = options['Use cross derivatives']

#        ostiumOptions = options['Ureter']
#        ostiumDefaultOptions = ostiumOptions.getScaffoldSettings()

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
            nodetemplate = nodes.createNodetemplate()
            nodetemplate.defineField(coordinates)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_VALUE, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS1, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D_DS2, 1)
            nodetemplate.setValueNumberOfVersions(coordinates, -1, Node.VALUE_LABEL_D2_DS1DS2, 1)

#        mesh = fm.findMeshByDimension(3)
#        eftfactory = eftfactory_bicubichermitelinear(mesh, useCrossDerivatives)
#        eft = eftfactory.createEftBasic()

#        elementtemplate = mesh.createElementtemplate()
#        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_CUBE)
#        elementtemplate.defineField(coordinates, -1, eft)

        mesh = fm.findMeshByDimension(2)
        bicubicHermiteBasis = fm.createElementbasis(2, Elementbasis.FUNCTION_TYPE_CUBIC_HERMITE)
        eft = mesh.createElementfieldtemplate(bicubicHermiteBasis)
        if not useCrossDerivatives:
            for n in range(4):
                eft.setFunctionNumberOfTerms(n*4 + 4, 0)
        elementtemplate = mesh.createElementtemplate()
        elementtemplate.setElementShapeType(Element.SHAPE_TYPE_SQUARE)
        result = elementtemplate.defineField(coordinates, -1, eft)

        cache = fm.createFieldcache()

#        neckGroup = AnnotationGroup(region, 'neck of bladder', FMANumber='unknown', lyphID='unknown')
#        bodyGroup = AnnotationGroup(region, 'body of bladder', FMANumber='unknown', lyphID='unknown')
#        annotationGroups = [neckGroup, bodyGroup]

#        neckMeshGroup = neckGroup.getMeshGroup(mesh)
#        bodyMeshGroup = bodyGroup.getMeshGroup(mesh)

        # create nodes
        nodeIdentifier = 1
        elementIdentifier = 1
        radiansPerElementAround = 2.0*math.pi/elementsCountAround
        radiansPerElementAlong = (math.pi/4)/elementsCountAlong


        x1List = []
        d1List = []
        radius1list = []

        nodeIdentifier, elementIdentifier  = \
            TubeBifurcationMesh(region, x1List,d1List,radius1list,
                                elementsCountAround,elementsCountAlong,
                                nodeIdentifier, elementIdentifier)

        fm.endChange()
#        return annotationGroups

