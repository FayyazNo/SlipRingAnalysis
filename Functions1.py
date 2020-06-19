
# This Code is Written By Fayyaz Nosouhi
# Tehran University November 2017
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import numpy






#####Parameters:
##nH=12
##n=12
##L1=170
##L2=139
##L3=10
##L4=6
##L5=20
##L6=25
##L7=12  #Redundant data L7=(L1-n*L8-2*L12)/(n-1)+L8
##L8=4
##L9=27
##L10=10
##L11=16
##L12=3.5
##
##L16=1.2
##L17=150
##L18=2.5
##L19=25
##
##L21=175.5
##R1=120
##R2=150
##R3=172.5
##R4=191
##R5=211
##R6=170
##R7=60
##R8=2
##R9=119
##R10=16.9
##D1=13
##D2=13
##D3=18 # diameter of the bolts whach has distance L21
##M=10
## 
##
###----------------------------------------------------------------
##Eshaft=208000
##vshaft=0.3
##Sysh=[[750., 0], [790, .03], [815, .05],[840, .08],[850, 0.1], ]
##Den_sh=7850e-12
##Exp_sh=1.21e-5
##k_sh=36
##
##Eslp=208000
##vslp=0.3
##Syslp=[[325, 0], [420, .04], [470, .07],[520, .11],[550, 0.15], ]
##Den_slp=7850e-12
##Exp_slp=1.2e-5
##k_slp=50
##
##Eins=24000
##vins=0.33
##Syins=[[1000, 0], [1100, .1], [1300, .2],[1500, .3],]
##Den_ins=1900e-12
##Exp_ins=1e-5
##k_ins=0.3
###---------------------------------------------------------------------
##
###---------------------------------------------------------------------
##FC=0.15 #Friction Coeffecient
##prsscon=[[50/1000, 3.5], [91, 50], [200, 100], ] #Pressure vs. thermal conductance
##Mass=1
##
##Tshaft=40
##teta1=73
##teta2=78
##teta3=80
##teta4=82
##teta5=76
##teta6=73
##teta7=71
##teta8=68
##T=[]
##for i in range (4*n+1):
##    T.append(80+i*.05)
##
##    
##w=314  # Rotational Speed
###----------------------------------------------------------------------

execfile('Var.py')
xL17= 1.5 #raio of the length of the shaft to L17
Rx=R3-(1.0*L5/L3)*(R3-R2);
L7=(L1-n*L8-2*L12)/(n-1)+L8
L13=0
L14=0
L15=0
L20=19
stptime_static=1
stptime_couple=1

#singpntCriteria='Linear'
#singpntCriteria='Quadratic'
##linDistance =[5, 15]
##quadDistance=[4, 8, 12]

#----------------------------------------------------------------------
Model='Model1'
mdb.Model(name=Model, modelType=STANDARD_EXPLICIT)


import os
os.chdir(directory)


ang2="global"
ang3="global"
ang4="global"
ang5="global"
ang6="global"
ang7="global"

ang2=asin(L21/(2.*R6-L9))
ang3=0
ang4=ang3+2.*ang2
ang5=atan(D3/(R6-L9/2.))
ang6=atan((D3/2.)/(R6-L9/2.))
ang7=ang4+ang5
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def myPart (nH, n, L1, L2, L3, L4, L5, L6, L7,L8 ,L9, L10, L11, L12, L13, L14,L15,\
L16, L17, L18, L19, L20, L21, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, D1, D2, D3, xL17, Model):
    
    Rx=R3-(1.0*L5/L3)*(R3-R2);      tet=atan(L5/(R3-Rx)); Lhole=(R4-R2)/sin(tet);           r_L17=1.5;
    # Slip-Ring==========================================================================
    s1 = mdb.models[Model].ConstrainedSketch(name='__profile__', sheetSize=1000.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.Line(point1=(R1, L5), point2=(Rx, L5))
    s1.Line(point1=(Rx, L5), point2=(R3, 0.0))
    s1.Line(point1=(R3, 0.0), point2=(R5, 0.0))
    s1.Line(point1=(R5, 0.0), point2=(R5, L12))
    s1.Line(point1=(R5, L12), point2=(R5-(L11-R8), L12))
    s1.Arc3Points(point1=(R5-(L11-R8), L12), point2=(R5-(L11-R8), L12+L8), point3=(R5-(L11-R8)-R8, L12+L8/2.))
    s1.Line(point1=(R5-(L11-R8), L12+L8), point2=(R5, L12+L8))
    LL=L12+L8
    for i in range(n-1):
        l1=L1-2*L12
        lt=(l1-n*L8)/(n-1)
        Dt=L11-R8
        s1.Line(point1=(R5, LL), point2=(R5, LL+lt))
        s1.Line(point1=(R5, LL+lt), point2=(R5-Dt,  LL+lt))
        s1.Arc3Points(point1=(R5-Dt,  LL+lt), point2=(R5-Dt,  LL+lt+L8), point3=(R5-Dt-R8,  LL+lt+L8/2.))
        s1.Line(point1=(R5-Dt,  LL+lt+L8), point2=(R5,  LL+lt+L8))    
        LL=(i+1)*(lt+L8)+L12+L8
        print(LL)
         
    s1.Line(point1=(R5,  LL), point2=(R5, LL+L12))
    s1.Line(point1=(R5, LL+L12), point2=(R6, LL+L12))
    s1.Line(point1=(R6, LL+L12), point2=(R6, LL+L12-L4))
    s1.Line(point1=(R6, LL+L12-L4), point2=(R6-L9, LL+L12-L4))
    s1.Line(point1=(R6-L9, LL+L12-L4), point2=(R1+L10, LL+L12-L6))
    s1.Line(point1=(R1+L10, LL+L12-L6), point2=(R1, LL+L12-L6))
    s1.Line(point1=(R1, LL+L12-L6), point2=(R1, L5))
    p = mdb.models[Model].Part(name='Part-1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models[Model].sketches['__profile__']
    #mdb.models[Model].parts['Part-1'].setValues(geometryRefinement=EXTRA_FINE)
    # Prat 1 attributs----------------------------
    p11 = mdb.models[Model].parts['Part-1']
    f11 = p11.faces
    e11 = p11.edges
    d11 = p11.datums
    c11 = p11.cells
    v11 = p11.vertices
    #---------------------------------------------
    #two bolt hole----------------------------------
    E1=e11.findAt((R6-L9/2., L1-L4, 0.),)
    F1=f11.findAt((R6-L9/2.,L1-L4, 0.),)
    t = p.MakeSketchTransform(sketchPlane=F1, sketchUpEdge=E1, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, L1-L4,  0.0))
    s = mdb.models[Model].ConstrainedSketch(name='__profile__', sheetSize=918.59, gridSpacing=22.96, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p11.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    #--------------------------------------------------------------------------------------------------
##    s.CircleByCenterPerimeter(center=(0.0, R6-L9/2), point1=(0, R6-L9/2+M/2))
##    g1=g.findAt((0, R6-L9/2+M/2),)
##    s.radialPattern(geomList=(g1, ), vertexList=(), number=nH-3, totalAngle=360.0-4*360/nH, centerPoint=(0.0, 0.0))
##    s.radialPattern(geomList=(g1, ), vertexList=(), number=2, totalAngle=360.0-2*360/nH, centerPoint=(0.0, 0.0))
    #----------------------------------------------------------------------------------------------------
    X1=(R6-L9/2.)*sin(ang3)
    Z1=(R6-L9/2.)*cos(ang3)
    X2=( R6-L9/2.+D3/2.)*sin(ang3)
    Z2=( R6-L9/2.+D3/2.)*cos(ang3)
    s.CircleByCenterPerimeter(center=(X1, Z1), point1=(X2, Z2))
    #------------------------------------
    X1=(R6-L9/2)*sin(ang4)
    Z1=(R6-L9/2)*cos(ang4)
    X2=( R6-L9/2+D3/2)*sin(ang4)
    Z2=( R6-L9/2+D3/2)*cos(ang4)
    s.CircleByCenterPerimeter(center=(X1, Z1), point1=(X2, Z2))
    E1=e11.findAt((R6-L9/2., L1-L4, 0.),)
    F1=f11.findAt((R6-L9/2.,L1-L4, 0.),)
    p11.CutExtrude(sketchPlane=F1, sketchUpEdge=E1, sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, depth=L19, 
        flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__']
    #---------------------------------------------
    E1=e11.findAt((R5-R5/10000.,0,0),)
    F1=f11.findAt((R5-R5/10000.,0,0),)
    #Longtitude holes created by pattern---------------------------------------------
    t = p.MakeSketchTransform(sketchPlane=F1, sketchUpEdge=E1, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
    s = mdb.models[Model].ConstrainedSketch(name='__profile__', sheetSize=989.94, gridSpacing=24.74, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p11.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    s.CircleByCenterPerimeter(center=(0.0, -R4), point1=(0, -R4-D1/2.))
    gg=g.findAt((0, -R4-D1/2. ),)
    s.radialPattern(geomList=(gg, ), vertexList=(), number=nH, totalAngle=360.0, centerPoint=(0.0, 0.0))
    p11.CutExtrude(sketchPlane=F1, sketchUpEdge=E1, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, depth=L2, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__']
    #Inclined Holes-----------------------------------------------------------
    for i in range (0,nH):
        i=i+1
        RR=R4-D1/2.
        X=RR*cos(i*2.*pi/nH)
        Z=RR*sin(i*2.*pi/nH)
        X1=RR*cos((i-.5)*2.*pi/nH)
        Z1=RR*sin((i-.5)*2.*pi/nH)
        E1=e11.findAt((X, 0, Z),)
        highlight(E1)
        hh=p11.DatumPlaneByLinePoint(line=d11[1], point=p11.InterestingPoint(edge=E1, rule=CENTER)).id
        pickedCells = c11.findAt((X1, 0, Z1),)
        highlight(pickedCells)
        print (i)
        t = p11.MakeSketchTransform(sketchPlane=d11[hh], sketchUpEdge=d11[1], sketchPlaneSide=SIDE1, origin=(0.0, 0, 0.0))
        s = mdb.models[Model].ConstrainedSketch(name='__profile__', sheetSize=1434.22, gridSpacing=35.85, transform=t)
        s.setPrimaryObject(option=SUPERIMPOSE)
        p11.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
        s.rectangle(point1=(-2*R5, 10*L1), point2=(0.0, -10*L1))
        p11.PartitionCellBySketch(sketchPlane=d11[hh], sketchUpEdge=d11[1], cells=pickedCells, sketch=s)
        s.unsetPrimaryObject()
        del mdb.models[Model].sketches['__profile__']
#--------------------------------------------------------------------
    for i in range (nH):
        Z1=(R1+R1/1000.)*sin(i*2.*pi/nH)
        X1=(R1+R1/1000.)*cos(i*2.*pi/nH)
        E1=e11.findAt((X1,L5,Z1),)
        Z2=R1*sin(i*2.*pi/nH)
        X2=R1*cos(i*2.*pi/nH)
        E2=e11.findAt((X2,L1/2.,Z2),)
        #highlight(e1)
        #highlight(e2)
        ss=p.DatumAxisByRotation(line=E1, axis=E2, angle=90.0).id
        Z3=Rx*sin(i*2.*pi/nH)
        X3=Rx*cos(i*2.*pi/nH)
        V1=v11.findAt((X3, L5, Z3))
        #highlight(v1)
        ss1=p11.DatumAxisByParToEdge(edge=d11[ss], point=V1).id
        Z4=R3*sin(i*2.*pi/nH)
        X4=R3*cos(i*2.*pi/nH)
        V2=v11.findAt((X4,0,Z4),)
        #highlight(v2)
        ss2=p11.DatumPlaneByLinePoint(line=d11[ss1], point=V2).id
        d = p11.datums
        X=R2*cos(i*2.*pi/nH)
        Z=R2*sin(i*2.*pi/nH)
        t = p11.MakeSketchTransform(sketchPlane=d11[ss2], sketchUpEdge=d11[1], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(X, L3, Z))
        s = mdb.models[Model].ConstrainedSketch(name='__profile__',  sheetSize=1379.71, gridSpacing=34.49, transform=t)
        s.setPrimaryObject(option=SUPERIMPOSE)
        s.CircleByCenterPerimeter(center=(0, 0), point1=(0, D2/2.))
        p11.CutExtrude(sketchPlane=d11[ss2], sketchUpEdge=d11[1], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, depth=Lhole, flipExtrudeDirection=ON)
        p11.CutExtrude(sketchPlane=d11[ss2], sketchUpEdge=d11[1], sketchPlaneSide=SIDE1, 
        sketchOrientation=RIGHT, sketch=s, depth=200, flipExtrudeDirection=OFF)
    s.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__']
    #--------------------------------------------------------------------------------------
    ## two partition for preload in surface of ring
    #
    sketchPlane1=f11.findAt(((R6-L9/2.)*cos(-ang5), L1-L4, (R6-L9/2.)*sin(-ang5)), )
    sketchUpEdge1=e11.findAt((1.01*R6, L1, 0), )
    t = p11.MakeSketchTransform(sketchPlane=sketchPlane1,\
                    sketchUpEdge=sketchUpEdge1, sketchPlaneSide=SIDE1, origin=(0.0, L1-L4, 0.0))
    s1 = mdb.models[Model].ConstrainedSketch(name='__profile__',  sheetSize=744.37, gridSpacing=18.6, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    #
    s1.Line(point1=(0, 0), point2=(1.2*R6*sin(-ang5), 1.2*R6*cos(-ang5)))
    s1.Line(point1=(0, 0), point2=(1.2*R6*sin(ang4+ang5), 1.2*R6*cos(ang4+ang5)))
    #
    pickedFaces_1 = f11.findAt(((R6-L9/2.)*cos(-ang5),L1-L4 ,(R6-L9/2.)*sin(-ang5)), )
    pickedFaces_2 = f11.findAt(((R6-L9/2.)*cos(ang7+.01),L1-L4 ,(R6-L9/2.)*sin(ang7+.01)), )
    Faces=(pickedFaces_1, pickedFaces_2, )
    sketchUpEdge2=e11.findAt((.999*R6, L1-L4, 0), )
    p11.PartitionFaceBySketch(sketchUpEdge=sketchUpEdge2, faces=Faces, sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__']
    #
    pickedFaces_1 = f11.findAt(((R6-L9/2.)*cos(-ang5),L1-L4 ,(R6-L9/2.)*sin(-ang5)), )
    pickedFaces_2 = f11.findAt(((R6-L9/2.)*cos(ang7+.01),L1-L4 ,(R6-L9/2.)*sin(ang7+.01)), )
    Faces=(pickedFaces_1, pickedFaces_2, )
    ## Datum plane for creating section in depth L20 for preload
    face_hole1 = f11.findAt(((R6-L9/2.),(L1-L4)*.99 ,(R6-L9/2.)*tan(-ang6)), )
    face_hole2 = f11.findAt(((R6-L9/2.),(L1-L4)*.99 ,(R6-L9/2.)*tan(+ang6)), )
    face_hole3 = f11.findAt(((R6-L9/2.)*cos(ang4)-D3/2.,(L1-L4)*.99 ,(R6-L9/2.)*sin(ang4)), )
    face_hole4 = f11.findAt(((R6-L9/2.)*cos(ang4)+D3/2.,(L1-L4)*.99 ,(R6-L9/2.)*sin(ang4)), )
    #
    if (face_hole3==face_hole4):
        allHoleFaces_part=(face_hole1, face_hole2, face_hole3, )
    else:
        allHoleFaces_part=(face_hole1, face_hole2, face_hole3, face_hole4, )
    #
    face_ref=f11.findAt(((R6-L9/2.)*cos(ang5*2),(L1-L4) ,(R6-L9/2.)*sin(ang5*2)), )  
    datumPlaneForHole=p11.DatumPlaneByOffset(plane=face_ref, flip=SIDE2, offset=L20)
    p11.PartitionFaceByDatumPlane(datumPlane=p11.datums[datumPlaneForHole.id], faces=allHoleFaces_part)    
    # Shaft=====================================================================
    s1 = mdb.models[Model].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.Line(point1=(R9, L18), point2=(R9+L16, L18))
    s1.Line(point1=(R9+L16, L18), point2=(R9+L16, L18-((xL17-1)/2.)*L17))
    s1.Line(point1=(R9+L16, L18-((xL17-1)/2.)*L17), point2=(R7, L18-((xL17-1)/2.)*L17))
    s1.Line(point1=(R7, L18-((xL17-1)/2.)*L17), point2=(R7, L18+((xL17+1)/2.)*L17))
    s1.Line(point1=(R7, L18+((xL17+1)/2.)*L17), point2=(R9+L16, L18+((xL17+1)/2.)*L17))
    s1.Line(point1=(R9+L16, L18+((xL17+1)/2.)*L17), point2=(R9+L16, L18+L17))
    s1.Line(point1=(R9+L16, L18+L17), point2=(R9, L18+L17))
    s1.Line(point1=(R9, L18+L17), point2=(R9, L18))
    p22 = mdb.models[Model].Part(name='Shaft', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p22 = mdb.models[Model].parts['Shaft']
    p22.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__']
    #--------------------------------------------------------------
    # Prat 2 attributes----------------------------
    p22 = mdb.models[Model].parts['Shaft']
    f22 = p22.faces
    e22 = p22.edges
    d22 = p22.datums
    c22 = p22.cells
    LL=((1-xL17)*L17)/2.+2.5
    t = p22.MakeSketchTransform(sketchPlane=f22.findAt((R7+.01*R7,LL, 0), ),\
                                sketchUpEdge=e22.findAt((R7*cos(.01), LL, R7*sin(.01)), ), sketchPlaneSide=SIDE1, origin=(0.0, LL, 0.0))
    #-----------------------------------------------------
    s1 = mdb.models[Model].ConstrainedSketch(name='__profile__',  sheetSize=744.37, gridSpacing=18.6, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p22.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    #------------------------------------------------------
    for i in range (nH):
        s1.Line(point1=(.2*R7*cos(i*2.*pi/nH), .2*R7*sin(i*2.*pi/nH)), point2=(1.2*R9*cos(i*2.*pi/nH), 1.2*R9*sin(i*2.*pi/nH)))
    #------------------------------------------------------   
    pickedFaces = f22.findAt((R7+.01*R7,LL, 0), )
    p22.PartitionFaceBySketch(sketchUpEdge=e22.findAt((R7*cos(.01), LL, R7*sin(.01)), ), faces=pickedFaces, sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__'] 
    for i in range (1, nH+1):
        pickedCells = c22.findAt((R7*cos(i*2.*pi/nH), L17/2, R7*sin(i*2.*pi/nH) ), )
        pickedEdges = e22.findAt(((R7+L16/2.)*cos(i*2.*pi/nH), LL, (R7+L16/2.)*sin(i*2.*pi/nH) ),)
        highlight(pickedEdges)
        p22.PartitionCellByExtrudeEdge(line=e22.findAt((R7, L17/2, 0),), cells=pickedCells, edges=pickedEdges, 
            sense=FORWARD)
    # Prat 2 attributes----------------------------
    p22 = mdb.models[Model].parts['Shaft']
    f22 = p22.faces
    e22 = p22.edges
    d22 = p22.datums
    c22 = p22.cells
    v22 = p22.vertices
    #---------------------------------------------
    # Insulation=====================================================================
    s1 = mdb.models[Model].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.Line(point1=(R9, L18), point2=(R9+L16, L18))
    s1.Line(point1=(R9+L16, L18), point2=(R9+L16, L18+L17))
    s1.Line(point1=(R9+L16, L18+L17), point2=(R9, L18+L17))
    s1.Line(point1=(R9, L18+L17), point2=(R9, L18))
    p33 = mdb.models[Model].Part(name='Insulation', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p33 = mdb.models[Model].parts['Insulation']
    p33.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models[Model].sketches['__profile__']
    # Prat 3 attributes----------------------------
    p33 = mdb.models[Model].parts['Insulation']
    f33 = p33.faces
    e33 = p33.edges
    d33 = p33.datums
    c33 = p33.cells
    v33 = p33.vertices
    t = p33.MakeSketchTransform(sketchPlane=f33.findAt((R9+.001, L17+L18, 0), ),\
                                sketchUpEdge=e33.findAt((R9*cos(.01), L17+L18, R9*sin(.01)), ), sketchPlaneSide=SIDE1, origin=(0.0, L17+L18, 0.0))
    #-----------------------------------------------------
    s1 = mdb.models[Model].ConstrainedSketch(name='__profile__',  sheetSize=744.37, gridSpacing=18.6, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=SUPERIMPOSE)
    p33.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    #------------------------------------------------------
    for i in range (nH+1):
        s1.Line(point1=(.2*R9*cos(i*2.*pi/nH), .2*R9*sin(i*2.*pi/nH)), point2=(1.2*R9*cos(i*2.*pi/nH), 1.2*R9*sin(i*2.*pi/nH)))
    #------------------------------------------------------  
    pickedFaces = f33.findAt((R9+.001, L17+L18, 0), )
    p33.PartitionFaceBySketch(sketchUpEdge=e33.findAt((R9*cos(.01), L17+L18, R9*sin(.01)), ), faces=pickedFaces, sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models[Model].sketches['__profile__'] 
    #------------------------------------------------------
    for i in range (nH):
        pickedCells = c33.findAt((R9*cos(i*2.*pi/nH), L17/2, R9*sin(i*2.*pi/nH) ), )
        pickedEdges = e33.findAt(((R9+L16/2.)*cos(i*2.*pi/nH), L17+L18, (R9+L16/2.)*sin(i*2.*pi/nH) ),)
        highlight(pickedEdges)
        p33.PartitionCellByExtrudeEdge(line=e33.findAt((R9, L17/2, 0),), cells=pickedCells, edges=pickedEdges, 
            sense=REVERSE)
    # Prat 3 attributes----------------------------
    p33 = mdb.models[Model].parts['Insulation']
    f33 = p33.faces
    e33 = p33.edges
    d33 = p33.datums
    c33 = p33.cells
    v33 = p33.vertices
    #
    #Partitioning around the start and end of the contact on the all 3 parts
    for i in range (nH):
         X1=(R9+L16)*cos(i*2.*pi/nH)
         Z1=(R9+L16)*sin(i*2.*pi/nH)
         X2=(R9)*cos(i*2.*pi/nH)
         Z2=(R9)*sin(i*2.*pi/nH)
         X3=(R7)*cos(i*2.*pi/nH)
         Z3=(R7)*sin(i*2.*pi/nH)
         X4=(R1)*cos(i*2.*pi/nH)
         Z4=(R1)*sin(i*2.*pi/nH)     
         #Partitioning side 1 & 2  Insulation on upper surface
         edge1INS=e33.findAt((X1, L1/2 , Z1),)
         p33.PartitionEdgeByParam(edges=edge1INS, parameter=1-(L5-L18)/L17)    
         edge1INS=e33.findAt((X1, L1/2 , Z1),)  
         p33.PartitionEdgeByParam(edges=edge1INS, parameter=1-(L17+L18-L1+L6)/(L17+L18-L5))
         #Partitioning side 1 & 2  Insulation on lower surface
         edge2INS=e33.findAt((X2, L1/2 , Z2),)
         p33.PartitionEdgeByParam(edges=edge2INS, parameter=1-(L5-L18)/L17)
         edge2INS=e33.findAt((X2, L1/2 , Z2),)
         p33.PartitionEdgeByParam(edges=edge2INS, parameter=1-(L17+L18-L1+L6)/(L17+L18-L5))
         #
         #Partitioning side 1 & 2  SLR on inner surface
         edge1SLR=e11.findAt((X4, L1/2 , Z4),)
         p11.PartitionEdgeByParam(edges=edge1SLR, parameter=1-(L5-L18)/L17)
         edge1SLR=e11.findAt((X4, L1/2 , Z4),)
         p11.PartitionEdgeByParam(edges=edge1SLR, parameter=1-(L17+L18-L1+L6)/(L17+L18-L5))
         #
         #Partitioning side 1 & 2  Shaft on upper surface
         edge1SHFT=e22.findAt((X2, L1/2 , Z2),)
         p22.PartitionEdgeByParam(edges=edge1SHFT, parameter=1-(L5-L18)/L17)    
         edge1SHFT=e22.findAt((X2, L1/2 , Z2),)  
         p22.PartitionEdgeByParam(edges=edge1SHFT, parameter=1-(L17+L18-L1+L6)/(L17+L18-L5))
         #Partitioning side 1 & 2  Shaft on lower surface
         edge2SHFT=e22.findAt((X3, L1/2 , Z3),)
         p22.PartitionEdgeByParam(edges=edge2SHFT, parameter=1-(L5-L18)/L17)
         edge2SHFT=e22.findAt((X3, L1/2 , Z3),)
         p22.PartitionEdgeByParam(edges=edge2SHFT, parameter=1-(L17+L18-L1+L6)/(L17+L18-L5))
#
    #Creat sets For singular points and for internal surface of SlipRing
    set1=[]
    #
    for i in range (nH):
        X=R1*cos(i*2.*pi/nH+0.0001)
        Z=R1*sin(i*2.*pi/nH+0.0001)
        faces = f11.findAt(((X, L1/2, Z),))
        set1.append(faces)
        highlight(faces)
    #
    setSurfSLR=p11.Set(faces=set1, name='Set_contact')
    tet=atan(L5/(R3-Rx))
    Lx= R4-D1/2.-R2
    Lxx=Lx-D2/(2.*cos(tet))
    LsingPnt=L3+Lxx/tan(tet)
    set2=[]
    # Creat set for treatment of singular points stress
    for i in range (nH):
        X=(R4-D1/2.)*cos(i*2.*pi/nH)
        Z=(R4-D1/2.)*sin(i*2.*pi/nH)
        vertices1 = v11.findAt(((X, LsingPnt, Z),))
        vertices2 = v11.findAt(((X, 0, Z),))
        set2.append(vertices1)
        set2.append(vertices2)
    # set of all nodes in intersection of the two holes on slip ring+ ponts in begining of hole (each to points create a path
    #from beginig of the hole to the singular points)
    p11.Set(vertices=set2, name='singHoleIntersect')
    #
    set2=[]
    for i in range (nH):
        X=(R1)*cos(i*2.*pi/nH)
        Z=(R1)*sin(i*2.*pi/nH)
        vertices1 = v11.findAt(((X, L5, Z),))
        vertices2 = v11.findAt(((X, L1-L6, Z),))
        set2.append(vertices1)
        set2.append(vertices2)
    # set of nodes in the begining and the end of the contact in slip ring. 
    p11.Set(vertices=set2, name='singContactSLR')
    #
    set2=[]
    for i in range (nH):
        X=(R9+L16)*cos(i*2.*pi/nH)
        Z=(R9+L16)*sin(i*2.*pi/nH)
        vertices1 = v33.findAt(((X, L5, Z),))
        vertices2 = v33.findAt(((X, L1-L6, Z),))
        set2.append(vertices1)
        set2.append(vertices2)
    # set of nodes in the begining and the end of the contact in insulation. 
    p33.Set(vertices=set2, name='singContactINS')
    #
    return p11, f11, d11, c11, e11, v11, p22, f22, d22, c22, e22, v22,  p33, f33, d33, c33, e33, v33


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def myMaterial(Eshaft, vshaft, Sysh, Den_sh, Exp_sh, k_sh,Eins, vins, Den_ins, Exp_ins, k_ins, Eslp, vslp, Syslp,\
               Den_slp, Exp_slp, k_slp, R7, L1, R9,nH, L2, Model, c11 ,c22, c33, p11, p22, p33):
    #Material=
    mdb.models[Model].Material(name='S_26NiCrMoV14-5')
    mdb.models[Model].materials['S_26NiCrMoV14-5'].Elastic(table=((Eshaft, vshaft), ))
    mdb.models[Model].materials['S_26NiCrMoV14-5'].Plastic(table=((Sysh)))
    mdb.models[Model].materials['S_26NiCrMoV14-5'].Density(table=((Den_sh, ), ))
    mdb.models[Model].materials['S_26NiCrMoV14-5'].Expansion(table=((Exp_sh, ), ))
    mdb.models[Model].materials['S_26NiCrMoV14-5'].Conductivity(table=((k_sh, ), ))
    mdb.models[Model].HomogeneousSolidSection(name='shaft', material='S_26NiCrMoV14-5', thickness=None)
#----------------------------------------------------------
    mdb.models[Model].Material(name='EPGC 203')
    mdb.models[Model].materials['EPGC 203'].Elastic(table=((Eins, vins), ))
        #mdb.models[Model].materials['EPGC 203'].Plastic(table=((Syins)))
    mdb.models[Model].materials['EPGC 203'].Density(table=((Den_ins, ), ))
    mdb.models[Model].materials['EPGC 203'].Expansion(table=((Exp_ins, ), ))
    mdb.models[Model].materials['EPGC 203'].Conductivity(table=((k_ins, ), ))
    mdb.models[Model].HomogeneousSolidSection(name='insulation', material='EPGC 203', thickness=None)
#-----------------------------------------------------------
    mdb.models[Model].Material(name='S355J2N')
    mdb.models[Model].materials['S355J2N'].Elastic(table=((Eslp, vslp), ))
    mdb.models[Model].materials['S355J2N'].Plastic(table=((Syslp)))
    mdb.models[Model].materials['S355J2N'].Density(table=((Den_slp, ), ))
    mdb.models[Model].materials['S355J2N'].Expansion(table=((Exp_slp, ), ))
    mdb.models[Model].materials['S355J2N'].Conductivity(table=((k_slp, ), ))
    mdb.models[Model].HomogeneousSolidSection(name='slipring', material='S355J2N', thickness=None)
#---------------------------------------------------------------
    p22.SectionAssignment(region=(c22 ,), sectionName='shaft', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    p33.SectionAssignment(region=(c33 , ), sectionName='insulation', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    cells1=[]
    for i in range (nH):
        Z1=(R1+0.01)*sin(i*2*pi/nH+0.1)
        X1=(R1+0.01)*cos(i*2*pi/nH+.1)
        cells1.append(c11.findAt((X1, L2/2., Z1),))
        highlight(c11.findAt((X1, L2/2., Z1),))
    cells1=(cells1, )
    p11.SectionAssignment(region=cells1, sectionName='slipring', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def myAssembly (p11, p22, p33, Model):
    rootAssembly = mdb.models[Model].rootAssembly
    rootAssembly.DatumCsysByDefault(CARTESIAN)
    rootAssembly.Instance(name='Part-1-1', part=p11, dependent=OFF)
    rootAssembly.Instance(name='Part-2-1', part=p22, dependent=OFF)
    rootAssembly.Instance(name='Part-3-1', part=p33, dependent=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(datumAxes=OFF, datumPlanes=OFF)
    session.viewports['Viewport: 1'].enableMultipleColors()
    cmap = session.viewports['Viewport: 1'].colorMappings['Part instance']
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].disableMultipleColors() 
    return rootAssembly

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def myStep(stptime_static, Model, no_step):
    if 'Step-1' in mdb.models[Model].steps:
        del mdb.models[Model].steps['Step-1']
    #
    if 'Step-2' in mdb.models[Model].steps:
        del mdb.models[Model].steps['Step-2']
    #
    if 'Step-3' in mdb.models[Model].steps:
        del mdb.models[Model].steps['Step-3']
    #
    mdb.models[Model].StaticStep(name='Step-1', previous='Initial',timePeriod=stptime_static)
    mdb.models[Model].CoupledTempDisplacementStep(name='Step-2', previous='Step-1', response=STEADY_STATE, deltmx=None, cetol=None, creepIntegration=None, amplitude=RAMP)
    mdb.models[Model].steps['Step-2'].setValues(stabilizationMagnitude=0.0005, stabilizationMethod=DISSIPATED_ENERGY_FRACTION, continueDampingFactors=False, adaptiveDampingRatio=0.08)
    mdb.models[Model].FieldOutputRequest(name='F-Output-1', 
        createStepName='Step-1', variables=('S', 'MISES', 'U', 'PEEQ', 'CSTATUS'))
    mdb.models[Model].FieldOutputRequest(name='F-Output-2', 
        createStepName='Step-2', variables=('S', 'MISES', 'U', 'PEEQ', 'NT', 'CSTATUS'))
    #
    if no_step==3:
        mdb.models[Model].CoupledTempDisplacementStep(name='Step-3', previous='Step-2', response=STEADY_STATE, deltmx=None, cetol=None, creepIntegration=None, amplitude=RAMP)
        mdb.models[Model].steps['Step-3'].setValues(stabilizationMagnitude=0.0005, stabilizationMethod=DISSIPATED_ENERGY_FRACTION, continueDampingFactors=False, adaptiveDampingRatio=0.08)   
        mdb.models[Model].FieldOutputRequest(name='F-Output-3', 
        createStepName='Step-3', variables=('S', 'MISES', 'U', 'PEEQ', 'NT', 'CSTATUS'))
        mdb.models[Model].FieldOutputRequest(name='F-Output-3', 
             createStepName='Step-3', variables=('S', 'MISES', 'U', 'PEEQ', 'NT', 'CSTATUS'))



        
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def myInteraction (R9, L17,R6, L9, L1, L4, FC, presscon, RAssmb, Model ):
    s1 = RAssmb.instances['Part-3-1'].faces
    s_tie2=[]
    for i in range (nH):
        highlight(s1.findAt((R9*cos(i*2.*pi/nH+0.0005), L17/2, R9*sin(i*2.*pi/nH+0.0005)),))
        s_tie2.append(s1.findAt(((R9*cos(i*2.*pi/nH+0.0005), L17/2, R9*sin(i*2.*pi/nH+0.0005)),)))
    #----------------------------------------------------
    region2=RAssmb.Surface(side1Faces=s_tie2, name='s_tie2')
    #----------------------------------------------------
    s1 = RAssmb.instances['Part-2-1'].faces
    s_tie1=[]
    for i in range (nH):
        highlight(s1.findAt((R9*cos(i*2.*pi/nH+0.0005), L17/2, R9*sin(i*2.*pi/nH+0.0005)),))
        s_tie1.append(s1.findAt(((R9*cos(i*2.*pi/nH+0.0005), L17/2, R9*sin(i*2.*pi/nH+0.0005)),)))
    #----------------------------------------------------
    region1=RAssmb.Surface(side1Faces=s_tie1, name='s_tie1')
    #----------------------------------------------------
    mdb.models[Model].Tie(name='Constraint-2', master=region1, slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    #-----------------------------------------------------
    mdb.models[Model].ContactProperty('IntProp-1')
    mdb.models[Model].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((FC, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mdb.models[Model].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    mdb.models[Model].interactionProperties['IntProp-1'].ThermalConductance(
        definition=TABULAR, clearanceDependency=OFF, pressureDependency=ON, 
        temperatureDependencyP=OFF, massFlowRateDependencyP=OFF, 
        dependenciesP=0, pressureDepTable=prsscon)
    #-------------------------------------------------------------------
    s1 = RAssmb.instances['Part-3-1'].faces
    s_con2=[]
    for i in range (nH):
        highlight(s1.findAt(((R9+L16)*cos(i*2.*pi/nH+0.0005), L17/2, (R9+L16)*sin(i*2.*pi/nH+0.0005)),))
        s_con2.append(s1.findAt((((R9+L16)*cos(i*2.*pi/nH+0.0005), L17/2, (R9+L16)*sin(i*2.*pi/nH+0.0005)),)))
    #----------------------------------------------------
    region2=RAssmb.Surface(side1Faces=s_con2, name='s_con2')
    #----------------------------------------------------
    s1 = RAssmb.instances['Part-1-1'].faces
    s_con1=[]
    for i in range (nH):
        Z1=(R1)*sin(i*2.*pi/nH+0.0005)
        X1=(R1)*cos(i*2.*pi/nH+0.0005)
        s_con1.append(s1.findAt(((X1,L1/2.,Z1),)))
    #-------------------------------------------------
    region1=RAssmb.Surface(side1Faces=s_con1, name='s_con1')
    #---------------------------------------------------------------------
    mdb.models[Model].SurfaceToSurfaceContactStd(name='Int-1', 
        createStepName='Initial', master=region1, slave=region2, 
        sliding=FINITE, thickness=ON, interactionProperty='IntProp-1', 
        adjustMethod=NONE, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)
    #----------------------------------------------------
    #Partioning shaft
    s1 = RAssmb.instances['Part-2-1'].faces
    c1 = RAssmb.instances['Part-2-1'].cells
    for i in range (nH):
        pickedCells = c1.findAt((R9*sin(i*2.*pi/nH+0.000001), L17/2., R9*cos(i*2.*pi/nH+0.000001)), )
        pickedFaces1 = s1.findAt(((R9+R9/100000.)*sin(i*2.*pi/nH+0.000001), L17+L18, (R9+R9/100000.)*cos(i*2.*pi/nH+0.000001)), )
        RAssmb.PartitionCellByExtendFace(extendFace=pickedFaces1, cells=pickedCells)
        pickedFaces2 = s1.findAt(((R9+R9/100000.)*sin(i*2.*pi/nH+0.000001), L18, (R9+R9/100000.)*cos(i*2.*pi/nH+0.000001)), )
        pickedCells = c1.findAt((R9*sin(i*2.*pi/nH+0.000001), L17/2., R9*cos(i*2.*pi/nH+0.000001)), )
        RAssmb.PartitionCellByExtendFace(extendFace=pickedFaces2, cells=pickedCells)
    #------------------------------------------------------------
    # define reference points for conections mass
    X1=(R6-L9/2.)*sin(ang3)
    Z1=(R6-L9/2.)*cos(ang3)
    X2=(R6-L9/2)*sin(ang4)
    Z2=(R6-L9/2)*cos(ang4)
    #------------------------------------
    rp1=RAssmb.ReferencePoint(point=(Z1, L1-L4, X1)).id
    rp2=RAssmb.ReferencePoint(point=(Z2, L1-L4, X2)).id
    r1 = RAssmb.referencePoints
    refPoints1=(r1[rp1], )
    refPoints2=(r1[rp2], )
    #--------------------------------------
    RAssmb.engineeringFeatures.PointMassInertia(
        name='Inertia-1', region=refPoints1, mass=Mass, i11=i11, i22=i22, i33=i33, 
        alpha=0.0, composite=0.0)
    RAssmb.engineeringFeatures.PointMassInertia(
        name='Inertia-2', region=refPoints2, mass=Mass, i11=i11, i22=i22, i33=i33, 
        alpha=0.0, composite=0.0)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4


def myBoundaryCondition (nH, n, R1, R3, Rx, R5,  R6, R7, R9, R10, L1, L2,  L4, L5, L6, L8, L9, L10, L11,L12,L17, L18, xL17,D3, teta1,\
                         teta2, teta3, teta4, teta5,teta6, teta7, teta8,Tshaft,  T ,w, w_Tol, no_step , Model, RAssmb):
    
    f1 = RAssmb.instances['Part-1-1'].faces
    k2=[]
    for i in range (nH):
        Z1=(R1+.1)*sin(i*2*pi/nH+0.0005)
        X1=(R1+.1)*cos(i*2*pi/nH+0.0005)
        k2.append(f1.findAt((X1, L5, Z1),))
        #highlight(f1.findAt((X1, L5, Z1),))
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta1', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta1, amplitude=UNSET)
#--------------------------------------
    angle=atan(L5/float(R3-Rx))
    k2=[]
    for i in range (nH):
        Z1=(R3-1)*sin(2*i*pi/nH+0.001)
        X1=(R3-1)*cos(2*i*pi/nH+.001)
        k2.append(f1.findAt((X1, tan(angle), Z1),))
        #highlight(f1.findAt((X1, tan(angle), Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta2', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta2, amplitude=UNSET)
#-----------------------------------
    k2=[]
    for i in range (nH):
        Z1=(R5-.001)*sin(i*2*pi/nH+0.001)
        X1=(R5-.001)*cos(i*2*pi/nH+.001)
        k2.append(f1.findAt((X1, 0, Z1),))
        #highlight(f1.findAt((X1, 0, Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta3', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta3, amplitude=UNSET)
#-----------------------------------
    k2=[]
    for i in range (nH):
        Z1=(R5-.001)*sin(i*2*pi/nH+0.001)
        X1=(R5-.001)*cos(i*2*pi/nH+.001)
        k2.append(f1.findAt((X1, L1, Z1),))
        #highlight(f1.findAt((X1, L1, Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta4', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta4, amplitude=UNSET)
#-----------------------------------
    k2=[]
    for i in range (nH):
        Z1=(R6)*sin(i*2*pi/nH+0.001)
        X1=(R6)*cos(i*2*pi/nH+.001)
        k2.append(f1.findAt((X1, L1-L4/2., Z1),))
        #highlight(f1.findAt((X1, L1-L4/2., Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta5', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta5, amplitude=UNSET)
#-----------------------------------
    k2=[]
    k2.append(f1.findAt(((R6-L9/2.)*cos(-ang5),L1-L4 ,(R6-L9/2.)*sin(-ang5)), ))
    k2.append(f1.findAt(((R6-L9/2.)*cos(ang5+ang4+0.01),L1-L4 ,(R6-L9/2.)*sin(ang5+ang4+0.01)), ))
    for i in range (nH):
        Z1=(R6-.01)*sin(i*2*pi/nH+pi/nH)
        X1=(R6-.01)*cos(i*2*pi/nH+pi/nH )
        k2.append(f1.findAt((X1, L1-L4, Z1),))
        #highlight(f1.findAt((X1, L1-L4, Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta6', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta6, amplitude=UNSET)
#-----------------------------------
    angle=atan((L6-L4)/float(R6-R1-L9-L10))
    k2=[]
    for i in range (nH):
        Z1=(R6-L9-1)*sin(i*2*pi/nH+0.001)
        X1=(R6-L9-1)*cos(i*2*pi/nH+.001)
        k2.append(f1.findAt((X1, L1-L4-tan(angle), Z1),))
        #highlight(f1.findAt((X1, L1-L4-tan(angle), Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta7', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta7, amplitude=UNSET)
#-----------------------------------
    k2=[]
    for i in range (nH):
        Z1=(R1+L10/2.)*sin(i*2*pi/nH+0.001)
        X1=(R1+L10/2.)*cos(i*2*pi/nH+.001)
        k2.append(f1.findAt((X1, L1-L6, Z1),))
        #highlight(f1.findAt((X1, L1-L6, Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='teta8', createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=teta8, amplitude=UNSET)
#-----------------------------------
    k2=[]
    for i in range (nH):
        Z1=(R5)*sin(i*2*pi/nH+0.001)
        X1=(R5)*cos(i*2*pi/nH+.001)
        k2.append(f1.findAt((X1, 0.001*L1, Z1),))
        highlight(f1.findAt((X1, 0.001*L1, Z1),))
#-----------------------------------
    region=(k2,)
    mdb.models[Model].TemperatureBC(name='T'+str(0), createStepName='Step-2', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=T[0], amplitude=UNSET)
#-----------------------------------
    L77=(L1-n*L8-2*L12)/(n-1)+L8
    for j in range (n):
        k3=[]
        k4=[]
        k5=[]
        k6=[]
        for i in range (nH):
            Z1=(R5-L11/2.)*sin(i*2*pi/nH+0.001)
            X1=(R5-L11/2.)*cos(i*2*pi/nH+.001)
            k3.append(f1.findAt((X1, L12+j*L77, Z1),))
            #highlight(f1.findAt((X1, L12+j*L77, Z1),))
            Z1=(R5-L11)*sin(i*2*pi/nH+0.1)
            X1=(R5-L11)*cos(i*2*pi/nH+.1)
            k4.append(f1.findAt((X1, L12+L8/2.+j*L77, Z1),))
            #highlight(f1.findAt((X1, L12+L8/2.+j*L77, Z1),))
            Z1=(R5-L11/2.)*sin(i*2*pi/nH+0.1)
            X1=(R5-L11/2.)*cos(i*2*pi/nH+.1)
            k5.append(f1.findAt((X1, L12+L8+j*L77, Z1),))
            #highlight(f1.findAt((X1, L12+L8+j*L77, Z1),))
            Z1=(R5)*sin(i*2*pi/nH+0.1)
            X1=(R5)*cos(i*2*pi/nH+.1)
            k6.append(f1.findAt((X1, L12+L8+j*L77+.001*L1, Z1),))
            #highlight(f1.findAt((X1, L12+L8+j*L77+.001, Z1),))      
        region1=(k3,)
        region2=(k4,)
        region3=(k5,)
        region4=(k6,)
        mdb.models[Model].TemperatureBC(name='T'+str(j*4+1), createStepName='Step-2', region=region1, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=T[4*j+1], amplitude=UNSET)
        mdb.models[Model].TemperatureBC(name='T'+str(j*4+2), createStepName='Step-2', region=region2, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=T[4*j+2], amplitude=UNSET)
        mdb.models[Model].TemperatureBC(name='T'+str(j*4+3), createStepName='Step-2', region=region3, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=T[4*j+3], amplitude=UNSET)
        mdb.models[Model].TemperatureBC(name='T'+str(j*4+4), createStepName='Step-2', region=region4, fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=T[4*j+4], amplitude=UNSET)
#-----------------------------------
    c1 = RAssmb.instances['Part-2-1'].cells
    mdb.models[Model].TemperatureBC(name='Tshaft', createStepName='Step-2', region=(c1, ), fixed=OFF, distributionType=UNIFORM, fieldName='', magnitude=Tshaft, amplitude=UNSET)
#-----------------------------------
    cells1=[]
    c1 = RAssmb.instances['Part-1-1'].cells
    for i in range (nH):
        Z1=(R1+0.01)*sin(i*2*pi/nH+0.1)
        X1=(R1+0.01)*cos(i*2*pi/nH+.1)
        cells1.append(c1.findAt((X1, L2/2., Z1),))
        highlight(c1.findAt((X1, L2/2., Z1),))
#-----------------------------------
    LL=((1-xL17)*L17)/2.+2.5
    LL1=((1-xL17)*L17)/2.+2.5+L17*xL17
    s1 = RAssmb.instances['Part-2-1'].faces
    v1 = RAssmb.instances['Part-2-1'].vertices
    k1=[]
    k1.append(v1.findAt((R7,LL,0),))
    k1.append(v1.findAt((R7*cos(2.*pi/nH*floor(1*nH/4.)),LL,R7*sin(2.*pi/nH*floor(1*nH/4.))),))
    k1.append(v1.findAt((R7*cos(2.*pi/nH*floor(2*nH/4.)),LL,R7*sin(2.*pi/nH*floor(2*nH/4.))),))
    k1.append(v1.findAt((R7*cos(2.*pi/nH*floor(3*nH/4.)),LL,R7*sin(2.*pi/nH*floor(3*nH/4.))),))
    k2=[]
    for i in range (nH):
        t2=s1.findAt(((R7+.01*R7)*cos(i*2*pi/nH+.001),LL1, (R7+.01*R7)*sin(i*2*pi/nH+.001)), )
        highlight(t2)

        k2.append(t2)
#bolt pre load---------------------------
    f1 = RAssmb.instances['Part-1-1'].faces
    f_preload=[]
    f_preload.append(f1.findAt(((1.01*(R6-L9)*cos(-ang5+0.01),(L1-L4) ,1.01*(R6-L9)*sin(-ang5+0.01)), )))
    for i in range (nH):
        f_preload.append(f1.findAt(((1.01*(R6-L9)*cos(i*2*pi/nH+.01),(L1-L4) ,1.01*(R6-L9)*sin(i*2*pi/nH+0.01)), )))
        if ((i+1)*2*pi/nH>ang7):
            break
    #
    f_preload_assem_surf = RAssmb.Surface(side1Faces=f_preload, name='f_preload_assem_surf')
    pressure_bolt=2*boltLoad/((R6**2-(R6-L9)**2)*(ang5+ang7)/2.)
    mdb.models['Model1'].Pressure(name='pressure-bolt', createStepName='Step-1', 
        region=f_preload_assem_surf, distributionType=UNIFORM, field='', magnitude=pressure_bolt, amplitude=UNSET)
    #
    allHoleFaces_assem=[]
    allHoleFaces_assem.append(f1.findAt((((R6-L9/2.),(L1-L4)*.99 ,(R6-L9/2.)*tan(-ang6)), )))
    allHoleFaces_assem.append(f1.findAt((((R6-L9/2.),(L1-L4)*.99 ,(R6-L9/2.)*tan(+ang6)), )))
    allHoleFaces_assem.append(f1.findAt((((R6-L9/2.)*cos(ang4)-D3/2.,(L1-L4)*.99 ,(R6-L9/2.)*sin(ang4)), )))
    face_hole3 = f11.findAt(((R6-L9/2.)*cos(ang4)-D3/2.,(L1-L4)*.99 ,(R6-L9/2.)*sin(ang4)), )
    face_hole4 = f11.findAt(((R6-L9/2.)*cos(ang4)+D3/2.,(L1-L4)*.99 ,(R6-L9/2.)*sin(ang4)), )
    if (face_hole3!=face_hole4):
        allHoleFaces_assem.append(f1.findAt((((R6-L9/2.)*cos(ang4)+D3/2.,(L1-L4)*.99 ,(R6-L9/2.)*sin(ang4)), )))
    #
    allHoleFaces_assem_surf = RAssmb.Surface(side1Faces=allHoleFaces_assem, name='allHoleFaces_assem_surf')
    traction_bolt=2*boltLoad/(2*pi*D3*L20)
    v1 = RAssmb.instances['Part-1-1'].vertices
    ver1=v1.findAt((R1, L5, 0),)
    ver2=v1.findAt((R1, L1-L6, 0),)
    mdb.models['Model1'].SurfaceTraction(name='traction_bolt', createStepName='Step-1', 
        region=allHoleFaces_assem_surf, magnitude=traction_bolt, directionVector=(ver1, ver2), 
        distributionType=UNIFORM, field='', localCsys=None)
    #-------------------------------------------------------
    X=RAssmb.DatumCsysByThreePoints(name='Datum csys-cylinderical', coordSysType=CYLINDRICAL, 
        origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, 1.0)).id
    datum = mdb.models[Model].rootAssembly.datums[X]
    mdb.models[Model].YsymmBC(name='BC_1', createStepName='Initial', 
        region=(k1, ), localCsys=datum)
    mdb.models[Model].YsymmBC(name='BC_2', createStepName='Initial', 
        region=(k2, ), localCsys=None)
    cell1 = RAssmb.instances['Part-1-1'].cells
    cell2 = RAssmb.instances['Part-2-1'].cells
    cell3 = RAssmb.instances['Part-3-1'].cells
    region = (cell1, )+(cell2, )+(cell3, )
    mdb.models[Model].RotationalBodyForce(name='Load-1', createStepName='Step-2', region=region,\
                    magnitude=w, centrifugal=ON, rotaryAcceleration=OFF, point1=(0.0, 0.0, 0.0), point2=(0.0, 8.0, 0.0))
#
    if no_step==3:
            mdb.models[Model].loads['Load-1'].setValuesInStep(stepName='Step-3', 
                 magnitude=w_Tol)   


   #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

def Mymesh (R1, R6, R7, R9, L1,L2, L4, L9, L21, D3, shaftMeshSize, slipMeshSize,shrinkMeshSize, RAssmb, Model):
    e1 = RAssmb.instances['Part-1-1'].edges
    pickedEdges = []
    for i in range (nH):
        pickedEdges.append(e1.findAt(((R1*cos(i*2*pi/nH+0.1),L5,R1*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt(((R1*cos(i*2*pi/nH+0.1),L1-L6,R1*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt(((R1*cos(i*2*pi/nH),L1/2.,R1*sin(i*2*pi/nH+0)))))
        highlight( pickedEdges[3*i])
        highlight( pickedEdges[3*i+1])
        highlight( pickedEdges[3*i+2])
    #--------------------------------------------------------------------
    RAssmb.seedEdgeBySize(edges=pickedEdges, size=shrinkMeshSize, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
    #--------------------------------------------------------------------
    e1 = RAssmb.instances['Part-2-1'].edges
    pickedEdges = []
    for i in range (nH):
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH+0.1),L18,R9*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt((((R9+L16)*cos(i*2*pi/nH+0.1),L18,(R9+L16)*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH),L1/2.,R9*sin(i*2*pi/nH+0)))))
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH+0.1),L18+L17,R9*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt((((R9+L16)*cos(i*2*pi/nH+0.1),L18+L17,(R9+L16)*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH),L1/2,R9*sin(i*2*pi/nH)))))   
    #--------------------------------------------------------------------    
    RAssmb.seedEdgeBySize(edges=pickedEdges, size=shrinkMeshSize, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
    #--------------------------------------------------------------------  
    e1 =  RAssmb.instances['Part-3-1'].edges
    pickedEdges = []
    for i in range (nH):
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH+0.1),L18,R9*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt((((R9+L16)*cos(i*2*pi/nH+0.1),L18,(R9+L16)*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH+0.1),L18+L17,R9*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt((((R9+L16)*cos(i*2*pi/nH+0.1),L18+L17,(R9+L16)*sin(i*2*pi/nH+0.1)))))
        pickedEdges.append(e1.findAt(((R9*cos(i*2*pi/nH),L1/2.,R9*sin(i*2*pi/nH)))))
        pickedEdges.append(e1.findAt((((R9+L16)*cos(i*2*pi/nH),L1/2.,(R9+L16)*sin(i*2*pi/nH)))))
    #----------------------------------------------------------------------
    RAssmb.seedEdgeBySize(edges=pickedEdges, size=shrinkMeshSize, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)
    cells1=[]
    c1 = RAssmb.instances['Part-1-1'].cells
    for i in range (nH):
        Z1=(R1+0.01)*sin(i*2*pi/nH+0.1)
        X1=(R1+0.01)*cos(i*2*pi/nH+.1)
        cells1.append(c1.findAt((X1, L2/2., Z1),))
#-----------------------------------
    c2 = RAssmb.instances['Part-2-1'].cells
    cells2=c2.findAt((R7+0.001*R7, L2/2., 0),)
    c3 = RAssmb.instances['Part-3-1'].cells
    cells3=c3.findAt((R9+0.001*R9, L2/2., 0),)
    RAssmb  = mdb.models[Model].rootAssembly
    partInstances =(RAssmb.instances['Part-1-1'], )
    RAssmb.seedPartInstance(regions=partInstances, size=slipMeshSize, deviationFactor=0.1, minSizeFactor=0.1)
    RAssmb.setMeshControls(regions=cells1, elemShape=TET, technique=FREE, sizeGrowthRate=1.0)
    RAssmb.generateMesh(regions=partInstances)
#-------------------------------------
    partInstances =( RAssmb.instances['Part-2-1'], )
    RAssmb.seedPartInstance(regions=partInstances, size=shaftMeshSize, deviationFactor=0.1, minSizeFactor=0.1)
    RAssmb.setMeshControls(regions=c2, elemShape=HEX, technique=STRUCTURED, sizeGrowthRate=1.0)
    RAssmb.generateMesh(regions=partInstances)
#---------------------------------------
    partInstances =( RAssmb.instances['Part-3-1'], )
    RAssmb.seedPartInstance(regions=partInstances, size=shrinkMeshSize, deviationFactor=0.1, minSizeFactor=0.1)
    RAssmb.setMeshControls(regions=c3, elemShape=HEX, technique=STRUCTURED, sizeGrowthRate=1.0)
    RAssmb.generateMesh(regions=partInstances)
#---------------------------------------    
    elemType1 = mesh.ElemType(elemCode=C3D8T, elemLibrary=STANDARD, secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=STANDARD)
    pickedRegions =((c3), (c2), (cells1), )
    RAssmb.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
    #mesh around 2 holes    
##    ang1=2.*pi/nH
##    ang2=asin(L21/(2.*R6-L9))
##    ang3=ang1
##    ang4=ang1+2.*ang2
##    ang5=atan(L9/(2.*R6-L9))
##    X1=(R6-L9/2.)*cos(ang3)
##    Z1=(R6-L9/2.)*sin(ang3)+(D3/2.)
##    Z2=(R6-L9/2.)*sin(ang3)-(D3/2.)
##    Z3=(R6-L9/2.)*sin(ang3)+(L9/2.)
##    Z4=(R6-L9/2.)*sin(ang3)-(L9/2.)
##    X2=(R6-L9/2.)*cos(ang4)
##    Z5=(R6-L9/2.)*sin(ang4)+(D3/2.)
##    Z6=(R6-L9/2.)*sin(ang4)-(D3/2.)
##    Z7=(R6-L9/2.)*sin(ang4)+(L9/2.)
##    Z8=(R6-L9/2.)*sin(ang4)-(L9/2.)
##    pickedEdges=[]
##    e1 = RAssmb.instances['Part-1-1'].edges
##    pickedEdges.append(e1.findAt((X1, L1-L4, Z1),))
##    pickedEdges.append(e1.findAt((X1, L1-L4, Z2),))
##    pickedEdges.append(e1.findAt((X1, L1-L4, Z3),))
##    pickedEdges.append(e1.findAt((X1, L1-L4, Z4),))
##    pickedEdges.append(e1.findAt((X2, L1-L4, Z5),))
##    pickedEdges.append(e1.findAt((X2, L1-L4, Z6),))
##    pickedEdges.append(e1.findAt((X2, L1-L4, Z7),))
##    pickedEdges.append(e1.findAt((X2, L1-L4, Z8),))
##    highlight(pickedEdges[2])
##    RAssmb.seedEdgeBySize(edges=pickedEdges, size=R8, deviationFactor=0.1, minSizeFactor=0.1, constraint=FINER)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
def myJob(jobName, Model, numCPU):
    mdb.Job(name=jobName, model=Model, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', multiprocessingMode=DEFAULT, numCpus=numCPU, numDomains=2*numCPU, 
        numGPUs=0)
    
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Interpolation of the S_N curve, This function interpolates S base on the N.   

def interpolateS_N(S_N_data, N):
    N1=log10(S_N_data[0][0])
    N2=log10(S_N_data[0][1])
    N3=log10(S_N_data[0][2])                                   
    N4=log10(S_N_data[0][3])
    S1=log10(S_N_data[1][0])
    S2=log10(S_N_data[1][1])
    S3=log10(S_N_data[1][2])                                    
    S4=log10(S_N_data[1][3])
    N=log10(N)
    if (N>=N1 and N<N2):
        S=S2+(S1-S2)*((N-N2)*1.0)/((N1-N2)*1.0)
    elif (N>=N2 and N<N3):
        S=S3+(S2-S3)*((N-N3)*1.0)/((N2-N3)*1.0)
    elif (N>=N3 and N<=N4):
        S=S4+(S3-S4)*((N-N4)*1.0)/((N3-N4)*1.0)
    elif (N<N1 or N>N4):
        S=1e+38
    return (10.**S)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Interpolation of xy data  

def intrplt(xy, x):
    for i in range (1,len(xy)):
        if  x>=xy[i-1][0] and x<xy[i][0]:
            y=xy[i-1][1]+((x-xy[i-1][0])/(xy[i][0]-xy[i-1][0]))*(xy[i][1]-xy[i-1][1])
            break
    return y





#******************************************************************************       

                  #M O D E L   C R E A T E
               
#******************************************************************************
p11, f11, d11, c11, e11, v11, p22, f22, d22, c22, e22, v22,  p33, f33, d33, c33, e33, v33 = \
myPart (nH, n, L1, L2, L3, L4, L5, L6, L7,L8 ,L9, L10, L11, L12, L13, L14,L15,\
L16, L17, L18, L19, L20, L21, R1, R2, R3, R4, R5,R6, R7, R8,R9, R10, D1, D2, D3, xL17, Model)

myMaterial(Eshaft, vshaft, Sysh, Den_sh, Exp_sh, k_sh,Eins, vins, Den_ins, Exp_ins, k_ins, Eslp,\
           vslp, Syslp, Den_slp, Exp_slp, k_slp, R7, L1, R9,nH, L2, Model, c11 ,c22, c33, p11, p22, p33)
RAssmb=myAssembly (p11, p22, p33, Model)

myStep(stptime_static, Model,2)

myInteraction (R9, L17,R6, L9, L1, L4, FC, prsscon, RAssmb, Model )



#******************************************************************************       
                  #M E S H  D E P E N D E N C Y  C H E C K
                  #Mesh Convergence for Regular Points

#******************************************************************************
if (meshConvergence==1):
    i=0
    wallSteps=10
    cd=0
    # two argument after w are w_Tol and no_step respectiveley
    myBoundaryCondition (nH, n, R1, R3, Rx, R5,  R6, R7, R9, R10,  L1, L2,  L4, L5, L6, L8, L9, L10, L11,L12,L17, L18, xL17,D3, teta1,\
                          teta2, teta3, teta4, teta5,teta6, teta7, teta8,Tshaft,  T ,w, 0, 2, Model, RAssmb)
    f = open('MeshConvergenceResults.txt', 'wb')
    f.write("***************************************************\n\n")
    f.write("Mesh Convergence Analyis\n\n")
    f.write("***************************************************\n")
    f.write("Run= %i\n" %1)
    f.write("Mesh Size of SlipRing= %5.4f\n" %slipMeshSize)
    f.write("Mesh Size of Shaft= %5.4f\n" %slipMeshSize)
    f.write("Mesh Size of insulation= %5.4f\n" %slipMeshSize)
    f.write("--------------------------------------------------\n")
    X1=R3*cos(-1*2*pi/nH)
    Z1=R3*sin(-1*2*pi/nH)
    X2=R3*cos(0*2*pi/nH)
    Z2=R3*sin(0*2*pi/nH)
    X3=R3*cos(1*2*pi/nH)
    Z3=R3*sin(1*2*pi/nH)
    SS=[]
    difrrenceSSLR=1e9
    while (True):
        i=i+1
        Mymesh (R1, R6, R7, R9, L1,L2, L4, L9, L21, D3, shaftMeshSize, slipMeshSize,shrinkMeshSize, RAssmb, Model)
        jobName2= 'MeshConvergenceRun'+str(i)
        directory2=directory+'/'+jobName2
        myJob(jobName2, Model, numCPU)
        mdb.saveAs(directoryname+'/'+jobName2)
        mdb.jobs[jobName2].submit(consistencyChecking=OFF)
        mdb.jobs[jobName2].waitForCompletion()
        modb = session.openOdb(name=directory2+'.odb')
        odb = session.odbs[directory2+'.odb']
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=-1)
        session.Path(name='Path-1', type=CIRCUMFERENTIAL, expression=((X1, 0, Z1), (X2,0, Z2), (X3, 0, Z3)),\
                     circleDefinition=POINT_ARC, numSegments=50, startAngle=0, endAngle=2*360./nH, radius=CIRCLE_RADIUS)
        pth = session.paths['Path-1']
        SS.append(session.XYDataFromPath(name='S_'+jobName2, path=pth, includeIntersections=False, 
            pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED, 
            labelType=TRUE_DISTANCE))
        #
        if (i!=1):
            x=abs(SS[i-1]-SS[i-2])
            s=0
            for h in range (len (x)):
                s=s+x[h][1]  
            #
            difrrenceSSLR=s/(len (x))
             #
            f.write("Run= %i\n" %(i))
            f.write("Mesh Size of SlipRing= %5.4f\n" %slipMeshSize)
            f.write("Mesh Size of Shaft= %5.4f\n" %shaftMeshSize)
            f.write("Mesh Size of insulation= %5.4f\n" %shrinkMeshSize)
            f.write("Difrence of Mises stress in a region of slip ring in two successive steps= %5.4f\n" % difrrenceSSLR)  
            f.write("-------------------------------------------------\n")
            shaftMeshSize=shaftMeshSize*mesh_ref_fac
            slipMeshSize=slipMeshSize*mesh_ref_fac
            shrinkMeshSize=shrinkMeshSize*mesh_ref_fac
            #------------------------- 
            if  (abs(difrrenceSSLR) <1) :
                cd=1
        else:
            shaftMeshSize=shaftMeshSize*mesh_ref_fac
            slipMeshSize=slipMeshSize*mesh_ref_fac
            shrinkMeshSize=shrinkMeshSize*mesh_ref_fac
        #-------------------------
        if (cd==1 or i==wallSteps):
            break
        #--------------------------                                        
    f.close()                                             


##
###******************************************************************************       
##
##                  #S I M P L E  A N A L Y S I S
##              #Extrapolation of Stress in singular points
##                 
###******************************************************************************
if meshConvergence==1:
    modb = session.openOdb(name=directory2+'.odb')
    odb = session.odbs[directory2+'.odb']
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)
else:
    myBoundaryCondition (nH, n, R1, R3, Rx, R5,  R6, R7, R9, R10,  L1, L2,  L4, L5, L6, L8, L9, L10, L11,L12,L17, L18, xL17,D3, teta1,\
                             teta2, teta3, teta4, teta5,teta6, teta7, teta8,Tshaft,  T ,w,0,2,  Model, RAssmb)
    Mymesh (R1, R6, R7, R9, L1,L2, L4, L9, L21, D3, shaftMeshSize, slipMeshSize,shrinkMeshSize, RAssmb, Model)
    jobName0= 'SimpleAnalysis'
    directory0=directory+'/'+jobName0
    myJob(jobName0, Model, numCPU)
    mdb.saveAs(directoryname+'/'+jobName0)
    mdb.jobs[jobName0].submit(consistencyChecking=OFF)
    mdb.jobs[jobName0].waitForCompletion()
    modb = session.openOdb(name=directory0+'.odb')
    odb = session.odbs[directory0+'.odb']
    session.viewports['Viewport: 1'].setValues(displayedObject=odb)


quadDistance=[14, 19, 24]
linDistance =[11, 16]
tet=atan(L5/(R3-Rx))
Lx= R4-D1/2.-R2
Lxx=Lx-D2/(2.*cos(tet))
LsingPnt=L3+Lxx/tan(tet)
setSingHole=modb.rootAssembly.instances['PART-1-1'].nodeSets['SINGHOLEINTERSECT'].nodes


f0 = open('SingularPointsAssessment.txt', 'wb')
f0.write('Evaluating stress values in singular points\n')
f0.write('-------------------------------------------\n\n')
f0.write('Singular points located on two hole intersection\n\n')
f0.write('Node Label No.  Mises    Quad. Interp.  Lin. Interp.\n\n')
#
instancesName=['PART-1-1', 'PART-2-1', 'PART-3-1']
leaf = dgo.LeafFromPartInstance(partInstanceName=('PART-1-1', 'PART-2-1', 'PART-3-1' ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)
leaf = dgo.LeafFromPartInstance(partInstanceName=(instancesName[0], ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.add(leaf=leaf)
for i in range (len(setSingHole)):
    for j in range (i+1,len(setSingHole)):
        if (setSingHole[i].coordinates[0]-setSingHole[j].coordinates[0])==0 and \
           (setSingHole[i].coordinates[2]-setSingHole[j].coordinates[2])==0 :
            label1=setSingHole[i].label
            label2=setSingHole[j].label
            #
            if setSingHole[i].coordinates[1]>setSingHole[j].coordinates[1]:
                lb=[label2, label1]
            else:
                lb=[label1, label2]
            #
            session.Path(name='Path-3', type=NODE_LIST, expression=(('PART-1-1', lb), ))
            pth = session.paths['Path-3']
            session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=-1)
            xy = session.XYDataFromPath(name='xyhole'+str(i), path=pth, includeIntersections=True, 
                      pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
            # Quadratic Interpolation
            d1=intrplt(xy, LsingPnt-quadDistance[2])
            d2=intrplt(xy, LsingPnt-quadDistance[1])
            d3=intrplt(xy, LsingPnt-quadDistance[0])
            xx=[LsingPnt-quadDistance[2], LsingPnt-quadDistance[1], LsingPnt-quadDistance[0]]
            z1 = numpy.polyfit(xx, [d1, d2, d3], 2)
            pp1 = numpy.poly1d(z1)
            stressSingPointQuad=pp1(LsingPnt)
            stressSingPointQuad
            # Linear Interpolation
            d1=intrplt(xy, LsingPnt-linDistance[1])
            d2=intrplt(xy, LsingPnt-linDistance[0])
            xx=[LsingPnt-linDistance[1], LsingPnt-linDistance[0]]
            z1 = numpy.polyfit(xx, [d1, d2], 1)
            pp1 = numpy.poly1d(z1)
            stressSingPointLin=pp1(LsingPnt)
            stressSingPointLin
            f0.write("%6d        %5.3f       %5.3f       %5.3f\n " %(lb[0], xy[len(xy)-1][1], stressSingPointQuad, stressSingPointLin))
            break

#
f0.write('-------------------------------------------\n\n')
f0.write('Singular points located on two hole intersection\n\n')
f0.write('Node Label No.  Mises    Quad. Interp.  Lin. Interp.\n\n')
setSingContact=modb.rootAssembly.instances['PART-1-1'].nodeSets['SINGCONTACTSLR'].nodes

leaf = dgo.LeafFromPartInstance(partInstanceName=('PART-1-1', 'PART-2-1', 'PART-3-1' ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)
leaf = dgo.LeafFromPartInstance(partInstanceName=(instancesName[0], ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.add(leaf=leaf)
for i in range (len(setSingContact)):
    for j in range (i+1,len(setSingContact)):
        if (setSingContact[i].coordinates[0]-setSingContact[j].coordinates[0])==0 and \
           (setSingContact[i].coordinates[2]-setSingContact[j].coordinates[2])==0 :
            label1=setSingContact[i].label
            label2=setSingContact[j].label
            #
            if setSingContact[i].coordinates[1]>setSingContact[j].coordinates[1]:
                lb=[label2, label1]
            else:
                lb=[label1, label2]
            #
            session.Path(name='Path-4', type=NODE_LIST, expression=(('PART-1-1', lb), ))
            pth = session.paths['Path-4']
            session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=-1)
            xy = session.XYDataFromPath(name='xyconatctSLR'+str(i), path=pth, includeIntersections=True, 
                      pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
            # Quadratic Interpolation
            d1=intrplt(xy, quadDistance[0])
            d2=intrplt(xy, quadDistance[1])
            d3=intrplt(xy, quadDistance[2])
            xx=[quadDistance[0], quadDistance[1], quadDistance[2]]
            z1 = numpy.polyfit(xx, [d1, d2, d3], 2)
            pp1 = numpy.poly1d(z1)
            stressSingPointQuad=pp1(0)
            stressSingPointQuad
            # Linear Interpolation
            d1=intrplt(xy, linDistance[0])
            d2=intrplt(xy, linDistance[1])
            xx=[linDistance[0], linDistance[1]]
            z1 = numpy.polyfit(xx, [d1, d2], 1)
            pp1 = numpy.poly1d(z1)
            stressSingPointLin=pp1(0)
            stressSingPointLin
            f0.write("%6d        %5.3f       %5.3f       %5.3f\n " %(lb[0], xy[0][1], stressSingPointQuad, stressSingPointLin))
            # Quadratic Interpolation
            d1=intrplt(xy, L1-L5-L6-quadDistance[2])
            d2=intrplt(xy, L1-L5-L6-quadDistance[1])
            d3=intrplt(xy, L1-L5-L6-quadDistance[0])
            xx=[L1-L5-L6-quadDistance[2], L1-L5-L6-quadDistance[1], L1-L5-L6-quadDistance[0]]
            z1 = numpy.polyfit(xx, [d1, d2, d3], 2)
            pp1 = numpy.poly1d(z1)
            stressSingPointQuad=pp1(xy[len(xy)-1][0])
            stressSingPointQuad
            # Linear Interpolation
            d1=intrplt(xy, L1-L5-L6-linDistance[1])
            d2=intrplt(xy, L1-L5-L6-linDistance[0])
            xx=[L1-L5-L6-linDistance[1], L1-L5-L6-linDistance[0]]
            z1 = numpy.polyfit(xx, [d1, d2], 1)
            pp1 = numpy.poly1d(z1)
            stressSingPointLin=pp1(xy[len(xy)-1][0])
            stressSingPointLin
            f0.write("%6d        %5.3f       %5.3f       %5.3f\n " %(lb[1], xy[len(xy)-1][1], stressSingPointQuad, stressSingPointLin))
            break

#
f0.write('-------------------------------------------\n\n')
f0.write('Singular points located on the begining and end of the contact in SInsulation\n')
f0.write('Node Label No.  Mises    Quad. Interp.  Lin. Interp.\n\n')
setSingContactINS=modb.rootAssembly.instances['PART-3-1'].nodeSets['SINGCONTACTINS'].nodes

leaf = dgo.LeafFromPartInstance(partInstanceName=('PART-1-1', 'PART-2-1', 'PART-3-1' ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.remove(leaf=leaf)
leaf = dgo.LeafFromPartInstance(partInstanceName=(instancesName[2], ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.add(leaf=leaf)
for i in range (len(setSingContactINS)):
    for j in range (i+1,len(setSingContact)):
        if (setSingContactINS[i].coordinates[0]-setSingContactINS[j].coordinates[0])==0 and \
           (setSingContactINS[i].coordinates[2]-setSingContactINS[j].coordinates[2])==0 :
            label1=setSingContactINS[i].label
            label2=setSingContactINS[j].label
            #
            if setSingContactINS[i].coordinates[1]>setSingContactINS[j].coordinates[1]:
                lb=[label2, label1]
            else:
                lb=[label1, label2]
            #
            session.Path(name='Path-5', type=NODE_LIST, expression=(('PART-3-1', lb), ))
            pth = session.paths['Path-5']
            session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=-1)
            xy = session.XYDataFromPath(name='xyconatctINS'+str(i), path=pth, includeIntersections=True, 
                      pathStyle=PATH_POINTS, numIntervals=10, shape=UNDEFORMED, labelType=TRUE_DISTANCE)
            # Quadratic Interpolation
            d1=intrplt(xy, quadDistance[0])
            d2=intrplt(xy, quadDistance[1])
            d3=intrplt(xy, quadDistance[2])
            xx=[quadDistance[0], quadDistance[1], quadDistance[2]]
            z1 = numpy.polyfit(xx, [d1, d2, d3], 2)
            pp1 = numpy.poly1d(z1)
            stressSingPointQuad=pp1(0)
            stressSingPointQuad
            # Linear Interpolation
            d1=intrplt(xy, linDistance[0])
            d2=intrplt(xy, linDistance[1])
            xx=[linDistance[0], linDistance[1]]
            z1 = numpy.polyfit(xx, [d1, d2], 1)
            pp1 = numpy.poly1d(z1)
            stressSingPointLin=pp1(0)
            stressSingPointLin
            f0.write("%6d        %5.3f       %5.3f       %5.3f\n " %(lb[0], xy[0][1], stressSingPointQuad, stressSingPointLin))
            # Quadratic Interpolation
            d1=intrplt(xy, L1-L5-L6-quadDistance[2])
            d2=intrplt(xy, L1-L5-L6-quadDistance[1])
            d3=intrplt(xy, L1-L5-L6-quadDistance[0])
            xx=[L1-L5-L6-quadDistance[2], L1-L5-L6-quadDistance[1], L1-L5-L6-quadDistance[0]]
            z1 = numpy.polyfit(xx, [d1, d2, d3], 2)
            pp1 = numpy.poly1d(z1)
            stressSingPointQuad=pp1(xy[len(xy)-1][0])
            stressSingPointQuad
            # Linear Interpolation
            d1=intrplt(xy, L1-L5-L6-linDistance[1])
            d2=intrplt(xy, L1-L5-L6-linDistance[0])
            xx=[L1-L5-L6-linDistance[1], L1-L5-L6-linDistance[0]]
            z1 = numpy.polyfit(xx, [d1, d2], 1)
            pp1 = numpy.poly1d(z1)
            stressSingPointLin=pp1(xy[len(xy)-1][0])
            stressSingPointLin
            f0.write("%6d        %5.3f       %5.3f       %5.3f\n " %(lb[1], xy[len(xy)-1][1], stressSingPointQuad, stressSingPointLin))
            break

#
f0.close()



#*******************************************************************************************
                      #S E P E R A T I O N   I N   O V E R  S P E E D
# Seperation test at over speed. This section calculate percent of seperated element
# in the conact region of the slipRing and insulation at high speed which called w_overSpeed
                                        
#********************************************************************************************                                
if (seperationCheck==1):
    myStep(stptime_static, Model,2)
    # two argument after w are w_Tol and no_step respectiveley
    myBoundaryCondition (nH, n, R1, R3, Rx, R5,  R6, R7, R9, R10,  L1, L2,  L4, L5, L6, L8, L9, L10, L11,L12,L17, L18, xL17,D3, teta1,\
                          teta2, teta3, teta4, teta5,teta6, teta7, teta8,Tshaft,  T ,w_overSpeed, 0, 2, Model, RAssmb)
    Mymesh (R1, R6, R7, R9, L1,L2, L4, L9, L21, D3, shaftMeshSize, slipMeshSize,shrinkMeshSize, RAssmb, Model)
    jobName3= 'SeperationCheck'
    directory3=directory+'/'+jobName3
    myJob(jobName3, Model, numCPU)
    mdb.saveAs(directoryname+'/'+jobName3)
    mdb.jobs[jobName3].submit(consistencyChecking=OFF)
    mdb.jobs[jobName3].waitForCompletion()
    modb = session.openOdb(name=directory3+'.odb')
    # SET-contact should be defined in PART modules of the Abaqus as a "set" wich
    # include inner face of the slipRing and has contact with insulation
    setContact=modb.rootAssembly.instances['PART-1-1'].nodeSets['SET_CONTACT']        
    # CSTATUS values of 1, 2 and -3.8e+38 stand for closed(sticking), 
    # closed(slipping) and open states respectivelye  in cantact for each node.
    STICK=0
    OPEN=0
    SLIP=0
    CSTATUS=[]
    field=modb.steps['Step-1'].frames[-1].fieldOutputs['CSTATUS'].getSubset(region=setContact)
    for i in range (len(setContact.nodes)):
        CSTATUS.append(field.values[i].data)
#                                        
    for i in range (len(setContact.nodes)):
        if (CSTATUS[i]==2):
            STICK=STICK+1
        elif (CSTATUS[i]==1):
            SLIP=SLIP+1
        else:
            OPEN=OPEN+1
#                                       
    sepNodes=OPEN/(len(setContact.nodes)*1.)*100
    stickNodes= STICK/(len(setContact.nodes)*1.)*100
    slipNodes=SLIP/(len(setContact.nodes)*1.)*100
    f3 = open('SeperationCheckResults.txt', 'wb')
    f3.write("*****************************************************\n")
    f3.write("Seperation of the Slip Ring and Shadt in Over Speed\n")
    f3.write("This analysis is acomplised based of the elements of the interior surface of the slip ring.\n")
    f3.write("*****************************************************\n")
    f3.write("Percent of the Seperated Elements= %5.4f\n" %sepNodes)
    f3.write("Percent of the Sticked Elements= %5.4f\n"   %stickNodes)
    f3.write("Percent of the Slipped Elements= %5.4f\n"   %slipNodes)
    f3.write("-----------------------------------------------------\n")                                     
#                                             
    f3.close()
                                         

#************************************************************************************************
                            #L O W   C Y C L E   F A T I G U E
#Low cycle fatigue whcih include two load cases: first, last frame of step 1 after shrink process
# and second, last frame of step 2 where rotation forces and thermal loads are applied to.                        
                                    
#************************************************************************************************

if (lowCycleAnalysis==1):                                   
    myStep(stptime_static, Model,2)
    # two argument after w are w_Tol and no_step respectiveley
    myBoundaryCondition (nH, n, R1, R3, Rx, R5,  R6, R7, R9, R10,  L1, L2,  L4, L5, L6, L8, L9, L10, L11,L12,L17,   L18, xL17,D3, teta1,\
                          teta2, teta3, teta4, teta5,teta6, teta7, teta8,Tshaft,  T ,w, 0, 2, Model, RAssmb)
    Mymesh (R1, R6, R7, R9, L1,L2, L4, L9, L21, D3, shaftMeshSize, slipMeshSize,shrinkMeshSize, RAssmb, Model)
    jobName4= 'LowCycleFatigueAnalysis'
    directory4=directory+'/'+jobName4
    myJob(jobName4, Model, numCPU)
    mdb.saveAs(directoryname+'/'+jobName4)
    mdb.jobs[jobName4].submit(consistencyChecking=OFF)
    mdb.jobs[jobName4].waitForCompletion()
    modb = session.openOdb(name=directory4+'.odb')
    #
    S_LowCycleSLR=interpolateS_N(S_N_SLR, N_LowCycle)
    S_LowCycleSHFT=interpolateS_N(S_N_SHFT, N_LowCycle)
    S_LowCycleINS=interpolateS_N(S_N_INS, N_LowCycle)
    n1=modb.rootAssembly.instances['PART-1-1']
    n2=modb.rootAssembly.instances['PART-2-1']
    n3=modb.rootAssembly.instances['PART-3-1']
    #SlipRing Analysis
    misesStep1SLR=[]
    misesStep2SLR=[]
    misesAmpSLR=[]
    notPassedElSLR=[]
    field1=modb.steps['Step-1'].frames[-1].fieldOutputs['S'].getSubset(region=n1)
    field2=modb.steps['Step-2'].frames[-1].fieldOutputs['S'].getSubset(region=n1)
    misesStep1SLR=([g.mises for g in field1.values])
    misesStep2SLR=([g.mises for g in field2.values])
    misesSLR_a=[]
    misesSLR_m=[]
    for i in range(len(n1.elements)):
        misesSLR_a.append(abs(misesStep1SLR[i]-misesStep2SLR[i])/2.)
        misesSLR_m.append(abs(misesStep1SLR[i]+misesStep2SLR[i])/2.)
    #----------------------------------------
    sf_SLR=[]
    for i in range(len(n1.elements)):
        aa=(misesSLR_a[i]/S_LowCycleSLR)+(misesSLR_m[i]/S_U_SLR)
        sf_SLR.append((1./aa))
        if  (aa>1):
            notPassedElSLR.append(n1.elements[i].label)
    #Shaft Analysis
    misesStep1SHFT=[]
    misesStep2SHFT=[]
    misesAmpSHFT=[]
    notPassedElSHFT=[]
    field1=modb.steps['Step-1'].frames[-1].fieldOutputs['S'].getSubset(region=n2)
    field2=modb.steps['Step-2'].frames[-1].fieldOutputs['S'].getSubset(region=n2)
    misesStep1SHFT=([g.mises for g in field1.values])
    misesStep2SHFT=([g.mises for g in field2.values])
    misesSHFT_a=[]
    misesSHFT_m=[]
    for i in range(len(misesStep1SHFT)):
         misesSHFT_a.append(abs(misesStep1SHFT[i]-misesStep2SHFT[i])/2.)
         misesSHFT_m.append(abs(misesStep1SHFT[i]+misesStep2SHFT[i])/2.)
    #----------------------------------------
    sf_SHFT=[]     
    for i in range(len(misesStep1SHFT)):
        aa=(misesSHFT_a[i]/S_LowCycleSHFT)+(misesSHFT_m[i]/S_U_SHFT)
        sf_SHFT.append((1./aa))
        if  (aa >1):
            notPassedElSHFT.append(i)
    #Insulation Analysis
    misesStep1INS=[]
    misesStep2INS=[]
    misesAmpINS=[]
    notPassedElINS=[]
    field1=modb.steps['Step-1'].frames[-1].fieldOutputs['S'].getSubset(region=n3)
    field2=modb.steps['Step-2'].frames[-1].fieldOutputs['S'].getSubset(region=n3)
    misesStep1INS=([g.mises for g in field1.values])
    misesStep2INS=([g.mises for g in field2.values])
    misesINS_a=[]
    misesINS_m=[]
    for i in range(len(misesStep1INS)):
        misesINS_a.append(abs(misesStep1INS[i]-misesStep2INS[i])/2.)
        misesINS_m.append(abs(misesStep1INS[i]+misesStep2INS[i])/2.)
    #----------------------------------------
    sf_INS=[]                              
    for i in range(len(misesStep1INS)):
        aa=(misesINS_a[i]/S_LowCycleINS)+(misesINS_m[i]/S_U_INS)
        sf_INS.append((1./aa))
        if  (aa>1):
            notPassedElINS.append(n3.elements[i/8].label)
    #--------------------------
    notPassSLPpercent = 100.*len(notPassedElSLR) /len( misesSLR_a)
    notPassSHFTpercent= 100.*len(notPassedElSHFT)/len( misesSHFT_a)
    notPassINSpercent = 100.*len(notPassedElINS)/len( misesINS_a)
    #Low cycle fatigue whcih include two load cases: first, last frame of step 1 after shrink process
    # and second, last frame of step 2 where rotation forces and thermal loads are applied to. 
    f4 = open('LowCycleFatigue.txt', 'wb')
    f4.write("*****************************************************\n")
    f4.write("Low cycle fatigue whcih include two load cases: \n")
    f4.write("First: last frame of step 1 after shrink process\n")
    f4.write("And second: last frame of step 2 where rotation forces\n")
    f4.write("And thermal loads are applied to.\n")
    f4.write("*****************************************************\n")
    f4.write(" No. Cycle= %5.4f \n" %  (N_LowCycle))
    f4.write(" S of Slip Ring in the given cycle number= %5.4f \n" %  (S_LowCycleSLR ))
    f4.write(" S of Shaft in the given cycle numbe= %5.4f \n"     %  (S_LowCycleSHFT))
    f4.write(" S of Insulation in the given cycle number= %5.4f \n"%  (S_LowCycleINS))
    f4.write(" Percent of SlipRing elements which not passed high cycle fatigue analysis= %5.4f\n" %  notPassSLPpercent)
    f4.write(" Minimum safty factor in SlipRing= %5.4f\n" %  min(sf_SLR))
    f4.write(" Percent of Shaft elements which not passed high cycle fatigue analysis= %5.4f\n" %  notPassSHFTpercent)
    f4.write(" Minimum safty factor in Sahft= %5.4f\n" %  min(sf_SHFT))
    f4.write(" Percent of Insulation elements which not passed high cycle fatigue analysis= %5.4f\n" %  notPassINSpercent)
    f4.write(" Minimum safty factorin Insulation= %5.4f\n" %  min(sf_INS))
    f4.write("-----------------------------------------------------\n")
    f4.write("SlipRing\n\n")
    f4.write("El.    SF\n")                       
    for i in range (len(misesSLR_a)):
           f4.write("%i       %5.4f\n" % (n1.elements[i].label,  sf_SLR[i]))              
    #
    f4.write("-----------------------------------------------------\n")
    f4.write("Shaft\n\n")
    f4.write("El.   Int. Point    SF\n")                       
    for i in range (len(misesSHFT_a)):
           f4.write("%i      %i          %5.4f\n" % (n2.elements[i/8].label, i%8 ,sf_SHFT[i]))
    #
    f4.write("-----------------------------------------------------\n")
    f4.write("SlipRing\n\n")
    f4.write("El.   Int. Point    SF\n")                         
    for i in range (len(misesINS_a)):
           f4.write("%i      %i          %5.4f\n" % (n3.elements[i/8].label, i%8 ,sf_INS[i]))
   #                                     
    f4.close()     




#************************************************************************************************
                            #HIGH   C Y C L E   F A T I G U E
#High cycle fatigue which assume fluctuation in rotational speed as load cases.
                                    
#************************************************************************************************
if (highCycleAnalysis==1):
    S_HighCycleSLR=interpolateS_N(S_N_SLR, N_HighCycle)
    S_HighCycleSHFT=interpolateS_N(S_N_SHFT, N_HighCycle)
    S_HighCycleINS=interpolateS_N(S_N_INS, N_HighCycle)    
    #
    wFluct=w_tol
    highSpeed=w+wFluct
    lowSpeed=w-wFluct
    #
    myStep(stptime_static, Model,3)
    # two argument after w are w_Tol and no_step respectiveley
    myBoundaryCondition (nH, n, R1, R3, Rx, R5,  R6, R7, R9, R10,  L1, L2,  L4, L5, L6, L8, L9, L10, L11,L12,L17, L18, xL17,D3, teta1,\
                          teta2, teta3, teta4, teta5,teta6, teta7, teta8,Tshaft,  T ,highSpeed, lowSpeed, 3, Model, RAssmb)
    Mymesh (R1, R6, R7, R9, L1,L2, L4, L9, L21, D3, shaftMeshSize, slipMeshSize,shrinkMeshSize, RAssmb, Model)
    jobName6= 'HighCycleFatigueAnalysis'
    directory6=directory+'/'+jobName6
    myJob(jobName6, Model, numCPU)
    mdb.saveAs(directoryname+'/'+jobName6)
    mdb.jobs[jobName6].submit(consistencyChecking=OFF)
    mdb.jobs[jobName6].waitForCompletion()
    modb = session.openOdb(name=directory6+'.odb')                                     
    #
    n1=modb.rootAssembly.instances['PART-1-1']
    n2=modb.rootAssembly.instances['PART-2-1']
    n3=modb.rootAssembly.instances['PART-3-1']
    #SlipRing Analysis
    misesHighSpeedSLR=[]
    misesLowSpeedSLR=[]
    misesAmpSLR=[]
    notPassedElSLR=[]
    field1=modb.steps['Step-2'].frames[-1].fieldOutputs['S'].getSubset(region=n1)
    field2=modb.steps['Step-3'].frames[-1].fieldOutputs['S'].getSubset(region=n1)
    misesHighSpeedSLR=([g.mises for g in field1.values])
    misesLowSpeedSLR= ([g.mises for g in field2.values])
    misesSLR_a=[]
    misesSLR_m=[]
    for i in range(len(n1.elements)):
        misesSLR_a.append(abs(misesHighSpeedSLR[i]-misesLowSpeedSLR[i])/2.)
        misesSLR_m.append(abs(misesHighSpeedSLR[i]+misesLowSpeedSLR[i])/2.)
    #----------------------------------------                                 
    sf_SLR=[]
    for i in range(len(misesHighSpeedSLR)):                         
        aa=(misesSLR_a[i]/S_HighCycleSLR)+(misesSLR_m[i]/S_U_SLR)
        sf_SLR.append((1./aa))
        if  (aa>1):
             notPassedElSLR.append(n1.elements[i].label)
    #Shaft Analysis
    misesHighSpeedSHFT=[]
    misesLowSpeedSHFT=[]
    misesAmpSHFT=[]
    notPassedElSHFT=[]
    field1=modb.steps['Step-2'].frames[-1].fieldOutputs['S'].getSubset(region=n2)
    field2=modb.steps['Step-3'].frames[-1].fieldOutputs['S'].getSubset(region=n2)
    misesHighSpeedSHFT=([g.mises for g in field1.values])
    misesLowSpeedSHFT =([g.mises for g in field2.values])
    misesSHFT_a=[]
    misesSHFT_m=[]
    for i in range(len(misesLowSpeedSHFT)):
        misesSHFT_a.append(abs(misesHighSpeedSHFT[i]-misesLowSpeedSHFT[i])/2.)
        misesSHFT_m.append(abs(misesHighSpeedSHFT[i]+misesLowSpeedSHFT[i])/2.)
    #----------------------------------------                                      
    sf_SHFT=[]
    for i in range(len(misesLowSpeedSHFT)):                         
        aa=(misesSHFT_a[i]/S_HighCycleSHFT)+(misesSHFT_m[i]/S_U_SHFT)
        sf_SHFT.append((1./aa))
        if  (aa>1):
            notPassedElSHFT.append(n2.elements[i/8].label)
    #Insulation Analysis
    misesHighSpeedINS=[]
    misesLowSpeedINS=[]
    misesAmpINS=[]
    notPassedElINS=[]
    field1=modb.steps['Step-2'].frames[-1].fieldOutputs['S'].getSubset(region=n3)
    field2=modb.steps['Step-3'].frames[-1].fieldOutputs['S'].getSubset(region=n3)
    misesHighSpeedINS=([g.mises for g in field1.values])
    misesLowSpeedINS=([g.mises for g in field2.values])
    misesINS_a=[]
    misesINS_m=[]
    for i in range(len(misesLowSpeedINS)):
        misesINS_a.append(abs(misesLowSpeedINS[i]-misesHighSpeedINS[i])/2)
        misesINS_m.append(abs(misesLowSpeedINS[i]+misesHighSpeedINS[i])/2)
    #----------------------------------------                                 
    sf_INS=[]
    for i in range(len(misesLowSpeedINS)):
        aa=(misesINS_a[i]/S_HighCycleINS)+(misesINS_m[i]/S_U_INS)
        sf_INS.append((1./aa))
        if  (aa >1):
            notPassedElINS.append(n3.elements[i/8].label)
    #
    notPassSLPpercent = 100.*len(notPassedElSLR) /len( misesSLR_a)
    notPassSHFTpercent= 100.*len(notPassedElSHFT)/len( misesSHFT_a)
    notPassINSpercent = 100.*len(notPassedElINS)/len( misesINS_a)
    #High cycle fatigue which assume fluctuation in rotational speed as load cases. 
    f5 = open('HighCycleFatigue.txt', 'wb')
    f5.write("*****************************************************\n")
    f5.write("High cycle fatigue whcih include two load cases: \n")
    f5.write("First: rotational speed W+W_flactuation\n")
    f5.write("And second: rotational speed W-W_flactuation\n")
    f5.write("*****************************************************\n")
    f5.write(" No. Cycle= %5.4f= \n" %  (N_HighCycle))
    f5.write(" S of Slip Ring in the given cycle number %5.4f=\n" %  (S_HighCycleSLR ))
    f5.write(" S of Shaft in the given cycle number%5.4f= \n"     %  (S_HighCycleSHFT))
    f5.write(" S of Insulation in the given cycle number %5.4f=\n"%  (S_HighCycleINS))
    f5.write(" Percent of SlipRing elements which not passed high cycle fatigue analysis= %5.4f\n" %  notPassSLPpercent)
    f5.write(" Minimum safty factor in SlipRing= %5.4f\n" %  min(sf_SLR))
    f5.write(" Percent of Shaft elements which not passed high cycle fatigue analysis= %5.4f\n" %  notPassSHFTpercent)
    f5.write(" Minimum safty factor in Sahft= %5.4f\n" %  min(sf_SHFT))
    f5.write(" Percent of Insulation elements which not passed high cycle fatigue analysis= %5.4f\n" %  notPassINSpercent)
    f5.write(" Minimum safty factorin Insulation= %5.4f\n" %  min(sf_INS))
    f5.write("-----------------------------------------------------\n")
    f5.write("SlipRing\n\n")
    f5.write("El.    SF\n")                       
    for i in range (len(misesSLR_a)):
           f5.write("%i       %5.4f\n" % (n1.elements[i].label,  sf_SLR[i]))              
    #
    f5.write("-----------------------------------------------------\n")
    f5.write("Shaft\n\n")
    f5.write("El.   Int. Point    SF\n")                       
    for i in range (len(misesSHFT_a)):
           f5.write("%i      %i          %5.4f\n" % (n2.elements[i/8].label, i%8 ,sf_SHFT[i]))
    #
    f5.write("-----------------------------------------------------\n")
    f5.write("SlipRing\n\n")
    f5.write("El.   Int. Point    SF\n")                         
    for i in range (len(misesINS_a)):
           f5.write("%i      %i          %5.4f\n" % (n3.elements[i/8].label, i%8 ,sf_INS[i]))
   #                                     
    f5.close()     

     
leaf = dgo.LeafFromPartInstance(partInstanceName=('PART-1-1', 'PART-2-1', 'PART-3-1' ))
session.viewports['Viewport: 1'].odbDisplay.displayGroup.add(leaf=leaf)

     
