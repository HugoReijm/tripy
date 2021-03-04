import math
import time
import numpy as np
import numpy.linalg as npla
from point import point
from edge import edge
from triangle import triangle
import random
import matplotlib.pyplot as plt
        
def quickhull(P,errtol=1e-12):
    if len(P)<=1:
        print("Error: Unable to make a convex hull due to lack of sufficient points")
        return []
    elif len(P)==2:
        return [edge(P[0],P[1])]
    else:
        e0=edge(P[0],P[1])
        colinearbool=True
        for i in range(2,len(P)):
            if abs(e0.cross(edge(P[0],P[i])))>abs(errtol):
                colinearbool=False
                break
        if colinearbool:
            print("Error: point set is colinear, unable to make convex hull")
            return []
       
    anchor1=P[0]
    for i in range(1,len(P)):
        if P[i].x<anchor1.x:
            anchor1=P[i]
        elif P[i].x==anchor1.x and P[i].y<anchor1.y:
            anchor1=P[i]
    
    anchor2=P[0]
    for i in range(1,len(P)):
        if P[i].x>anchor2.x:
            anchor2=P[i]
        elif P[i].x==anchor2.x and P[i].y>anchor2.y:
            anchor2=P[i]
    
    anchorU=anchor1
    heightU=0
    anchorL=anchor2
    heightL=0
    for p in P:
        res=-(anchor2.y-anchor1.y)*p.x+(anchor2.x-anchor1.x)*p.y-anchor2.x*anchor1.y+anchor2.y*anchor1.x
        if (anchor2.x-anchor1.x)*res>abs(errtol) and res>heightU+abs(errtol):
            anchorU=p
            heightU=res
        elif (anchor2.x-anchor1.x)*res<-abs(errtol) and res+abs(errtol)<heightL:
            anchorL=p
            heightL=res
    
    e0=edge(anchor1,anchor2)
    hullU=[edge(anchor1,anchorU),edge(anchorU,anchor2)]
    hullL=[edge(anchor1,anchorL),edge(anchorL,anchor2)]
    
    exteriorU=[]
    exteriorL=[]
    triU=triangle(anchor1,anchorU,anchor2)
    triU.edges=[e0]+hullU
    triL=triangle(anchor1,anchorL,anchor2)
    triL.edges=[e0]+hullL
    for p in P:
        if not p.is_point(anchor1) and not p.is_point(anchorU) and not p.is_point(anchorL) and not p.is_point(anchor2):
            res=-(anchor2.y-anchor1.y)*p.x+(anchor2.x-anchor1.x)*p.y-anchor2.x*anchor1.y+anchor2.y*anchor1.x
            if (anchor2.x-anchor1.x)*res>abs(errtol) and not triU.point_triangle_intersect(p,includeboundary=False):
                exteriorU.append(p)
            elif (anchor2.x-anchor1.x)*res<-abs(errtol) and not triL.point_triangle_intersect(p,includeboundary=False):
                exteriorL.append(p)
    
    counter=0
    while len(exteriorU)>0 and counter<len(P):
        for i in range(len(hullU)):
            res1=hullU[i].points[1].x-hullU[i].points[0].x
            if res1!=0:
                anchor=hullU[i].points[0]
                height=0
                for p in exteriorU:
                    res2=-(hullU[i].points[1].y-hullU[i].points[0].y)*p.x+(hullU[i].points[1].x-hullU[i].points[0].x)*p.y-hullU[i].points[1].x*hullU[i].points[0].y+hullU[i].points[1].y*hullU[i].points[0].x
                    if res1*res2>abs(errtol) and res2>height+abs(errtol):
                        anchor=p
                        height=res2

                if height>abs(errtol) or not anchor.is_point(hullU[i].points[0]):
                    e1=edge(hullU[i].points[0],anchor)
                    e2=edge(anchor,hullU[i].points[1])
                    hullU+=[e1,e2]
                    
                    if height>abs(errtol):
                        tri=triangle(hullU[i].points[0],hullU[i].points[1],anchor)
                        for j in range(len(exteriorU)-1,-1,-1):
                            if tri.point_triangle_intersect(exteriorU[j],includeboundary=False):
                                del exteriorU[j]    
                    
                    exteriorU.remove(anchor)
                    del hullU[i]
                    break
        counter+=1

    counter=0
    while len(exteriorL)>0 and counter<len(P):
        for i in range(len(hullL)):
            res1=hullL[i].points[1].x-hullL[i].points[0].x
            if res1!=0:
                anchor=hullL[i].points[0]
                height=0
                for p in exteriorL:
                    res2=-(hullL[i].points[1].y-hullL[i].points[0].y)*p.x+(hullL[i].points[1].x-hullL[i].points[0].x)*p.y-hullL[i].points[1].x*hullL[i].points[0].y+hullL[i].points[1].y*hullL[i].points[0].x
                    if res1*res2<-abs(errtol) and res2+abs(errtol)<height:
                        anchor=p
                        height=res2

                if height<-abs(errtol) or not anchor.is_point(hullL[i].points[0]):
                    e1=edge(hullL[i].points[0],anchor)
                    e2=edge(anchor,hullL[i].points[1])
                    hullL+=[e1,e2]
                    
                    if height<-abs(errtol):
                        tri=triangle(hullL[i].points[0],hullL[i].points[1],anchor)
                        for j in range(len(exteriorL)-1,-1,-1):
                            if tri.point_triangle_intersect(exteriorL[j],includeboundary=False):
                                del exteriorL[j]    
                    
                    exteriorL.remove(anchor)
                    del hullL[i]
                    break
        counter+=1
        
    return hullU+hullL
    
def scan_triangulation(P,errtol=1e-12):
    N=len(P)
    P.sort(key=lambda p:p.x)
    E=[]
    T=[]
    
    if N<2:
        return E,T
    elif N==2:
        e=edge(P[0],P[1])
        e.update()
        E.append(e)
        return E,T
    elif all([abs((P[1].x-P[0].x)*(P[i].y-P[0].y)-(P[i].x-P[0].x)*(P[1].y-P[0].y))<=abs(errtol) for i in range(2,N)]):
        print("Error: all points are colinear")
        E=[edge(P[i],P[i+1]) for i in range(N-1)]
        return E,T
    
    startindex=-1
    for i in range(2,N):
        t0=triangle(P[0],P[1],P[i])
        if t0.area()>abs(errtol):
            T.append(t0)
            e1=edge(P[0],P[1])
            e1.update()
            e2=edge(P[0],P[i])
            e2.update()
            e3=edge(P[1],P[i])
            e3.update()
            t0.edges=[e1,e2,e3]
            t0.update()
            E=[e1,e2,e3]
            startindex=i
            break
        else:
            continue
    
    if startindex==-1:
        print("All points are colinear")
    else:
        for i in range(2,N):
            if i!=startindex:
                newE=[]
                newEbool=[]
                for e in E:
                    if not e.enclosed:
                        e1=edge(e.points[0],P[i])
                        match1bool=False
                        for j in range(len(newE)):
                            if e1.is_edge(newE[j]):
                                e1=newE[j]
                                e1bool=newEbool[j]
                                match1bool=True
                                break
                        if not match1bool:
                            e1bool=any([elem.edge_edge_intersect(e1,includeboundary=False) for elem in E])
                        
                        e2=edge(e.points[1],P[i])
                        match2bool=False
                        for j in range(len(newE)):
                            if e2.is_edge(newE[j]):
                                e2=newE[j]
                                e2bool=newEbool[j]
                                match2bool=True
                                break
                        if not match2bool:
                            e2bool=any([elem.edge_edge_intersect(e2,includeboundary=False) for elem in E])
                        
                        if not e1bool and not e2bool:
                            t=triangle(e.points[0],e.points[1],P[i])
                            t.edges=[e,e1,e2]
                            if all([not t.triangle_triangle_intersect(tri,includeboundary=False) for tri in T]):
                                e.enclosed=True
                                if not match1bool:
                                    e1.update()
                                    newE.append(e1)
                                    newEbool.append(e1bool)
                                else:
                                    e1.enclosed=True
                                if not match2bool:
                                    e2.update()
                                    newE.append(e2)
                                    newEbool.append(e2bool)
                                else:
                                    e2.enclosed=True
                                t.update()
                                T.append(t)
                E=E+newE
    return P,E,T

def delaunize(P,E,T,errtol=1e-12):
    delaunaybool=False
    counter=0
    N=len(E)
    while not delaunaybool and counter<N*(N-1)/2+1:
        delaunaybool=True
        for e in E:
            if e.enclosed:
                p1=e.triangles[0].points[0]
                if p1 in e.points:
                    p1=e.triangles[0].points[1]
                    if p1 in e.points:
                        p1=e.triangles[0].points[2]
                
                p2=e.triangles[1].points[0]
                if p2 in e.points:
                    p2=e.triangles[1].points[1]
                    if p2 in e.points:
                        p2=e.triangles[1].points[2]
                
                if e.triangles[0].inCircumcircle(p2):
                    if all([abs(edge(p1,e.points[i]).cross(edge(p2,e.points[i])))>abs(errtol) for i in range(2)]):
                        t1=e.triangles[0]
                        t2=e.triangles[1]
                                              
                        for elem in t1.edges:
                            elem.triangles.remove(t1)
                        for elem in t1.points:
                            elem.triangles.remove(t1)
                        for elem in t2.edges:
                            elem.triangles.remove(t2)
                        for elem in t2.points:
                            elem.triangles.remove(t2)
                        for elem in e.points:
                            elem.edges.remove(e)
                        
                        temppoints=[elem for elem in e.points]
                        e.points=[p1,p2]
                        
                        tempt1edges=[elem for elem in t1.edges if not e.is_edge(elem)]
                        tempt2edges=[elem for elem in t2.edges if not e.is_edge(elem)]
                        if tempt1edges[0].edge_edge_intersect(tempt2edges[0],includeboundary=True):
                            t1.edges=[tempt1edges[0],tempt2edges[0],e]
                            t2.edges=[tempt1edges[1],tempt2edges[1],e]
                        else:
                            t1.edges=[tempt1edges[0],tempt2edges[1],e]
                            t2.edges=[tempt1edges[1],tempt2edges[0],e]
                        
                        if temppoints[0] in tempt1edges[0].points:
                            t1.points=[p1,p2,temppoints[0]]
                            t2.points=[p1,p2,temppoints[1]]
                        else:
                            t1.points=[p1,p2,temppoints[1]]
                            t2.points=[p1,p2,temppoints[0]]
                            
                        t1.update()
                        t2.update()
                        e.update()
                        t1.triagarea=None
                        t2.triagarea=None
                        e.edgelength=None
                        delaunaybool=False
            counter+=1

def insertVertexDelaunay(point,P,E,T):
    polyedges=[]
    staredges=[]
    constraintpolypoints=[]
    for i in range(len(T)-1,-1,-1):
        if T[i].inCircumcircle(point):
            if T[i].point_triangle_intersect(point):
                for j in range(len(T[i].edges)):
                    if T[i].edges[j].point_edge_intersect(point,includeboundary=False):
                        if T[i].edges[j].constraint:
                            constraintpolypoints+=T[i].edges[j].points
                        break
                for e in T[i].edges:
                    if e not in polyedges:
                        polyedges.append(e)
                for p in T[i].points:
                    starE=edge(point,p)
                    if all([not starE.is_edge(elem) for elem in staredges]):
                        staredges.append(starE)
                
                for p in T[i].points:
                    p.triangles.remove(T[i])
                for e in T[i].edges:
                    e.triangles.remove(T[i])
                del T[i]
            else:
                intersectbool=False
                for j,e in enumerate(T[i].edges):
                    for p in T[i].points:
                        if not p.is_point(e.points[0]) and not p.is_point(e.points[1]):
                            starE=edge(point,p)
                            if starE.edge_edge_intersect(e,includeboundary=False):
                                #if all([not starE.edge_edge_intersect(E[k],includeboundary=False) for k in constraintIndex]):
                                if all([not starE.edge_edge_intersect(elem,includeboundary=False) for elem in E if elem.constraint]):
                                    for k in range(j):
                                        if T[i].edges[k] not in polyedges:
                                            polyedges.append(T[i].edges[k])
                                    for k in range(j+1,3):
                                        if T[i].edges[k] not in polyedges:
                                            polyedges.append(T[i].edges[k])
                                    for elemP in T[i].points:    
                                        starE=edge(point,elemP)
                                        if all([not starE.is_edge(elem) for elem in staredges]):
                                            staredges.append(starE)
                                    intersectbool=True
                                    for elemP in T[i].points:
                                        elemP.triangles.remove(T[i])
                                    for elemE in T[i].edges:
                                        elemE.triangles.remove(T[i])
                                    del T[i]
                            break
                    if intersectbool:
                        break
    
    for i in range(len(polyedges)-1,-1,-1):
        if any([polyedges[i].edge_edge_intersect(starE,includeboundary=False) for starE in staredges]):
            polyedges[i].points[0].edges.remove(polyedges[i])
            polyedges[i].points[1].edges.remove(polyedges[i])
            E.remove(polyedges[i])
            del polyedges[i]
            #indexE=E.index(polyedges[i])
            #del polyedges[i]
            #del E[indexE]
            #for j in range(len(constraintIndex)-1,-1,-1):
            #    if constraintIndex[j]>indexE:
            #        constraintIndex[j]-=1
            #    elif constraintIndex[j]==indexE:
            #        del constraintIndex[j]
    
    for e in polyedges:
        e1=edge(e.points[0],point)
        e2=edge(e.points[1],point)
        
        t=triangle(e.points[0],e.points[1],point)
        t.edges=[e]
        matchbool=False
        for elem in E:
            if e1.is_edge(elem):
                matchbool=True
                t.edges.append(elem)
                break
        if not matchbool:
            if e.points[0] in constraintpolypoints:
                e1.constraint=True
                #constraintIndex.append(len(E))
            E.append(e1)
            t.edges.append(e1)
            e1.update()
            e1.enclosed=True
        matchbool=False
        for elem in E:
            if e2.is_edge(elem):
                matchbool=True
                t.edges.append(elem)
                break
        if not matchbool:
            if e.points[1] in constraintpolypoints:
                e2.constraint=True
                #constraintIndex.append(len(E))
            E.append(e2)
            t.edges.append(e2)
            e2.update()
            e2.enclosed=True
        t.update()
        T.append(t)

def Bowyer_Watson(P,constraints=[]):    
    maxX=P[0].x
    minX=P[0].x
    maxY=P[0].y
    minY=P[0].y
    for i in range(1,len(P)):
        if P[i].x>maxX:
            maxX=P[i].x
        elif P[i].x<minX:
            minX=P[i].x
        if P[i].y>maxY:
            maxY=P[i].y
        elif P[i].y<minY:
            minY=P[i].y
    cx=(maxX+minX)/2
    cy=(maxY+minY)/2
    r=math.sqrt((maxX-minX)**2+(maxY-minY)**2)/2+1
    p01=point(cx,cy+2*r)
    p02=point(cx+math.sqrt(3)*r,cy-r)
    p03=point(cx-math.sqrt(3)*r,cy-r)
    P.insert(0,p03)
    P.insert(0,p02)
    P.insert(0,p01)
    
    e01=edge(p01,p02)
    e01.update()
    e02=edge(p01,p03)
    e02.update()
    e03=edge(p02,p03)
    e03.update()
    E=[e01,e02,e03]
    
    t0=triangle(p01,p02,p03)
    t0.edges=[e01,e02,e03]
    t0.update()
    T=[t0]
    
    for p in P[3:]:
        insertVertexDelaunay(p,P,E,T)
        
    for i in range(len(T)-1,-1,-1):
        if (p01 in T[i].points) or (p02 in T[i].points) or (p03 in T[i].points):
            for elem in T[i].points:
                elem.triangles.remove(T[i])
            for elem in T[i].edges:
                elem.triangles.remove(T[i])
                elem.enclosed=False
            del T[i]
    
    for i in range(len(E)-1,-1,-1):
        if (p01 in E[i].points) or (p02 in E[i].points) or (p03 in E[i].points):
            for elem in E[i].points:
                elem.edges.remove(E[i])
            del E[i]
    
    del P[0]
    del P[0]
    del P[0]
    
    if constraints is not None:
        #constraintIndex=[]
        if isinstance(constraints,list) or isinstance(constraints,tuple) or isinstance(constraints,np.ndarray):
            if all([isinstance(constr,edge) for constr in constraints]):
                for i in range(len(constraints)-1,-1,-1):
                    matchbool=False
                    for j in range(len(E)):
                        if constraints[i].is_edge(E[j]):
                            matchbool=True
                            E[j].constraint=True
                            #constraintIndex.append(j)
                            break
                    if matchbool:
                        #print("Constraint "+constraints[i].to_string()+" is already included in the triangulation")
                        del constraints[i]
                for i in range(len(constraints)-1):
                    for j in range(len(constraints)-1,i,-1):
                        if constraints[i].is_edge(constraints[j]):
                            print("Constraint "+constraints[j].to_string()+" is already included in constraints set")
                            del constraints[j]
                        elif constraints[i].edge_edge_intersect(constraints[j],includeboundary=False):
                            print("Constraint "+constraints[j].to_string()+" intersects with constraint "+constraints[i].to_string()+"; the former has been deleted")
                            del constraints[j]
                for i in range(len(constraints)-1,-1,-1):
                    p1bool=False
                    p2bool=False
                    for p in P:
                        if constraints[i].points[0].is_point(p):
                            p1bool=True
                        elif constraints[i].points[1].is_point(p):
                            p2bool=True
                        if p1bool and p2bool:
                            break
                    if not p1bool or not p2bool:
                        print("One of the points in Constraint "+constraints[i].to_string()+" is not included in the triangulation")
                        del constraints[i]
                
                for constr in constraints:
                    constr.constraint=True
                
                #constrained(P,E,T,constraints,constraintIndex)
                constrained(P,E,T,constraints)
    return P,E,T

def constrained(P,E,T,constraints):
    def patch(constr,anchor,cavityPoints,involvedE):
        nonlocal E
        nonlocal T
        while not anchor.points[1].is_point(constr.points[0]):
            paths=[]
            for e in anchor.points[1].edges:
                if not e.enclosed and not e.is_edge(anchor):
                    if e.points[1].is_point(anchor.points[1]):
                        e.swap()
                    paths.append(e)    
            if len(paths)==0:
                anchor.swap()
                cavityPoints.append(anchor.points[1])
            elif len(paths)==1:
                anchor=paths[0]
                cavityPoints.append(anchor.points[1])
            else:
                cwIndex=0
                thetaprime=math.atan2(anchor.points[0].y-anchor.points[1].y,
                                  anchor.points[0].x-anchor.points[1].x)
                if thetaprime<0:
                    thetaprime+=2*math.pi
                minAngle=math.atan2(paths[0].points[1].y-paths[0].points[0].y,paths[0].points[1].x-paths[0].points[0].x)-thetaprime
                while minAngle<=0:
                    minAngle+=2*math.pi
                for i in range(1,len(paths)):
                    res=math.atan2(paths[i].points[1].y-paths[i].points[0].y,paths[i].points[1].x-paths[i].points[0].x)-thetaprime
                    while res<=0:
                        res+=2*math.pi
                    if res<minAngle:
                        minAngle=res
                        cwIndex=i
                anchor=paths[cwIndex]
                cavityPoints.append(anchor.points[1])
        
        segments=[cavityPoints]
        segmentE=[constr]
        while len(segments)>0:
            for i in range(len(segments)-1,-1,-1):
                if len(segments[i])==3:
                    e1=segments[i][1].edges[0]
                    e2=segments[i][1].edges[1]
                    ref1=edge(segments[i][0],segments[i][1])
                    ref2=edge(segments[i][-1],segments[i][1])
                    for j in range(len(segments[i][1].edges)):
                        if segments[i][1].edges[j].is_edge(ref1):
                            e1=segments[i][1].edges[j]
                        elif segments[i][1].edges[j].is_edge(ref2):
                            e2=segments[i][1].edges[j]
                    involvedE+=[e1,e2]
                    tri=triangle(segments[i][0],segments[i][1],segments[i][-1])
                    tri.edges=[e1,e2,segmentE[i]]
                    tri.update()
                    T.append(tri)
                    del segments[i]
                    del segmentE[i]
                else:
                    for j in range(1,len(segments[i])-1):
                        tri=triangle(segments[i][0],segments[i][j],segments[i][-1])
                        validbool=True
                        for k in range(1,len(segments[i])-1):
                            if k!=j and tri.inCircumcircle(segments[i][k],includeboundary=False):
                                validbool=False
                                break
                        if validbool:
                            if j==1:
                                e1=segments[i][1].edges[0]
                                ref1=edge(segments[i][0],segments[i][1])
                                for k in range(1,len(segments[i][1].edges)):
                                    if segments[i][1].edges[k].is_edge(ref1):
                                        e1=segments[i][1].edges[k]
                                        break
                                e2=edge(segments[i][1],segments[i][-1])
                                e2.update()
                                E.append(e2)
                                involvedE+=[e1,e2]
                                tri.edges=[e1,e2,segmentE[i]]
                                tri.update()
                                T.append(tri)
                                segments[i]=segments[i][1:]
                                segmentE[i]=e2
                            elif j==len(segments[i])-2:
                                e1=edge(segments[i][0],segments[i][-2])
                                e1.update()
                                e2=segments[i][-2].edges[0]
                                ref2=edge(segments[i][-2],segments[i][-1])
                                for k in range(1,len(segments[i][-2].edges)):
                                    if segments[i][-2].edges[k].is_edge(ref2):
                                        e2=segments[i][-2].edges[k]
                                        break
                                E.append(e1)
                                involvedE+=[e1,e2]
                                tri.edges=[e1,e2,segmentE[i]]
                                tri.update()
                                T.append(tri)
                                segments[i]=segments[i][0:-1]
                                segmentE[i]=e1
                            else:
                                e1=edge(segments[i][0],segments[i][j])
                                e1.update()
                                e2=edge(segments[i][j],segments[i][-1])
                                e2.update()
                                E+=[e1,e2]
                                involvedE+=[e1,e2]
                                tri.edges=[segmentE[i],e1,e2]
                                tri.update()
                                T.append(tri)
                                segments.append(segments[i][j:])
                                segmentE.append(e2)
                                segments[i]=segments[i][0:j+1]
                                segmentE[i]=e1
                            break
    
    for constr in constraints:
        for i in range(len(T)-1,-1,-1):
            if T[i].edge_triangle_intersect(constr,includeboundary=False) or any([elem.edge_edge_intersect(constr,includeboundary=False) for elem in T[i].edges]):
                for elem in T[i].points:
                    elem.triangles.remove(T[i])
                for elem in T[i].edges:
                    if elem.edge_edge_intersect(constr,includeboundary=False):
                        try:
                            elem.points[0].edges.remove(elem)
                            elem.points[1].edges.remove(elem)
                            E.remove(elem)
                            #indexE=E.index(elem)
                            #del E[indexE]
                            #for j in range(len(constraintIndex)):
                            #    if constraintIndex[j]>indexE:
                            #        constraintIndex[j]-=1
                            #    elif constraintIndex[j]==indexE:
                            #        del constraintIndex[j]
                        except:
                            pass
                    else:
                        elem.triangles.remove(T[i])
                        elem.enclosed=False
                del T[i]
        
        #constr.enclosed=True
        constr.update()
        #constraintIndex.append(len(E))
        E.append(constr)
        
        involvedE=[constr]
        
        anchor=constr
        cavityPoints=[constr.points[1]]
        patch(constr,anchor,cavityPoints,involvedE)
        
        anchor=constr
        anchor.swap()
        cavityPoints=[constr.points[1]]
        patch(constr,anchor,cavityPoints,involvedE)
        
        for e in involvedE:
            if len(e.triangles)<2:
                e.enclosed=False
            else:
                e.enclosed=True
                
def Ruppert(P,E,T,theta=20,areatol=None,errtol=1e-12):
        theta*=math.pi/180
        count=0
        while theta<0 and count<100:
            theta+=2*math.pi
            count+=1
        if count>=100:
            print("Parameter theta too small. Resorting to theta=20 degrees")
            theta=math.pi/9
        count=0
        while theta>2*math.pi and count<100:
            theta-=2*math.pi
            count+=1
        if count>=100:
            print("Parameter theta too large. Resorting to theta=20 degrees")
            theta=math.pi/9
            
        center=[]
        for p in P:
            constraintTotal=sum([1 for e in p.edges if e.constraint])
            if constraintTotal>1:
                cascadebool=False
                minEdgeLength=p.edges[0].length()
                minAbsCosTheta=1
                for i in range(len(p.edges)-1):
                    if p.edges[i].constraint:
                        for j in range(1,len(p.edges)):
                            if p.edges[j].constraint:
                                a=max(p.edges[i].length(),p.edges[j].length())
                                b=min(p.edges[i].length(),p.edges[j].length())
                                if b<minEdgeLength:
                                    minEdgeLength=b
                                if p.edges[i].points[0].is_point(p) and p.edges[j].points[0].is_point(p):
                                    c2=(p.edges[i].points[1].x-p.edges[j].points[1].x)**2+(p.edges[i].points[1].y-p.edges[j].points[1].y)**2
                                elif p.edges[i].points[0].is_point(p) and p.edges[j].points[1].is_point(p):
                                    c2=(p.edges[i].points[1].x-p.edges[j].points[0].x)**2+(p.edges[i].points[1].y-p.edges[j].points[0].y)**2
                                elif p.edges[i].points[1].is_point(p) and p.edges[j].points[0].is_point(p):
                                    c2=(p.edges[i].points[0].x-p.edges[j].points[1].x)**2+(p.edges[i].points[0].y-p.edges[j].points[1].y)**2
                                else:
                                    c2=(p.edges[i].points[0].x-p.edges[j].points[0].x)**2+(p.edges[i].points[0].y-p.edges[j].points[0].y)**2
                                costheta=(a**2+b**2-c2)/(2*a*b)
                                if abs(costheta)<minAbsCosTheta:
                                    minAbsCosTheta=abs(costheta)
                                if costheta>abs(1e-12):
                                    n=math.floor(math.log(a,2)+math.log(costheta,2)-math.log(b,2))
                                    if b*2**n<a*costheta-abs(errtol) and a<b*costheta*2**(n+1)-abs(errtol):
                                        cascadebool=True
                    if cascadebool:
                        break
                if cascadebool:
                    sintheta=math.sqrt(1-minAbsCosTheta**2)
                    if areatol is not None and abs(areatol)<minEdgeLength**2*sintheta/2-abs(errtol):
                        minEdgeLength=math.sqrt(2*abs(areatol)/sintheta)
                    for e in p.edges:
                        if e.constraint:
                            if e.points[0].is_point(p):
                                center.append(point(p.x+0.5*(e.points[1].x-p.x)*minEdgeLength/e.length(),
                                                     p.y+0.5*(e.points[1].y-p.y)*minEdgeLength/e.length()))
                            else:
                                center.append(point(p.x+0.5*(e.points[0].x-p.x)*minEdgeLength/e.length(),
                                                     p.y+0.5*(e.points[0].y-p.y)*minEdgeLength/e.length()))
        for c in center:
            insertVertexDelaunay(c,P,E,T)
    
        
        donebool=False
        count=0
        while not donebool and count<500:
            donebool=True
            for e in E:
                if e.constraint:
                    encroachedbool=False
                    for elem in e.points[0].edges:
                        if (elem.points[0].is_point(e.points[0]) and e.inCircumcircle(elem.points[1],includeboundary=False)) or (elem.points[1].is_point(e.points[0]) and e.inCircumcircle(elem.points[0],includeboundary=False)):
                            center=point((e.points[1].x+e.points[0].x)/2,(e.points[1].y+e.points[0].y)/2)
                            insertVertexDelaunay(center,P,E,T)
                            encroachedbool=True
                            donebool=False
                            break
                    if not encroachedbool:
                        for elem in e.points[1].edges:
                            if (elem.points[0].is_point(e.points[1]) and e.inCircumcircle(elem.points[1],includeboundary=False)) or (elem.points[1].is_point(e.points[1]) and e.inCircumcircle(elem.points[0],includeboundary=False)):
                                center=point((e.points[1].x+e.points[0].x)/2,(e.points[1].y+e.points[0].y)/2)
                                insertVertexDelaunay(center,P,E,T)
                                encroachedbool=True
                                donebool=False
                                break
                if not donebool:
                    break
            count+=1
        
        donebool=False
        count=0
        while not donebool and count<500:
            donebool=True
            for t in T:
                a=t.edges[0].length()
                b=t.edges[1].length()
                c=t.edges[2].length()
                minLengthIndex=0
                if a<=b and a<=c:
                    minAngle=math.acos((b**2+c**2-a**2)/(2*b*c))
                elif b<=a and b<=c:
                    minLengthIndex=1
                    minAngle=math.acos((a**2+c**2-b**2)/(2*a*c))
                else:
                    minLengthIndex=2
                    minAngle=math.acos((a**2+b**2-c**2)/(2*a*b))
                if minAngle<theta-abs(errtol) or (areatol is not None and t.area()>abs(areatol)+abs(errtol)):
                    if sum([1 for e in t.edges if e.constraint])<2:
                        center,rad=t.circumcircle()
                        if 2*math.asin(t.edges[minLengthIndex].length()/(2*rad))<theta-abs(errtol):
                            edgeCenter=point((t.edges[minLengthIndex].points[0].x+t.edges[minLengthIndex].points[1].x)/2,
                                             (t.edges[minLengthIndex].points[0].y+t.edges[minLengthIndex].points[1].y)/2)
                            dist=center.distance(edgeCenter)
                            distprime=t.edges[minLengthIndex].length()/(2*math.tan(theta/2))
                            center.x=edgeCenter.x+0.99*(center.x-edgeCenter.x)*distprime/dist
                            center.y=edgeCenter.y+0.99*(center.y-edgeCenter.y)*distprime/dist
                        encroachedbool=False
                        for e in E:
                            if e.constraint and e.inCircumcircle(center,includeboundary=False):
                                insertVertexDelaunay(point((e.points[1].x+e.points[0].x)/2,(e.points[1].y+e.points[0].y)/2),P,E,T)
                                encroachedbool=True
                                break
                        if not encroachedbool:
                            insertVertexDelaunay(center,P,E,T)
                        donebool=False
                if not donebool:
                    break
            count+=1