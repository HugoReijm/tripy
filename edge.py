import math
import numpy.linalg as npla
from point import point

class edge:
    def __init__(self,p1,p2):
        self.points=[p1,p2]
        self.triangles=[]
        self.edgelength=None
        self.enclosed=False
        self.constraint=False
        self.index=-1
        self.circumcenter=None
        self.circumradius=None
        
    def is_edge(self,e,errtol=1e-12):
        return isinstance(e,edge) and ((self.points[0].is_point(e.points[0],errtol=errtol) and self.points[1].is_point(e.points[1],errtol=errtol)) or (self.points[0].is_point(e.points[1],errtol=errtol) and self.points[1].is_point(e.points[0],errtol=errtol)))
        
    def update(self):
        for p in self.points:
            p.edges.append(self)
        
    def length(self):
        if self.edgelength is None:
            self.edgelength=math.sqrt((self.points[0].x-self.points[1].x)**2+(self.points[0].y-self.points[1].y)**2)
        return self.edgelength
    
    def cross(self,edge):
        return (self.points[1].x-self.points[0].x)*(edge.points[1].y-edge.points[0].y)-(self.points[1].y-self.points[0].y)*(edge.points[1].x-edge.points[0].x)
    
    def point_edge_intersect(self,pointVar,includeboundary=True,errtol=1e-12):
        if includeboundary:
            if abs((self.points[1].x-self.points[0].x)*(pointVar.y-self.points[0].y)-(pointVar.x-self.points[0].x)*(self.points[1].y-self.points[0].y))<abs(errtol):
                if ((self.points[1].x-self.points[0].x)**2+(self.points[1].y-self.points[0].y)**2-(pointVar.x-self.points[0].x)**2-(pointVar.y-self.points[0].y)**2>=-abs(errtol)
                    and (self.points[1].x-self.points[0].x)**2+(self.points[1].y-self.points[0].y)**2-(pointVar.x-self.points[1].x)**2-(pointVar.y-self.points[1].y)**2>=-abs(errtol)):
                    return True
        elif abs(pointVar.x-self.points[0].x)>abs(errtol) or abs(pointVar.y-self.points[0].y)>abs(errtol):
            if abs((self.points[1].x-self.points[0].x)*(pointVar.y-self.points[0].y)-(pointVar.x-self.points[0].x)*(self.points[1].y-self.points[0].y))<abs(errtol):
                if ((self.points[1].x-self.points[0].x)**2+(self.points[1].y-self.points[0].y)**2-(pointVar.x-self.points[0].x)**2-(pointVar.y-self.points[0].y)**2>abs(errtol)
                    and (self.points[1].x-self.points[0].x)**2+(self.points[1].y-self.points[0].y)**2-(pointVar.x-self.points[1].x)**2-(pointVar.y-self.points[1].y)**2>abs(errtol)):
                    return True
        return False

    def edge_edge_intersect(self,edge,includeboundary=True,errtol=1e-12):
        if self.is_edge(edge):
            return True
        
        orient1=(self.points[1].y-self.points[0].y)*(edge.points[0].x-self.points[1].x)-(self.points[1].x-self.points[0].x)*(edge.points[0].y-self.points[1].y)
        orient2=(self.points[1].y-self.points[0].y)*(edge.points[1].x-self.points[1].x)-(self.points[1].x-self.points[0].x)*(edge.points[1].y-self.points[1].y)
        orient3=(edge.points[1].y-edge.points[0].y)*(self.points[0].x-edge.points[1].x)-(edge.points[1].x-edge.points[0].x)*(self.points[0].y-edge.points[1].y)
        orient4=(edge.points[1].y-edge.points[0].y)*(self.points[1].x-edge.points[1].x)-(edge.points[1].x-edge.points[0].x)*(self.points[1].y-edge.points[1].y)
        
        if includeboundary:
            if (orient1*orient2<-abs(errtol)) and (orient3*orient4<-abs(errtol)):
                return True
            elif (orient1*orient2<-abs(errtol)) and (abs(orient3*orient4)<=abs(errtol)):
                return True
            elif (abs(orient1*orient2)<=abs(errtol)) and (orient3*orient4<-abs(errtol)):
                return True
            elif (abs(orient1*orient2)<=abs(errtol)) and (abs(orient3*orient4)<=abs(errtol)):
                if any([self.point_edge_intersect(p) for p in edge.points]):
                    return True
                elif any([edge.point_edge_intersect(p) for p in self.points]):
                    return True
        else:
            if (orient1*orient2<-abs(errtol)) and (orient3*orient4<-abs(errtol)):
                return True
            elif (abs(orient1*orient2)<=abs(errtol)) and (abs(orient3*orient4)<=abs(errtol)):
                if any([self.point_edge_intersect(p,includeboundary=False) for p in edge.points]):
                    return True
                elif any([edge.point_edge_intersect(p,includeboundary=False) for p in self.points]):
                    return True
        return False
    
    def circumcircle(self):
        if self.circumcenter is None or self.circumradius is None:
            self.circumcenter=point((self.points[0].x+self.points[1].x)/2,(self.points[0].y+self.points[1].y)/2)
            self.circumradius=math.sqrt((self.points[1].x-self.points[0].x)**2+(self.points[1].y-self.points[0].y)**2)/2
        return self.circumcenter,self.circumradius
    
    def inCircumcircle(self,pointVar,includeboundary=True,errtol=1e-12):
        pointd=point((self.points[0].x+self.points[1].x)/2,(self.points[0].y+self.points[1].y)/2)
        pointc=point(self.points[0].y-pointd.y+pointd.x,pointd.x-self.points[0].x+pointd.y)
        if (self.points[1].x-self.points[0].x)*(pointc.y-self.points[1].y)-(self.points[1].y-self.points[0].y)*(pointc.x-self.points[1].x)>0:
            res=npla.det([[self.points[0].x-pointVar.x,self.points[1].x-pointVar.x,pointc.x-pointVar.x],
                          [self.points[0].y-pointVar.y,self.points[1].y-pointVar.y,pointc.y-pointVar.y],
                          [(self.points[0].x-pointVar.x)**2+(self.points[0].y-pointVar.y)**2,(self.points[1].x-pointVar.x)**2+(self.points[1].y-pointVar.y)**2,(pointc.x-pointVar.x)**2+(pointc.y-pointVar.y)**2]])
        else:
            res=npla.det([[self.points[0].x-pointVar.x,pointc.x-pointVar.x,self.points[1].x-pointVar.x],
                          [self.points[0].y-pointVar.y,pointc.y-pointVar.y,self.points[1].y-pointVar.y],
                          [(self.points[0].x-pointVar.x)**2+(self.points[0].y-pointVar.y)**2,(pointc.x-pointVar.x)**2+(pointc.y-pointVar.y)**2,(self.points[1].x-pointVar.x)**2+(self.points[1].y-pointVar.y)**2]])
        if includeboundary:
            return res>=-abs(errtol)
        else:
            return res>abs(errtol)
    
    def swap(self):
        tempP=self.points[0]
        self.points[0]=self.points[1]
        self.points[1]=tempP
    
    def to_string(self):
        return "Edge <"+self.points[0].to_string()+", "+self.points[1].to_string()+">"
    
    def copy(self):
        e=edge(self.points[0].copy(),self.points[1].copy())
        if self.edgelength is not None:
            e.edgelength=self.edgelength
        return e
    
    def draw(self,plotaxis,color="black",alpha=1):
        plotaxis.plot([self.points[0].x,self.points[1].x],
                      [self.points[0].y,self.points[1].y],color=color,alpha=alpha)
        
    def kill(self,E=None):
        for p in self.points:
            try:
                p.edges.remove(self)
            except ValueError:
                pass
        for t in self.triangles:
            try:
                t.edges.remove(self)
            except ValueError:
                pass
        if isinstance(E,list):
            try:
                E.remove(self)
            except ValueError:
                pass
        del self