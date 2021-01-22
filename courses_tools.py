import numpy as np
from scipy.optimize import minimize

def altitude2density(h): 
    p0 = 101325 #pa
    T0=288.15 #k
    g=9.8 #m/s-2
    L=0.0065 # temperature lapse rate: k/m
    R = 8.31447 #J/(mol.K)
    M = 0.0289644 #kg/mol

    T = T0-L*h
    p = p0 * (1.- L*h/T0)**(g*M/(R*L))

    rho = p*M/(R*T)
    return rho


def weibull_distribution(v,A,k):
    v_adim = 1.*v/A
    fv = k * v_adim**(k-1.) * np.exp(- (v_adim)**k )/A
    return fv

data_weibull = np.array([[  0.00000000e+00,   5.00000000e+00],
       [  1.00000000e+00,   5.70000000e+00],
       [  2.00000000e+00,   7.50000000e+00],
       [  3.00000000e+00,   8.40000000e+00],
       [  4.00000000e+00,   1.00000000e+01],
       [  5.00000000e+00,   1.03000000e+01],
       [  6.00000000e+00,   9.50000000e+00],
       [  7.00000000e+00,   8.00000000e+00],
       [  8.00000000e+00,   6.10000000e+00],
       [  9.00000000e+00,   5.00000000e+00],
       [  1.00000000e+01,   4.80000000e+00],
       [  1.10000000e+01,   4.50000000e+00],
       [  1.20000000e+01,   2.70000000e+00],
       [  1.30000000e+01,   2.70000000e+00],
       [  1.40000000e+01,   1.80000000e+00],
       [  1.50000000e+01,   1.50000000e+00],
       [  1.60000000e+01,   1.35000000e+00],
       [  1.70000000e+01,   1.20000000e+00],
       [  1.80000000e+01,   1.10000000e+00],
       [  1.90000000e+01,   1.00000000e+00],
       [  2.00000000e+01,   7.50000000e-01],
       [  2.10000000e+01,   5.00000000e-01],
       [  2.20000000e+01,   3.00000000e-01],
       [  2.30000000e+01,   2.00000000e-01],
       [  2.40000000e+01,   1.00000000e-01],
       [  2.50000000e+01,   0.00000000e-02]])

def J_weibull(param):
    A=param[0]
    k=param[1]
    n=0.
    for data in data_weibull:
        v=data[0]
        d=data[1]/100.
        n+=np.abs(weibull_distribution(v,A,k) - d)
    return n

results_weibull = minimize(J_weibull,[6,2]).x
A_opt = results_weibull[0]
k_opt = results_weibull[1]


class triangle:
    def __init__(self,A=False,B=False,C=False,a=False,b=False,c=False,angleA=False,angleB=False,angleC=False,degrees=True):
        self.A=A
        self.B=B
        self.C=C
        self.angleA=angleA
        self.angleB=angleB
        self.angleC=angleC
        self.a=a
        self.b=b
        self.c=c
        self.degrees = degrees
        self.set_r(self.degrees)

    def set_r(self,degrees) :
        if degrees is True :
            self.r = 2.*np.pi /360.
        else :
            self.r = 1.

    def ssl(self,side,side_known):
        return self.side_via_sine_law(side,side_known)

    def side_via_sine_law(self,side,side_known):
        if side.lower() == side_known.lower():
            print('sides have to be different')
            return -1
        if side.lower()=='a':
            S = self.a
            angleS = self.angleA
        if side.lower()=='b':
            S = self.b
            angleS = self.angleB
        if side.lower()=='c':
            S = self.c
            angleS = self.angleC
        if side_known.lower()=='a':
            Sk = self.a
            angleSk = self.angleA
        if side_known.lower()=='b':
            Sk = self.b
            angleSk = self.angleB
        if side_known.lower()=='c':
            Sk = self.c
            angleSk = self.angleC
        if angleS is False :
            print('angle opposite of the side ' + side.lower() + ' has to be known') 
            return -1
        if angleSk is False :
            print('angle opposite of the side ' + side_known.lower() + ' has to be known') 
            return -1
        if Sk is False :
            print('side ' + side_known.lower() + ' has to be known') 
            return -1

        S = Sk * np.sin(self.r*angleS)/np.sin(self.r*angleSk)

        return S

    def asl(self,side,side_known):
        return self.angle_via_sine_law(side,side_known)
    def angle_via_sine_law(self,side,side_known):
        if side.lower() == side_known.lower():
            print('sides have to be different')
            return -1
        if side.lower()=='a':
            S = self.a
            angleS = self.angleA
        if side.lower()=='b':
            S = self.b
            angleS = self.angleB
        if side.lower()=='c':
            S = self.c
            angleS = self.angleC
        if side_known.lower()=='a':
            Sk = self.a
            angleSk = self.angleA
        if side_known.lower()=='b':
            Sk = self.b
            angleSk = self.angleB
        if side_known.lower()=='c':
            Sk = self.c
            angleSk = self.angleC
        if S is False :
            print('The side ' + side.lower() + ' has to be known') 
            return -1
        if angleSk is False :
            print('angle opposite of the side ' + side_known.lower() + ' has to be known') 
            return -1
        if Sk is False :
            print('side ' + side_known.lower() + ' has to be known') 
            return -1

        angleS = (1./self.r) * np.arcsin( S * np.sin(self.r*angleSk)/Sk) 
        return angleS

    def scl(self,side):
        return self.side_via_cos_law(side)
    
    def side_via_cos_law(self,side):
        if side.lower() == side_known.lower():
            print('sides have to be different')
            return -1
        if side.lower()=='a':
            S = self.a
            angleS = self.angleA
            S2 = self.b
            S3 = self.c
        if side.lower()=='b':
            S = self.b
            angleS = self.angleB
            S2 = self.a
            S3 = self.c
        if side.lower()=='c':
            S = self.c
            angleS = self.angleC
            S2 = self.a
            S3 = self.b
        
        if angleS is False :
            print('angle opposite of the side ' + side.lower() + ' has to be known') 
            return -1
        if (S2 is False) or (S3 is False):
            print('side adjacent to' + side.lower() + ' has to be known') 
            return -1

        S = np.sqrt(S2**2 + S3*2 -2*S2*S3 * np.cos(self.r*angleS) )
        return S


        return S


    def acl(self,side):
        return self.angle_via_cos_law(side)
    def angle_via_cos_law(self,side):
        if side.lower()=='a':
            S = self.a
            angleS = self.angleA
            S2 = self.b
            S3 = self.c
        if side.lower()=='b':
            S = self.b
            angleS = self.angleB
            S2 = self.a
            S3 = self.c
        if side.lower()=='c':
            S = self.c
            angleS = self.angleC
            S2 = self.a
            S3 = self.b
        
        if S is False :
            print('angle opposite of the angle ' + side.lower() + ' has to be known') 
            return -1
        if (S2 is False) or (S3 is False):
            print('side adjacent to the angle ' + side.lower() + ' has to be known') 
            return -1

        angleS = (1./self.r)*np.arccos( (S**2 - S2**2 - S3**2)/(2.*S2*S3)  )
        return angleS

    def unknown_sides(self):
        sides = ''
        if self.a is False :
            sides += 'a'
        if self.b is False :
            sides += 'b'
        if self.c is False :
            sides += 'c'
        return sides
    def unknown_angles(self):
        angles = ''
        if self.angleA is False :
            angles += 'a'
        if self.angleB is False :
            angles += 'b'
        if self.angleC is False :
            angles += 'c'
        return angles
    def set_side(self,s,n):
        if s == 'a' :
            self.a = n
        if s == 'b' :
            self.b = n
        if s == 'c' :
            self.c = n
    def set_angle(self,a,n):
        if a == 'a' :
            self.angleA = n
        if a == 'b' :
            self.angleB = n
        if a == 'c' :
            self.angleC = n

    def complete_angle(self):

        if self.degrees == True :
            angle_triangle = 180
        else :
            angle_triangle = np.pi
        unknown_angles =self.unknown_angles()
        if len(unknown_angles)==1 :
            if unknown_angles == 'a' :
                self.angleA = angle_triangle - self.angleB - self.angleC
            if unknown_angles == 'b' :
                self.angleB = angle_triangle - self.angleA - self.angleC
            if unknown_angles == 'c' :
                self.angleC = angle_triangle - self.angleA - self.angleB
        else :
            if len(self.unknown_sides()) == 0 :
                for a in unknown_angles:
                    self.set_angle(a,self.acl(a))
            else : pass

    def solve(self):
        what2use = None
        what2find = None
        unknown_sides = self.unknown_sides()
        unknown_angles = self.unknown_angles()
        if len(unknown_angles) == 1 :
            self.complete_angle()
            self.solve()
            unknown_angles = self.unknown_angles()
        if len(unknown_sides) == 3 :
            print('at least one side is needed')
            return -1
        if len(unknown_angles) == 0 :
            for known_side in 'abc':
                if known_side in unknown_sides : pass
                else : 
                    for side in 'abc':
                        if side in unknown_sides : pass

                
        if len(unknown_sides) == 0 :
            what2use = 'cos_law'
        

    
    
