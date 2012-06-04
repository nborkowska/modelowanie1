#! /usr/bin/env python

import sys
import numpy as np
import math


class Atom(object):
    
    def __init__(self, position, mass=1.0, energy = 0.0):
        self.position = position
        self.mass = mass
        self.energy = energy                     #energia kinetyczna
        self.force = np.zeros(Atoms.dim)
        self.velocity = np.zeros(Atoms.dim)

class RandomSample(object):

    @staticmethod 
    def getSample(mu, sigma, dim):                  #losowanie polozen poczatkowych
        return 1.5*np.random.normal(mu, sigma, dim) #to juz zwraca wektor
    
    @staticmethod
    def getSampleSet(mu, sigma, size, dim):
        return [RandomSample.getSample(mu, sigma, dim) for i in range(size)]

class Atoms(object):            # zbior czasteczek 
    
    dim = 1                    #na sile, bo mozemy miec kilka pudelek z atomami ale musimy rozwac ten sam wymiar ;]

    def __init__(self, size, atoms=[]):            #size - l.czasteczek, dim-wymiar
        self.size = size
        self.atoms = [Atom(vector) for vector in \
                RandomSample.getSampleSet(0, 10, self.size, self.dim)]

class ForceField(object):
""" 'Abstract' class"""

    def singleForce(self, atom):
        raise NotImplementedError
    
    def singleEnergy(self, atom):
        raise NotImplementedError
    
    def pairForce(self, atom):
        raise NotImplementedError
    
    def pairEnergy(self, atom1, atom2):
        raise NotImplementedError

class SoftWalls(ForceField):
    
    f, L = 0.2, 5                               #10, 10 nie wychodzi, nie potrzeba atrybutow instancji, wystarcza atrybuty klasy

    def singleForce(self, atom):         
        distance = np.linalg.norm(atom.position)       #tu we wzorze sa 2 rozne r_i!
        if distance < self.L:
            atom.force = np.zeros(Atoms.dim) 
        else:
            atom.force = self.f*(self.L-distance)*atom.position/distance # bo tu jest juz rozniczka, a minus ze wzoru na sile
            #jesli podzielimy wektor przez jego dlugosc, to dostaniemy wektor unormowany - czyli o dlugosci 1

    def singleEnergy(self, atom):
        distance = np.linalg.norm(atom.position)
        if distance < self.L:
            atom.energy = 0
        else:
            atom.energy = 0.5*self.f*(self.L - distance)**2
    
    def pairForce(self, atom):
        pass
    
    def pairEnergy(self, atom1, atom2):
        pass


class MBM(ForceField):
    
    a, b, c, d = 5, 10, 3, 0.02     #parametry ze skryptu
    
    def singleForce(self, atom):
        x = atom.position
        atom.force = self.a*math.exp(-self.b*(x-1)**2)*2*self.b*(x-1) \
                + self.c*math.exp(-(x+1)**2)*2*(1+x) \
                + 4*self.d*x**3
    
    def singleEnergy(self, atom):
        x = atom.position[0]
        atom.energy = self.a*math.exp(-self.b*(x-1)**2) \
                - self.c*math.exp(-(x+1)**2) + self.d*x**4

    def pairForce(self, atom):
        pass
    
    def pairEnergy(self, atom1, atom2):
        pass


class LenardJones(ForceField):
    
    R, e = 1.0, 1.0

    def singleForce(self, atom):
        pass
    
    def singleEnergy(self, atom):
        pass
    
    def setParams(self, atom1, atom2):
        direction = atom1.position-atom2.position # wektor kierunkowy
        distance = np.linalg.norm(atom1.position-atom2.position) #odleglosc euklidesowa
        a = (self.R/distance)**6
        print 'distance: ', distance, '\na: ', a
        return direction/distance, a
    
    def pairForce(self, atom1, atom2):
        normalized, a = self.setParams(atom1, atom2)
        F = -12.0*self.e*a*(a-1)*normalized
        atom1.force += F
        atom2.force -= F

    def pairEnergy(self, atom1, atom2):
        normalized, a = self.setParams(atom1, atom2)
        E = self.e*a*(a-2)
        #print 'ss',a**2 - 2*a co robic z ujemna energia :(?
        atom1.energy, atom2.energy = E/2, E/2

def help():

    info = "Usage: \n ./symulacje.py p i n s t \nwhere \n \
            p - used potential, available options: \n \
                0 - soft walls, 1 - mbm, 2 - Lenard-Jones \n \n \
            i - intergation function, available options: \n \
                0 - basic Verlet, 1 - Velocity Verlet, 2 - Leapfrog\n \n \
            n - number of molecules \n \n \
            s - number of simulation steps \n \n \
            t - step length \n \n \
            Example: ./symulacje.py 1 0 10 25 1 \n"
    
    return info

def main(*args):
    
    if len(args) < 6:
        print 'Too few parameters.\n', help()
        return 0
    else:
        potential, integration, no_molecules, no_steps, steps_length = args[1:]
        a=Atoms(5)
        b = SoftWalls()
        c = MBM()
        d = LenardJones()
        for i in a.atoms:
            print i.position, b.singleForce(i), c.singleForce(i)
        print a.atoms[0].force
        c.singleEnergy(a.atoms[0])
        b.singleEnergy(a.atoms[3])
        d.pairEnergy(a.atoms[1], a.atoms[2])
        print a.atoms[0].energy, a.atoms[1].energy, a.atoms[3].energy
                                         
if __name__ == '__main__':
    sys.exit(main(*sys.argv))

"""
a=Atoms(5)
for i in a.atoms:
    print i.position.array()
wersja ustalona jest chyba taka, ze liczymy sily tak ze jak 5 atomow i miekkie scianki
to licze miekkie scianki tylko, a jak lj to robie pary i licze tylko dla par"""
"""Klasa wektor zostala zastapiona numpy.array - nie trzeba kombinowac z dzieleniem etc"""
