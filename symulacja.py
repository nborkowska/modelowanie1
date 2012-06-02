#! /usr/bin/env python

import sys
import numpy as np
import math


class Atom(object):
    
    def __init__(self, position, mass=1):
        self.position = position
        self.mass = mass

class RandomSample(object):

    @staticmethod 
    def getSample(mu, sigma, dim):                  #losowanie polozen poczatkowych
        return 1.5*np.random.normal(mu, sigma, dim) #to juz zwraca wektor
    
    @staticmethod
    def getSampleSet(mu, sigma, size, dim):
        return [RandomSample.getSample(mu, sigma, dim) for i in range(size)]

class Atoms(object):            # zbior czasteczek 
    
    dim = 1                    #na sile, bo mozemy miec kilka pudelek z tomami ale musimy rozwac ten sam wymiar ;]

    def __init__(self, size, atoms=[]):            #size - l.czasteczek, dim-wymiar
        self.size = size
        self.atoms = [Atom(vector) for vector in \
                RandomSample.getSampleSet(0, 10, self.size, Atoms.dim)]

class ForceField(object):
    
    def singleForce(self, atom):           #albo nazwac to od razu liczeniem sil, ale u kazdego bedzie co innego, chyba
        #ze wprowadze liczenie pochodnej jako klase w ogole
        raise NotImplementedError        #bo jest zdefiniowane dopiero w podklasach

class SoftWalls(ForceField):
    
    f, L = 0.2, 5                               #10, 10 nie wychodzi, nie potrzeba atrybutow instancji, wystarcza atrybuty klasy

    def singleForce(self, atom):         
        distance = np.linalg.norm(atom.position)       #tu we wzorze sa 2 rozne r_i!
        if distance < self.L:
            return np.zeros(Atoms.dim) 
        else:
            return self.f*(self.L-distance)*atom.position/distance # bo tu jest juz rozniczka, a minus ze wzoru na sile
            #jesli podzielimy wektor przez jego dlugosc, to dostaniemy wektor unormowany - czyli o dlugosci 1

class MBM(ForceField):
    
    a, b, c, d = 5, 10, 3, 0.02     #parametry ze skryptu
    
    def singeForce(self, atom):
        x = atom.position.x
        return None


class LenardJones(ForceField):
    pass

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
        for i in a.atoms:
            print i.position, b.singleForce(i)
                                         
if __name__ == '__main__':
    sys.exit(main(*sys.argv))

"""
a=Atoms(5)
for i in a.atoms:
    print i.position.array()
wersja ustalona jest chyba taka, ze liczymy sily tak ze jak 5 atomow i miekkie scianki
to licze miekkie scianki tylko, a jak lj to robie pary i licze tylko dla par"""
"""Klasa wektor zostala zastapiona numpy.array - nie trzeba kombinowac z dzieleniem etc"""
