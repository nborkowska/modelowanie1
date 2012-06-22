#! /usr/bin/env python

import sys
import numpy as np
import math
import matplotlib.pyplot as plt

class Atom(object):
    
    def __init__(self, position, mass=1.0, energy = 0.0):
        self.position = position
        self.mass = mass
        self.energy = energy                     #energia potencjalna, kinetyczna jest wyliczana
        self.force = np.zeros(Atoms.dim)
        self.velocity = np.zeros(Atoms.dim)
    
    def getKinEnergy(self):
        return 0.5*self.mass*np.linalg.norm(self.velocity)**2

class RandomSample(object):

    @staticmethod 
    def getSample(mu, sigma, dim):                  #losowanie polozen poczatkowych
        return 1.5*np.random.normal(mu, sigma, dim)#to juz zwraca wektor
    
    @staticmethod
    def getSampleSet(mu, sigma, size, dim):
        return [RandomSample.getSample(mu, sigma, dim) for i in range(size)]

class Atoms(object):            # zbior czasteczek 
    
    dim = 1                    

    def __init__(self, size, atoms=[]):            #size - l.czasteczek, dim-wymiar
        self.size = size
        """self.atoms = [Atom(vector) for vector in \
                RandomSample.getSampleSet(0, 10, self.size, self.dim)]"""
        self.atoms = [Atom(np.array([i])) for i in range(-self.size/2,self.size/2)]
        #self.atoms = [Atom(np.array([i/10.0])) for i in range(-self.size, self.size, 2)] # do testu wykresu mbm

    def resetFAndE(self):
        for atom in self.atoms:
            self.force = np.zeros(self.dim)
            self.energy = 0


class ForceField(object):
    """ 'Abstract' class"""

    def singleForce(self, atom):
        raise NotImplementedError
    
    def singleEnergy(self, atom):
        raise NotImplementedError
    
    def pairForce(self, atom1, atom2):
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
            

    def singleEnergy(self, atom):
        distance = np.linalg.norm(atom.position)
        if distance < self.L:
            atom.energy = 0
        else:
            atom.energy = 0.5*self.f*(self.L - distance)**2
    
    def pairForce(self, atom1, atom2):
        pass
    
    def pairEnergy(self, atom1, atom2):
        pass


class MBM(ForceField):
    
    a, b, c, d = 5.0, 10.0, 3.0, 0.02     #parametry ze skryptu
    
    def singleForce(self, atom):
        x = atom.position
        atom.force = -self.a*math.exp(-self.b*(x-1)**2)*2*self.b*(x-1) \
                - self.c*math.exp(-(x+1)**2)*2*(1+x) \
                - (x**3)*4*self.d
    
    def singleEnergy(self, atom):
        x = atom.position[0]
        atom.energy = -self.a*math.exp(-self.b*(x-1)**2) \
                - self.c*math.exp(-(x+1)**2) + self.d*x**4

    def pairForce(self, atom1, atom2):
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
        return direction/distance, a
    
    def pairForce(self, atom1, atom2):
        normalized, a = self.setParams(atom1, atom2)
        F = -12.0*self.e*a*(a-1)*normalized
        atom1.force += F
        atom2.force -= F

    def pairEnergy(self, atom1, atom2):
        normalized, a = self.setParams(atom1, atom2)
        E = self.e*a*(a-2)
        atom1.energy, atom2.energy = E/2, E/2


class Simulation(object):

    def __init__(self, potential, integration, no_molecules, noSteps, stepSize):
        self.potential = potential
        self.integration = integration
        self.no_molecules = int(no_molecules)
        self.noSteps = int(noSteps)
        self.stepSize = float(stepSize)

    def start(self):                     # w tym bedzie wmieszane juz pisanie do pliku (funkcja do podzielenia na mniejsze) 
        energy = open('energy.csv', 'w') 
        trajectory = open('trajectory.xyz', 'w') 
        
        system = Atoms(self.no_molecules)

        #######################################
        """test wykresu potencjalu mbm z ukladaniem atomow na prostej co 0.2 tylko na potrzeby wykresu,
        w celu takiego ulozenia nalezy zmienic self.atoms w klasie Atoms """
        energies, a, system2 = [], MBM(), Atoms(40)
        positions = [i.position for i in system2.atoms]
        for i in system2.atoms:
            a.singleEnergy(i)
            energies.append(i.energy)
        plt.plot(positions, energies)
        plt.ylim([-6,4])
        #plt.savefig("potencjalmbm.svg")
        plt.close()

        ############################################
        
        prevPos = [x.position for x in system.atoms]
        prevVel = prevFor = np.zeros((self.no_molecules,Atoms.dim))
        previous = zip(prevPos, prevVel, prevFor)        #dla pojedynczego atomu previous[index] to 3-elementowa krotka
        verlet = {'0': BaseVerlet(), '1': VelocityVerlet(), '2': LeapFrog()}.get(self.integration)
        potential = {'0': SoftWalls(), '1': MBM(), '2': LenardJones()}.get(self.potential)
        energy_result=[]
        
        for step in range(self.noSteps):
            system.resetFAndE()
            trajectory.write(str(self.noSteps)+'\nkomentarz\n')
            totalPotEnergy = totalKinEnergy = 0
            for i in range(self.no_molecules):
                potential.singleForce(system.atoms[i])
                potential.singleEnergy(system.atoms[i])
                for atom in system.atoms[i+1:]:
                    potential.pairForce(system.atoms[i], atom)
                    potential.pairEnergy(system.atoms[i], atom)
                
                totalPotEnergy += system.atoms[i].energy
                previous[i] = verlet.step(system.atoms[i], previous[i], self.stepSize) #od razu sie ustawia nowe previous dla tego atomu
                totalKinEnergy += system.atoms[i].getKinEnergy()
                trajectory.write(str(i)+'\t'+str(system.atoms[i].position[0])+'\t0.000\t0.000\n')
            avPotEnergy = 1.0*totalPotEnergy/self.no_molecules
            avKinEnergy = 1.0*totalKinEnergy/self.no_molecules
            energy_result.append(avPotEnergy+avKinEnergy)
            energy.write(str(avPotEnergy)+'\t'+str(avKinEnergy)+'\t'+str(avPotEnergy+avKinEnergy)+'\n')
        
        plt.plot(energy_result)
        plt.savefig("energia_calkowita.svg")
        plt.close()
        energy.close()
        trajectory.close()



class Integration(object):
    """ Abstract """
    
    def step(self, atom, previous, stepSize):                      #nastepny krok dla pojedynczego atomu i pojedynczego kroku 
        raise NotImplementedError


class BaseVerlet(Integration):
    
    def step(self, atom, previous, stepSize):
        current = (atom.position, atom.velocity, atom.force)
        atom.position = 2*atom.position-previous[0]+(stepSize**2)*atom.force/atom.mass 
        atom.velocity = (atom.position - current[0])/stepSize
        return current
    
    
class VelocityVerlet(Integration):
    
    def step(self, atom, previous, stepSize):
        current = (atom.position, atom.velocity, atom.force)
        atom.velocity = atom.velocity+stepSize*0.5*(previous[2]/atom.mass - atom.force/atom.mass)
        atom.position = atom.position+atom.velocity*stepSize+(stepSize**2)*0.5*atom.force/atom.mass
        return current


class LeapFrog(Integration):
    
    def step(self, atom, previous, stepSize):
        current = (atom.position, atom.velocity, atom.force)
        atom.velocity = atom.velocity+stepSize*atom.force/atom.mass
        atom.position = atom.position+atom.velocity*stepSize
        return current
    

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
        potential, integration, no_molecules, noSteps, stepSize = args[1:] #slaaabe 
        available = range(3)
        if int(potential) not in available or int(integration) not in available:
            print 'Incorrect potential/integration. \n', help()
            return 0
        else:
            simulation = Simulation(potential, integration, no_molecules, noSteps, stepSize)              #i tak korzystam w verlecie z globalnej wartosci stepSize... czyli tej co powyzej
            simulation.start()
                                                 
if __name__ == '__main__':
    sys.exit(main(*sys.argv))
"""
a=Atoms(5)
for i in a.atoms:
    print i.position.array()
wersja ustalona jest chyba taka, ze liczymy sily tak ze jak 5 atomow i miekkie scianki
to licze miekkie scianki tylko, a jak lj to robie pary i licze tylko dla par"""
"""Klasa wektor zostala zastapiona numpy.array - nie trzeba kombinowac z dzieleniem etc"""
